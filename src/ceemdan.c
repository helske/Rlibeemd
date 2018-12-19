/* Copyright 2013 Perttu Luukko

 * This file is part of libeemd.

 * libeemd is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.

 * libeemd is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with libeemd.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "ceemdan.h"

// Main CEEMDAN decomposition routine definition
libeemd_error_code ceemdan(double const* restrict input, size_t N,
		double* restrict output, size_t M,
		unsigned int ensemble_size, double noise_strength, unsigned int
		S_number, unsigned int num_siftings, unsigned long int rng_seed, int threads) {
	gsl_set_error_handler_off();
	// Validate parameters
	libeemd_error_code validation_result = validate_eemd_parameters(ensemble_size, noise_strength, S_number, num_siftings);
	if (validation_result != EMD_SUCCESS) {
		return validation_result;
	}
	// For empty data we have nothing to do
	if (N == 0) {
		return EMD_SUCCESS;
	}
	// For M == 1 the only "IMF" is the residual
	if (M == 1) {
		memcpy(output, input, N*sizeof(double));
		return EMD_SUCCESS;
	}
	if (M == 0) {
		M = emd_num_imfs(N);
	}
	const double one_per_ensemble_size = 1.0/ensemble_size;
	// Initialize output data to zero
	memset(output, 0x00, M*N*sizeof(double));
	// Each thread gets a separate workspace if we are using OpenMP
	eemd_workspace** ws = NULL;
	// All threads need to write to the same row of the output matrix
	// so we need only one shared lock
	lock* output_lock = malloc(sizeof(lock));
	init_lock(output_lock);
	// The threads also share the same precomputed noise
	double* noises = malloc(ensemble_size*N*sizeof(double));
	// Since we need to decompose this noise by EMD, we also need arrays for storing
	// the residuals
	double* noise_residuals = malloc(ensemble_size*N*sizeof(double));
	// Don't start unnecessary threads if the ensemble is small
	#ifdef _OPENMP
	int old_maxthreads = 1;
	if (threads>0) {
	  old_maxthreads = omp_get_max_threads();
	  omp_set_num_threads(threads);    
	}
	if (omp_get_num_threads() > (int)ensemble_size) {
	  omp_set_num_threads((int)ensemble_size);
	}
	#endif
	size_t num_threads;
	// The following section is executed in parallel
	#pragma omp parallel
	{
		#ifdef _OPENMP
	  num_threads = (size_t)omp_get_num_threads();
	  const size_t thread_id = (size_t)omp_get_thread_num();
		#if EEMD_DEBUG >= 1
		#pragma omp single
		REprintf("Using %d thread(s) with OpenMP.\n", num_threads);
		#endif
		#else
		num_threads = 1;
		const size_t thread_id = 0;
		#endif
		#pragma omp single
		{
			ws = malloc(num_threads*sizeof(eemd_workspace*));
		}
		// Each thread allocates its own workspace
		ws[thread_id] = allocate_eemd_workspace(N);
		eemd_workspace* w = ws[thread_id];
		// Precompute and store white noise, since for each mode of the data we
		// need the same mode of the corresponding realization of noise
		#pragma omp for
		for (size_t en_i=0; en_i<ensemble_size; en_i++) {
			// set rng seed based on ensemble member to ensure
			// reproducibility even in a multithreaded case
			set_rng_seed(w, rng_seed+en_i);
			for (size_t j=0; j<N; j++) {
				noises[N*en_i+j] = gsl_ran_gaussian(w->r, 1.0);
			}
		}
	} // Return to sequental mode
	// Allocate memory for the residual shared among all threads
	double* restrict res = malloc(N*sizeof(double));
	// For the first iteration the residual is the input signal
	array_copy(input, N, res);
	// Each mode is extracted sequentially, but we use parallelization in the inner loop
	// to loop over ensemble members
	for (size_t imf_i=0; imf_i<M; imf_i++) {
		// Provide a pointer to the output vector where this IMF will be stored
		double* const imf = &output[imf_i*N];
		// Then we go parallel to compute the different ensemble members
		libeemd_error_code sift_err = EMD_SUCCESS;
		#pragma omp parallel
		{
			#ifdef _OPENMP
			const int thread_id = omp_get_thread_num();
			#else
			const int thread_id = 0;
			#endif
			eemd_workspace* w = ws[thread_id];
			unsigned int sift_counter = 0;
			#pragma omp for
			for (size_t en_i=0; en_i<ensemble_size; en_i++) {
				// Check if an error has occured in other threads
				#pragma omp flush(sift_err)
				if (sift_err != EMD_SUCCESS) {
					continue;
				}
				// Provide a pointer to the noise vector and noise residual used by
				// this ensemble member
				double* const noise = &noises[N*en_i];
				double* const noise_residual = &noise_residuals[N*en_i];
				// Initialize input signal as data + noise.
				// The noise standard deviation is noise_strength times the
				// standard deviation of input data divided by the standard
				// deviation of the noise. This is used to fix the SNR at each
				// stage.
				const double noise_sd = gsl_stats_sd(noise, 1, N);
				const double noise_sigma = (noise_sd != 0)? noise_strength*gsl_stats_sd(res, 1, N)/noise_sd : 0;
				array_addmul_to(res, noise, noise_sigma, N, w->x);
				// Sift to extract first EMD mode
				sift_err = _sift(w->x, w->emd_w->sift_w, S_number, num_siftings, &sift_counter);
				#pragma omp flush(sift_err)
				// Sum to output vector
				get_lock(output_lock);
				array_add(w->x, N, imf);
				release_lock(output_lock);
				// Extract next EMD mode of the noise. This is used as the noise for
				// the next mode extracted from the data
				if (imf_i == 0) {
					array_copy(noise, N, noise_residual);
				}
				else {
					array_copy(noise_residual, N, noise);
				}
				sift_err = _sift(noise, w->emd_w->sift_w, S_number, num_siftings, &sift_counter);
				#pragma omp flush(sift_err)
				array_sub(noise, N, noise_residual);
			}
		} // Parallel section ends
		if (sift_err != EMD_SUCCESS) {
			return sift_err;
		}
		// Divide with ensemble size to get the average
		array_mult(imf, N, one_per_ensemble_size);
		// Subtract this IMF from the previous residual to form the new one
		array_sub(imf, N, res);
	}
	// Save final residual
	get_lock(output_lock);
	array_add(res, N, output+N*(M-1));
	release_lock(output_lock);
	// Free global resources
	for (size_t thread_id=0; thread_id<num_threads; thread_id++) {
		free_eemd_workspace(ws[thread_id]);
	}
	free(ws); ws = NULL;
	free(res); res = NULL;
	free(noise_residuals); noise_residuals = NULL;
	free(noises); noises = NULL;
	destroy_lock(output_lock);
	free(output_lock); output_lock = NULL;
	
#ifdef _OPENMP
	if (threads>0) {
	  omp_set_num_threads(old_maxthreads);    
	}
#endif
	return EMD_SUCCESS;
}
