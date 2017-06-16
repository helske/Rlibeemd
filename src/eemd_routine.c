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

#include "eemd_routine.h"

// Main EEMD decomposition routine definition
libeemd_error_code eemd(double const* restrict input, size_t N,
		double* restrict output, size_t M,
		unsigned int ensemble_size, double noise_strength, unsigned int
		S_number, unsigned int num_siftings, unsigned long int rng_seed) {
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
	if (M == 0) {
		M = emd_num_imfs(N);
	}
	// The noise standard deviation is noise_strength times the standard deviation of input data
	const double noise_sigma = (noise_strength != 0)? gsl_stats_sd(input, 1, N)*noise_strength : 0;
	// Initialize output data to zero
	memset(output, 0x00, M*N*sizeof(double));
	// Each thread gets a separate workspace if we are using OpenMP
	eemd_workspace** ws = NULL;
	// The locks are shared among all threads
	lock** locks;
	// Don't start unnecessary threads if the ensemble is small
	#ifdef _OPENMP
	if (omp_get_num_threads() > (int)ensemble_size) {
	  omp_set_num_threads((int)ensemble_size);
	}
	#endif
	unsigned int ensemble_counter = 0;
	// The following section is executed in parallel
	libeemd_error_code emd_err = EMD_SUCCESS;
	#pragma omp parallel
	{
		#ifdef _OPENMP
	  const size_t num_threads = (size_t)omp_get_num_threads();
	  const size_t thread_id = (size_t)omp_get_thread_num();
		#if EEMD_DEBUG >= 1
		#pragma omp single
		REprintf("Using %d thread(s) with OpenMP.\n", num_threads);
		#endif
		#else
		const size_t num_threads = 1;
		const size_t thread_id = 0;
		#endif
		#pragma omp single
		{
			ws = malloc(num_threads*sizeof(eemd_workspace*));
			locks = malloc(M*sizeof(lock*));
			for (size_t i=0; i<M; i++) {
				locks[i] = malloc(sizeof(lock));
				init_lock(locks[i]);
			}
		}
		// Each thread allocates its own workspace
		ws[thread_id] = allocate_eemd_workspace(N);
		eemd_workspace* w = ws[thread_id];
		// All threads share the same array of locks
		w->emd_w->locks = locks;
		// Loop over all ensemble members, dividing them among the threads
		#pragma omp for
		for (size_t en_i=0; en_i<ensemble_size; en_i++) {
			// Check if an error has occured in other threads
			#pragma omp flush(emd_err)
			if (emd_err != EMD_SUCCESS) {
				continue;
			}
			// Initialize ensemble member as input data + noise
			if (noise_strength == 0.0) {
				array_copy(input, N, w->x);
			}
			else {
				// set rng seed based on ensemble member to ensure
				// reproducibility even in a multithreaded case
				set_rng_seed(w, rng_seed+en_i);
				for (size_t i=0; i<N; i++) {
					w->x[i] = input[i] + gsl_ran_gaussian(w->r, noise_sigma);
				}
			}
			// Extract IMFs with EMD
			emd_err = _emd(w->x, w->emd_w, output, M, S_number, num_siftings);
			#pragma omp flush(emd_err)
			#pragma omp atomic
			ensemble_counter++;
			#if EEMD_DEBUG >= 1
			REprintf("Ensemble iteration %u/%u done.\n", ensemble_counter, ensemble_size);
			#endif
		}
		// Free resources
		free_eemd_workspace(w);
		#pragma omp single
		{
			free(ws); ws = NULL;
			for (size_t i=0; i<M; i++) {
				destroy_lock(locks[i]);
				free(locks[i]);
			}
			free(locks); locks = NULL;
		}
	} // End of parallel block
	if (emd_err != EMD_SUCCESS) {
		return emd_err;
	}
	// Divide output data by the ensemble size to get the average
	if (ensemble_size != 1) {
		const double one_per_ensemble_size = 1.0/ensemble_size;
		array_mult(output, N*M, one_per_ensemble_size);
	}
	return EMD_SUCCESS;
}
