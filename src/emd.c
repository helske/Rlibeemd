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

#include "emd.h"

libeemd_error_code _sift(double* restrict input, sifting_workspace*
		restrict w, unsigned int S_number, unsigned int num_siftings,
		unsigned int* sift_counter) {
	const size_t N = w->N;
	// Provide some shorthands to avoid excessive '->' operators
	double* const maxx = w->maxx;
	double* const maxy = w->maxy;
	double* const minx = w->minx;
	double* const miny = w->miny;
	// Initialize counters that keep track of the number of siftings
	// and the S number
	*sift_counter = 0;
	unsigned int S_counter = 0;
	// Numbers of extrema and zero crossings are initialized to dummy values
	size_t num_max = (size_t)(-1);
	size_t num_min = (size_t)(-1);
	size_t num_zc = (size_t)(-1);
	size_t prev_num_max = (size_t)(-1);
	size_t prev_num_min = (size_t)(-1);
	size_t prev_num_zc = (size_t)(-1);
	while (num_siftings == 0 || *sift_counter < num_siftings) {
		(*sift_counter)++;
		#if EEMD_DEBUG >= 1
		if (*sift_counter == 10000) {
		  REprintf("Something is probably wrong. Sift counter has reached 10000.\n");
		}
		#endif
		prev_num_max = num_max;
		prev_num_min = num_min;
		prev_num_zc = num_zc;
		// Find extrema and count zero crossings
		emd_find_extrema(input, N, maxx, maxy, &num_max, minx, miny, &num_min, &num_zc);
		// Check if we are finished based on the S-number criteria
		if (S_number != 0) {
			const int max_diff = (int)num_max - (int)prev_num_max;
			const int min_diff = (int)num_min - (int)prev_num_min;
			const int zc_diff = (int)num_zc - (int)prev_num_zc;
			if (abs(max_diff)+abs(min_diff)+abs(zc_diff) <= 1) {
				S_counter++;
				if (S_counter >= S_number) {
					const int num_diff = (int)num_min + (int)num_max - 4 - (int)num_zc;
					if (abs(num_diff) <= 1) {
						// Number of extrema has been stable for S_number steps
						// and the number of *interior* extrema and zero
						// crossings differ by at most one -- we are converged
						// according to the S-number criterion
						break;
					}
				}
			}
			else {
				S_counter = 0;
			}
		}
		// Fit splines, choose order of spline based on the number of extrema
		libeemd_error_code max_errcode = emd_evaluate_spline(maxx, maxy, num_max, w->maxspline, w->spline_workspace);
		if (max_errcode != EMD_SUCCESS) {
			return max_errcode;
		}
		libeemd_error_code min_errcode = emd_evaluate_spline(minx, miny, num_min, w->minspline, w->spline_workspace);
		if (min_errcode != EMD_SUCCESS) {
			return min_errcode;
		}
		// Subtract envelope mean from the data
		for (size_t i=0; i<N; i++) {
			input[i] -= 0.5*(w->maxspline[i] + w->minspline[i]);
		}
	}
	return EMD_SUCCESS;
}

libeemd_error_code _emd(double* restrict input, emd_workspace* restrict w,
		double* restrict output, size_t M,
		unsigned int S_number, unsigned int num_siftings) {
	// Provide some shorthands to avoid excessive '->' operators
	const size_t N = w->N;
	double* const res = w->res;
	lock** locks = w->locks;
	if (M == 0) {
		M = emd_num_imfs(N);
	}
	// We need to store a copy of the original signal so that once it is
	// reduced to an IMF we have something to subtract the IMF from to form
	// the residual for the next iteration
	array_copy(input, N, res);
	// Loop over all IMFs to be separated from input
	unsigned int sift_counter;
	for (size_t imf_i=0; imf_i<M-1; imf_i++) {
		if (imf_i != 0) {
			// Except for the first iteration, restore the previous residual
			// and use it as an input
			array_copy(res, N, input);
		}
		// Perform siftings on input until it is an IMF
		libeemd_error_code sift_err = _sift(input, w->sift_w, S_number, num_siftings, &sift_counter);
		if (sift_err != EMD_SUCCESS) {
			return sift_err;
		}
		// Subtract this IMF from the saved copy to form the residual for
		// the next round
		array_sub(input, N, res);
		// Add the discovered IMF to the output matrix. Use locks to ensure
		// other threads are not writing to the same row of the output matrix
		// at the same time
		get_lock(locks[imf_i]);
		array_add(input, N, output+N*imf_i);
		release_lock(locks[imf_i]);
		#if EEMD_DEBUG >= 2
		REprintf("IMF %zd saved after %u siftings.\n", imf_i+1, sift_counter);
		#endif
	}
	// Save final residual
	get_lock(locks[M-1]);
	array_add(res, N, output+N*(M-1));
	release_lock(locks[M-1]);
	return EMD_SUCCESS;
}

size_t emd_num_imfs(size_t N) {
	if (N == 0) {
		return 0;
	}
	if (N <= 3) {
		return 1;
	}
	return (size_t)(log2(N));
}
