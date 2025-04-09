// Changes for Rlibeemd:
// Include extras.h
// move typedef for libeemd_error_code to extras.h
// Removed unnecessary functions
//   emd_report_if_error
//   emd_report_to_file_if_error
// Moved bemd to bemd.h 

#include "extras.h"

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

#ifndef _EEMD_H_
#define _EEMD_H_

#ifndef EEMD_DEBUG
#define EEMD_DEBUG 0
#endif

#if EEMD_DEBUG == 0
#ifndef NDEBUG
#define NDEBUG
#endif
#endif

#include <stddef.h>
//#include <complex.h>
#include <stdbool.h>
// No need for this in Rlibeemd
//#include <stdio.h>
// ditto
// extern const char* libeemd_version;


//*** Removed in Rlibeemd ***//
// // Helper functions to print an error message if an error occured
// void emd_report_if_error(libeemd_error_code err);
// void emd_report_to_file_if_error(FILE* file, libeemd_error_code err);

// Main EEMD decomposition routine as described in:
//   Z. Wu and N. Huang,
//   Ensemble Empirical Mode Decomposition: A Noise-Assisted Data Analysis
//   Method, Advances in Adaptive Data Analysis,
//   Vol. 1, No. 1 (2009) 1–41
//
// Parameters 'input' and 'N' denote the input data and its length,
// respectively. Output from the routine is written to array 'output', which
// needs to be able to store at least N*M doubles, where M is the number of
// Intrinsic Mode Functions (IMFs) to compute. If M is set to zero, a value of
// M = emd_num_imfs(N) will be used, which corresponds to a maximal number of
// IMFs. Note that the final residual is also counted as an IMF in this
// respect, so you most likely want at least num_imfs=2. The following
// parameters are the ensemble size and the relative noise standard deviation,
// respectively. These are followed by the parameters for the stopping
// criterion. The stopping parameter can be defined by a S-number (see the
// article for details) or a fixed number of siftings. If both are specified,
// the sifting ends when either criterion is fulfilled. The final parameter is
// the seed given to the random number generator. A value of zero denotes a
// RNG-specific default value.
//
// To compute the original EMD decomposition you can use this function with
// ensemble_size = 1 and noise_strength = 0.
libeemd_error_code eemd(double const* __restrict input, size_t N,
		double* __restrict output, size_t M,
		unsigned int ensemble_size, double noise_strength, unsigned int
		S_number, unsigned int num_siftings, unsigned long int rng_seed, int threads);

// A complete variant of EEMD as described in:
//   M. Torres et al,
//   A Complete Ensemble Empirical Mode Decomposition with Adaptive Noise
//   IEEE Int. Conf. on Acoust., Speech and Signal Proc. ICASSP-11,
//   (2011) 4144-4147
//
// Parameters are identical to routine eemd
libeemd_error_code ceemdan(double const* __restrict input, size_t N,
		double* __restrict output, size_t M,
		unsigned int ensemble_size, double noise_strength, unsigned int
		S_number, unsigned int num_siftings, unsigned long int rng_seed, int threads);

// A method for finding the local minima and maxima from input data specified
// with parameters x and N. The memory for storing the coordinates of the
// extrema and their number are passed as the rest of the parameters. The
// arrays for the coordinates must be at least size N. The method also checks whether
// found minima are negative and maxima are positive, and returns this as boolean value.
bool emd_find_extrema(double const* __restrict x, size_t N,
		double* __restrict maxx, double* __restrict maxy, size_t* num_max_ptr,
		double* __restrict minx, double* __restrict miny, size_t* num_min_ptr);

// Return the number of IMFs that can be extracted from input data of length N,
// including the final residual.
size_t emd_num_imfs(size_t N);

// This routine evaluates a cubic spline with nodes defined by the arrays x and
// y, each of length N. The spline is evaluated using the not-a-node end point
// conditions (same as Matlab). The y values of the spline curve will be
// evaluated at integer points from 0 to x[N-1], and these y values will be
// written to the array spline_y. The endpoint x[N-1] is assumed to be an
// integer, and the x values are assumed to be in ascending order, with x[0]
// equal to 0. The workspace required is 5*N-10 doubles, except that N==2
// requires no extra memory. For N<=3 the routine falls back to polynomial
// interpolation, same as Matlab.
//
// This routine is mainly exported so that it can be tested separately to
// produce identical results to the Matlab routine 'spline'.
libeemd_error_code emd_evaluate_spline(double const* __restrict x, double const* __restrict y,
		size_t N, double* __restrict spline_y, double* spline_workspace);

#endif // _EEMD_H_
