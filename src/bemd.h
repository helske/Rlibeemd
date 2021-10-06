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

#include <stddef.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <gsl/gsl_errno.h>

#include "array.h"
#include "lock.h"
#include "extrema.h"
#include "spline.h"
#include "eemd.h"

// Bivariate EMD described as scheme 2 in:
//   G. Rilling, P. Flandrin, P. Goncalves and J. M. Lilly,
//   Bivariate Empirical Mode Decomposition
//   IEEE Signal Processing Letters, vol. 14, no. 12, pp. 936-939, Dec. 2007.
//
// Parameters 'directions' and 'num_directions' define a vector of directions (phi_k in
// the article) used for the decomposition. 
libeemd_error_code bemd(double _Complex const* restrict input, size_t N,
  double const* restrict directions, size_t num_directions,
  double _Complex* restrict output, size_t M,
  unsigned int num_siftings);

// For BEMD sifting we need arrays for storing the found maxima of the signal,
// memory required to form the spline envelopes, and a shared lock to compute
// different directions in parallel.
typedef struct {
  // Number of samples in the signal
  size_t N;
  // Input signal projected to a particular direction in the complex plane
  double* projected_signal;
  // Found maxima
  double* restrict maxx;
  double* restrict maxy;
  size_t num_max;
  // Upper and lower envelope spline values
  double* restrict maxspline;
  // Extra memory required for spline evaluation
  double* restrict spline_workspace;
  // Lock
  lock* output_lock;
} bemd_sifting_workspace;

bemd_sifting_workspace* allocate_bemd_sifting_workspace(size_t N, lock* output_lock);
void free_bemd_sifting_workspace(bemd_sifting_workspace* w);

