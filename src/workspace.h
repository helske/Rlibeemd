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

#ifndef _EEMD_WORKSPACE_H_
#define _EEMD_WORKSPACE_H_

#include <stddef.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>

#include "lock.h"

// Necessary workspace memory structures for various EMD operations

// For sifting we need arrays for storing the found extrema of the signal, and memory required
// to form the spline envelopes
typedef struct {
	// Number of samples in the signal
	size_t N;
	// Found extrema
	double* restrict maxx;
	double* restrict maxy;
	double* restrict minx;
	double* restrict miny;
	// Upper and lower envelope spline values
	double* restrict maxspline;
	double* restrict minspline;
	// Extra memory required for spline evaluation
	double* restrict spline_workspace;
} sifting_workspace;

sifting_workspace* allocate_sifting_workspace(size_t N);
void free_sifting_workspace(sifting_workspace* w);

// For EMD we need space to do the sifting and somewhere to save the residual from the previous run.
// We also leave room for an array of locks to protect multi-threaded EMD.
typedef struct {
	size_t N;
	// Previous residual for EMD
	double* restrict res;
	// What is needed for sifting
	sifting_workspace* restrict sift_w;
	// A pointer for shared locks. These locks are used to make EMD thread-safe
	// even when several threads run EMD with the same output matrix (we'll do
	// this in EEMD).
	lock** locks;
} emd_workspace;

emd_workspace* allocate_emd_workspace(size_t N);
void free_emd_workspace(emd_workspace* w);

// EEMD needs a random number generator in addition to emd_workspace. We also need a place to store
// the member of the ensemble (input signal + realization of noise) to be worked on.
typedef struct {
	size_t N;
	// The random number generator
	gsl_rng* r;
	// The ensemble member signal
	double* restrict x;
	// What is needed for running EMD
	emd_workspace* restrict emd_w;
} eemd_workspace;

eemd_workspace* allocate_eemd_workspace(size_t N);
void set_rng_seed(eemd_workspace* w, unsigned long int rng_seed);
void free_eemd_workspace(eemd_workspace* w);

#endif // _EEMD_WORKSPACE_H_
