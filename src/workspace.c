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

#include "workspace.h"

// sifting_workspace

sifting_workspace* allocate_sifting_workspace(size_t N) {
	sifting_workspace* w = malloc(sizeof(sifting_workspace));
	w->N = N;
	w->maxx = malloc(N*sizeof(double));
	w->maxy = malloc(N*sizeof(double));
	w->minx = malloc(N*sizeof(double));
	w->miny = malloc(N*sizeof(double));
	w->maxspline = malloc(N*sizeof(double));
	w->minspline = malloc(N*sizeof(double));
	// Spline evaluation requires 5*m-10 doubles where m is the number of
	// extrema. The worst case scenario is that every point is an extrema, so
	// use m=N to be safe.
	const size_t spline_workspace_size = (N > 2)? 5*N-10 : 0;
	w->spline_workspace = malloc(spline_workspace_size*sizeof(double));
	return w;
}

void free_sifting_workspace(sifting_workspace* w) {
	free(w->spline_workspace); w->spline_workspace = NULL;
	free(w->minspline); w->minspline = NULL;
	free(w->maxspline); w->maxspline = NULL;
	free(w->miny); w->miny = NULL;
	free(w->minx); w->minx = NULL;
	free(w->maxy); w->maxy = NULL;
	free(w->maxx); w->maxx = NULL;
	free(w); w = NULL;
}

// emd_workspace

emd_workspace* allocate_emd_workspace(size_t N) {
	emd_workspace* w = malloc(sizeof(emd_workspace));
	w->N = N;
	w->res = malloc(N*sizeof(double));
	w->sift_w = allocate_sifting_workspace(N);
	w->locks = NULL; // The locks are assumed to be allocated and freed independently
	return w;
}

void free_emd_workspace(emd_workspace* w) {
	free_sifting_workspace(w->sift_w);
	free(w->res); w->res = NULL;
	free(w); w = NULL;
}

// emd_workspace

eemd_workspace* allocate_eemd_workspace(size_t N) {
	eemd_workspace* w = malloc(sizeof(eemd_workspace));
	w->N = N;
	w->r = gsl_rng_alloc(gsl_rng_mt19937);
	w->x = malloc(N*sizeof(double));
	w->emd_w = allocate_emd_workspace(N);
	return w;
}

void set_rng_seed(eemd_workspace* w, unsigned long int rng_seed) {
	gsl_rng_set(w->r, rng_seed);
}

void free_eemd_workspace(eemd_workspace* w) {
	free_emd_workspace(w->emd_w);
	free(w->x); w->x = NULL;
	gsl_rng_free(w->r); w->r = NULL;
	free(w); w = NULL;
}
