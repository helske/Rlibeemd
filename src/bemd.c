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

#include "bemd.h"

bemd_sifting_workspace* allocate_bemd_sifting_workspace(size_t N, lock* output_lock) {
  bemd_sifting_workspace* w = malloc(sizeof(bemd_sifting_workspace));
  w->N = N;
  w->projected_signal = malloc(N*sizeof(double));
  w->maxx = malloc(N*sizeof(double));
  w->maxy = malloc(N*sizeof(double));
  w->maxspline = malloc(N*sizeof(double));
  // Spline evaluation requires 5*m-10 doubles where m is the number of
  // extrema. The worst case scenario is that every point is an extrema, so
  // use m=N to be safe.
  const size_t spline_workspace_size = (N > 2)? 5*N-10 : 0;
  w->spline_workspace = malloc(spline_workspace_size*sizeof(double));
  w->output_lock = output_lock;
  return w;
}

void free_bemd_sifting_workspace(bemd_sifting_workspace* w) {
  free(w->projected_signal); w->projected_signal = NULL;
  free(w->maxx); w->maxx = NULL;
  free(w->maxy); w->maxy = NULL;
  free(w->maxspline); w->maxspline = NULL;
  free(w->spline_workspace); w->spline_workspace = NULL;
  free(w); w = NULL;
}

// double complex* __restrict x -> Rcomplex* x
static libeemd_error_code _bemd_sift_once(Rcomplex* x, size_t N, double const* __restrict directions, size_t num_directions, bemd_sifting_workspace* w) {
  libeemd_error_code errcode = EMD_SUCCESS;
  //double complex* __restrict m = calloc(N, sizeof(double complex));
  Rcomplex* m = calloc(N, sizeof(Rcomplex));
  double* const px = w->projected_signal;
  // TODO: handle different directions in parallel
  for (size_t direction_i=0; direction_i<num_directions; direction_i++) {
    const double phi = directions[direction_i];
    const double sin_phi = sin(phi);
    const double cos_phi = cos(phi);
    // Project signal
    for (size_t i=0; i<N; i++) {
      // const double a = creal(x[i]);
      // const double b = cimag(x[i]);
      const double a = x[i].r; // Rcomplex
      const double b = x[i].i;
      px[i] = a*cos_phi + b*sin_phi;
    }
    // Find maxima
    emd_find_maxima(px, N, w->maxx, w->maxy, &(w->num_max));
    // Fit spline
    errcode = emd_evaluate_spline(w->maxx, w->maxy, w->num_max, w->maxspline, w->spline_workspace);
    if (errcode != EMD_SUCCESS) {
      return errcode;
    }
    // Add to m
    for (size_t i=0; i<N; i++) {
      //m[i] += cexp(phi*I) * (w->maxspline)[i];
      // Rcomplex
      double amp = w->maxspline[i];
      double re_exp = cos(phi);
      double im_exp = sin(phi);
      m[i].r += amp * re_exp;
      m[i].i += amp * im_exp;
    }
  }
  // Scale m
  complex_array_mult(m, N, 2.0/(double)num_directions);
  // Subtract mean from input
  complex_array_sub(m, N, x);
  // Done
  free(m); m = NULL;
  return errcode;
}

//double _Complex const* __restrict input -> Rcomplex* input
//double _Complex* __restrict output -> const Rcomplex* output
libeemd_error_code bemd(const Rcomplex* input, size_t N,
  double const* __restrict directions, size_t num_directions,
  Rcomplex* output, size_t M,
  unsigned int num_siftings) {
  gsl_set_error_handler_off();
  if (M == 0) {
    M = emd_num_imfs(N);
  }
  libeemd_error_code bemd_err = EMD_SUCCESS;
  // Create a read-write copy of input data
  //double complex* const x = malloc(N*sizeof(double complex));
  Rcomplex* const x = malloc(N*sizeof(Rcomplex));
  complex_array_copy(input, N, x);
  //double complex* const res = malloc(N*sizeof(double complex));
  Rcomplex* const res = malloc(N*sizeof(Rcomplex));
  // For the first iteration, the residual is the original input data
  complex_array_copy(input, N, res);
  bemd_sifting_workspace* w = allocate_bemd_sifting_workspace(N, NULL);
  // Loop over all IMFs to be separated from input
  for (size_t imf_i=0; imf_i<M-1; imf_i++) {
    if (imf_i != 0) {
      // Except for the first iteration, restore the previous residual
      // and use it as an input
      complex_array_copy(res, N, x);
    }
    // Perform siftings on x until it is an IMF
    for (unsigned int sift_counter=0; sift_counter<num_siftings; sift_counter++) {
      bemd_err = _bemd_sift_once(x, N, directions, num_directions, w);
      if (bemd_err != EMD_SUCCESS) {
        return bemd_err;
      }
    }
    // Subtract this IMF from the saved copy to form the residual for
    // the next round
    complex_array_sub(x, N, res);
    // Write the discovered IMF to the output matrix
    complex_array_copy(x, N, output+N*imf_i);
  }
  // Save final residual
  complex_array_copy(res, N, output+N*(M-1));
  free_bemd_sifting_workspace(w);
  free(res);
  free(x);
  return bemd_err;
}
