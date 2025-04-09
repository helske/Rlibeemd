// Changes for Rlibeemd:
// Moved complex versions to array_complex.h
// defined restrict as __restrict

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

#ifndef _EEMD_ARRAY_COMPLEX_H_
#define _EEMD_ARRAY_COMPLEX_H_

#include <string.h>
#include <complex.h>
#include <stddef.h>
#include <R_ext/Complex.h>

// Versions for complex-valued arrays
// Original:
// static inline void complex_array_copy(double _Complex const* __restrict src, size_t n, double _Complex* __restrict dest) {
//   memcpy(dest, src, n*sizeof(double _Complex));
// }
// 
// static inline void complex_array_sub(double _Complex const* src, size_t n, double _Complex* dest) {
//   for (size_t i=0; i<n; i++)
//     dest[i] -= src[i];
// }
// 
// static inline void complex_array_mult(double _Complex* dest, size_t n, double val) {
//   for (size_t i=0; i<n; i++)
//     dest[i] *= val;
// }
// Using Rcomplex:
static inline void complex_array_copy(const Rcomplex* src, size_t n, Rcomplex* dest) {
  memcpy(dest, src, n * sizeof(Rcomplex));
}

static inline void complex_array_sub(const Rcomplex* src, size_t n, Rcomplex* dest) {
  for (size_t i = 0; i < n; i++) {
    dest[i].r -= src[i].r;
    dest[i].i -= src[i].i;
  }
}

static inline void complex_array_mult(Rcomplex* dest, size_t n, double val) {
  for (size_t i = 0; i < n; i++) {
    dest[i].r *= val;
    dest[i].i *= val;
  }
}
#endif // _EEMD_ARRAY_COMPLEX_H_
