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

#ifndef _EEMD_ARRAY_H_
#define _EEMD_ARRAY_H_

#include <string.h>

// Helper functions for working with data arrays
inline void array_copy(double const* restrict src, size_t n, double* restrict dest) {
	memcpy(dest, src, n*sizeof(double));
}

inline void array_add(double const* src, size_t n, double* dest) {
	for (size_t i=0; i<n; i++)
		dest[i] += src[i];
}

inline void array_add_to(double const* src1, double const* src2, size_t n, double* dest) {
	for (size_t i=0; i<n; i++)
		dest[i] = src1[i] + src2[i];
}

inline void array_addmul_to(double const* src1, double const* src2, double val, size_t n, double* dest) {
	for (size_t i=0; i<n; i++)
		dest[i] = src1[i] + val*src2[i];
}

inline void array_sub(double const* src, size_t n, double* dest) {
	for (size_t i=0; i<n; i++)
		dest[i] -= src[i];
}

inline void array_mult(double* dest, size_t n, double val) {
	for (size_t i=0; i<n; i++)
		dest[i] *= val;
}

#endif // _EEMD_ARRAY_H_
