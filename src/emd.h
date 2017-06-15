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

#ifndef _EEMD_EMD_H_
#define _EEMD_EMD_H_

#include <string.h>
#include <stddef.h>
#include <math.h>

#include "array.h"
#include "error.h"
#include "workspace.h"
#include "eemd.h"

// This file contains helper functions for doing simple EMD. They are then used
// to build more complicated routines like EEMD and CEEMDAN. Note that these
// functions use internal data structures and are not intended to be used
// outside libeemd. If you need to compute the ordinary EMD, use the public
// eemd() routine.

// Helper function for applying the sifting procedure to input until it is
// reduced to an IMF according to the stopping criteria given by S_number and
// num_siftings. The required number of siftings is saved to sift_counter.
libeemd_error_code _sift(double* restrict input, sifting_workspace*
		restrict w, unsigned int S_number, unsigned int num_siftings,
		unsigned int* sift_counter);

// Helper function for extracting all IMFs from input using the sifting
// procedure defined by _sift. The contents of the input array are destroyed in
// the process.
libeemd_error_code _emd(double* restrict input, emd_workspace* restrict w,
		double* restrict output, size_t M,
		unsigned int S_number, unsigned int num_siftings);

#endif // _EEMD_EMD_H_
