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

#ifndef _EEMD_EXTREMA_H_
#define _EEMD_EXTREMA_H_

#include <stddef.h>
#include <assert.h>

#include "eemd.h"

// Helper function for extrapolating data at the ends. For a line passing
// through (x0, y0), (x1, y1), and (x, y), return y for a given x.
inline double linear_extrapolate(double x0, double y0,
		double x1, double y1, double x) {
	assert(x1 != x0);
	return y0 + (y1-y0)*(x-x0)/(x1-x0);
}

#endif // _EEMD_EXTREMA_H_
