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

#include "extrema.h"

bool emd_find_extrema(double const* restrict x, size_t N,
  double* restrict maxx, double* restrict maxy, size_t* nmax,
  double* restrict minx, double* restrict miny, size_t* nmin) {
  // Set the number of extrema to zero initially
  *nmax = 0;
  *nmin = 0;
  // Handle empty array as a special case
  if (N == 0) {
    return true;
  }
  // Add the ends of the data as both local minima and maxima. These
  // might be changed later by linear extrapolation.
  maxx[0] = 0;
  maxy[0] = x[0];
  (*nmax)++;
  minx[0] = 0;
  miny[0] = x[0];
  (*nmin)++;
  // If we had only one data point this is it
  if (N == 1) {
    return true;
  }
  // Now starts the main extrema-finding loop. The loop detects points where
  // the slope of the data changes sign. In the case of flat regions at the
  // extrema, the center point of the flat region will be considered the
  // extremal point. While detecting extrema, the loop also counts the number
  // of zero crossings that occur.
  bool all_extrema_good = true;
  enum slope { UP, DOWN, NONE };
  enum slope previous_slope = NONE;
  int flat_counter = 0;
  for (size_t i=0; i<N-1; i++) {
    if (x[i+1] > x[i]) { // Going up
      if (previous_slope == DOWN) {
        // Was going down before -> local minimum found
        minx[*nmin] = (double)(i)-(double)(flat_counter)/2;
        miny[*nmin] = x[i];
        (*nmin)++;
        if (x[i] >= 0) { // minima need to be negative
          all_extrema_good = false;
        }
      }
      previous_slope = UP;
      flat_counter = 0;
    }
    
    else if (x[i+1] < x[i]) { // Going down
      if (previous_slope == UP) {
        // Was going up before -> local maximum found
        maxx[*nmax] = (double)(i)-(double)(flat_counter)/2;
        maxy[*nmax] = x[i];
        (*nmax)++;
        if (x[i] <= 0) { // maxima need to be positive
          all_extrema_good = false;
        }
      }
      previous_slope = DOWN;
      flat_counter = 0;
    }
    else { // Staying flat
      flat_counter++;
#if EEMD_DEBUG >= 3
      REprintf("Warning: a flat slope found in data. The results will differ from the reference EEMD implementation.\n");
#endif
    }
  }
  // Add the other end of the data as extrema as well.
  maxx[*nmax] = (double)(N-1);
  maxy[*nmax] = x[N-1];
  (*nmax)++;
  minx[*nmin] = (double)(N-1);
  miny[*nmin] = x[N-1];
  (*nmin)++;
  // If we have at least two interior extrema, test if linear extrapolation provides
  // a more extremal value.
  if (*nmax >= 4) {
    const double max_el = linear_extrapolate(maxx[1], maxy[1],
      maxx[2], maxy[2], 0);
    if (max_el > maxy[0])
      maxy[0] = max_el;
    const double max_er = linear_extrapolate(maxx[*nmax-3], maxy[*nmax-3],
      maxx[*nmax-2], maxy[*nmax-2], (double)(N-1));
    if (max_er > maxy[*nmax-1])
      maxy[*nmax-1] = max_er;
  }
  if (*nmin >= 4) {
    const double min_el = linear_extrapolate(minx[1], miny[1],
      minx[2], miny[2], 0);
    if (min_el < miny[0])
      miny[0] = min_el;
    const double min_er = linear_extrapolate(minx[*nmin-3], miny[*nmin-3],
      minx[*nmin-2], miny[*nmin-2], (double)(N-1));
    if (min_er < miny[*nmin-1])
      miny[*nmin-1] = min_er;
  }
  return all_extrema_good;
}

void emd_find_maxima(double const* restrict x, size_t N, double* restrict maxx, double* restrict maxy, size_t* nmax) {
  // Set the number of maxima to zero initially
  *nmax = 0;
  // Handle empty array as a special case
  if (N == 0) {
    return;
  }
  // Add the ends of the data as maxima. These might be changed later by
  // linear extrapolation.
  maxx[0] = 0;
  maxy[0] = x[0];
  (*nmax)++;
  // If we had only one data point this is it
  if (N == 1) {
    return;
  }
  // Now starts the main extrema-finding loop. The loop detects points where
  // the slope of the data changes sign. In the case of flat regions at the
  // extrema, the center point of the flat region will be considered the
  // extremal point. While detecting extrema, the loop also counts the number
  // of zero crossings that occur.
  enum slope { UP, DOWN, NONE };
  enum slope previous_slope = NONE;
  int flat_counter = 0;
  for (size_t i=0; i<N-1; i++) {
    if (x[i+1] > x[i]) { // Going up
      previous_slope = UP;
      flat_counter = 0;
    }
    else if (x[i+1] < x[i]) { // Going down
      if (previous_slope == UP) {
        // Was going up before -> local maximum found
        maxx[*nmax] = (double)(i)-(double)(flat_counter)/2;
        maxy[*nmax] = x[i];
        (*nmax)++;
      }
      previous_slope = DOWN;
      flat_counter = 0;
    }
    else { // Staying flat
      flat_counter++;
#if EEMD_DEBUG >= 3
      REprintf("Warning: a flat slope found in data. The results will differ from the reference EEMD implementation.\n");
#endif
    }
  }
  // Add the other end of the data as extrema as well.
  maxx[*nmax] = (double)(N-1);
  maxy[*nmax] = x[N-1];
  (*nmax)++;
  // If we have at least two interior extrema, test if linear extrapolation provides
  // a more extremal value.
  if (*nmax >= 4) {
    const double max_el = linear_extrapolate(maxx[1], maxy[1],
      maxx[2], maxy[2], 0);
    if (max_el > maxy[0])
      maxy[0] = max_el;
    const double max_er = linear_extrapolate(maxx[*nmax-3], maxy[*nmax-3],
      maxx[*nmax-2], maxy[*nmax-2], (double)(N-1));
    if (max_er > maxy[*nmax-1])
      maxy[*nmax-1] = max_er;
  }
  return;
}
