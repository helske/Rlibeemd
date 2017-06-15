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

#ifndef _EEMD_LOCK_H_
#define _EEMD_LOCK_H_

// If we are using OpenMP for parallel computation, we need locks to ensure
// that the same output data is not written by several threads at the same
// time.
#ifdef _OPENMP
#include <omp.h>
typedef omp_lock_t lock;
inline void init_lock(lock* l) { omp_init_lock(l); }
inline void destroy_lock(lock* l) { omp_destroy_lock(l); }
inline void get_lock(lock* l) { omp_set_lock(l); }
inline void release_lock(lock* l) { omp_unset_lock(l); }
#else
// If we don't use OpenMP, we provide a dummy lock that does nothing. This
// avoids littering the code with too many #ifdefs for _OPENMP.
typedef char lock;
inline void init_lock(__attribute__((unused)) lock* l) {}
inline void destroy_lock(__attribute__((unused)) lock* l) {}
inline void get_lock(__attribute__((unused)) lock* l) {}
inline void release_lock(__attribute__((unused)) lock* l) {}
#endif

#endif // _EEMD_LOCK_H_
