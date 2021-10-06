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

#include "spline.h"

libeemd_error_code emd_evaluate_spline(double const* restrict x, double const* restrict y,
		size_t N, double* restrict spline_y, double* restrict spline_workspace) {
	gsl_set_error_handler_off();
	const size_t n = N-1;
	const size_t max_j = (size_t)x[n];
	if (N <= 1) {
		return EMD_NOT_ENOUGH_POINTS_FOR_SPLINE;
	}
	// perform more assertions only if EEMD_DEBUG is on,
	// as this function is meant only for internal use
	#if EEMD_DEBUG >= 1
	if (x[0] != 0) {
		return EMD_INVALID_SPLINE_POINTS;
	}
	for (size_t i=1; i<N; i++) {
		if (x[i] <= x[i-1]) {
			return EMD_INVALID_SPLINE_POINTS;
		}
	}
	#endif
	// Fall back to linear interpolation (for N==2) or polynomial interpolation
	// (for N==3)
	if (N <= 3) {
		int gsl_status = gsl_poly_dd_init(spline_workspace, x, y, N);
		if (gsl_status != GSL_SUCCESS) {
			REprintf("Error reported by gsl_poly_dd_init: %s\n",
				gsl_strerror(gsl_status));
			return EMD_GSL_ERROR;
		}
		for (size_t j=0; j<=max_j; j++) {
		  spline_y[j] = gsl_poly_dd_eval(spline_workspace, x, N, (double)j);
		}
		return EMD_SUCCESS;
	}
	// For N >= 4, interpolate by using cubic splines with not-a-node end conditions.
	// This algorithm is described in "Numerical Algorithms with C" by
	// G. Engeln-Müllges and F. Uhlig, page 257.
	//
	// Extra homework assignment for anyone reading this: Implement this
	// algorithm in GSL, so that next time someone needs these end conditions
	// they can just use GSL.
	const size_t sys_size = N-2;
	double* const c = spline_workspace;
	double* const diag = c+N;
	double* const supdiag = diag + sys_size;
	double* const subdiag = supdiag + (sys_size-1);
	double* const g = subdiag + (sys_size-1);
	// Define some constants for easier comparison with Engeln-Mullges & Uhlig
	// and let the compiler optimize them away.
	const double h_0 = x[1]-x[0];
	const double h_1 = x[2]-x[1];
	const double h_nm1 = x[n]-x[n-1];
	const double h_nm2 = x[n-1]-x[n-2];
	// Describe the (N-2)x(N-2) linear system Ac=g with the tridiagonal
	// matrix A defined by subdiag, diag and supdiag
	// first row
	diag[0] = h_0 + 2*h_1;
	supdiag[0] = h_1 - h_0;
	g[0] = 3.0/(h_0 + h_1)*((y[2]-y[1]) - (h_1/h_0)*(y[1]-y[0]));
	// rows 2 to n-2
	for (size_t i=2; i<=n-2; i++) {
		const double h_i = x[i+1] - x[i];
		const double h_im1 = x[i] - x[i-1];

		subdiag[i-2] = h_im1;
		diag[i-1] = 2*(h_im1 + h_i);
		supdiag[i-1] = h_i;
		g[i-1] = 3.0*((y[i+1]-y[i])/h_i - (y[i]-y[i-1])/h_im1);
	}
	// final row
	subdiag[n-3] = h_nm2 - h_nm1;
	diag[n-2] = 2*h_nm2 + h_nm1;
	g[n-2] = 3.0/(h_nm1 + h_nm2)*((h_nm2/h_nm1)*(y[n]-y[n-1]) - (y[n-1]-y[n-2]));
	// Solve to get c_1 ... c_{n-1}
	gsl_vector_view diag_vec = gsl_vector_view_array(diag, n-1);
	gsl_vector_view supdiag_vec = gsl_vector_view_array(supdiag, n-2);
	gsl_vector_view subdiag_vec = gsl_vector_view_array(subdiag, n-2);
	gsl_vector_view g_vec = gsl_vector_view_array(g, n-1);
	gsl_vector_view solution_vec = gsl_vector_view_array(c+1, n-1);
	int gsl_status = gsl_linalg_solve_tridiag(&diag_vec.vector,
			                                    &supdiag_vec.vector,
												&subdiag_vec.vector,
												&g_vec.vector,
												&solution_vec.vector);
	if (gsl_status != GSL_SUCCESS) {
	  REprintf("Error reported by gsl_linalg_solve_tridiag: %s\n",
				gsl_strerror(gsl_status));
		return EMD_GSL_ERROR;
	}
	// Compute c[0] and c[n]
	c[0] = c[1] + (h_0/h_1)*(c[1]-c[2]);
	c[n] = c[n-1] + (h_nm1/h_nm2)*(c[n-1]-c[n-2]);
	// The coefficients b_i and d_i are computed from the c_i's, so just
	// evaluate the spline at the required points. In this case it is easy to
	// find the required interval for spline evaluation, since the evaluation
	// points j just increase monotonically from 0 to max_j.
	size_t i = 0;
	for (size_t j=0; j<=max_j; j++) {
		if (j > x[i+1]) {
			i++;
			assert(i < n);
		}
		const double dx = (double)j-x[i];
		if (dx == 0) {
			spline_y[j] = y[i];
			continue;
		}
		// Compute coefficients b_i and d_i
		const double h_i = x[i+1] - x[i];
		const double a_i = y[i];
		const double b_i = (y[i+1]-y[i])/h_i - (h_i/3.0)*(c[i+1]+2*c[i]);
		const double c_i = c[i];
		const double d_i = (c[i+1]-c[i])/(3.0*h_i);
		// evaluate spline at x=j using the Horner scheme
		spline_y[j] = a_i + dx*(b_i + dx*(c_i + dx*d_i));
	}
	return EMD_SUCCESS;
}
