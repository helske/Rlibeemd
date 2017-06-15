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

// Common error reporting and validation routines

#include "error.h"

libeemd_error_code validate_eemd_parameters(unsigned int ensemble_size, double noise_strength, unsigned int S_number, unsigned int num_siftings) {
	if (ensemble_size < 1) {
		return EMD_INVALID_ENSEMBLE_SIZE;
	}
	if (noise_strength < 0) {
		return EMD_INVALID_NOISE_STRENGTH;
	}
	if (ensemble_size == 1 && noise_strength > 0) {
		return EMD_NOISE_ADDED_TO_EMD;
	}
	if (ensemble_size > 1 && noise_strength == 0) {
		return EMD_NO_NOISE_ADDED_TO_EEMD;
	}
	if (S_number == 0 && num_siftings == 0) {
		return EMD_NO_CONVERGENCE_POSSIBLE;
	}
	return EMD_SUCCESS;
}

//*** Removed in Rlibeemd ***//

/*
// Helper functions for printing what error codes mean
void emd_report_to_file_if_error(FILE* file, libeemd_error_code err) {
	if (err == EMD_SUCCESS) {
		return;
	}
	fprintf(file, "libeemd error: ");
	switch (err) {
		case EMD_INVALID_ENSEMBLE_SIZE :
			fprintf(file, "Invalid ensemble size (zero or negative)\n");
			break;
		case EMD_INVALID_NOISE_STRENGTH :
			fprintf(file, "Invalid noise strength (negative)\n");
			break;
		case EMD_NOISE_ADDED_TO_EMD :
			fprintf(file, "Positive noise strength but ensemble size is one (regular EMD)\n");
			break;
		case EMD_NO_NOISE_ADDED_TO_EEMD :
			fprintf(file, "Ensemble size is more than one (EEMD) but noise strength is zero\n");
			break;
		case EMD_NO_CONVERGENCE_POSSIBLE :
			fprintf(file, "Stopping criteria invalid: would never converge\n");
			break;
		case EMD_NOT_ENOUGH_POINTS_FOR_SPLINE :
			fprintf(file, "Spline evaluation tried with insufficient points\n");
			break;
		case EMD_INVALID_SPLINE_POINTS :
			fprintf(file, "Spline evaluation points invalid\n");
			break;
		case EMD_GSL_ERROR :
			fprintf(file, "Error reported by GSL library\n");
			break;
		default :
			fprintf(file, "Error code with unknown meaning. Please file a bug!\n");
	}
}
void emd_report_if_error(libeemd_error_code err) {
	emd_report_to_file_if_error(stderr, err);
}
*/