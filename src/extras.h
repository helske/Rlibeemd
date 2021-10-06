/*
 ** Stuff needed for R/C++ compatibility:
 ** Added #include <R_ext/Print.h>
 ** Changed calls fprintf(stderr,...) to R compatible REprintf(...)
 ** Removed unnecessary functions
 **  emd_report_if_error
 **  emd_report_to_file_if_error
 **  
 **  
 */

#ifndef _EXTRAS_H_
#define _EXTRAS_H_

#include <R_ext/Print.h>

#ifndef restrict
#define restrict // nothing
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

// Possible error codes returned by functions eemd, ceemdan and
// emd_evaluate_spline
typedef enum {
  EMD_SUCCESS = 0,
  // Errors from invalid parameters
  EMD_INVALID_ENSEMBLE_SIZE = 1,
  EMD_INVALID_NOISE_STRENGTH = 2,
  EMD_NOISE_ADDED_TO_EMD = 3,
  EMD_NO_NOISE_ADDED_TO_EEMD = 4,
  EMD_NO_CONVERGENCE_POSSIBLE = 5,
  EMD_NOT_ENOUGH_POINTS_FOR_SPLINE = 6,
  EMD_INVALID_SPLINE_POINTS = 7,
  // Other errors
  EMD_GSL_ERROR = 8,
  EMD_NO_CONVERGENCE_IN_SIFTING = 9
} libeemd_error_code;


void printError(libeemd_error_code err);

#endif