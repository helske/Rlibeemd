#include <Rcpp.h>

extern "C"
{
  #include "eemd.h"
}

using namespace Rcpp;

void printError(libeemd_error_code err){
switch (err) {
  	case EMD_INVALID_ENSEMBLE_SIZE :
			stop("Invalid ensemble size (zero or negative)");
		case EMD_INVALID_NOISE_STRENGTH :
			stop("Invalid noise strength (negative)");
		case EMD_NOISE_ADDED_TO_EMD :
			stop("Positive noise strength but ensemble size is one (regular EMD)");
		case EMD_NO_NOISE_ADDED_TO_EEMD :
			stop("Ensemble size is more than one (EEMD) but noise strength is zero");
		case EMD_NO_CONVERGENCE_POSSIBLE :
			stop("Stopping criteria invalid: would never converge");
		case EMD_NOT_ENOUGH_POINTS_FOR_SPLINE :
			stop("Spline evaluation tried with insufficient points");
		case EMD_INVALID_SPLINE_POINTS :
			stop("Spline evaluation points invalid");
		case EMD_GSL_ERROR :
			stop("Error reported by GSL library");
    case EMD_NO_CONVERGENCE_IN_SIFTING :
      stop("Convergence not reached after sifting 10000 times");
		default :
			stop("Error code with unknown meaning. Please file a bug!");
	}
}
