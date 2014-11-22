#include <Rcpp.h>
#ifndef restrict
#define restrict // nothing
#endif
extern "C"
{
  #include "eemd.h"
}

using namespace Rcpp;

void printError(libeemd_error_code err){
switch (err) {
  	case EMD_INVALID_ENSEMBLE_SIZE :
			stop("Invalid ensemble size (zero or negative)");
			break;
		case EMD_INVALID_NOISE_STRENGTH :
			stop("Invalid noise strength (negative)");
			break;
		case EMD_NOISE_ADDED_TO_EMD :
			stop("Positive noise strength but ensemble size is one (regular EMD)");
			break;
		case EMD_NO_NOISE_ADDED_TO_EEMD :
			stop("Ensemble size is more than one (EEMD) but noise strength is zero");
			break;
		case EMD_NO_CONVERGENCE_POSSIBLE :
			stop("Stopping criteria invalid: would never converge");
			break;
		case EMD_NOT_ENOUGH_POINTS_FOR_SPLINE :
			stop("Spline evaluation tried with insufficient points");
			break;
		case EMD_INVALID_SPLINE_POINTS :
			stop("Spline evaluation points invalid");
			break;
		case EMD_GSL_ERROR :
			stop("Error reported by GSL library");
			break;
		default :
			stop("Error code with unknown meaning. Please file a bug!");
	}
}