#include <Rcpp.h>
#ifndef restrict
#define restrict // nothing
#endif
extern "C"
{
  #include "eemd.h"
}

using namespace Rcpp;
//' CEEMDAN decomposition
//' 
//' Decompose input data to Intrinsic Mode Functions (IMFs) with the
//' Complete Ensemble Empirical Mode Decomposition with Adaptive Noise (CEEMDAN)
//' algorithm, a variant of EEMD.
//'
//' The size of the ensemble and the relative magnitude of the added noise are
//' given by parameters \code{ensemble_size} and \code{noise_strength}, respectively.  The
//' stopping criterion for the decomposition is given by either a S-number or
//' an absolute number of siftings. In the case that both are positive numbers,
//' the sifting ends when either of the conditions is fulfilled.
//'
//' @export
//' @name ceemdan
//' @inheritParams eemd
//' @return Matrix of size MxN, where M = \code{emd_num_imfs}(N). The rows of the array are the
//'        IMFs of the input signal, with the last row being the final residual.
//' @references M. Torres et al, A Complete Ensemble Empirical Mode Decomposition with Adaptive Noise
//'   IEEE Int. Conf. on Acoust., Speech and Signal Proc. ICASSP-11,
//'   (2011) 4144-4147
//' @seealso \code{\link{eemd}} 
// [[Rcpp::export(ceemdan)]]
NumericMatrix ceemdanR(NumericVector input, unsigned int ensemble_size=250, 
double noise_strength=0.2, unsigned int S_number=4, unsigned int num_siftings=0, 
unsigned long int rng_seed=0){
  
  size_t N = input.size();
  NumericMatrix output(N,emd_num_imfs(N));
  ceemdan(input.begin(), N, output.begin(), ensemble_size, noise_strength, S_number, num_siftings, rng_seed);
  return output;
}
