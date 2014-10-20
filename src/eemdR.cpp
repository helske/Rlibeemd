#include <Rcpp.h>
#ifndef restrict
#define restrict // nothing
#endif
extern "C"
{
  #include "eemd.h"
}

using namespace Rcpp;

//' EEMD Decomposition
//' 
//' Decompose input data to Intrinsic Mode Functions (IMFs) with the
//' Ensemble Empirical Mode Decomposition algorithm [1].
//'
//' The size of the ensemble and the relative magnitude of the added noise are
//' given by parameters \code{ensemble_size} and \code{noise_strength}, respectively.  The
//' stopping criterion for the decomposition is given by either a S-number or
//' an absolute number of siftings. In the case that both are positive numbers,
//' the sifting ends when either of the conditions is fulfilled.
//'
//' @export
//' @name eemd
//' @param input Vector of length N. he input signal to decompose.
//' @param ensemble_size Number of copies of the input signal to use as the ensemble.
//' @param noise_strength Standard deviation of the Gaussian random numbers used as additional
//'         noise. **This value is relative** to the standard deviation of the input signal.
//' @param S_number Integer. Use the S-number stopping criterion [2] for the EMD procedure with the given values of `S`.
//'        That is, iterate until the number of extrema and zero crossings in the
//'        signal differ at most by one, and stay the same for S consecutive
//'        iterations. Typical values are in the range 3--8. If \code{S_number} is
//'        zero, this stopping criterion is ignored.      
//' @param num_siftings Use a maximum number of siftings as a stopping criterion. If
//'        `num_siftings` is zero, this stopping criterion is ignored.
//' @param rng_seed A seed for the random number generator. A value of zero denotes
//'        an implementation-defined default value.
//' @return Matrix of size MxN, where M = \code{emd_num_imfs}(N). The rows of the array are the
//'        IMFs of the input signal, with the last row being the final residual.
//' 
//' @references
//'       [1] Z. Wu and N. Huang, "Ensemble Empirical Mode Decomposition: A 
//'       Noise-Assisted Data Analysis Method", Advances in Adaptive Data Analysis,
//'       Vol. 1 (2009) 1–41 \cr
//'       [2] N. E. Huang, Z. Shen and S. R. Long, "A new view of nonlinear water
//'       waves: The Hilbert spectrum", Annual Review of Fluid Mechanics, Vol. 31
//'       (1999) 417–457
//' @seealso \code{\link{emd}}
// [[Rcpp::export(eemd)]]
NumericMatrix eemdR(NumericVector input, unsigned int ensemble_size=250, 
double noise_strength=0.2, unsigned int S_number=4, unsigned int num_siftings=0, 
unsigned long int rng_seed=0){
  
  size_t N = input.size();
  NumericMatrix output(N,emd_num_imfs(N));
  eemd(input.begin(), N, output.begin(), ensemble_size, noise_strength, S_number, num_siftings, rng_seed);
  return output;
}