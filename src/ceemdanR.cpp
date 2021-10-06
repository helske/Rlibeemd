#include <Rcpp.h>

extern "C"
{
  #include "eemd.h"
}

using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix ceemdanR(NumericVector input, double num_imfs=0, unsigned int ensemble_size=250, 
double noise_strength=0.2, unsigned int S_number=4, unsigned int num_siftings=50, 
unsigned long int rng_seed=0, int threads=0){ 
  
  size_t N = input.size();
  size_t M = 0;
  if(num_imfs==0){
    M = emd_num_imfs(N);
  } else {
    M = (size_t)num_imfs;
  }
  NumericMatrix output(N,M);
  libeemd_error_code err = ceemdan(input.begin(), N, output.begin(), M, ensemble_size, 
    noise_strength, S_number, num_siftings, rng_seed, threads);
  

  
  if(err!=EMD_SUCCESS){
    printError(err);
  }
  return output;
}
