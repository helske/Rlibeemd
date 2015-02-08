#include <Rcpp.h>
#ifndef restrict
#define restrict // nothing
#endif
extern "C"
{
  #include "eemd.h"
}
#include "extras.h"

using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix eemdR(NumericVector input, double num_imfs=0, unsigned int ensemble_size=250, 
double noise_strength=0.2, unsigned int S_number=4, unsigned int num_siftings=50, 
unsigned long int rng_seed=0, int threads=0){
  
  #ifdef _OPENMP
  int old_maxthreads = 1;
  if (threads>0) {
    old_maxthreads = omp_get_max_threads();
    omp_set_num_threads(threads);    
  }
  #endif
  
  size_t N = input.size();
  size_t M = 0;
  if(num_imfs==0){
    M = emd_num_imfs(N);
  } else {
    M = (size_t)num_imfs;
  }
  NumericMatrix output(N,M);
  libeemd_error_code err = eemd(input.begin(), N, output.begin(), M, ensemble_size, noise_strength, S_number, num_siftings, rng_seed);
  
  #ifdef _OPENMP
  if (threads>0) {
    omp_set_num_threads(old_maxthreads);    
  }
  #endif
  
  if(err!=EMD_SUCCESS){
    printError(err);
  }
  return output;
}