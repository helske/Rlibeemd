#include <Rcpp.h>
extern "C"
{
  #include "bemd.h"
}
using namespace Rcpp;

// [[Rcpp::export]]
ComplexMatrix bemdR(ComplexVector input, NumericVector directions,
  double num_imfs = 0, unsigned int num_siftings = 50){
  
  size_t N = input.size();
  size_t M = 0;
  if(num_imfs==0){
    M = emd_num_imfs(N);
  } else {
    M = (size_t)num_imfs;
  }
  size_t D = directions.size();
  
  ComplexMatrix output(N,M);
  
  libeemd_error_code err = bemd(reinterpret_cast<double _Complex const*>(input.begin()), N, 
    directions.begin(), D, reinterpret_cast<double _Complex*>(output.begin()), M, num_siftings);
 
  if(err!=EMD_SUCCESS){
    printError(err);
  }
  return output;
}
