#include <Rcpp.h>
#ifndef restrict
#define restrict // nothing
#endif
extern "C"
{
  #include "eemd.h"
}

using namespace Rcpp;
//' Number of IMFs
//' 
//' Return the number of IMFs extracted from input data of length N, including
//' the final residual. This is just 1+[log_2(N)] for N>3.
//' @export
//' @name nIMFs
//' @param N An integer defining the length of input data.
//' @return The number of IMFs which would be extracted from input data of length N, including
//' the final residual.
// [[Rcpp::export(emd_num_imfs)]]
int emd_num_imfsR(double N) {
return emd_num_imfs((size_t)N);
}
