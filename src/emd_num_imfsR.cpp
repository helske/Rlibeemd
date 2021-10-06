#include <Rcpp.h>

extern "C"
{
  #include "eemd.h"
}

using namespace Rcpp;

// [[Rcpp::export]]
int emd_num_imfsR(double N) {
return emd_num_imfs((size_t)N);
}
