#include <Rcpp.h>

extern "C"
{
  #include "eemd.h"
}

using namespace Rcpp;
// [[Rcpp::export]]
List extremaR(NumericVector x){
  
  size_t N = x.size();
  NumericVector maxx(x.size());
  NumericVector maxy(x.size());
  NumericVector minx(x.size());
  NumericVector miny(x.size());
  size_t nmax;
  size_t nmin;
  emd_find_extrema(x.begin(), N, maxx.begin(), maxy.begin(), &nmax, minx.begin(), miny.begin(), &nmin);
  
  return List::create(Named("x_max") = head(maxx,nmax),Named("y_max") = head(maxy,nmax),
  Named("x_min") = head(minx,nmin), Named("y_min") = head(miny,nmin));
  
}
