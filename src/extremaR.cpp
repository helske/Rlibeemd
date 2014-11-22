#include <Rcpp.h>
#ifndef restrict
#define restrict // nothing
#endif
extern "C"
{
  #include "eemd.h"
}

using namespace Rcpp;
//' Find the local extrema of the input data
//' 
//' Find the local minima and maxima from input data \code{x}. This includes the
//' artificial extrema added to the ends of the data as specified in the
//' original EEMD article.
//' @export
//' @name extrema
//' @param x Vector of length x.
//' @return a list with components ...
//'
//' @references Z. Wu and N. Huang, "Ensemble Empirical Mode Decomposition: A
//'       Noise-Assisted Data Analysis Method", Advances in Adaptive Data Analysis,
//'       Vol. 1 (2009) 1-41
// [[Rcpp::export]]
List extrema(NumericVector x){
  
  size_t N = x.size();
  NumericVector maxx(x.size());
  NumericVector maxy(x.size());
  NumericVector minx(x.size());
  NumericVector miny(x.size());
  size_t nmax = 0;
  size_t nmin = 0;
  bool all_extrema_good = emd_find_extrema(x.begin(), N, maxx.begin(), maxy.begin(), &nmax, minx.begin(), miny.begin(), &nmin);
  
  return List::create(Named("max_x") = head(maxx,nmax),Named("max_y") = head(maxy,nmax),
  Named("min_x") = head(minx,nmin), Named("min_y") = head(miny,nmin), Named("all_extrema_good") = wrap(all_extrema_good));
  
}
