#include <Rcpp.h>
#ifndef restrict
#define restrict // nothing
#endif
extern "C"
{
  #include "eemd.h"
}

using namespace Rcpp;

//' Cubic Spline Interpolation
//' 
//' Perform cubic spline interpolation with nodes defined by the vectors x and
//' y, each of length N. The spline is evaluated using the not-a-node end point
//' conditions (same as Matlab). The y values of the spline curve will be
//' evaluated at integer points from 0 to x[N-1]. The endpoint x[N-1] is assumed to be an
//' integer, and the x values are assumed to be in ascending order, with x[0]
//' equal to 0.
//' @export
//' @name cspline
//' @param x,y Vectors giving the coordinates of the nodes.
//' @return The cubic spline curve.
// [[Rcpp::export]]
NumericVector cspline(NumericVector x, NumericVector y){
  size_t N = x.size();
  NumericVector spline(x.size());
  NumericVector spline_workspace(5*N-10);
  emd_evaluate_spline(x.begin(), y.begin(), N, spline.begin(), spline_workspace.begin());
  return spline;
}
