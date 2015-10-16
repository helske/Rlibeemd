#include <Rcpp.h>
#include <gsl/gsl_errno.h>
using namespace Rcpp;

// [[Rcpp::export]]
void gslErrorHandlerOff() {
  gsl_set_error_handler_off();
}
