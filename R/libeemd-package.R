#' libeemd: Ensemble empirical mode decomposition (EEMD) and its complete variant (CEEMDAN)
#'
#' Package libeemd contains functions for the ensemble empirical mode decomposition (EEMD), 
#' its complete variant (CEEMDAN) or the regular empirical mode decomposition (EMD).
#' 
#' Package is based on the libeemd C library.
#' @references P. Luukko, J. Helske and E. Räsänen, 
#' "Introducing libeemd: a C library and a Python module for performing the
#' ensemble empirical mode decomposition", Submitted.
#'       Z. Wu and N. Huang, "Ensemble Empirical Mode Decomposition: A 
#'       Noise-Assisted Data Analysis Method", Advances in Adaptive Data Analysis,
#'       Vol. 1 (2009) 1–41. \cr
#'        N. E. Huang, Z. Shen and S. R. Long, "A new view of nonlinear water
#'       waves: The Hilbert spectrum", Annual Review of Fluid Mechanics, Vol. 31
#'       (1999) 417–457M. \cr
#'       Torres et al, A Complete Ensemble Empirical Mode Decomposition with Adaptive Noise
#'   IEEE Int. Conf. on Acoust., Speech and Signal Proc. ICASSP-11,
#'   (2011) 4144-4147.
#' @docType package
#' @name libeemd
#' @aliases libeemd
#' @useDynLib libeemd
#' @importFrom Rcpp evalCpp
#' @examples
# x = 0:511
# y = cos(2*pi/30*x) + cos(2*pi/34*x)
# out<-eemd(y, S_number=4,ensemble_size=250, noise_strength=0.2,
#                       num_siftings=0, rng_seed=0)
# 
# ts.plot(out)
