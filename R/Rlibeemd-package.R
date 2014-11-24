#' Rlibeemd: Ensemble empirical mode decomposition (EEMD) and its complete variant (CEEMDAN)
#'
#' Package Rlibeemd contains functions for the ensemble empirical mode decomposition (EEMD), 
#' its complete variant (CEEMDAN) or the regular empirical mode decomposition (EMD).
#' 
#' Package is based on the libeemd C library.
#' @references
#' \itemize{
##'  \item{P. Luukko, J. Helske and E. Rasanen, 
#' "Introducing libeemd: a C library and a Python module for performing the
#' ensemble empirical mode decomposition", Submitted.}
#' \item{ Z. Wu and N. Huang, "Ensemble Empirical Mode Decomposition: A 
#'       Noise-Assisted Data Analysis Method", Advances in Adaptive Data Analysis,
#'       Vol. 1 (2009) 1--41.}
#' \item{ N. E. Huang, Z. Shen and S. R. Long, "A new view of nonlinear water
#'       waves: The Hilbert spectrum", Annual Review of Fluid Mechanics, Vol. 31
#'       (1999) 417--457.}
#' \item{Torres et al, A Complete Ensemble Empirical Mode Decomposition with Adaptive Noise
#'   IEEE Int. Conf. on Acoust., Speech and Signal Proc. ICASSP-11,
#'   (2011) 4144--4147.}
#'   }
#' @docType package
#' @name Rlibeemd
#' @aliases Rlibeemd
#' @useDynLib Rlibeemd
#' @importFrom Rcpp evalCpp
NULL