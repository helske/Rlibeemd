#' CEEMDAN decomposition
#' 
#' Decompose input data to Intrinsic Mode Functions (IMFs) with the
#' Complete Ensemble Empirical Mode Decomposition with Adaptive Noise (CEEMDAN)
#' algorithm [1], a variant of EEMD.
#'
#' The size of the ensemble and the relative magnitude of the added noise are
#' given by parameters \code{ensemble_size} and \code{noise_strength}, respectively.  The
#' stopping criterion for the decomposition is given by either a S-number [2] or
#' an absolute number of siftings. In the case that both are positive numbers,
#' the sifting ends when either of the conditions is fulfilled.
#'
#' @export
#' @name ceemdan
#' @inheritParams eemd
#' @return Time series object of class \code{"mts"} where series corresponds to
#'        IMFs of the input signal, with the last series being the final residual.
#' @references
#' \enumerate{ 
#'  \item{M. Torres et al, A Complete Ensemble Empirical Mode Decomposition with Adaptive Noise
#'   IEEE Int. Conf. on Acoust., Speech and Signal Proc. ICASSP-11,
#'   (2011) 4144--4147}
#'  \item{N. E. Huang, Z. Shen and S. R. Long, "A new view of nonlinear water
#'       waves: The Hilbert spectrum", Annual Review of Fluid Mechanics, Vol. 31
#'       (1999) 417--457}
#'       }
#' @seealso \code{\link{eemd}} 
#' @examples
#' imfs <- ceemdan(UKgas)
#' # trend extraction
#' ts.plot(UKgas, imfs[,ncol(imfs)], col = 1:2, main = "Quarterly UK gas consumption", ylab = "Million therms")
#' 
#' # CEEMDAN for logarithmic demand, note that increasing ensemble size 
#' # from default will produce smoother results
#' imfs <- ceemdan(log(UKgas))
#' plot(ts.union(Seasonal = imfs[,1], Irregular = rowSums(imfs[,2:5]), Trend = imfs[,6]), 
#'      main = "Quarterly UK gas consumption")
ceemdan <- function(input, num_imfs = 0, ensemble_size = 250L, noise_strength = 0.2, S_number = 4L, 
                    num_siftings = 50L, rng_seed = 0L) {
  output<-.Call('Rlibeemd_ceemdanR', PACKAGE = 'Rlibeemd', input, num_imfs, ensemble_size, 
             noise_strength, S_number, num_siftings, rng_seed)
  if(inherits(input,"ts")){
    tsp(output)<-tsp(input)
  } else tsp(output)<-c(1,nrow(output),1)
  if(ncol(output)>1){
  class(output)<-c("mts","ts","matrix")
  colnames(output)<-c(paste("IMF",1:(ncol(output)-1)),"Residual")
  } else class(output)<-"ts"
  output
}