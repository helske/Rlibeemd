#' Bivariate EMD decomposition
#' 
#' Function \code{bemd} implements the Bivariate EMD (Scheme 2 in the cited article).
#'
#' @export
#' @name bemd
#' @param input Complex vector of length N. The input signal to decompose.
#' @param directions Vector of directions (phi_k in the article) used for decomposition.
#' @param num_imfs Number of Intrinsic Mode Functions (IMFs) to compute. If num_imfs is set to zero, a value of
#'        num_imfs = emd_num_imfs(N) will be used, which corresponds to a maximal number of
#'        IMFs. Note that the final residual is also counted as an IMF in this
#'        respect, so you most likely want at least num_imfs=2.
#' @param num_siftings Use a maximum number of siftings as a stopping criterion. If
#'        \code{num_siftings} is zero, this stopping criterion is ignored. Default is 50.
#' @return Time series object of class \code{"mts"} where series corresponds to
#'        IMFs of the input signal, with the last series being the final residual.
#'  @references
#' \enumerate{
#'       \item{G. Rilling, P. Flandrin, P. Goncalves and J. M. Lilly,
#'   "Bivariate Empirical Mode Decomposition",
#'   IEEE Signal Processing Letters, Vol. 14 (2007) 936--939}
#'       }
#' @examples 
#' N <- 512 
#' t <- 2 * pi * (0:(N-1))/N
#' input <- cos(0.3 * t) * exp(2i * t) + 0.3 * abs(sin(2.3 * t)) * exp(17i * t)
#' 
#' # Use evenly spaced angles as directions
#' num_directions <- 64
#' directions <- 2 * pi * 1:num_directions / num_directions
#' output <- bemd(input, directions, num_imfs = 4, num_siftings = 10)
#' plot(Re(input), Im(input), xlim = c(-1, 2))
#' for(i in 1:4)
#'   points(Re(output[,i]), Im(output[,i]), col = 1 + i)
#' legend("bottomright", col = 1:5, legend = c("signal", paste0("IMF ",1:4)), pch = 1)
#' 
bemd <- function(input, directions, num_imfs = 0, num_siftings = 50L) {
  if (!all(is.finite(input))) 
    stop("'input' must contain finite values only.")
  if (num_imfs < 0)
    stop("Argument 'num_imfs' must be non-negative integer.")
  if (num_siftings < 0)
    stop("Argument 'num_siftings' must be non-negative integer.")
  if (!all(is.finite(directions))) 
    stop("'input' must contain finite values only.")
  
  output <- bemdR(input, directions,num_imfs, num_siftings)
  if (inherits(input, "ts")) {
    tsp(output) <- tsp(input)
  } else tsp(output) <- c(1, nrow(output), 1)
  if (ncol(output) > 1) {
    class(output) <- c("mts", "ts", "matrix")
    colnames(output) <- c(paste("IMF", 1:(ncol(output) - 1)), "Residual")
  } else class(output) <- "ts"
  output
}