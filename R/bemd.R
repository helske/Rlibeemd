#' Bivariate EMD decomposition
#' 
#' Function \code{bemd} implements the Bivariate EMD (Scheme 2 in the cited article).
#'
#' @export
#' @name bemd
#' @param input Complex vector of length N. The input signal to decompose.
#' @param directions Vector of directional angles (in radians) to use for the decomposition, 
#' or an integer defining the number of equally spaced angles to use.
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
#' imfs <- bemd(input, directions, num_imfs = 4, num_siftings = 10)
#' 
#' # plot the data
#' plot(Re(input), Im(input), xlim = c(-1, 2))
#' # plot signal and the imfs
#' for(i in 1:4)
#'   points(Re(imfs[,i]), Im(imfs[,i]), col = 1 + i)
#' legend("bottomright", col = 1:5, legend = c("signal", paste0("IMF ",1:4)), pch = 1)
#' 
#' data("float")
#' plot(float, type = "l")
#' signal <- float[, 1] + float[, 2] * 1i
#' imfs <- bemd(signal, num_siftings = 10, num_imfs = 4)
#' 
#' # plot the data and the imfs
#' oldpar <- par()
#' par(mfrow = c(5, 1), mar = c(0.5, 4.5, 0.5, 0.5), oma = c(4, 0, 2, 0))
#' ts.plot(float, col = 1:2, lty = 1:2, ylab = "signal", gpars = list(xaxt = "n"))
#' for(i in 1:4) {
#'   ts.plot(Re(imfs[, i]), Im(imfs[, i]), col = 1:2, lty = 1:2, 
#'     ylab = if(i < 4) paste("IMF", i) else "residual", gpars = list(xaxt = "n"))
#'  }
#' axis(1)
#' title(xlab = "Time (days)", main = "Bivariate EMD decomposition", outer = TRUE)
#' par(oldpar)
bemd <- function(input, directions = 64L, num_imfs = 0L, num_siftings = 50L) {
  
  if (!all(is.finite(input))) 
    stop("'input' must contain finite values only.")
  if (num_imfs < 0)
    stop("Argument 'num_imfs' must be non-negative integer.")
  if (num_siftings < 0)
    stop("Argument 'num_siftings' must be non-negative integer.")
  if (!all(is.finite(directions))) 
    stop("'input' must contain finite values only.")
  if (!is.complex(input)) 
    stop("Argument 'input' must be a complex vector. ")
  
  if (length(directions) == 1) {
    if(directions <= 0) stop("Argument 'directions' must be a numeric vector of positive integer. ")
    directions <- 2 * pi * 0:(directions - 1) / directions
  }
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