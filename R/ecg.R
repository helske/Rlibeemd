#'
#'  Electrocardiogram Data
#'  
#'  Example ECG data from MIT-BIH Normal Sinus Rhythm Database, ECG1 of record 16265,
#'  first 2049 observations (0 to 16 seconds with sampling interval of 0.0078125 seconds) 
#'  
#' @name ECG  
#' @docType data
#' @format A time series object.
#' @source MIT-BIH Normal Sinus Rhythm Database, PhysioBank ATM, \url{http://www.physionet.org/cgi-bin/atm/ATM}
#' @keywords datasets
#' @examples
#' data(ECG)
#' plot(ECG)
NULL
#'
#'  Float Data
#'  
#'  The data are a position record from an acoustically tracked subsurface oceanographic float, 
#'  used as an example data in Rilling et al (2007).
#'  
#' @name float
#' @docType data
#' @format A time series object.
#' @source  http://wfdac.whoi.edu
#' @keywords datasets
#'  @references
#' \enumerate{
#'       \item{G. Rilling, P. Flandrin, P. Goncalves and J. M. Lilly,
#'   "Bivariate Empirical Mode Decomposition",
#'   IEEE Signal Processing Letters, Vol. 14 (2007) 936--939}
#'       }
#' @examples
#' data(float)
#' plot(float, type = "l"))
NULL
