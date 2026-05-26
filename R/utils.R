#' Triangular distribution random number generator
#'
#' Generate random numbers following a triangular distribution,
#' useful for simulating continuous variables whose median is not
#' at the center of the range.
#'
#' @param n Sample size
#' @param min_val Minimum value
#' @param mode_val Mode (peak location)
#' @param max_val Maximum value
#' @return A numeric vector of length n
#' @export
#'
#' @examples
#' rtriang(100, 28, 64, 90)
rtriang <- function(n, min_val, mode_val, max_val) {
  u <- stats::runif(n)
  f <- (mode_val - min_val) / (max_val - min_val)
  ifelse(u < f,
         min_val + sqrt(u * f) * (max_val - min_val),
         max_val - sqrt((1 - u) * (1 - f)) * (max_val - min_val))
}
