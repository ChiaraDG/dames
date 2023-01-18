#' Print the estimated coefficients
#'
#' @param object Object of class TwoPhaseMM
#' @param ... other inputs
#'
#' @return Estimated coefficients from a TwoPhaseMM object
#' @export
#'
coef.TwoPhaseMM <- function (object, ...) {
  list('beta'=object$beta, 'alpha'=object$alpha)
}
