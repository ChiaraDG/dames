#' Compute the inverse of logit function
#'
#' @param x numeric vector
#'
#' @keywords internal
#' @return The inverse logit of the input x
#' @export
#'
expit <- function(x){
  if(!is.numeric(x)){stop("x needs to be numeric")}
  exp(x)/(1+exp(x))
  }
