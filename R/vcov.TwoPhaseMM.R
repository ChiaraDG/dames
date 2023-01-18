#' Print Covariance Matrix
#'
#' @param object Object of class TwoPhaseMM
#' @param ... additional input
#' @return Return a list with two elements: the covariance matrix of the coefficients in the
#' marginal mean model and the covariance matrix of the coefficients in the dependence model
#'
#' @export
#'
vcov.TwoPhaseMM <- function (object, ...) {
  n_beta     = length(object$beta)
  n_alpha    = length(object$alpha)
  beta_ix    = seq(n_beta)
  alpha_ix   = seq(n_beta+1, n_beta+n_alpha)
  object_vc  = object$mod.cov
  if(object$control['useRobCov'] == 1) object_vc  = object$rob.cov
  list('beta'=object_vc[beta_ix,beta_ix], 'alpha'=object_vc[alpha_ix,alpha_ix])
}
