#' Compute the Log-Likelihood and the Score for Each Subject
#'
#' @param params  Vector of parameters from the marginal mean and the dependence model parameters
#' @param subjectData List of lists with subject-specific data
#' @param ParamLengths Vector of three elements: the number of parameters in the marginal mean model, the number of parameters for the transition component of the dependence model, the number of paramters for the dependence component of the dependence model
#' @param CondLike TRUE/FALSE. Set to TRUE from ascertainment corrected likelihood
#' @param EmpiricalCheeseCalc TRUE/FALSE. Set to TRUE for robust standard errors estimation
#' @param Q A scalar to denote the number of quadrature points used for GH numerical integration. Default is 10.
#' @param W A scalar resulting from GH numerical integration
#' @param Z A scalar resulting from GH numerical integration
#' @param AdaptiveQuad TRUE/FALSE. Set to TRUE for adaptive quadrature
#'
#' @keywords internal
#' @return Log-likelihood and ascertainment correction (if CondLike is set to be FALSE) computed for each subject in the data
#' @export
#'

LogLScoreCalc <-
  function( params, subjectData, CondLike = FALSE, ParamLengths,
            EmpiricalCheeseCalc = FALSE, Q, W, Z, AdaptiveQuad=FALSE){

    betaM <- params[1:ParamLengths[1]]
    gamma <- params[ (ParamLengths[1]+1): (ParamLengths[1]+ParamLengths[2]) ]
    sigma <- exp(params[(ParamLengths[1]+ParamLengths[2]+1): (ParamLengths[1]+ParamLengths[2]+ParamLengths[3])])

    .Call("LogLScoreCalc_CALL", betaM, gamma, sigma, subjectData, CondLike,
          EmpiricalCheeseCalc, Q, W, Z,  AdaptiveQuad)
  }
