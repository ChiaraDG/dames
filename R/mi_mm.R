#' Multiple Imputation for Two-Phase ODS Designs and Marginalized Transition and Latent Variable Models
#'
#' Estimates regression parameters from a marginalized transition and latent variable model for scenarios where
#' data are collected using a two-phase ODS design and the expensive exposure of interest is binary.
#'
#' @param mean.formula  Formula for the marginal mean model in which a binary variable is regressed on a set of covariates. This is a required argument
#' @param lv.formula Formula for the latent part in the dependence model. This is optional and of the form ~ covariates. Default is set to NULL.
#' @param t.formula Formula for the transition part in the dependence model. This is optional and of the form ~ covariates. Default is set to NULL.
#' @param id Column name for the id variable as recorded in the dataset. This is a required argument
#' @param data Full cohort data. This dataset should not have missing values for variables other than X. It should
#' contain columns named \code{id}, \code{X}, and \code{sampled}
#' @param sampled Column name for the sampling indicator variable as recorded in the dataset. This is a required argument
#' @param X Column name for the expensive exposure variable as recorded in the dataset. This is a required argument
#' @param samp.probs Vector of 3 numbers indicating the probability that a subject is sampled from the none
#'                   some and all stratum respectively. This is a required argument when \code{method} is set equal to indirect. It is not
#'                   used when \code{method} is set equal to direct, Default is set to NULL
#' @param maxiter A scalar indicating the maximum number of iteration in the maximization algorithm. Default is 100
#' @param Q A scalar to denote the number of quadrature points used for GH numerical integration. Default is 10.
#' @param M Number of imputations. Default is 5.
#' @param m A scalar giving the number of iterations within each imputation. The default is 5. The argument is only when \code{method} is set to direct.
#' @param marg.exp.formula Marginal exposure model of the form: X ~ Z where Z are the inexpensive covariates in \code{data}. This is a required argument.
#' @param method Character indication whether the indirect MI or the direct MI should be performed. Default is indirect. Note that indirect MI can only be used under an NSA design
#' @param verbose An indicator of whether the steps of the maximization algorithm should be printed
#'
#' @importFrom dplyr n
#'
#' @return This function returns marginal mean (beta) and dependence parameters (alpha) along with the associated model and empirical covariance matrices.
#' @export
#'
#' @examples
#' \donttest{
#' data(exampledat)
#' # perform the NSA design
#' design       <- nsa_design(Y = "Y", id = "id", n.sample = c(85, 100, 15), data = exampledat)
#' # Assume X is not collected for subjects not sampled
#' exampledat$X <- ifelse(exampledat$id %in% design$sampled.id, exampledat$X, NA)
#' # estimate the parameters using Direct MI
#' mod.mi2 <- mi_mm(mean.formula = Y ~ X + time + Z + X:time, t.formula = ~1,
#'                  lv.formula = NULL, id = "id", X = "X",
#'                  data = exampledat, marg.exp.formula = X ~ Z, M = 5,
#'                 method = "direct")
#' # print the results
#' summary(mod.mi2)
#' }
#'
mi_mm <- function(mean.formula, lv.formula = NULL, t.formula = NULL, id,
                  data, sampled = NULL, X,
                  samp.probs = NULL, maxiter=100, Q=10, M=5, m = 5, marg.exp.formula,
                  method = "indirect", verbose=FALSE){

  if(is.null(samp.probs) & method == "indirect"){
    stop("If method is set to indirect, sampling probabilities need to be specified")}
  if(is.null(sampled) & method == "indirect"){
    stop("If method is set to indirect, a variable name indicating the sampling indicator in the data.frame need to be specified")}
  if(is.null(lv.formula) & is.null(t.formula)){
    stop("At least one of t.formula and lv.formula needs to be specified")}
  if(!is.null(samp.probs) & method == "direct"){
    warning("Sampling probabilities are not used when method is set to direct")}
  if(!is.character(X) | !is.character(id)){
    stop("X and/or id need to be character vectors")}

  id       = data[, id]
  Sampled  = data[, sampled]
  ## Keep those who were sampled to create the ODS dataset
  odsdat   = data[Sampled==1,]

  if(method == "indirect"){
    FIT  = acml_mm(mean.formula = mean.formula, lv.formula=lv.formula, t.formula=t.formula,
              data=odsdat, id=id,
              samp.probs=samp.probs, maxiter=maxiter, verbose=FALSE, Q=Q, return_args=TRUE)
    imps = ImputeDataInd(FIT=FIT, DAT=data, M=M, ROB=FALSE, MARG.EXP.FORMULA = marg.exp.formula,
                          verbose=verbose, EXPVAR = X, Sampled = Sampled)
    final= FitImputedData(FIT, ImputedData = imps, id, iter.lim = maxiter, q = Q)
  }

  if(method == "direct"){

    final =  DirectImputation(mean.formula =mean.formula,
                            t.formula = t.formula,
                            lv.formula = lv.formula, id = id,
                            data = data, Q = Q, m = m, M = M,
                            X = X, marg.exp.formula = marg.exp.formula, verbose = verbose)

  }

  # set up results for printing
  maxtime   = max(data %>% group_by(id) %>% summarise(mtime = n()) %>% pull())
  N         = length(unique(id))
  nms       = list()
  nms$beta  = colnames(model.matrix(mean.formula, data))
  nms$alpha = c(if(!is.null(t.formula)){paste('gamma',colnames(model.matrix(t.formula, data)),sep=':')},
                if(!is.null(lv.formula)){paste('log(sigma)',colnames(model.matrix(lv.formula, data)),sep=':')})

  out = list(call    = match.call(),
             alpha   = final$coefficients[(length(nms$beta) + 1):length(final$coefficients)],
             beta    = final$coefficients[1:length(nms$beta)],
             mod.cov = final$variance,
             control = c(1, N, 0, 0, maxtime, M))
  names(out$control)  = c("convergence_code", "n_subj", "se_comp", "useRobCov", 'max_n_visit',
                          'num_imps')
  names(out$beta)     = nms$beta
  names(out$alpha)    = nms$alpha


  class(out) = "TwoPhaseMM"


  out
}
