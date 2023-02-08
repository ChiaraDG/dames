#' Ascertainment Corrected Maximum Likelihood Estimator for the NSA design
#' The function estimates the parameters of interest for the NSA design using the ascertainment correct maximum likelihood (ACML)
#'
#' @importFrom  stats glm
#' @importFrom  stats model.response
#' @importFrom  stats na.omit
#'
#' @param mean.formula  Formula for the marginal mean model in which a binary variable is regressed on a set of covariates. This is a required argument
#' @param lv.formula Formula for the latent part in the dependence model. This is optional and of the form ~ covariates
#' @param t.formula Formula for the transition part in the dependence model. This is optional and of the form ~ covariates
#' @param id Vector with \code{id} variable. This is a required argument
#' @param data  Dataframe. This is a required argument
#' @param inits An optional list of length 3 containing initial values for marginal mean parameters and all dependence parameters. The format of the list should be: (1) estimates of the mean parameters, (2) estimates of the transition parameters (or NULL if only fitting a mLV model) and (3) estimates of the latent variable parameters (or NULL if only fitting a mT model). If NULL, initial values will be automatically generated. Default is set to NULL
#' @param samp.probs Vector of 3 numbers indicating the probability that a subject is sampled from the none
#'                   some and all stratum respectively. If not specified, the sampling probabilities are set to 1
#' @param offset an optional offset term.
#' @param Q A scalar to denote the number of quadrature points used for GH numerical integration. Default is 10.
#' @param step.max A scalar indicating the step maximum in the maximization algorithm. Default is 1
#' @param step.tol A scalar indicating the step tolerance in the maximization algorithm. Default is 1e-6
#' @param hess.eps A scalar indicating the perturbation value in the Hessian computation. Default is 1e-7
#' @param adapt.quad An indicator if adaptive quadrature is to be used.  NOT CURRENTLY IMPLEMENTED.
#' @param verbose An indicator of whether the steps of the maximization algorithm should be printed
#' @param maxiter A scalar indicating the maximum number of iteration in the maximization algorithm. Default is 100
#' @param return_args If TRUE returns the parameters used for the likelihood computation. Default is FALSE
#'
#' @return This function returns marginal mean (beta) and dependence parameters (alpha) along with the associated model and empirical covariance matrices.
#'
#' @export
#'
#' @examples
#' \donttest{
#' # perform the NSA design
#' data("exampledat")
#' design       <- nsa_design(Y = "Y", id = "id", n.sample = c(85, 100, 15), data = exampledat)
#' # Assume X is not collected for subjects not sampled
#' exampledat$X <- ifelse(exampledat$id %in% design$sampled.id, exampledat$X, NA)
#' # estimate the parameters using ACML
#' mod.acml    <- acml_mm(mean.formula = Y ~ X + time + Z + X:time, t.formula = ~1, id = id,
#'                        samp.probs = design$sample.probs, data = exampledat)
#' # print the results
#' summary(mod.acml)
#'
#' }


acml_mm <- function(mean.formula, lv.formula = NULL, t.formula = NULL, id, data, inits = NULL,
                    samp.probs = c(1, 1, 1), offset = NULL, Q = 10,
                    step.max = 1, step.tol = 1e-06, hess.eps = 1e-07,
                    adapt.quad = FALSE, verbose = FALSE, maxiter = 100, return_args=FALSE) {

  if(is.null(lv.formula) & is.null(t.formula)) {
    stop('Specify association model (both lv.formula and t.formula arguments cannot be NULL.')}
  if(!is.data.frame(data)) {
    data <- as.data.frame(data)
    warning('data converted to data.frame.')
  }

  terms = unique( c(all.vars(mean.formula), all.vars(lv.formula), all.vars(t.formula),
                    as.character(substitute(id))) )
  data  = data[,terms]
  if(any(is.na(data))) data = na.omit(data)
  id0   =  as.character(substitute(id))
  id    = data$id    = data[ , id0 ]

  if(Q<=1) Q <- 2
  if(is.null(lv.formula)) q = 1

  mean.f = model.frame(mean.formula, data)
  mean.t = attr(mean.f, "terms")
  y  = model.response(mean.f,'numeric')
  uy = unique(y)
  x  = model.matrix(mean.formula,mean.f)

  x.t = x.lv = matrix(0, ncol=1, nrow=length(y))
  if(!is.null(t.formula))   x.t  = model.matrix(t.formula,model.frame(t.formula, data))
  if(!is.null(lv.formula))  x.lv = model.matrix(lv.formula, model.frame(lv.formula, data))

  # transform samp.probs in matrix form
  samp.probs = matrix(samp.probs, nrow = length(y), ncol = 3, byrow = TRUE)


  # weights need to be set to 1
  samp.probi = matrix(1,nrow=length(y),ncol=1)
  useROBCOV <- FALSE


  if (!is.null(inits)) {
    inits <- unlist(inits)
  }

  if(is.null(inits)) {
    inits = c(glm(mean.formula,family='binomial',data=data)$coef, rep(1, ncol(x.t) + ncol(x.lv)))
    if(any(is.na(inits))) {
      omit_dup_col = which(is.na(inits))
      x            = x[,-c(omit_dup_col)]
      inits        = inits[-c(omit_dup_col)]
    }
  }

  if(is.null(offset)) {
    offset <- rep(0, length(y))
  }

  mm.fit = twophase_cond(params=inits, id=id, X=x, Y=y, Xgam=x.t, Xsig=x.lv, Q=Q, condlike=TRUE,
                    sampprobs=samp.probs, sampprobi=samp.probi, offset=offset,
                    stepmax=step.max, steptol=step.tol, hess.eps=hess.eps,
                    AdaptiveQuad=adapt.quad, verbose=verbose,iterlim=maxiter)

  nms = list()
  nms$beta = colnames(x)
  nms$alpha = c( if(!is.null(t.formula)){paste('gamma',colnames(x.t),sep=':')},
                 if(!is.null(lv.formula)){paste('log(sigma)',colnames(x.lv),sep=':')})

  out = NULL
  out$call    = match.call()
  out$logLik = mm.fit$logL
  out$beta  = mm.fit$beta
  out$alpha = mm.fit$alpha
  out$mod.cov = mm.fit$modelcov
  out$rob.cov = mm.fit$empiricalcov
  names(out$beta) = nms$beta
  names(out$alpha) = nms$alpha
  colnames(out$mod.cov) = rownames(out$mod.cov) = colnames(out$rob.cov) = rownames(out$rob.cov) = unlist(nms)
  out$control = with(mm.fit, c(condlike, AdaptiveQuad, code, niter, length(table(id)), max(table(id)),
                               useROBCOV, 0))
  names(out$control) <- c('cond.like','adaptive.quad','convergence_code','n_iter','n_subj','max_n_visit',
                          'useRobCov', 'se_comp')

  aic = function(l=mm.fit$logL,k=nrow(mm.fit$modelcov)) 2*k-2*l
  bic = function(l=mm.fit$logL,k=nrow(mm.fit$modelcov),n=length(table(id))) -2*l + k *log(n)
  deviance = function(l=mm.fit$logL) -2*l
  out$info_stats = c(aic(),bic(),mm.fit$logL,deviance())
  names(out$info_stats) = c('AIC','BIC','logLik','Deviance')
  out$LogLikeSubj = mm.fit$LogLikeSubj
  out$ObsInfoSubj = mm.fit$ObsInfoSubj
  out$ACSubj      = mm.fit$ACSubj
  out$LLSC_args   = mm.fit$LLSC_args

  if(return_args) attr(out,'args') <- list('mean.formula'=mean.formula, 't.formula'=t.formula, 'lv.formula'=lv.formula,
                                           'id'=id0, 'samp.probs'=samp.probs,'samp.probi'=samp.probi,'offset'=offset)




  class(out) = "TwoPhaseMM"

  return(out)
}
