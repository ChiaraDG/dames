#' Create Subject-Specific Data to Use in the SMLE function
#'
#' @importFrom stats setNames
#'
#' @param Y Character indicating the name of the outcome variable in the dataset
#' @param X Character indicating the name of the expensive exposure in the dataset
#' @param mean.formula Formula for the marginal mean model in which a binary variable is regressed on a set of covariates.
#' @param t.formula Formula for the transition part in the dependence model. This is of the form ~ covariates
#' @param lv.formula Formula for the latent part in the dependence model. This is of the form ~ covariates
#' @param id Column name with the \code{id} variable.
#' @param m Number of distinct values of the expensive exposure observed in phase two
#' @param weights Vector of weights for the M-step in the EM-algorithm
#' @param d dataset
#'
#' @keywords internal
#' @return Returns a list of data per each subject to be used in the likelihood computation
#' @export
#'

ptheta.input <- function(Y, X, mean.formula, t.formula, lv.formula, id, m, weights = NULL, d){

  if(is.null(weights)){weights = rep(1, length(d[,Y]))}
  # change the id names in the dataset as id.
  names(d)[names(d) == id]     = "id"

  n.var =  ncol(model.matrix(mean.formula, d))

  # create subject specific data
  x.transition = x.latent <- matrix(0, ncol = 1, nrow = length(d[,Y]))
  if(!is.null(t.formula)){
    x.transition   = model.matrix(t.formula, model.frame(t.formula, d))
    x.transition   = as.matrix(x.transition)
  }
  if(!is.null(lv.formula)){
    x.latent       = model.matrix(lv.formula, model.frame(lv.formula, d))
    x.latent       = as.matrix(x.latent)
  }
  Xgam.col         = ncol(x.transition)
  Xsig.col         = ncol(x.latent)
  paramslength     = c(n.var, Xgam.col, Xsig.col)

  # Lagged response
  Ylag             = rep(0, length(d[,Y]))
  for(k in 2:length(d[,Y])) {if(d$id[k] == d$id[k-1]) Ylag[k] = d[,Y][k-1]}

  sampProbi        = weights
  sampprobs        = matrix(1, ncol = 3, nrow = length(d[,Y]))
  offset           = rep(0, length(d[,Y]))
  mean.f           = model.frame(mean.formula, d)
  mean.pred        = model.matrix(mean.formula, mean.f)

  if(is.null(X)){
    subjectData    = list()
    uid            = as.character(unique(d$id))

    id.tmp         = split(d$id, d$id)
    id.tmp         = lapply(id.tmp, function(x)as.character(unique(x)))
    X.tmp          = split(mean.pred,  d$id)
    X.tmp          = lapply(X.tmp, function(x)matrix(x, ncol=ncol(mean.pred)))
    Y.tmp          = split(d[,Y],  d$id)
    Y.tmp          = lapply(Y.tmp, as.double)
    Xgam.tmp       = split(x.transition, d$id)
    Xgam.tmp       = lapply(Xgam.tmp, function(x)matrix(x, ncol=Xgam.col))
    Xsig.tmp       = split(x.latent, d$id)
    Xsig.tmp       = lapply(Xsig.tmp, function(x)matrix(x, ncol=Xsig.col))
    Ylag.tmp       = split(Ylag, d$id)
    Ylag.tmp       = lapply(Ylag.tmp, as.double)
    SampProbi.tmp  = split(sampProbi, d$id)
    SampProbi.tmp  = lapply(SampProbi.tmp, unique)
    SampProbs.tmp  = split(sampprobs, d$id)
    SampProbs.tmp  = lapply(SampProbs.tmp, function(x)matrix(x, ncol=ncol(sampprobs))[1,])
    offset.tmp     = split(offset, d$id)

    names(id.tmp) = names(X.tmp) = names(Y.tmp) = names(Xgam.tmp) = names(Xsig.tmp) = names(Ylag.tmp) =
      names(SampProbi.tmp) = names(SampProbs.tmp) = names(offset.tmp) = uid

    subjectData        = Map(list, id.tmp, X.tmp, Y.tmp, Ylag.tmp, Xgam.tmp, Xsig.tmp,
                             SampProbi.tmp, SampProbs.tmp, offset.tmp)
    subjectData        = unique(subjectData)
    names(subjectData) = uid
  }

  if(!is.null(X)){
    subjectData   = list()
    uid           = as.character(paste0(rep(unique(d$id), times = m), "_",
                                      rep(1:m, each = length(unique(d$id)))))

    id.tmp        = split(d$id, list(d$id, d[,X]))
    id.tmp        = lapply(id.tmp, function(x)as.character(unique(x)))
    X.tmp         = split(mean.pred,  list(d$id, d[,X]))
    X.tmp         = lapply(X.tmp, function(x)matrix(x, ncol=ncol(mean.pred)))
    Y.tmp         = split(d[,Y],  list(d$id, d[,X]))
    Y.tmp         = lapply(Y.tmp, as.double)
    Xgam.tmp      = split(x.transition, list(d$id, d[,X]))
    Xgam.tmp      = lapply(Xgam.tmp, function(x)matrix(x, ncol=Xgam.col))
    Xsig.tmp      = split(x.latent, list(d$id, d[,X]))
    Xsig.tmp      = lapply(Xsig.tmp, function(x)matrix(x, ncol=Xsig.col))
    Ylag.tmp      = split(Ylag, list(d$id, d[,X]))
    Ylag.tmp      = lapply(Ylag.tmp, as.double)
    SampProbi.tmp = split(sampProbi, list(d$id, d[,X]))
    SampProbi.tmp = lapply(SampProbi.tmp, unique)
    SampProbs.tmp = split(sampprobs, list(d$id, d[,X]))
    SampProbs.tmp = lapply(SampProbs.tmp, function(x)matrix(x, ncol=ncol(sampprobs))[1,])
    offset.tmp    = split(offset, list(d$id, d[,X]))

    names(id.tmp) = names(X.tmp) = names(Y.tmp) = names(Xgam.tmp) = names(Xsig.tmp) = names(Ylag.tmp) =
      names(SampProbi.tmp) = names(SampProbs.tmp) = names(offset.tmp) = uid

    subjectData        = Map(list, id.tmp, X.tmp, Y.tmp, Ylag.tmp, Xgam.tmp, Xsig.tmp,
                             SampProbi.tmp, SampProbs.tmp, offset.tmp)
    subjectData        = unique(subjectData)
    names(subjectData) = uid
  }

  # name each element of subjectData
  subjectData          = lapply(subjectData, function(x){setNames(x, c("id", "X", "Y", "Ylag",
                                                                       "Xgam", "Xsig", "SampProbi",
                                                                       "SampProbs", "offset"))})

  paramslength         =  c(n.var, Xgam.col, Xsig.col)

  out                  = list(subjectData = subjectData, p.length = paramslength)

  return(out)

}



#' Compute the Log-Likelihood for a marginalized transition and latent variable model
#'
#'
#' @param sData List of subject-specific data
#' @param params Vector of model parameters
#' @param p.length Length of the vector of model parameters
#' @param X Character indicating the name of the expensive exposure in the dataset
#' @param m Number of distinct values of the expensive exposure observed in phase two
#' @param Q Number of quadrature point to be used if there is a latent term in the dependence model
#' @param d Dataset
#'
#' @keywords internal
#' @return Returns the log-likelihood computed for each subject
#' @export
#'
ptheta.fun <- function(sData, params, p.length, X, m, Q, d){

  AdaptiveQuad = get.GH(Q)
  W            = AdaptiveQuad$w
  Z            = AdaptiveQuad$z

  lik          = LogLScoreCalc(params, subjectData = sData, CondLike = FALSE,
                               ParamLengths = p.length, EmpiricalCheeseCalc = FALSE, Q = Q,
                               W = W, Z = Z)
  ptheta       = attr(lik, "LogLikeSubj")
  u.id         = unique(d$id)

  if(!is.null(X)){
    ptheta.dat = data.frame(id = rep(u.id, m), ptheta = ptheta)
  }

  if(is.null(X)){
    ptheta.dat = data.frame(id = u.id, ptheta = ptheta)
  }

  return(ptheta.dat)
}


#' Compute the observed data log-likelihood (Equation 2 in Tao et al 2017)
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr select
#' @importFrom dplyr distinct
#'
#' @param X Character indicating the name of the expensive exposure in the dataset
#' @param Z Character indicating the name of the inexpensive covariates in the dataset
#' @param id Column name with the \code{id} variable.
#' @param Bspline Matrix of B-spline basis used for the sieve approximation of P(X|Z)
#' @param n_sieve Number of sieves used for the B-spline approximation of P(X|Z)
#' @param thetas Vector of parameters
#' @param sData.n1 List of subject-specific data (includes only subjects who were not sample for expensive exposure ascertainment)
#' @param sData.n2 List of subject-specific data (includes only subjects who were sample for expensive exposure ascertainment)
#' @param theta.length Length of the vector of parameters theta
#' @param n1 Number of subjects who were not sampled for expensive exposure ascertainment
#' @param m Number of distinct values of the expensive exposure observed in phase two
#' @param Q Number of quadrature point to be used if there is a latent term in the dependence model
#' @param p p_{ij} used to approximate P(X|Z)
#' @param data.n1 Dataset including only subjects not sampled for phase two
#' @param data.n2 Dataset including only subjects sampled for phase two
#'
#' @keywords internal
#'
#' @return Return the observed data log-likelihood for the sample
#' @export
#'
#'
obs.loglik <- function(X, Z = NULL, id, Bspline = NULL, n_sieve = NULL,
                       thetas, sData.n1, sData.n2, theta.length,
                       n1, m, Q, p, data.n1, data.n2){

  ## define subject-specific data for the those selected in phase 1 and those who were not
  dat.inex.n2 <- data.n2 %>% select(id, X, Z, colnames(Bspline), "k") %>% distinct()
  dat.inex.n1 <- data.n1 %>% select(id, X, Z, colnames(Bspline), "k") %>% distinct()

  ## Get the spline terms for the likelihood computation
  if(is.null(Z) | is.null(Bspline)){
    splines  <- rep(1, n1)
  } else {
    splines  <- dat.inex.n1[,colnames(dat.inex.n1) %in% colnames(Bspline)]
  }

  ## compute the part of the likelihood for those subjects sampled in phase 2
  # start with ptheta
  log.ptheta.n2      <- ptheta.fun(sData = sData.n2,
                                   X = NULL, params = thetas, p.length = theta.length,
                                   m = m, d = data.n2, Q = Q)

  log.ptheta.n2      <- sum(log.ptheta.n2$ptheta)
  # then get Blogp for all splines and for all plausible x
  logp               <-  log(p)
  logp[logp == -Inf] <- 0

  if(is.null(Z) | is.null(Bspline)){
    Blogp    <- matrix(table(dat.inex.n2[,X]), ncol = 1)*logp
  } else {
    Blogp              <- rowsum(x = dat.inex.n2[,colnames(Bspline)],
                                 group = dat.inex.n2$k, reorder = TRUE)*logp
  }

  loglik.n2          <- log.ptheta.n2 +  sum(Blogp)

  ## compute the part of the likelihood for those subjects not sampled in phase two
  ptheta.n1   <- ptheta.fun(sData = sData.n1,
                            X = NULL, params = thetas, p.length = theta.length,
                            m = m, d = data.n1, Q = Q)

  ptheta.n1   <- exp(ptheta.n1$ptheta)
  lik.n1      <- rowSums(rowsum(c(ptheta.n1*p[dat.inex.n1[, "k"], ])*splines,
                                group = rep(1:n1, times = m)))
  loglik.n1   <- log(lik.n1)
  loglik.n1[loglik.n1 == -Inf] <- 0

  ## sum the two terms
  loglik.all <- loglik.n2 + sum(loglik.n1)

  ## return log-likelihood
  return(loglik.all)
}


#' Compute the profile likelihood needed to estimate the coefficients' standard errors
#'
#' @importFrom dplyr select
#' @importFrom dplyr distinct
#' @importFrom magrittr %>%
#'
#' @param X Character indicating the name of the expensive exposure in the dataset
#' @param Z Character indicating the name of the inexpensive covariates in the dataset
#' @param Bspline Matrix of B-spline basis used for the sieve approximation of P(X|Z)
#' @param n_sieve Number of sieves used for the B-spline approximation of P(X|Z)
#' @param sData.n1 List of subject-specific data (includes only subjects who were not sample for expensive exposure ascertainment)
#' @param sData.n2 List of subject-specific data (includes only subjects who were sample for expensive exposure ascertainment)
#' @param theta.length Length of the vector of parameters theta
#' @param id Character indicating the name of the subjects id variable in the dataset
#' @param data.n1 Dataset including only subjects not sampled for phase two
#' @param data.n2 Dataset including only subjects sampled for phase two
#' @param thetas Vector of parameters
#' @param p.est p_{ijk} estimated from the EM algorithm
#' @param n1 number of subjects not sampled for expensive exposure ascertainment
#' @param m Number of distinct values of the expensive exposure observed in phase two
#' @param Q Number of quadrature points for the estimation of the random effect in the latent component of the dependence model.
#' @param maxiter Maximum number of iterations for the EM algorithm
#' @param tol Tolerance to decide on algorithm convergence
#'
#' @keywords internal
#' @return Return the profile log-likelihood
#' @export
#'

pl.loglik <- function(X, Z = NULL, Bspline = NULL, n_sieve = NULL,
                           sData.n1, sData.n2, theta.length, id,
                           data.n1, data.n2, thetas,
                           p.est, n1, m, Q,
                           maxiter = 100, tol = 1e-4){

  # define subject-specific data for the those selected in phase 1 and those who were not
  dat.inex.n2 <- data.n2 %>% dplyr::select(id, X, Z, colnames(Bspline), "k") %>% distinct() %>%
    as.data.frame()
  dat.inex.n1 <- data.n1 %>% dplyr::select(id, X, Z, colnames(Bspline), "k") %>% distinct() %>%
    as.data.frame()

  ## Define constant quantities that do not change in either the E-step or the M-step
  # to compute \psi_{kji} I need the b-spline basis of those not sampled
  # ordered by X and within X we have the ID
  if(is.null(Z) | is.null(Bspline)){
    splines  <- rep(1, n1)
    p.kj     <- matrix(table(dat.inex.n2[,X]), ncol = 1)
  } else {
    splines  <- dat.inex.n1[,colnames(dat.inex.n1) %in% colnames(Bspline)]
    # p_{kj} for those who have been sampled.
    p.kj     <- rowsum(x = dat.inex.n2[,colnames(Bspline)], group = dat.inex.n2$k, reorder = TRUE)
  }

  # first get p_theta for each value of k: this is constant through the computation
  ptheta     <- ptheta.fun(sData = sData.n1, X = X, params = thetas, p.length = theta.length,
                            m = m, Q = Q, d = data.n1)


  ptheta     <- exp(ptheta$ptheta)

  ## Initialise values for the algorithm
  p       <- p.est
  # convergence flag and number of steps
  CONVERGED <- FALSE
  nstep     <- 1


  ## Start the EM algorithm
  while(nstep <= maxiter & CONVERGED == FALSE){


    ######################################## E-step ############################################
    # compute \psi_{kji} for unsampled subjects
    # then get Bp
    p.mat   <- p[rep(seq_len(m), each =  n1), ]
    Bp      <- splines*p.mat

    # the numerator is the product of ptheta and pB
    psi.num <- ptheta*Bp    ## size (m*n1) x n.splines

    # get the denominator:
    # First, for each ID sum the quantity in the numerator over all possible values k = 1, ...m of X
    psi.den <- rowsum(psi.num, group = rep(seq(1, n1), times = m))
    # Then sum of the splines (this gives one quantity per person)
    psi.den <- rowSums(psi.den)
    # avoid division by 0
    psi.den[psi.den == 0] <- 1
    # divide
    psi <- psi.num/psi.den

    ######################################## M-step ############################################

    # update p_{kj}
    p.num <- p.kj + rowsum(psi, group = rep(seq(1, m), each = n1), reorder = TRUE)
    p.den <- matrix(rep(colSums(p.num), m), ncol = n_sieve, nrow = m, byrow = T)
    new.p <- p.num/p.den

    # check for convergence for each parameters of interest
    if(all(abs(p - new.p) < tol)){CONVERGED = TRUE}
    # update the parameters
    p       <- new.p

    # if not converge increase the number of steps
    nstep <- nstep + 1
  }

  if(nstep > maxiter){warning("Maximum number of iteration reached without convergence")}

  # set the coefficient to be used in the profile likelihood
  th   <- thetas
  np   <- new.p

  # compute the profile log-likelihood using the observed data likelihood
  p.lik  <- obs.loglik(X = X, Z = Z, id = id, Bspline = Bspline, n_sieve = n_sieve, thetas = th,
                       sData.n1 = sData.n1, sData.n2 = sData.n2, theta.length = theta.length,
                       m = m, Q = Q, p = np, n1 = n1,
                       data.n1 = data.n1, data.n2 = data.n2)

  return(p.lik)
}




#' Sieve Maximum Likelihood Estimator (SMLE) for Two-Phase Designs and Marginalized Transition and Latent Variable Models
#'
#' Estimate the parameters of interest
#'
#' @importFrom dplyr distinct
#' @importFrom dplyr group_by
#' @importFrom dplyr summarise
#' @importFrom dplyr pull
#'
#' @param Y Column name for outcome variable as recorded in the dataset. This is a required argument
#' @param X Column name for expensive exposure variable as recorded in the dataset. This is a required argument
#' @param Z Vector of columns name for inexpensive variables as recorded in the dataset. This is an optional argument
#' @param Time Column name for the time variable as recorded in the dataset. This is a required argument
#' @param Bspline Matrix with B-spline basis. This is a required argument.
#' @param n_sieve Number of sieves for the B-spline approximation
#' @param mean.formula Formula for the marginal mean model in which a binary variable is regressed on a set of covariates. This is a required argument
#' @param t.formula  Formula for the transition part in the dependence model. This is optional and of the form ~ covariates
#' @param lv.formula  Formula for the latent part in the dependence model. This is optional and of the form ~ covariates
#' @param id Column name for the \code{id} variable as recorded in the dataset. This is a required argument
#' @param data Dataframe containing \code{Y}, \code{id}, \code{Z} and \code{X}. This is a required argument
#' @param Q Number of quadrature points for the estimation of the random effect in the latent component of the dependence model. Default is set to 10
#' @param maxiter Maximum number of iterations for the SMLE algorithm. Maximum is set to 100. This is an optional argument
#' @param tol  Tolerance for the convergence criterion in the EM algorithm. The default value is \code{1E-4}. This is an optional argument
#' @param no_SE  If \code{TRUE}, the variances of the parameter estimators will not be estimated. The default value is \code{FALSE}. This is an optional argument
#' @param verbose If \code{TRUE}, show estimated coefficients at each step of the EM algorithm. The default value is \code{FALSE}. This is an optional argument
#' @param hnscale  Specifies the scale of the perturbation constant in the variance estimation. Default is set to 1. This is an optional argument
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
#' # set up P(X|Z)
#' Bspline_Z <- cbind(as.numeric(exampledat$Z == 0), as.numeric(exampledat$Z == 1))
#' n_sieve             <- ncol(Bspline_Z)
#' colnames(Bspline_Z) <- paste("bs", 1:n_sieve, sep="")
#' exampledat          <- cbind(exampledat, Bspline_Z)
#' # estimate the parameters using SMLE
#' mod.smle <- smle_mm(mean.formula = Y ~ X + time + Z + X:time, t.formula = ~1, lv.formula = NULL,
#' Y = "Y", X = "X", Z = "Z", Time = "time",
#' id = "id", Bspline = Bspline_Z, n_sieve = n_sieve, no_SE = FALSE,
#' data = exampledat)
#' summary(mod.smle)
#' }
smle_mm <- function(Y, X, Z = NULL, Time, Bspline = NULL, n_sieve = NULL,
                    mean.formula = NULL,
                    t.formula = NULL,
                    lv.formula = NULL, id, data, Q = 10,
                    maxiter = 100, tol = 1e-4, no_SE = TRUE, verbose = FALSE,
                    hnscale = 1){

  # check the inputs are given correctly
  if(is.null(mean.formula)){stop("The mean model needs to be specified")}
  if(is.null(t.formula) & is.null(lv.formula)){
    stop("Either transition or latent model needs to be specified")}

  if(is.null(X)){stop("The expensive covariate X needs to be provided")}
  if(is.null(Y)){stop("The outcome variable Y needs to be provided")}
  if(is.null(Time)){stop("The time variable needs to be provided")}
  if(is.null(id)){stop("The id variable needs to be provided")}

  if(!is.null(X) & !is.character(X)){stop("X needs to be a character")}
  if(!is.null(Y) & !is.character(Y)){stop("Y needs to be a character")}
  if(!is.null(Time) & !is.character(Time)){stop("Time needs to be a character")}
  if(!is.null(Z) & !is.character(Z)){stop("Z needs to be a character")}
  if(!is.null(id) & !is.character(id)){stop("id needs to be a character")}

  # change the id names in the dataset as id.
  names(data)[names(data) == id]     <- "id"

  # Define some constant values

  N        <- length(unique(data$id))                            ## phase 1 sample size
  n2       <- length(unique(data[!is.na(data[,X]),"id"]))        ## phase 2 sample size (sampled subjects)
  n1       <- N - n2                                             ## number of unsampled subjects

  n.var    <- ncol(model.matrix(mean.formula, data))               ## covariates in mean model
  if(!is.null(t.formula)){                                ## covariates in the transition model
    n.tr <- ncol(model.matrix(t.formula, data))
  } else {
    n.tr <- 0
  }
  if(!is.null(lv.formula)){                                    ## covariates in the latent model
    n.lt <- ncol(model.matrix(lv.formula, data))
  } else {
    n.lt <- 0
  }

  # total number of parameters
  n.par    <- n.var + n.tr + n.lt
  par.vec  <- c(n.var, n.tr, n.lt)

  # spline parameters
  if(is.null(Bspline)){
    n_sieve = 1
  }

  ## Reorder dataset so that the first rows are those who have been selected for phase two
  # no inexpensive covariates
  if(is.null(Z)){
    dat.all           <- data.frame(id = data$id, data[,Y], data[,X], data[,Time])
    colnames(dat.all) <- c("id", Y, X, Time)
  }
  # inexpensive covariates all independent of the X. Set the data so can get p(X)
  if(!is.null(Z) & is.null(Bspline)){
    dat.all           <- data.frame(id = data$id, data[,Y], data[,X], data[,Time],
                                    data[,Z])
    colnames(dat.all) <- c("id", Y, X, Time, Z)
  }
  # inexpensive covariates with some being correlated to X, Set the data so can get p(X|Z)
  if(!is.null(Z) & !is.null(Bspline)){
    dat.all           <- data.frame(id = data$id, data[,Y], data[,X], data[,Time],
                                    data[,Z], Bspline)
    colnames(dat.all) <- c("id", Y, X, Time, Z, colnames(Bspline))
  }

  dat.all  <- dat.all[order(dat.all[,X], dat.all$id),]

  ## Create two datasets one for those who have been sampled and one for those who have not been sampled
  # number of observations for whom we have the expensive covariate
  n2.obs   <- sum(!is.na(dat.all[,X]))
  # data set for those with the expensive covariate
  dat.n2   <- dat.all[1:n2.obs,]
  # data set for those without the expensive covariate
  dat.n1   <- dat.all[-c(1:n2.obs),]
  n1.obs   <- nrow(dat.n1)

  ## Summarise observed distinct values for X and create the augmented dataset we need in the M-step
  x.obs      <- sort(unique(na.omit(data[,X])))
  m          <- length(x.obs)
  # create all possible X for those subjects not selected for phase 2
  dat.n1     <- dat.n1[rep(seq_len(n1.obs), each = m), ]
  dat.n1[,X] <- x.obs
  # order by X first, within X order by id
  dat.n1     <- dat.n1[order(dat.n1[,X], dat.n1$id),]
  # for each person assign a value k depending on which value X takes
  dat.n1$k   <- rep(1:m, each = n1.obs)
  dat.n2$k   <- rep(1:m, table(dat.n2[,X]))
  # stack together data for those being sampled with data for those not being sampled
  dat.stack  <- rbind(dat.n2, dat.n1)
  # rename expensive covariate to make the code easier
  names(dat.stack)[names(dat.stack) == X]     <- "X"
  # create new ID by grouping ID and group (this is for mm to work correctly)
  #dat.stack  <- dat.stack %>% group_by(id, X) %>% mutate(id.new = cur_group_id()) %>%
  #  ungroup %>% as.data.frame()

  # define subject-specific data for the those selected in phase 1 and those who were not
  # Subjects specific data are needed for the computation of weights in the EM
  dat.inex.n2 <- dat.n2 %>% select(id, X, Z, colnames(Bspline), "k") %>% distinct() %>%
    as.data.frame()
  dat.inex.n1 <- dat.n1 %>% select(id, X, Z, colnames(Bspline), "k") %>% distinct() %>%
    as.data.frame()

  ## Define constant quantities that do not change in either the E-step or the M-step
  if(is.null(Z) | is.null(Bspline)){
    # set a spline vector to 1. This is just to ease the coding in the E-step
    splines  <- rep(1, n1)
    # p_{kj} for those who have been sampled. This is also the initial value of p.
    p.kj     <- matrix(table(dat.inex.n2[,X]), ncol = 1)
    p.kj.den <- sum(p.kj)
    # P(X) across those who have been sampled (to use id X and Z independent)
    p.samp   <- p.kj/p.kj.den
  } else {
    # spline terms for those who have not been sampled (will helps computing the weights in the E-step)
    splines  <- dat.inex.n1[,colnames(dat.inex.n1) %in% colnames(Bspline)]
    # p_{kj} for those who have been sampled. This is also the initial value of p.
    p.kj     <- rowsum(x = dat.inex.n2[,colnames(Bspline)], group = dat.inex.n2$k, reorder = TRUE)
    p.kj.den <- matrix(rep(colSums(p.kj), m), nrow = m, byrow = T)
    # P(X|z) across those who have been sampled
    p.samp   <- p.kj/p.kj.den
  }

  ## Initialise values for the algorithm
  #thetas  <- rep(0, n.var + length(n.tr) + length(n.lt))

  if(n.lt != 0 & n.tr == 0){
    thetas  <- rep(0, n.var + length(n.tr) + n.lt)
  }
  if(n.lt == 0 & n.tr != 0){
    thetas  <- rep(0, n.var + n.tr + length(n.lt))
  }
  if(n.lt != 0 & n.tr != 0){
    thetas  <- rep(0, n.var + n.tr + n.lt)
  }

  p       <- p.samp
  # convergence flag and number of steps
  CONVERGED  <- CONVERGED.THETA <- CONVERGED.P <- FALSE
  nstep      <- 1

  # get the element of Ptheta that do not change with the EM iterations
  # Subject data necessary to compute the likelihood
  # Unsampled subjects
  p.theta.init.n1 <- ptheta.input(Y = Y, X = X, id = id,
                                  mean.formula = mean.formula, t.formula = t.formula,
                                  lv.formula = lv.formula,
                                  m = m, d = dat.n1)
  subjectDat.n1    <- p.theta.init.n1$subjectData
  subjectDat.n1.weights <- p.theta.init.n1$subjectData
  n.subjectDat.n1  <- length(subjectDat.n1)
  p.length         <- p.theta.init.n1$p.length
  # sample subjects
  p.theta.init.n2 <- ptheta.input(Y = Y, X = NULL, id = id,
                                  mean.formula = mean.formula, t.formula = t.formula,
                                  lv.formula = lv.formula,
                                  m = m, d = dat.n2)
  subjectDat.n2   <- p.theta.init.n2$subjectData

  # define elements needed for the M-step
  if(Q <= 1){Q <- 2}
  if(is.null(lv.formula)){Q <- 1}
  AdaptiveQuad <- get.GH(Q)
  W.ada <- AdaptiveQuad$w
  Z.ada <- AdaptiveQuad$z

  ## Start the EM algorithm
  while(nstep <= maxiter & CONVERGED == FALSE){

    ######################################## E-step ############################################
    # compute \psi_{kji} for unsampled subjects
    # first get p_theta for each value of k
    ptheta     <- ptheta.fun(sData = subjectDat.n1, params = thetas, p.length = p.length,
                             X = X, m = m, Q = Q, d = dat.n1)

    ptheta[,X] <- rep(x.obs, each = n1)
    ptheta     <- exp(ptheta$ptheta)

    # then get Bp
    p.mat   <- p[rep(seq_len(m), each =  n1), ]
    Bp      <- splines*p.mat

    # the numerator is the product of ptheta and pB
    psi.num <- ptheta*Bp    ## size (m*n1) x n.splines

    # get the denominator:
    # First, for each ID sum the quantity in the numerator over all possible values k = 1, ...m of X
    psi.den <- rowsum(psi.num, group = rep(seq(1, n1), times = m))
    # Then sum of the splines (this gives one quantity per person)
    psi.den <- rowSums(psi.den)
    # avoid division by 0
    psi.den[psi.den == 0] <- 1

    # divide
    psi <- psi.num/psi.den
    # compute q for the unsampled subjects summing over all possible sieves
    if(is.null(Z) | is.null(Bspline)){
      q <- psi
    } else {
      q   <- rowSums(psi)
    }

    # save q for each id
    q.dat     <- data.frame(id = unique(dat.n1$id), q = q)
    q.dat[,X] <- rep(x.obs, each = n1)
    q.dat$weights  <- 1/q.dat$q
    weights        <- split(q.dat$weights, list(q.dat[,id], q.dat[,X]))
    for(ls in 1: length(subjectDat.n1.weights)){
      subjectDat.n1.weights[[ls]]$SampProbi <- weights[[ls]]
    }
    sData.all <- append(subjectDat.n1.weights, subjectDat.n2)

    ######################################## M-step ############################################

    # fit a weighted marginalised model to update the parameters
    mod.smle  <- nlm(LogLScoreCalc, p = thetas, subjectData = sData.all,
                       CondLike = FALSE, ParamLengths = p.length,
                       Q = Q, AdaptiveQuad = FALSE,
                       W = W.ada, Z = Z.ada,
                       print.level = 0, stepmax = 1, steptol = 1e-6,
                       hessian = FALSE, iterlim = 500)
    new.thetas <- mod.smle$estimate
    loglik.m   <- -mod.smle$minimum

    # update p_{kj}
    p.num <- p.kj + rowsum(psi, group = rep(seq(1, m), each = n1), reorder = TRUE)
    p.den <- matrix(rep(colSums(p.num), m), ncol = n_sieve, nrow = m, byrow = T)
    new.p <- p.num/p.den

    # check for convergence for each parameters of interest
    if(all(abs(thetas - new.thetas) < tol)){CONVERGED.THETA = TRUE}
    if(all(abs(p - new.p) < tol)){CONVERGED.P = TRUE}
    # check for overall convergence
    if(CONVERGED.THETA == TRUE & CONVERGED.P == TRUE){CONVERGED = TRUE}

    # update the parameters
    thetas  <- new.thetas
    p       <- new.p

    if(verbose== TRUE){print(list(step = nstep, params = thetas))}

    # if not converge increase the number of steps
    nstep <- nstep + 1
  }

  if(nstep > maxiter){warning("Maximum number of iteration reached without convergence")}

  if(no_SE == TRUE){
    # compute the maximum cluster size
    maxtime <- max(dat.all %>% group_by(id) %>% summarise(mtime = max(dat.all[,Time])) %>% pull())
    # name variables for printing
    nms       = list()
    nms$beta  = colnames(model.matrix(mean.formula, data))
    nms$alpha = c(if(!is.null(t.formula)){paste('gamma',colnames(model.matrix(t.formula, data)),sep=':')},
                  if(!is.null(lv.formula)){paste('log(sigma)',colnames(model.matrix(lv.formula, data)),sep=':')})
    out <- list(call       = match.call(), alpha = new.thetas[(n.var + 1):length(new.thetas)],
                beta       = new.thetas[1:n.var],
                mod.cov    = matrix(NA, nrow = length(new.thetas), ncol = length(new.thetas)),
                info_stats = loglik.m,
                control    = c(CONVERGED*1, nstep - 1, N, no_SE, 0, maxtime))
    names(out$beta)        = nms$beta
    names(out$alpha)       = nms$alpha
    names(out$control)     = c("convergence_code", "n_iter", "n_subj", "se_comp", "useRobCov", 'max_n_visit')
    names(out$info_stats)  = c('logLik')
  }

  if(no_SE == FALSE){
    # do not perturb parameters
    pl.no <- obs.loglik(X = X, Z = Z,
                        Bspline = Bspline, n_sieve = n_sieve,
                        id = id, m = m, thetas = thetas,
                        theta.length = p.length, n1 = n1,
                        sData.n1 = subjectDat.n1, sData.n2 = subjectDat.n2, Q = Q,
                        p = p, data.n1 = dat.n1, data.n2 = dat.n2)

    # perturb one parameter at the time
    betas      <- new.thetas[1:n.var]

    if(is.null(lv.formula)){
      gammas     <- new.thetas[(n.var + 1):(n.var + n.tr)]
      new.thetas <- c(betas, gammas)
    }

    if(is.null(t.formula)){
      sigmas     <- new.thetas[(n.var + 2):length(new.thetas)]
      new.thetas <- c(betas, sigmas)
    }

    if(!is.null(t.formula) & !is.null(lv.formula)){
      gammas     <- new.thetas[(n.var + 1):(n.var + n.tr)]
      sigmas     <- new.thetas[(n.var + n.tr + 1):length(new.thetas)]
      new.thetas <- c(betas, gammas, sigmas)
    }

    npar       <- length(new.thetas)
    hn         <- hnscale*(1/sqrt(N))
    eps.params <- diag(rep(hn, npar))

    pl.est <- c()
    for(k in 1:npar){
      t.pert     <- new.thetas + eps.params[k, ]

      if(is.null(lv.formula)){
        t.pert     <- c(t.pert, 0)
      }
      if(is.null(t.formula)){
        t.pert     <- append(t.pert, values = 0, after = n.var)
      }
      pl.est[k]  <- pl.loglik(X = X, Z = Z, Bspline = Bspline, n_sieve = n_sieve,
                                   sData.n1 = subjectDat.n1, sData.n2 = subjectDat.n2,
                                   data.n1 = dat.n1, data.n2 = dat.n2,
                                   theta.length = p.length, id = id,
                                   thetas = t.pert, p.est = p, Q = Q,
                                   n1 = n1, m = m, maxiter = 1000, tol = tol)
    }

    # perturb all possible combination of parameters
    pl.mat <- matrix(NA, ncol = npar, nrow = npar)
    for(w in 1:npar){
      for(k in 1:w){
        t.pert      <- new.thetas + eps.params[w, ] + eps.params[k, ]

        if(is.null(lv.formula)){
          t.pert     <- c(t.pert, 0)
        }
        if(is.null(t.formula)){
          t.pert     <- append(t.pert, values = 0, after = n.var)
        }

        pl.mat[w,k] <- pl.loglik(X = X, Z = Z, Bspline = Bspline, n_sieve = n_sieve,
                                      sData.n1 = subjectDat.n1, sData.n2 = subjectDat.n2,
                                      data.n1 = dat.n1, data.n2 = dat.n2,
                                      theta.length = p.length, id = id,
                                      thetas = t.pert, p.est = p, Q = Q,
                                      n1 = n1, m = m, maxiter = 1000, tol = tol)
        pl.mat[k,w] <- pl.mat[w,k]
      }
    }

    # compute the diagonal element of the Hessian
    H.diag  <- c()
    for(m in 1:npar){
      H.diag[m]  <- pl.mat[m,m] - 2*pl.est[m] + pl.no
    }
    Hes     <- diag(H.diag)

    # compute the off diagonal elements
    H.offdiag      <- mat.row <- NULL
    for(w in 1:(npar - 1)){
      for(k in (w + 1):npar){
        mat.el   <- pl.mat[w, k] - pl.est[w] - pl.est[k] + pl.no
        mat.row  <- c(mat.row, mat.el)
      }
    }
    Hes[lower.tri(Hes, diag = FALSE)] <- mat.row
    Hes[upper.tri(Hes)] <- t(Hes)[upper.tri(t(Hes))]
    Hes                 <- Hes/(hn^2)
    covariance          <- solve(-Hes)
    se                  <- sqrt(diag(covariance))

    # compute the maximum cluster size
    maxtime <- max(dat.all %>% group_by(id) %>% summarise(mtime = max(dat.all[,Time])) %>% pull())

    # name variables for printing
    nms       = list()
    nms$beta  = colnames(model.matrix(mean.formula, data))
    nms$alpha = c(if(!is.null(t.formula)){paste('gamma',colnames(model.matrix(t.formula, data)),sep=':')},
                  if(!is.null(lv.formula)){paste('log(sigma)',colnames(model.matrix(lv.formula, data)),sep=':')})
    out <- list(
      call       = match.call(),
      alpha      = new.thetas[(n.var + 1):length(new.thetas)],
      beta       = new.thetas[1:n.var],
      mod.cov    = covariance,
      info_stats = loglik.m,
      control    = c(CONVERGED*1, nstep - 1, N, no_SE, 0, maxtime))
    names(out$beta)        = nms$beta
    names(out$alpha)       = nms$alpha
    names(out$control)     = c("convergence_code", "n_iter", "n_subj", "se_comp", "useRobCov", 'max_n_visit')
    names(out$info_stats)  = c('logLik')

  }

  class(out)  = "TwoPhaseMM"
  return(out)
}

