#' Transform the odds to a probability
#'
#' @param odds Odds
#'
#' @keywords internal
#' @return Probability
#'
Odds2Prob      <- function(odds){
  odds/(1+odds)
}

#' Prepare the data for the indirect imputation approach
#'
#' @param mean.formula Formula for the marginal mean model in which a binary variable is regressed on a set of covariates.
#' @param lv.formula  Formula for the latent part in the dependence model. This is optional and of the form ~ covariates
#' @param t.formula  Formula for the transition part in the dependence model. This is optional and of the form ~ covariates
#' @param id Vector of \code{id} variable. This is a required argument
#' @param data Dataframe. This is a required argument
#' @param inits  An optional list of length 3 containing initial values for marginal mean parameters and all dependence parameters. The format of the list should be: (1) estimates of the mean parameters, (2) estimates of the transition parameters (or NULL if only fitting a mLV model) and (3) estimates of the latent variable parameters (or NULL if only fitting a mT model). If NULL, initial values will be automatically generated. Default is set to NULL
#' @param samp.probs Vector of 3 numbers indicating the probability that a subject is sampled from the none
#'                   some and all stratum respectively. If not specified, the sampling probabilities are set to 1
#' @param samp.probi Vector of sampling probabilities for each individual
#' @param offset an optional offset term.
#'
#' @keywords internal
#'
#'
prep_data      <- function(mean.formula, lv.formula, t.formula, id, data, inits,
                           samp.probs, samp.probi, offset) {

  mean.f = model.frame(mean.formula, data)
  mean.t = attr(mean.f, "terms")

  y  = model.response(mean.f,'numeric')
  uy = unique(y)
  x  = model.matrix(mean.formula,mean.f)

  x.t = x.lv = matrix(0, ncol=1, nrow=length(y))
  if(!is.null(t.formula))   x.t  = model.matrix(t.formula,model.frame(t.formula, data))
  if(!is.null(lv.formula))  x.lv = model.matrix(lv.formula, model.frame(lv.formula, data))

  if(!is.matrix(samp.probs)) samp.probs = matrix(samp.probs, nrow = length(y), ncol = 3, byrow = TRUE)
  if(is.null(samp.probi)) samp.probi = matrix(1,nrow=length(y),ncol=1)
  if(!is.null(inits)) inits = unlist(inits)


  if(is.null(inits)) {
    inits = c(glm(mean.formula, family='binomial', data=data)$coef, rep(1, ncol(x.t) + ncol(x.lv)))
    if(any(is.na(inits))) {
      omit_dup_col = which(is.na(inits))
      x            = x[,-c(omit_dup_col)]
      inits        = inits[-c(omit_dup_col)]
    }
  }


  if(is.null(offset)) offset = rep(0, length(y))

  x.t = cbind(x.t, NULL)
  x.lv = cbind(x.lv, NULL)

  if ( length(inits) != ncol(x) + ncol(x.t) + ncol(x.lv) ) {
    stop("Parameter length incongruous with X, Xgam, and Xsig")}
  paramlengths <- c(ncol(x), ncol(x.t), ncol(x.lv))

  Ylag = rep(0, length(y))
  for(i in 2:length(y)) { if(id[i]==id[i-1]) Ylag[i]<-y[i-1] }

  id.tmp        = split(id, id)
  X.tmp         = split(x,id)
  Y.tmp         = split(y,id)
  Ylag.tmp      = split(Ylag,id)
  Xgam.tmp      = split(x.t,id)
  Xsig.tmp      = split(x.lv,id)
  SampProbi.tmp = split(samp.probi,id)
  SampProbs.tmp = split(samp.probs,id)
  offset.tmp    = split(offset,id)

  subjectData <- vector('list', length=length(unique(id)))
  subjectData <- list()
  uid <- as.character(unique(id))
  for(j in seq(along=uid)){
    i <- uid[j]
    subjectData[[j]] <- list(id=as.character(unique(id.tmp[[i]])),
                             X=matrix(X.tmp[[i]], ncol=ncol(x)),
                             Y=as.double(Y.tmp[[i]]),
                             Ylag=as.double(Ylag.tmp[[i]]),
                             Xgam=matrix(Xgam.tmp[[i]], ncol=ncol(x.t)),
                             Xsig=matrix(Xsig.tmp[[i]], ncol=ncol(x.lv)),
                             SampProbi=unique(SampProbi.tmp[[i]]),
                             SampProbs=matrix(SampProbs.tmp[[i]], ncol=ncol(samp.probs))[1,],
                             Offset=as.double(offset.tmp[[i]]))
  }
  names(subjectData) <- uid

  beta_nms  <- colnames(x)
  alpha_nms <- c( if(!is.null(t.formula)){paste('gamma',colnames(x.t),sep=':')},
                  if(!is.null(lv.formula)){paste('log(sigma)',colnames(x.lv),sep=':')})

  list(params=inits, paramlengths=paramlengths, subjectData=subjectData, samp.probi=tapply(samp.probi, id, unique),
       nms=list(beta_nms, alpha_nms))
}



#' Create imputed dataset for the indirect imputation approach
#'
#' @param FIT TwoPhaseMM object
#' @param DAT complete dataset with sampled flag (Sampled)
#' @param M Total number of imputations. The default number is set to 10.
#' @param ROB Indicator for the variance estimation. If FALSE use model-based covariance, otherwise use sandwich estimator. Default is set to FALSE.
#' @param MARG.EXP.FORMULA marginal exposure model formula (X | Z). This is the imputation model without the offset
#' @param verbose Indicator. If yes print the number of imputed dataset. Default is set to FALSE.
#' @param EXPVAR exposure / expensive variable: Xe
#' @param Sampled should be set to 1 if included in the sample and 0 otherwise
#'
#' @importFrom MASS mvrnorm
#' @importFrom stats coef
#' @importFrom stats vcov
#' @importFrom stats offset
#'
#' @keywords internal
#' @return Return the imputed dataset
#' @export
#'

ImputeDataInd   <- function(FIT, DAT, M=10, ROB=FALSE, MARG.EXP.FORMULA,
                            verbose=FALSE,EXPVAR=EXPVAR, Sampled=Sampled) {

  imp_data <- vector('list', M)

  for(ix in seq(M)) { if (verbose) cat("\r","Imputing dataset number:", ix)
    # Step 0: Update date files
    # Create versions of DAT with updated Xe values
    DAT_Xe1 = DAT_Xe0 = DAT

    DAT_Xe1[,EXPVAR][is.na(DAT_Xe1[,EXPVAR])] = 1
    DAT_Xe0[,EXPVAR][is.na(DAT_Xe0[,EXPVAR])] = 0
    DAT_Xe1[,EXPVAR] = 1
    DAT_Xe0[,EXPVAR] = 0


    # Non-sampled subjects
    DAT_S0     = DAT[Sampled==0,]
    DAT_S0_Xe1 = DAT_Xe1[Sampled==0, ]
    DAT_S0_Xe0 = DAT_Xe0[Sampled==0, ]
    dup        = duplicated(DAT_S0$id)
    DAT_S0_v0  = DAT_S0[!dup, ]

    # Sampled subjects
    DAT_S1     = DAT[Sampled==1,]
    DAT_S1_Xe1 = DAT_Xe1[Sampled==1, ]
    DAT_S1_Xe0 = DAT_Xe0[Sampled==1, ]
    dup        = duplicated(DAT_S1$id)
    DAT_S1_v0  = DAT_S1[!dup,]

    # Step 1: extract estimates from FIT and draw m^th theta from MVN
    theta      = c(FIT$beta, FIT$alpha)
    if(ROB) { cov_theta = FIT$rob.cov
    } else { cov_theta = FIT$mod.cov}

    theta.m    = MASS::mvrnorm(n = 1, mu=theta, Sigma=cov_theta)
    prep_input = attr(FIT,'args')

    # when fitting either a transition or a latent term only in the dependence model, we need
    # to add zeros to theta.m
    if(is.null(prep_input$lv.formula)){
      theta.m = c(theta.m, 0)
    }
    if(is.null(prep_input$t.formula)){
      ncol.x  = ncol(model.matrix(prep_input$mean.formula, data = DAT))
      theta.m = c(theta.m[1:ncol.x], 0, theta.m[(ncol.x + 1):length(theta.m)])
    }


    # Step 2: Using the estimated theta=(beta,alpha), calculate likelihood contributions for sampled subjects
    # for both Xe=1 and Xe=0.
    tmp_S1_Xe1 = prep_data(mean.formula=prep_input$mean.formula, lv.formula=prep_input$lv.formula,
                            t.formula=prep_input$t.formula,
                            id=DAT_S1_Xe1[,prep_input$id], data=DAT_S1_Xe1, inits=theta.m,
                            samp.probs=prep_input$samp.probs, samp.probi=prep_input$samp.probi,
                            offset=prep_input$offset)$subjectData
    tmp_S1_Xe0 = prep_data(mean.formula=prep_input$mean.formula, lv.formula=prep_input$lv.formula,
                            t.formula=prep_input$t.formula,
                            id=DAT_S1_Xe0[,prep_input$id], data=DAT_S1_Xe0, inits=theta.m,
                            samp.probs=prep_input$samp.probs, samp.probi=prep_input$samp.probi,
                            offset=prep_input$offset)$subjectData

    # These functions calculate the conditional likelihood contributions by and the ascertainment correction
    # for sampled subjects.  The difference of ascertainment corrections (i.e. log of the ratio
    # log(pr(S=1 | Xe=1, Xo)/pr(S=1 | Xe=1, Xo))= log(pr(S=1 | Xe=1, Xo))-log(pr(S=1 | Xe=0, Xo))) provides
    # the offset for the marginal exposure model
    LLSC_1 = LogLScoreCalc(params = theta.m, subjectData = tmp_S1_Xe1,
                            Q = FIT$LLSC_args$Q, W = FIT$LLSC_args$W, Z = FIT$LLSC_args$Z,
                            ParamLengths = FIT$LLSC_args$ParamLengths,
                            CondLike = FIT$control['cond.like']==1, EmpiricalCheeseCalc = ROB)
    LLSC_0 = LogLScoreCalc(params = theta.m, subjectData = tmp_S1_Xe0,
                            Q = FIT$LLSC_args$Q, W = FIT$LLSC_args$W, Z =  FIT$LLSC_args$Z,
                            ParamLengths = FIT$LLSC_args$ParamLengths,
                            CondLike = FIT$control['cond.like']==1, EmpiricalCheeseCalc = ROB)

    # Note: The ascertainment corrections used for the marginal exposure model offset
    #       are already on the log scale, so log(a/b) = loga-logb
    logAC_1 = attr(LLSC_1, "ACSubj")
    logAC_0 = attr(LLSC_0, "ACSubj")
    DAT_S1_v0$offset = logAC_1-logAC_0

    # Perform offsetted logistic regression [Note, only using first visit info i.e. no time-varying information.
    fit.exp = glm(MARG.EXP.FORMULA, offset=offset, data=DAT_S1_v0, family = "binomial")

    # Extract alpha estimates and sample from the sampling distribution
    alpha      = coef(fit.exp)
    cov_alpha  = vcov(fit.exp)

    alpha.m = mvrnorm(n = 1, mu=alpha, Sigma=cov_alpha)

    # Step 5: Apply the model and theta to calculate the ascertainment corrections for unsampled people,

    nr = nrow(DAT_S0_Xe1)

    tmp_S0_Xe1 = prep_data(mean.formula=prep_input$mean.formula, lv.formula=prep_input$lv.formula,
                            t.formula=prep_input$t.formula,
                            id=DAT_S0_Xe1[,prep_input$id], data=DAT_S0_Xe1, inits=theta.m,
                            samp.probs=prep_input$samp.probs[rep(1,nr),],
                            samp.probi=prep_input$samp.probi[rep(1,nr)],
                            offset=prep_input$offset[rep(1,nr)])$subjectData
    tmp_S0_Xe0 = prep_data(mean.formula=prep_input$mean.formula, lv.formula=prep_input$lv.formula,
                            t.formula=prep_input$t.formula,
                            id=DAT_S0_Xe0[,prep_input$id], data=DAT_S0_Xe0, inits=theta.m,
                            samp.probs=prep_input$samp.probs[rep(1,nr),],
                            samp.probi=prep_input$samp.probi[rep(1,nr)],
                            offset=prep_input$offset[rep(1,nr)])$subjectData

    LLSC_1 = LogLScoreCalc(params = theta.m, subjectData = tmp_S0_Xe1,
                            Q = FIT$LLSC_args$Q, W = FIT$LLSC_args$W, Z = FIT$LLSC_args$Z,
                            ParamLengths = FIT$LLSC_args$ParamLengths,
                            CondLike =  FIT$control['cond.like']==1, EmpiricalCheeseCalc = ROB)
    LLSC_0 = LogLScoreCalc(params = theta.m, subjectData = tmp_S0_Xe0,
                            Q = FIT$LLSC_args$Q, W = FIT$LLSC_args$W, Z =  FIT$LLSC_args$Z,
                            ParamLengths = FIT$LLSC_args$ParamLengths,
                            CondLike =  FIT$control['cond.like']==1, EmpiricalCheeseCalc = ROB)

    ## Applied to the non-sampled subjects: log(pr(S=1 | Xe=1, Xo)/pr(S=1 | Xe=1, Xo)).
    logAC_1 = attr(LLSC_1, "ACSubj") # Ascertainment corrections already log transformed
    logAC_0 = attr(LLSC_0, "ACSubj")
    offset_log_S0 = logAC_1-logAC_0

    # Create a temporary dataframe for non-sampled subjects in order to create the
    # model.matrix that is used to calculate marginal exposure odds
    tmp          = DAT_S0_v0
    tmp[,EXPVAR] = sample(c(0,1),size=nrow(tmp),replace=TRUE)
    mm_tmp       = model.matrix(MARG.EXP.FORMULA, tmp)

    # Compute the conditional exposure odds in non-sampled subjects
    # Notice that all conditional odds are equal: [Xe | Y, Xo, S=1] = [Xe | Y, Xo, S=0] = [Xe | Y, Xo] since we sampled based on Y
    # So, we multiply the odds (conditional on being sampled), pr[Xe=1 | Xo, S=1]/pr[Xe=0 | Xo, S=1]
    # by the likelihood ratio (conditional on being sampled), pr[Y | Xe=1, Xo, S=1]/ pr[Y | Xe=1, Xo, S=1]
    # in the unsampled subjects in order to obtain conditional odds in unsampled: pr[Xe=1 | Y, Xo, S=1]/pr[Xe=0 | Y, Xo, S=1]
    marg.odds.Xe.Xo.S1   = exp(mm_tmp %*% alpha.m)*exp(offset_log_S0)
    LR.Y.Xe.Xo.S1        = exp( attr(LLSC_1, "LogLikeSubj")-attr(LLSC_0, "LogLikeSubj"))
    cond.odds.Xe.Y.Xo.S0 =  LR.Y.Xe.Xo.S1*marg.odds.Xe.Xo.S1

    # Convert odds to probability, log(odds) = log(1/[1-p]) ---> p = odds/(1+odds)
    # Impute Xe for non-sampled subjects
    prXe1     = Odds2Prob(cond.odds.Xe.Y.Xo.S0)
    XeS0      = rbinom(nrow(tmp), size=1, prob=prXe1)

    # Add XeS0 info to DAT_S0_v0
    DAT_S0_v0[,EXPVAR] = XeS0
    DAT_S0 = DAT_S0[,-which(colnames(DAT_S0)==EXPVAR)]
    DAT_S0 = merge(DAT_S0, DAT_S0_v0[,c('id',EXPVAR)],by='id')
    DAT_imp = rbind(DAT_S1, DAT_S0[,colnames(DAT_S1)])

    # Store imputation data
    imp_data[[ix]] = DAT_imp
  } # end of for ix loop

  imp_data
}

#' Fit imputed data using indirect imputation
#'
#' @param FIT Model fit for imputation
#' @param ImputedData List of datasets where the expensive exposure has been imputed
#' @param verbose An indicator of whether the steps of the maximization algorithm should be printed. Default is set to FALSE
#' @param iter.lim  A scalar indicating the maximum number of iteration in the maximization algorithm. Default is 100
#' @param q A scalar to denote the number of quadrature points used for GH numerical integration. Default is 10.
#'
#' @keywords internal
#' @importFrom mitools MIcombine
#' @return The function returns the estimated coefficients and standard error from the indirect imputation approach
#'
FitImputedData <- function(FIT, ImputedData, id, verbose=FALSE, iter.lim = 1000, q = 10){
  M <- length(ImputedData)
  prep_input <- attr(FIT,'args')
  imp.mods <- ests <- vcovs <- list()
  for (i in 1:M){
    imp.mods[[i]]        <- wee_mm(prep_input$mean.formula,
                               data = ImputedData[[i]],
                               id = id,
                               lv.formula = prep_input$lv.formula,
                               t.formula = prep_input$t.formula, weights = rep(1, nrow(ImputedData[[i]])),
                               maxiter = iter.lim, verbose=verbose, Q=q)
    ests[[i]]  <- c(imp.mods[[i]]$beta, imp.mods[[i]]$alpha)
    vcovs[[i]] <- imp.mods[[i]]$mod.cov
  }
  out <- mitools::MIcombine(ests, vcovs)
  out
}


#' Impute the missing data and fit the complete data using the direct imputation method
#'
#' @param mean.formula Formula for the marginal mean model in which a binary variable is regressed on a set of covariates. This is a required argument
#' @param lv.formula Formula for the latent part in the dependence model. This is optional and of the form ~ covariates
#' @param t.formula Formula for the transition part in the dependence model. This is optional and of the form ~ covariates
#' @param id Column name for the id variable as recorded in the dataset. This is a required argument
#' @param data Full cohort data. This dataset should not have missing values for variables other than X. It should
#' contain columns named \code{id}, \code{X}, and \code{sampled}
#' @param sampled Column name for the sampling indicator variable as recorded in the dataset. This is a required argument
#' @param X Column name for the expensive exposure variable as recorded in the dataset. This is a required argument
#' @param samp.probs Vector of 3 numbers indicating the probability that a subject is sampled from the none
#'                   some and all stratum respectively. This is a required argument when \code{method} is set equal to indirect. It is not
#'                   used when \code{method} is set equal to direct
#' @param Q A scalar to denote the number of quadrature points used for GH numerical integration. Default is 10.
#' @param M Number of imputations. Default is 5.
#' @param marg.exp.formula Marginal exposure model of the form: X ~ Z where Z are the inexpensive covariates in \code{data}
#' @param m A scalar giving the number of iterations within each imputation. The default is 5.
#' @param verbise An indicator of whether the steps of the maximization algorithm should be printed. Default is set to FALSE
#'
#'
#' @importFrom binaryMM mm
#' @importFrom MASS mvrnorm
#' @importFrom mitools MIcombine
#' @importFrom dplyr sample_n
#' @keywords internal
#' @return The function returns the estimated coefficients and standard error from the direct imputation approach

DirectImputation <- function(mean.formula, lv.formula, t.formula, id, data, sampled, X,
                             samp.probs, Q=10, M=5, marg.exp.formula, m = 20, verbose = FALSE){

  # Step 0: prepare for multiple imputation
  DAT_Xe1 = DAT_Xe0 = data
  DAT_Xe1[,X] = 1
  DAT_Xe0[,X]= 0

  if(Q <= 1){Q <- 2}
  if(is.null(lv.formula)){Q <- 1}
  AdaptiveQuad <- get.GH(Q)
  W.ada <- AdaptiveQuad$w
  Z.ada <- AdaptiveQuad$z

  mi     = c(table(id))

  # Step 1: on sampled subjects only fit mm
  mstart = binaryMM::mm(mean.formula = mean.formula,
                        t.formula = t.formula, lv.formula = lv.formula,
                        id = id, data = data, q = Q, verbose = verbose)
  thetas    = c(mstart$beta, mstart$alpha)
  cov.theta = mstart$mod.cov

  # store the results
  coefficients = variance = list()

  for(j in 1:M){

    for(i in 1:m){
      # Step 2: sample the coefficients to compute the offsets
      new.theta = MASS::mvrnorm(n = 1, mu = thetas, Sigma = cov.theta)

      # Step 3: compute the offset using the new theta
      if(is.null(t.formula)){new.theta = append(new.theta, 0, after = (length(mstart$beta)))}
      if(is.null(lv.formula)){new.theta = c(new.theta, 0)}
      tmp_S1_Xe1 = prep_data(mean.formula = mean.formula, lv.formula = lv.formula,
                             t.formula = t.formula,
                             id=id, data=DAT_Xe1, inits=new.theta,
                             samp.probs=c(1,1,1), samp.probi=NULL, offset = NULL)
      sData_Xe1  = tmp_S1_Xe1$subjectData
      npar       = tmp_S1_Xe1$paramlengths
      sData_Xe0  = prep_data(mean.formula = mean.formula, lv.formula = lv.formula,
                             t.formula = t.formula,
                             id=id, data=DAT_Xe0, inits=new.theta,
                             samp.probs=c(1,1,1), samp.probi=NULL, offset = NULL)$subjectData

      LLSC_1 = LogLScoreCalc(params = new.theta, subjectData = sData_Xe1,
                             Q = Q, W = W.ada, Z = Z.ada,
                             ParamLengths = npar,
                             CondLike = FALSE, EmpiricalCheeseCalc = FALSE)
      LLSC_0 = LogLScoreCalc(params = new.theta, subjectData = sData_Xe0,
                             Q = Q, W = W.ada, Z = Z.ada,
                             ParamLengths = npar,
                             CondLike = FALSE, EmpiricalCheeseCalc = FALSE)
      offset.all  = rep(attr(LLSC_1,"LogLikeSubj") - attr(LLSC_0,"LogLikeSubj"), mi)
      data$offset = offset.all
      dat.wide    = data %>% group_by(id) %>% sample_n(1)

      # Step 4: fit an imputed glm
      fit.exp = glm(marg.exp.formula, offset = offset, data=dat.wide, family = "binomial")
      # Extract alpha estimates and sample from the sampling distribution
      alpha      = coef(fit.exp)
      cov_alpha  = vcov(fit.exp)
      new.alpha  = MASS::mvrnorm(n = 1, mu=alpha, Sigma=cov_alpha)

      # Step 5: for the unsampled subject impute X
      # Create a temporary dataframe for non-sampled subjects in order to create the
      # model.matrix that is used to calculate marginal exposure odds
      DAT_Xe1$offset       = offset.all
      tmp                  = DAT_Xe1 %>% group_by(id) %>% sample_n(1)
      mm_tmp               = model.matrix(marg.exp.formula, tmp)
      pred.prob            = expit(mm_tmp%*%new.alpha + tmp$offset)
      data.imp             = dat.wide
      imp.grp              = rbinom(nrow(data.imp), 1, pred.prob)
      X.new                = ifelse(is.na(data.imp$X), imp.grp, data.imp$X)
      data.imp$X.new       = X.new

      data.imp             = subset(data.imp, select = c(id, X.new))
      dat.final            = merge(data.imp, data, by = "id")
      dat.final$X          = ifelse(is.na(dat.final$X), dat.final$X.new, dat.final$X)
      dat.final            = as.data.frame(dat.final)

      # Step 6: fit a MM to get a new estimate for theta
      mod.imp = binaryMM::mm(mean.formula = mean.formula,
                             t.formula = t.formula, lv.formula = lv.formula,
                             id = id, data = dat.final, q = Q, verbose = verbose)
      thetas    = c(mod.imp$beta, mod.imp$alpha)
      cov.theta = mod.imp$mod.cov
    }

    coefficients[[j]] = c(mod.imp$beta, mod.imp$alpha)
    variance[[j]]     = mod.imp$mod.cov
  }

  out <- mitools::MIcombine(coefficients, variance)
  out
}
