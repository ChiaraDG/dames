#' Summarise the results
#'
#' @param object Object of class TwoPhaseMM
#' @param ... Additional input
#'
#' @return The function results result summaries of a TwoPhaseMM object
#' @export
#'

summary.TwoPhaseMM <-
  function(object, ...) {
    # get the object
    oo      = object
    # divide in alpha and beta
    m.i     = seq(length(oo$beta))
    a.i     = seq(length(oo$alpha))+length(oo$beta)
    # mean model paramters
    mean.oo        = with(oo, data.frame(estimate = beta, mod.se = sqrt(diag(mod.cov)[m.i])))
    mean.oo$X2     = with(mean.oo, (estimate/mod.se)^2)
    mean.oo$pX2    = with(mean.oo, pchisq(X2, df=1,lower.tail=FALSE))
    names(mean.oo) = c('Estimate','Model SE', 'Chi Square','Pr(>Chi)')
    # depedence model parameters
    assoc.oo        = with(oo, data.frame(estimate = alpha, mod.se = sqrt(diag(mod.cov)[a.i])))
    assoc.oo$X2     = with(assoc.oo, (estimate/mod.se)^2)
    assoc.oo$pX2    = with(assoc.oo, pchisq(X2, df=1,lower.tail=FALSE))
    names(assoc.oo) = c('Estimate','Model SE','Chi Square','Pr(>Chi)')

    if(object$control["useRobCov"] == 1) {
      mean.oo        = with(oo, data.frame(estimate = beta, rob.se = sqrt(diag(rob.cov)[m.i])))
      mean.oo$X2     = with(mean.oo, (estimate/rob.se)^2)
      mean.oo$pX2    = with(mean.oo, pchisq(X2, df=1,lower.tail=FALSE))
      names(mean.oo) = c('Estimate','Robust SE', 'Chi Square','Pr(>Chi)')

      assoc.oo        = with(oo, data.frame(estimate = alpha, rob.se = sqrt(diag(rob.cov)[a.i])))
      assoc.oo$X2     = with(assoc.oo, (estimate/rob.se)^2)
      assoc.oo$pX2    = with(assoc.oo, pchisq(X2, df=1,lower.tail=FALSE))
      names(assoc.oo) = c('Estimate','Robust SE','Chi Square','Pr(>Chi)')
      warning('When performing a weighted likelihood analysis (by specifying the weights argument), robust standard errors are reported. Model based standard errors will not be correct and should not be used.')
    }

    if(object$control["se_comp"] == 1) {
      mean.oo        = with(oo, data.frame(estimate = beta, rob.se = sqrt(diag(rob.cov)[m.i])))
      mean.oo$X2     = with(mean.oo, (estimate/rob.se)^2)
      mean.oo$pX2    = with(mean.oo, pchisq(X2, df=1,lower.tail=FALSE))
      names(mean.oo) = c('Estimate','SE', 'Chi Square','Pr(>Chi)')

      assoc.oo        = with(oo, data.frame(estimate = alpha, rob.se = sqrt(diag(rob.cov)[a.i])))
      assoc.oo$X2     = with(assoc.oo, (estimate/rob.se)^2)
      assoc.oo$pX2    = with(assoc.oo, pchisq(X2, df=1,lower.tail=FALSE))
      names(assoc.oo) = c('Estimate','SE','Chi Square','Pr(>Chi)')
      warning('To estimate the standard errors set no_SE = TRUE')
    }

    out = list(class = class(object), call = deparse(object$call), control=object$control,
               info=object$info_stats, mean.table = mean.oo, assoc.table=assoc.oo)
    class(out) = 'summary.TwoPhaseMM'
    out
  }
