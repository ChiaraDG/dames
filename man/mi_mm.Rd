% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mi_mm.R
\name{mi_mm}
\alias{mi_mm}
\title{Multiple Imputation for Two-Phase ODS Designs and Marginalized Transition and Latent Variable Models}
\usage{
mi_mm(
  mean.formula,
  lv.formula = NULL,
  t.formula = NULL,
  id,
  data,
  sampled = NULL,
  X,
  samp.probs = NULL,
  maxiter = 100,
  Q = 10,
  M = 5,
  m = 5,
  marg.exp.formula,
  method = "indirect",
  stepmax = 1,
  verbose = FALSE
)
}
\arguments{
\item{mean.formula}{Formula for the marginal mean model in which a binary variable is regressed on a set of covariates. This is a required argument}

\item{lv.formula}{Formula for the latent part in the dependence model. This is optional and of the form ~ covariates. Default is set to NULL.}

\item{t.formula}{Formula for the transition part in the dependence model. This is optional and of the form ~ covariates. Default is set to NULL.}

\item{id}{Column name for the id variable as recorded in the dataset. This is a required argument}

\item{data}{Full cohort data. This dataset should not have missing values for variables other than X. It should
contain columns named \code{id}, \code{X}, and \code{sampled}}

\item{sampled}{Column name for the sampling indicator variable as recorded in the dataset. This is a required argument}

\item{X}{Column name for the expensive exposure variable as recorded in the dataset. This is a required argument}

\item{samp.probs}{Vector of 3 numbers indicating the probability that a subject is sampled from the none
some and all stratum respectively. This is a required argument when \code{method} is set equal to indirect. It is not
used when \code{method} is set equal to direct, Default is set to NULL}

\item{maxiter}{A scalar indicating the maximum number of iteration in the maximization algorithm. Default is 100}

\item{Q}{A scalar to denote the number of quadrature points used for GH numerical integration. Default is 10.}

\item{M}{Number of imputations. Default is 5.}

\item{m}{A scalar giving the number of iterations within each imputation. The default is 5. The argument is only when \code{method} is set to direct.}

\item{marg.exp.formula}{Marginal exposure model of the form: X ~ Z where Z are the inexpensive covariates in \code{data}. This is a required argument.}

\item{method}{Character indication whether the indirect MI or the direct MI should be performed. Default is indirect. Note that indirect MI can only be used under an NSA design}

\item{stepmax}{Maximum step size allowed for the minimization algorithm in nlm. Default is set to 1.}

\item{verbose}{An indicator of whether the steps of the maximization algorithm should be printed}
}
\value{
This function returns marginal mean (beta) and dependence parameters (alpha) along with the associated model and empirical covariance matrices.
}
\description{
Estimates regression parameters from a marginalized transition and latent variable model for scenarios where
data are collected using a two-phase ODS design and the expensive exposure of interest is binary.
}
\examples{
\donttest{
data(exampledat)
# perform the NSA design
design       <- nsa_design(Y = "Y", id = "id", n.sample = c(85, 100, 15), data = exampledat)
# Assume X is not collected for subjects not sampled
exampledat$X <- ifelse(exampledat$id \%in\% design$sampled.id, exampledat$X, NA)
# estimate the parameters using Direct MI
mod.mi2 <- mi_mm(mean.formula = Y ~ X + time + Z + X:time, t.formula = ~1,
                 lv.formula = NULL, id = "id", X = "X",
                 data = exampledat, marg.exp.formula = X ~ Z, M = 5,
                method = "direct")
# print the results
summary(mod.mi2)
}

}
