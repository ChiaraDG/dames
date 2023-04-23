#' Simulated data set
#'
#' A simulated dataset.
#' Synthetic data were created to resemble data from a subset of participants in the Lung Health Study
#'
#'
#' @docType data
#'
#' @usage data(lhsexample)
#'
#' @format An object of class \code{data.frame} with 5125 observations and 8 variables:
#' \describe{
#' \item{\code{id}}{Subject identifier}
#' \item{\code{snp}}{A time-fixed binary covariate, indicating presence/absence of T-allele at snp}
#' \item{\code{age}}{Age at baseline (years)}
#' \item{\code{cigs}}{Number of cigarettes smoked per day in the year before the enrollment in the study}
#' \item{\code{bmi}}{BMI at baseline}
#' \item{\code{male}}{A binary variable indicating sex (1 = male, 0 = female) at baseline}
#' \item{\code{Y}}{A longitudinal binary outcome}
#' }
#'
#' @keywords datasets
#'
#'
#' @examples
#' data(lhsexample)
#' str(lhsexample)
"lhsexample"
