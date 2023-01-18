#' Simulated data set
#'
#' A simulated dataset.
#' Data were created from a marginalised latent and transition model with a transition term only.
#'
#'
#' @docType data
#'
#' @usage data(exampledat)
#'
#' @format An object of class \code{data.frame} with 2267 observations and 5 variables:
#' \describe{
#' \item{\code{id}}{Subject identifier}
#' \item{\code{X}}{A time-fixed binary covariate}
#' \item{\code{Z}}{A time-fixed binary covariate}
#' \item{\code{time}}{A time-varying covariate indicating study visit. Each subject is observed three to six times}
#' \item{\code{Y}}{A longitudinal binary outcome}
#' }
#'
#' @keywords datasets
#'
#'
#' @examples
#' data(exampledat)
#' str(exampledat)
"exampledat"
