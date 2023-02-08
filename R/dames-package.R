#' The "dames" package
#' @description The package allows user to design different two-phase outcome dependent sampling studies for a binary longitudinal outcome.
#' Together with the design the package provides different analysis methods to estimates the parameters of interest. Note that in each method,
#' the analysis is carried out using a marginalized transition and latent variable model
#'
#'
#'
#' @docType package
#' @name dames
#' @useDynLib dames, .registration = TRUE
#' @importFrom stats rbinom
#' @importFrom  binaryMM mm
#' @importFrom  dplyr filter
#' @importFrom  dplyr ungroup
#' @importFrom  stats model.matrix
#' @importFrom  stats model.frame
#' @importFrom  stats var
#' @importFrom  magrittr %>%
#' @importFrom  stats nlm
#' @import fastGHQuad
#' @importFrom  stats glm
#' @importFrom  stats model.response
#' @importFrom  stats na.omit
#' @references
#' Add the list of reference here
NULL
#> NULL
