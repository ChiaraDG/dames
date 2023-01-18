#' Residual-Based Design
#'
#' Function to perform residual-based designs. The function requires specification of the column name with the outcome, the column name with the ID,
#' the marginal mean model, at least one between the transition component or the latent component of the dependence model,
#' the number of subjects the user wants to sample, and the dataframe.
#' @importFrom  binaryMM mm
#' @importFrom  dplyr filter
#' @importFrom  dplyr ungroup
#' @importFrom  stats model.matrix
#' @importFrom  stats model.frame
#' @importFrom  stats var
#' @importFrom  stats nlm
#' @importFrom  magrittr %>%
#' @param Y Column name for outcome variable as recorded in the dataset. This is a required argument
#' @param id Column name for the \code{id} variable as recorded in the dataset. This is a required argument
#' @param mean.formula Formula for the marginal mean model in which a binary variable is regressed on a set of covariates. This is a required argument
#' @param t.formula Formula for the transition part in the dependence model. This is optional and of the form ~ covariates
#' @param lv.formula Formula for the transition part in the dependence model. This is optional and of the form ~ covariates
#' @param n.sample Number of people to be sampled with the residual-based design. This is a required argument
#' @param statistic Name of the statistic to be used for the residual-based design. Current statistics supported are \code{"mean"}, \code{"abs.mean"}, \code{"Variance"}. Default is \code{"mean"}
#' @param data Dataframe containing \code{Y} and \code{id}. This is a required argument
#'
#' @return A list with the following components:
#' \describe{
#' \item{sampled.id}{List of IDs number sampled}
#' }
#' @export
#'
#' @examples
#' \donttest{
#' data(exampledat)
#' res_design(Y = "Y", id = "id", mean.formula = Y ~ Z + time, t.formula = ~1, lv.formula = NULL,
#' n.sample = 200, statistic = "abs.mean", data = exampledat)
#' }
res_design <- function(Y, id,
                       mean.formula, t.formula = NULL,
                       lv.formula = NULL, n.sample, statistic = "mean", data = NULL){

  if(is.null(data)){stop("Name of the dataframe need to be specified")}
  if(!is.data.frame(data)){
    data <- data.frame(data)
    warning("Data was tranformed in a dataframe object")}
  if(is.null(mean.formula)){stop("Mean model needs to be specified")}
  if(is.null(t.formula) & is.null(lv.formula)){stop("Either transition or latent model needs to be specified")}
  if(is.null(n.sample)){stop("n.sample needs to be specified")}
  if(is.null(Y)){stop("Y needs to be specified")}
  if(!is.null(Y) & !is.character(Y)){stop("Y needs to be a character")}
  if(is.null(id)){stop("id needs to be specified")}
  if(!is.null(id) & !is.character(id)){stop("id needs to be a character")}
  if(!(Y %in% names(data))){stop(paste("There is variable in the dataframe named", Y))}
  if(!(id %in% names(data))){stop(paste("There is variable in the dataframe named", id))}

  names(data)[names(data) == id] <- "id"
  names(data)[names(data) == Y]  <- "Y"
  N.tot                          <- length(unique(data$id))
  if(n.sample > N.tot){
    warning("Number of people to be sampled greater than total number of people.
            Return all ids")
    out <- list(sampled.id = unique(data$id))
    }

  mod.phase1    <- mm(mean.formula, t.formula = t.formula, lv.formula = lv.formula, id = id, data = data)
  x.design      <- model.matrix(mean.formula, model.frame(mean.formula, data))
  mu.m.ij       <- expit(x.design%*%mod.phase1$beta)
  res           <- data[,"Y"] - mu.m.ij
  data          <- data.frame(data, res = res)

  if(statistic == "mean"){
    dat.res     <- data
    dat.res     <- dat.res %>% dplyr::group_by(id) %>%
      dplyr::mutate(Ysum = sum(Y), res.mean = mean(res)) %>% dplyr::ungroup() %>%
      data.frame()
    dat.res.unique  <- dat.res[, c("id", "res.mean", "Ysum")] %>%
      dplyr::group_by(id) %>% dplyr::slice(1) %>% dplyr::ungroup() %>%
      data.frame()
    dat.res.unique  <- dat.res.unique[order(-dat.res.unique$res.mean),]
    tokeep          <- dat.res.unique[c(1:ceiling(n.sample/2),
                                        (N.tot - ceiling(n.sample/2) + 1):N.tot),"id"]
    if(n.sample%%2 == 1){"n.sample is an odd number. Note that two additional ids were sampled"}
    out <- list(sampled.id = as.character(tokeep))
  }

  if(statistic == "abs.mean"){
    dat.res     <- data
    dat.res     <- dat.res %>% dplyr::group_by(id) %>%
      dplyr::mutate(Ysum = sum(Y), res.mean = abs(mean(res))) %>%
      dplyr::ungroup() %>% data.frame()
    dat.res.unique  <- dat.res[, c("id", "res.mean", "Ysum")] %>%
      dplyr::group_by(id) %>% dplyr::slice(1) %>% ungroup() %>% data.frame()
    dat.res.unique  <- dat.res.unique[order(-dat.res.unique$res.mean),]
    tokeep          <- dat.res.unique[1:n.sample, "id"]
    out <- list(sampled.id = as.character(tokeep))

  }
  if(statistic == "variance"){
    dat.res     <- data
    dat.res     <- dat.res %>% dplyr::group_by(id) %>%
      dplyr::mutate(Ysum = sum(Y), res.var = var(res)) %>%
      dplyr::ungroup() %>% data.frame()
    dat.res.unique  <- dat.res[, c("id", "res.var", "Ysum")] %>%
      dplyr::group_by(id) %>% dplyr::slice(1) %>% dplyr::ungroup() %>% data.frame()
    dat.res.unique  <- dat.res.unique[order(-dat.res.unique$res.var),]
    tokeep          <- dat.res.unique[1:n.sample,"id"]
    out <- list(sampled.id = as.character(tokeep))
  }

  out
}
