#' Title
#' None-Some-All (NSA) Design
#'
#' Function to perform the NSA design. The function requires specification of the column name with the outcome, the column name with the ID,
#' the number of subjects the user wants to sample from the none, some and all strata and the dataframe. For each stratum, if the
#' user-specified number of subjects of being sampled is greater than the actual number of subjects, the function
#' will return every subjects from that stratum with a warning.
#'
#' @importFrom stats rbinom
#' @param Y Column name with the \code{Y} variable. This is a required argument
#' @param id Column name with the \code{id} variable. This is a required argument
#' @param n.sample Vector for three values representing the number of people to be sampled from the none, some and all stratum respectively. This is a required argument
#' @param data Dataframe containing \code{Y} and \code{id}. This is a required argument
#'
#' @return A list with the following components:
#' \describe{
#' \item{sampled.id}{List of IDs number sampled in the none, some, all strata}
#' \item{sample.probs}{Vector of sampling probabilities from the none, some, all strata to be used in the analysis}
#' \item{n.stratum}{Vector specifying the number of subjects in the none, some, all strata}
#' \item{membership}{Vector specifying in which strata a subject belongs to (None, Some, All)}
#' }
#' @export
#'
#' @examples
#' \donttest{
#' data(exampledat)
#' nsa_design(Y = "Y", id = "id", n.sample = c(85, 100, 15), data = exampledat)
#' }
#'
nsa_design <- function(Y, id, n.sample, data = NULL){

  # errors/warning messages
  if(is.null(data)){stop("Dataframe needs to be specified")}
  if(!is.data.frame(data)){
    data <- data.frame(data)
    warning("Data was tranformed in a dataframe object")}
  if(is.null(Y)){stop("Y needs to be specified")}
  if(!is.null(Y) & !is.character(Y)){stop("Y needs to be a character")}
  if(is.null(id)){stop("id needs to be specified")}
  if(!is.null(id) & !is.character(id)){stop("id needs to be a character")}
  if(!(Y %in% names(data))){stop(paste("There is variable in the dataframe named", Y))}
  if(!(id %in% names(data))){stop(paste("There is variable in the dataframe named", id))}
  if(!is.null(n.sample) & !is.vector(n.sample)){stop("n.sample needs to be a vector")}

  if(is.null(n.sample) | length(n.sample) != 3){
    warning("Number of people to be sampled was not provided or vector provided was not of length 3. Returning all the ids")
    out        = as.character(unique(data[,id]))
    prob.none  = prob.some <- prob.all  <- 1
    } else {
      id           = data[,id]
      Y            = data[, Y]
      if(!is.numeric(Y)){stop("The outcome needs to be numeric")}
      if(is.numeric(Y) & any(Y != 1 & Y != 0)){stop("The outcome can only take 0 or 1 values")}

      # set up the design
      n.strata.all = tapply(Y, INDEX = id, FUN = mean)
      none.id      = names(n.strata.all[n.strata.all == 0])
      some.id      = names(n.strata.all[n.strata.all != 1 & n.strata.all != 0])
      all.id       = names(n.strata.all[n.strata.all == 1])
      n.none.tot   = length(none.id)
      n.some.tot   = length(some.id)
      n.all.tot    = length(all.id)

      # save which subject is in each stratum
      stratum.member = ifelse(unique(id) %in% unique(none.id), "None",
                              ifelse(unique(id) %in% unique(some.id), "Some", "All"))

      # sample
      if(n.sample[1] > n.none.tot){warning("Number to be sampled greater than number in the none stratum. Sample
                                       with probability 1 from the none stratum")}
      if(n.sample[2] > n.some.tot){warning("Number to be sampled greater than number in the some stratum. Sample
                                       with probability 1 from the some stratum")}
      if(n.sample[3] > n.all.tot){warning("Number to be sampled greater than number in the all stratum. Sample
                                       with probability 1 from the all stratum")}

      prob.none <- min(n.sample[1], n.none.tot)/n.none.tot
      prob.some <- min(n.sample[2], n.some.tot)/n.some.tot
      prob.all  <- min(n.sample[3], n.all.tot)/n.all.tot

      none.samp <- rbinom(n = n.none.tot, size = 1, prob = prob.none)
      some.samp <- rbinom(n = n.some.tot, size = 1, prob = prob.some)
      all.samp  <- rbinom(n = n.all.tot,  size = 1, prob = prob.all)
      out       <- c(none.id[none.samp == 1], some.id[some.samp == 1], all.id[all.samp == 1])

  }


  return(list(sampled.id   = out,
              sample.probs = c(prob.none, prob.some, prob.all),
              n.stratum    = c(n.none.tot, n.some.tot, n.all.tot),
              membership   = stratum.member))

}
