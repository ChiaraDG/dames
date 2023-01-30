test_that("acml_mm returns a TwoPhaseMM object", {

  data("exampledat")
  design       <- nsa_design(Y = "Y", id = "id", n.sample = c(85, 100, 15),
                             data = exampledat)
  # Assume X is not collected for subjects not sampled
  exampledat$X <- ifelse(exampledat$id %in% design$sampled.id, exampledat$X, NA)
  mod.acml    <- acml_mm(mean.formula = Y ~ X + time + Z + X:time,
                         t.formula = ~1,
                         id = id, samp.probs = design$sample.probs,
                         data = exampledat)

  expect_equal(class(mod.acml), "TwoPhaseMM")

})

test_that("acml_mm returns an error for misspecified models", {

  data("exampledat")
  design       <- nsa_design(Y = "Y", id = "id", n.sample = c(85, 100, 15),
                             data = exampledat)
  # Assume X is not collected for subjects not sampled
  exampledat$X <- ifelse(exampledat$id %in% design$sampled.id, exampledat$X, NA)
  expect_error(acml_mm(mean.formula = Y ~ X + time + Z + X:time,
                         id = id, samp.probs = design$sample.probs,
                         data = exampledat))

})

test_that("wee_mm returns an error for misspecified models", {

  data("exampledat")
  design       <- nsa_design(Y = "Y", id = "id", n.sample = c(85, 100, 15),
                             data = exampledat)
  # Assume X is not collected for subjects not sampled
  exampledat$X <- ifelse(exampledat$id %in% design$sampled.id, exampledat$X, NA)

  # vector of sampling probabilities per subject
  samp.probs <- ifelse(design$membership == "None", design$sample.probs[1],
                       ifelse(design$membership == "Some", design$sample.probs[2],
                              design$sample.probs[3]))
  # number of times each subject is observed
  time_i     <- exampledat %>% group_by(id) %>% summarise(mi = n()) %>% pull()
  # repeat the sampling probability for each subject based on the number
  # of times a subject is observed
  exampledat$samp.probs <- rep(samp.probs, time_i)
  samp.probs            <- exampledat[!is.na(exampledat$X),"samp.probs"]

  expect_error(summary(wee_mm(mean.formula = Y ~ X + time + Z + X:time,
                                id = id, weights = 1/samp.probs,
                                data = exampledat)))

})



test_that("wee_mm returns a warning", {

  data("exampledat")
  design       <- nsa_design(Y = "Y", id = "id", n.sample = c(85, 100, 15),
                             data = exampledat)
  # Assume X is not collected for subjects not sampled
  exampledat$X <- ifelse(exampledat$id %in% design$sampled.id, exampledat$X, NA)

  # vector of sampling probabilities per subject
  samp.probs <- ifelse(design$membership == "None", design$sample.probs[1],
                       ifelse(design$membership == "Some", design$sample.probs[2],
                              design$sample.probs[3]))
  # number of times each subject is observed
  time_i     <- exampledat %>% group_by(id) %>% summarise(mi = n()) %>% pull()
  # repeat the sampling probability for each subject based on the number
  # of times a subject is observed
  exampledat$samp.probs <- rep(samp.probs, time_i)
  samp.probs            <- exampledat[!is.na(exampledat$X),"samp.probs"]

  expect_warning(summary(wee_mm(mean.formula = Y ~ X + time + Z + X:time,
                        t.formula = ~1,
                        id = id, weights = 1/samp.probs,
                        data = exampledat)))

})
