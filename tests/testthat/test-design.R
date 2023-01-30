test_that("nsa_design and res_design results in a list of appropriate length", {

  data("exampledat")
  nsa <- nsa_design(Y = "Y", id = "id",
                    n.sample = c(85, 100, 15), data = exampledat)

  res <- res_design(Y = "Y", id = "id",
                    mean.formula = Y ~ Z + time, t.formula = ~1,
                    lv.formula = NULL,
                    n.sample = 200, statistic = "abs.mean", data = exampledat)

  expect_equal(class(nsa), "list")
  expect_equal(class(res), "list")

  expect_equal(length(nsa), 4)
  expect_equal(length(res), 1)

})
