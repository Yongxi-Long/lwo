test_that("calculate_win_odds returns numeric vector of correct length", {
  set.seed(123)

  baseprobs <- c(0.3, 0.4, 0.3)
  visits <- c("v0","v1","v2")
  corMat <- gen_corMatrix(n_visits = length(visits), rho = 0.4, corstr = "ar1")

  out <- calculate_win_odds(
    N_approx = 2000,
    baseprobs = baseprobs,
    covs_effects = c(age = 0.05),
    time_effects = c(0, 0),
    trt_ratio = 1,
    time_trt_effects = c(0, 0),
    visits = visits,
    corMatrix = corMat
  )

  expect_type(out, "double")
  expect_length(out, length(visits))
  expect_named(out, visits)
  expect_true(all(out > 0))
})


test_that("baseline visit win odds = 1 without treatment effect", {
  set.seed(123)

  baseprobs <- c(0.2,0.5,0.3)
  visits <- c("week0","week4")
  corMat <- gen_corMatrix(n_visits = length(visits), rho = 0.2, corstr = "ar1")

  out <- calculate_win_odds(
    N_approx = 1500,
    baseprobs = baseprobs,
    covs_effects = c(age = 0.01),
    time_effects = c(0),
    trt_ratio = 1,
    time_trt_effects = c(0),
    visits = visits,
    corMatrix = corMat
  )
  expect_equal(as.numeric(out[1]), 1, tolerance = 0.001)
})
