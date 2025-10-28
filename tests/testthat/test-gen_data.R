test_that("gen_data basic functionality works", {
  skip_if_not_installed("simstudy")
  skip_if_not_installed("reshape2")
  skip_if_not_installed("dplyr")

  set.seed(123)

  # Define inputs
  baseprobs <- c(0.3, 0.4, 0.3)
  covs_effects <- c(age = 0.1, pre_diarrhea = -0.3)
  visits <- paste0("visit", 0:3)
  corMat <- diag(length(visits))
  corMat[lower.tri(corMat)] <- 0.5
  corMat[upper.tri(corMat)] <- 0.5

  # Run generator
  dat <- gen_data(
    N = 50,
    baseprobs = baseprobs,
    covs_effects = covs_effects,
    visits = visits,
    corMatrix = corMat,
    trt_ratio = 1
  )

  # Check output structure
  expect_type(dat, "list")
  expect_named(dat, c("wide_format", "long_format"))
  expect_true(is.data.frame(dat$wide_format))
  expect_true(is.data.frame(dat$long_format))

  # Check dimensions
  expect_equal(nrow(dat$wide_format), 50)
  expect_true(all(visits %in% names(dat$wide_format)))

  # Each subject has all visits in long format
  expect_equal(
    nrow(dat$long_format),
    50 * length(visits)
  )

  # Scores should be numeric and in range
  expect_true(is.numeric(dat$long_format$score))
  expect_true(all(dat$long_format$score >= 1))
})

test_that("gen_data input validation works", {
  baseprobs <- c(0.3, 0.4, 0.3)
  covs_effects <- c(age = 0.1)
  visits <- c("v0", "v1")
  corMat <- diag(2)

  # Invalid N
  expect_error(gen_data(N = -10, baseprobs, covs_effects, visits = visits, corMatrix = corMat),
               "must be a positive integer")

  # baseprobs not summing to 1
  expect_error(gen_data(N = 10, baseprobs = c(0.5, 0.5, 0.2), covs_effects, visits = visits, corMatrix = corMat),
               "summing to 1")

  # Missing covariate names
  expect_error(gen_data(N = 10, baseprobs, covs_effects = c(0.2), visits = visits, corMatrix = corMat),
               "named numeric vector")

  # Wrong corMatrix size
  expect_error(gen_data(N = 10, baseprobs, covs_effects, visits = c("v0", "v1", "v2"), corMatrix = diag(2)),
               "dimensions matching the number of visits")

  # Invalid trt_ratio
  expect_error(gen_data(N=10,baseprobs,trt_ratio = -1, covs_effects, visits = c("v0", "v1", "v2"), corMatrix = diag(3)),
               "positive numeric value")

  # no multiple visits
  expect_error(gen_data(N=10,baseprobs, covs_effects, visits = c("v0"), corMatrix = diag(3)),
               "at least two elements")
})

test_that("gen_data handles time effects and time-by-treatment interaction correctly", {
  baseprobs <- c(0.3, 0.4, 0.3)
  covs_effects <- c(age = 0.2)
  visits <- paste0("v", 0:3)
  corMat <- diag(4)

  # Mismatched time_effects length
  expect_error(
    gen_data(
      N = 20,
      baseprobs,
      covs_effects,
      time_effects = c(0.1, 0.2),
      visits = visits,
      corMatrix = corMat
    ),
    "follow-up visits"
  )

  # Mismatched time_trt_interaction length
  expect_error(
    gen_data(
      N = 20,
      baseprobs,
      covs_effects,
      time_trt_effects = c(0.1, 0.2),
      visits = visits,
      corMatrix = corMat
    ),
    "follow-up visits"
  )
})

test_that("gen_data handles optional arguments extractTime and dropZ", {
  baseprobs <- c(0.3, 0.4, 0.3)
  covs_effects <- c(age = 0.2)
  visits <- c("visit0", "visit1", "visit2")
  corMat <- diag(3)

  dat <- gen_data(
    N = 10,
    baseprobs = baseprobs,
    covs_effects = covs_effects,
    visits = visits,
    corMatrix = corMat,
    extractTime = TRUE,
    dropZ = FALSE
  )

  expect_true("time" %in% names(dat$long_format))
  expect_true("z" %in% names(dat$long_format))
  expect_true(is.numeric(dat$long_format$time))
})
