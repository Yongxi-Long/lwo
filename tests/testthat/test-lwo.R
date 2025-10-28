test_that("lwo() basic functionality works on simple simulated data", {
  skip_if_not_installed("geepack")
  set.seed(123)

  # --- Simulate minimal longitudinal data ---
  n_subjects <- 10
  n_visits <- 3
  df <- expand.grid(
    id = 1:n_subjects,
    visit = 1:n_visits
  )
  df$age <- rnorm(nrow(df), 50, 10)
  df$treat <- rbinom(nrow(df), 1, 0.5)
  # outcome: ordinal-like numeric
  df$outcome <- rbinom(nrow(df), 2, plogis(-1 + 0.3 * df$treat + 0.02 * df$age))

  # --- Fit model ---
  mod <- lwo(
    outcome ~ treat + age,
    data = df,
    id.var = "id",
    visit.var = "visit",
    corstr = "independence",
    larger = TRUE
  )

  # --- Basic structure checks ---
  expect_s3_class(mod, "lwo")
  expect_type(mod$coefficients, "double")
  expect_true(is.matrix(mod$varcov))
  expect_equal(ncol(mod$model.matrix), length(mod$coefficients))
  expect_true(all(diag(mod$varcov) > 0))

  # --- Check that fitted values look plausible ---
  expect_true(all(mod$fitted.values >= 0 & mod$fitted.values <= 1))
  expect_length(mod$fitted.values, nrow(mod$model.matrix))

  # --- Check correlation structure name ---
  expect_true(mod$corstr %in% c("independence", "exchangeable", "ar1", "unstructured"))

  # --- Check that variance-covariance matrix is symmetric ---
  expect_equal(mod$varcov, t(mod$varcov))
})

test_that("lwo() throws informative errors for bad inputs", {
  skip_if_not_installed("geepack")

  df <- data.frame(
    id = rep(1:4, each = 2),
    visit = rep(1:2, times = 4),
    y = rbinom(8, 1, 0.5),
    x = rnorm(8)
  )

  # Missing id.var
  expect_error(
    lwo(y ~ x, data = df, visit.var = "visit"),
    "id.var"
  )

  # Missing visit.var
  expect_error(
    lwo(y ~ x, data = df, id.var = "id"),
    "visit.var"
  )

  # Unsupported correlation structure
  expect_error(
    lwo(y ~ x, data = df, id.var = "id", visit.var = "visit", corstr = "nonsense")
  )

  # Missing zcor when corstr = "fixed"
  expect_error(
    lwo(y ~ x, data = df, id.var = "id", visit.var = "visit", corstr = "fixed"),
    "zcor"
  )
})

test_that("lwo() handles different correlation structures and std.err options", {
  skip_if_not_installed("geepack")
  set.seed(123)

  df <- expand.grid(id = 1:6, visit = 1:3)
  df$x <- rnorm(nrow(df))
  df$y <- rbinom(nrow(df), 1, plogis(0.3 * df$x))

  # independence
  mod1 <- lwo(y ~ x, data = df, id.var = "id", visit.var = "visit", corstr = "independence")
  expect_s3_class(mod1, "lwo")

  # exchangeable
  mod2 <- lwo(y ~ x, data = df, id.var = "id", visit.var = "visit", corstr = "exchangeable")
  expect_s3_class(mod2, "lwo")

  # ar1
  mod3 <- lwo(y ~ x, data = df, id.var = "id", visit.var = "visit", corstr = "ar1")
  expect_s3_class(mod3, "lwo")

  # std.err = "san.se.modified"
  mod4 <- lwo(y ~ x, data = df, id.var = "id", visit.var = "visit",
              corstr = "independence", std.err = "san.se.modified")
  expect_s3_class(mod4, "lwo")
})

