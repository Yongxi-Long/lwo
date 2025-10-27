#' Calculate Marginal True Win Odds by Visit
#'
#' @description
#' Use simulation to compute the true (population-level) win odds at each visit
#' marginalized with respect to the covariate distribution
#'
#' @param N_approx Integer. Approximate number of subjects to simulate for estimating covariate distribution (default = 1e4).
#' @param baseprobs Numeric vector of baseline category probabilities when all covariates are set to zero/reference level (e.g., \code{c(0.2, 0.5, 0.3)}). Must sum to 1.
#' @param covs_effects Named numeric vector specifying the baseline covariate effects on the ordinal outcome
#'   (e.g., \code{c(age = 0.1, pre_diarrhea = -0.3)}). Names must match variables defined in \code{covariate_def}.
#' @param time_effects Optional numeric vector giving log-odds shifts across follow-up visits
#'   (length = number of visits - 1). Defaults to 0 if not specified.
#' @param trt_ratio Numeric. Treatment-to-control randomization ratio (default = 1, i.e., 1:1 randomization).
#' @param time_trt_effects Optional numeric vector giving additional time-varying treatment effects
#'   (length = number of visits - 1). Defaults to 0 if not specified.
#' @param visits Character vector of visit/measurement names (e.g., \code{c("visit0", "visit1", "visit2")}).
#'   The first element corresponds to the baseline visit. The remaining corresponds to follow-up visits
#' @param corMatrix Numeric correlation matrix specifying within-subject correlation across visits.
#'   Must be symmetric and have dimensions equal to \code{length(visits)}. Users can use the \code{gen_corMatrix()}
#'   function from this package to prepare it.
#' @param covariate_def Optional \pkg{simstudy} definition object created by \code{simstudy::defData()}.
#'   If omitted, the function uses built-in default covariates:
#'   a normally distributed \code{age} (mean 60, SD 10) and a binary \code{pre_diarrhea} (probability 0.4).
#'
#' @return Numeric vector of marginal win odds at each visit.
#' @examples
#' corMatrix <- gen_corMatrix(n_visits = 3,rho=0.6,corstr = "ar1")
#' estimands <- calculate_win_odds(N_approx = 1e4,
#'                                baseprobs = rev(c(0.06,0.11,0.12,0.50,0.21)),
#'                                covs_effects = c("age"=-0.005,"pre_diarrhea"= 0.23),
#'                                time_effects = c(0,0),
#'                                trt_ratio = 1,
#'                                time_trt_effects = c(1,1),
#'                                visits = visits <- c("week0","week4","week8"),
#'                                corMatrix = corMatrix)
#'
#' @export
#'
calculate_win_odds <- function(N_approx = 1e4,
                               baseprobs,
                               covs_effects,
                               time_effects = NULL,
                               trt_ratio,
                               time_trt_effects = NULL,
                               visits,
                               corMatrix,
                               covariate_def=NULL) {
  # ---- (1) Simulate data to approximate covariate distribution ----
  dat_full <- gen_data(
    N = N_approx,
    baseprobs = baseprobs,
    covs_effects = covs_effects,
    time_effects = time_effects,
    trt_ratio = trt_ratio,
    time_trt_effects = time_trt_effects,
    visits = visits,
    corMatrix = corMatrix,
    covariate_def = covariate_def,
    dropZ = FALSE
  )[["long_format"]]

  # ---- (2) Collapse covariate linear predictor (z) into strata ----
  dat_full$z_rounded <- round(dat_full$z, 1)
  table_z <- as.data.frame(table(dat_full$z_rounded), stringsAsFactors = FALSE)
  colnames(table_z) <- c("conditional_effect", "count")
  table_z$conditional_effect <- as.numeric(table_z$conditional_effect)
  table_z$freq <- table_z$count / sum(table_z$count)

  # ---- (3) Vectorized conditional win-odds calculation ----
  # Precompute per-stratum conditional win odds
  wos_conditional <- vapply(
    table_z$conditional_effect,
    FUN = function(ceff) {
      calculate_win_odds_multiple_visits(
        baseprobs = baseprobs,
        conditional_effect = ceff,
        time_effects = time_effects,
        time_trt_effects = time_trt_effects,
        visits = visits
      )
    },
    FUN.VALUE = numeric(length(visits))
  )

  # ---- (4) Compute marginal (population-level) win odds ----
  # Weighted average of conditional win odds across strata
  wos_marginal <- as.vector(wos_conditional %*% table_z$freq)

  # Name the results by visit
  names(wos_marginal) <- visits

  return(wos_marginal)
}


# ---- Helper: Calculate win odds for all visits given a conditional effect ----
calculate_win_odds_multiple_visits <- function(baseprobs,
                                               conditional_effect,
                                               time_effects,
                                               time_trt_effects,
                                               visits) {
  n_visits <- length(visits)

  # Control group probabilities
  probs_control_allvisits <- rbind(
    baseprobs,
    t(vapply(
      time_effects,
      FUN = function(x) calculate_multinomial_probs(x + conditional_effect, baseprobs),
      FUN.VALUE = numeric(length(baseprobs))
    ))
  )

  # Treatment group probabilities
  probs_trt_allvisits <- rbind(
    baseprobs,
    t(vapply(
      time_effects + time_trt_effects,
      FUN = function(x) calculate_multinomial_probs(x + conditional_effect, baseprobs),
      FUN.VALUE = numeric(length(baseprobs))
    ))
  )

  # Compute per-visit win odds
  vapply(
    seq_len(n_visits),
    FUN = function(i) calculate_win_odds_one_visit(
      p0 = probs_control_allvisits[i, ],
      p1 = probs_trt_allvisits[i, ]
    )[["true win odds"]],
    FUN.VALUE = numeric(1)
  )
}


# ---- Helper: Calculate win odds for one visit ----
calculate_win_odds_one_visit <- function(p0, p1) {
  k <- length(p0)
  cumsum_p1 <- cumsum(p1)
  rev_cumsum_p1 <- rev(cumsum(rev(p1)))

  # Probabilities of each pairwise comparison
  P_gt <- sum(rev_cumsum_p1[-1] * p0[-k])     # P(Y1 > Y0)
  P_lt <- sum(cumsum_p1[-k] * p0[-1])         # P(Y1 < Y0)
  P_eq <- sum(p1 * p0)                        # P(Y1 = Y0)

  WP <- P_gt + 0.5 * P_eq
  WO <- WP / (P_lt + 0.5 * P_eq)

  c("true win probability" = WP, "true win odds" = WO)
}
