utils::globalVariables(c("trt", "id", "visit", "z"))
#' Simulate Longitudinal Ordinal Data Under a Proportional Odds Model
#'
#' This function generates correlated ordinal outcome data over multiple visits,
#' based on a proportional odds model with covariates effects, time effects, and
#' time-by-treatment interaction effects. It uses the \pkg{simstudy} framework to define covariate
#' distributions and simulate correlated ordinal responses.
#'
#' @param N Integer. Number of subjects to simulate.
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
#' @param extractTime Logical (default = \code{FALSE}). If \code{TRUE}, extracts a numeric time variable
#'   from visit names (e.g., "visit3" to 3) in the long-format output.
#' @param covariate_def Optional \pkg{simstudy} definition object created by \code{simstudy::defData()}.
#'   If omitted, the function uses built-in default covariates:
#'   a normally distributed \code{age} (mean 60, SD 10) and a binary \code{pre_diarrhea} (probability 0.4).
#' @param dropZ Logical (default = \code{TRUE}). z is the variable recording the amount of shift on the log odds scale compared to baseprobs.
#' If FALSE, it will be included in the output data.frame
#'
#' @return A list with the following components:
#' \describe{
#'   \item{\code{wide_format}}{A data.frame in wide format, with one row per subject and one column per visit.}
#'   \item{\code{long_format}}{A data.frame in long format, with one row per subject-visit combination.}
#'   }
#'
#' @examples
#' set.seed(123)
#' def <- simstudy::defData(varname = "sex", dist = "binary", formula = 0.5)
#' def <- simstudy::defData(def, varname = "age50", dist = "normal", formula = 0, variance = 100)
#' baseprobs <- c(0.3,0.4,0.3)
#' corMat <- gen_corMatrix(n_visits = 4,rho=0.6,corstr = "cs")
#' dat <- gen_data(
#'  N = 200,
#'  baseprobs = baseprobs,
#'  covs_effects = c(sex = 0.5, age50 = -0.02),
#'  time_effects = c(0.1, 0.3, 0.5),
#'  trt_ratio = 1,
#'  time_trt_effects = c(0.05, 0.1, 0.15),
#'  visits = paste0("v", 1:4),
#'  corMatrix = corMat,
#'  covariate_def = def
#' )
#'
#' @export
gen_data <- function(N,
                     baseprobs,
                     covs_effects,
                     time_effects=NULL,
                     trt_ratio = 1,
                     time_trt_effects = NULL,
                     visits,
                     corMatrix,
                     extractTime = FALSE,
                     covariate_def = NULL,
                     dropZ = TRUE)
{
  # ---- (0) General sanity checks ----

  # 0.1 Check N
  if (!is.numeric(N) || length(N) != 1 || N <= 0 || N %% 1 != 0) {
    stop("`N` must be a positive integer.")
  }

  # 0.2 Check baseprobs
  if (!is.numeric(baseprobs) || any(baseprobs < 0) || abs(sum(baseprobs) - 1) > 1e-6) {
    stop("`baseprobs` must be a numeric vector of nonnegative probabilities summing to 1.")
  }

  # 0.3 Check covs_effects
  if (!is.numeric(covs_effects) || is.null(names(covs_effects))) {
    stop("`covs_effects` must be a *named numeric vector* (e.g., c(age = 0.2, pre_diarrhea = -0.5)).")
  }

  # 0.4 Check visits
  if (is.null(visits) || !is.character(visits) || length(visits) < 2) {
    stop("`visits` must be a character vector with at least two elements (e.g., c('visit0', 'visit1', 'visit2')).")
  }

  # 0.5 Check correlation matrix
  if (!is.matrix(corMatrix)) {
    stop("`corMatrix` must be a numeric matrix.")
  }
  if (nrow(corMatrix) != length(visits) || ncol(corMatrix) != length(visits)) {
    stop("`corMatrix` must be a square matrix with dimensions matching the number of visits.")
  }
  if (!all(abs(corMatrix - t(corMatrix)) < 1e-8)) {
    stop("`corMatrix` must be symmetric.")
  }

  # 0.6 Check trt_ratio
  if (!is.numeric(trt_ratio) || trt_ratio <= 0) {
    stop("`trt_ratio` must be a positive numeric value (e.g., 1 for 1:1 randomization).")
  }

  # 0.7 check optional parameters time_effects and time_trt_effects

  n_visits <- length(visits)
  n_followup <- n_visits - 1

  if (!is.null(time_effects) && length(time_effects) != n_followup) {
    stop("`time_effects` length should be number of follow-up visits (visits - 1).")
  }
  if (!is.null(time_trt_effects) && length(time_trt_effects) != n_followup) {
    stop("`time_trt_effects` length should be number of follow-up visits (visits - 1).")
  }

  # Default missing inputs to zeros
  if (is.null(time_effects)) {
    time_effects <- rep(0, n_followup)
  }
  if (is.null(time_trt_effects)) {
    time_trt_effects <- rep(0, n_followup)
  }

  # ---- (1) Marginal probabilities by time/treatment ----
  probs_control_allvisits <- t(sapply(time_effects, function(x)
    calculate_multinomial_probs(x, baseprobs)))
  probs_control_allvisits <- rbind(baseprobs, probs_control_allvisits)

  probs_trt_allvisits <- t(sapply(time_effects + time_trt_effects, function(x)
    calculate_multinomial_probs(x, baseprobs)))
  probs_trt_allvisits <- rbind(baseprobs, probs_trt_allvisits)


  # ---- (2) Define covariate distributions ----
  if (is.null(covariate_def)) {
    message("Covariate distribution not supplied!\nWill use build-in distributions for (continuous) age variable and (binary) preceding diarrhea variable.")
    covariate_def <- simstudy::defData(varname = "age", dist = "normal", formula = 60, variance = 10^2)
    covariate_def <- simstudy::defData(covariate_def, "pre_diarrhea", dist = "binomial", formula = 0.4, variance = 1)
  }

  # Random treatment assignment, default (1:1)
  covariate_def <- simstudy::defData(covariate_def, varname = "trt",
                           dist = "trtAssign",formula = paste0("1;",trt_ratio))

  # ---- (3) Define adjustment variable ----
  cov_names <- names(covs_effects)
  adj_formula <- paste(
    sapply(cov_names, function(v) {
      eff <- covs_effects[v]
      paste0(ifelse(eff >= 0, "+", ""), eff, "*", v)
    }),
    collapse = ""
  )

  covariate_def <- simstudy::defData(covariate_def,
                           varname = "z",
                           formula = adj_formula,
                           dist = "nonrandom")

  # ---- (4) Generate covariates and outcomes ----
  data_cov <- simstudy::genData(N, covariate_def)
  data_trt <- subset(data_cov, trt == 1)
  data_ctl <- subset(data_cov, trt == 0)

  data_trt_ord <- simstudy::genOrdCat(data_trt,
                                      adjVar = "z",
                                      baseprobs = probs_trt_allvisits,
                                      corMatrix = corMatrix,
                                      asFactor = FALSE)

  data_ctl_ord <- simstudy::genOrdCat(data_ctl,
                                      adjVar = "z",
                                      baseprobs = probs_control_allvisits,
                                      corMatrix = corMatrix,
                                      asFactor = FALSE)

  # ---- (5) Combine and reshape ----
  data_wide <- dplyr::bind_rows(data_trt_ord, data_ctl_ord) |>
    dplyr::rename_with(~visits, dplyr::starts_with("grp")) |>
    dplyr::arrange(id)

  data_long <- reshape2::melt(data_wide,
                              measure.vars = visits,
                              variable.name = "visit",
                              value.name = "score") |>
    dplyr::arrange(id)

  if (extractTime) {
    data_long <- data_long |>
      dplyr::mutate(time = as.numeric(stringr::str_extract(visit, "\\d+$")))
  }

  if(dropZ)
  {
    data_wide <- data_wide |>
      dplyr::select(-z)
    data_long <- data_long |>
      dplyr::select(-z)
  }
  return(list(wide_format = data_wide,
              long_format = data_long))

}
