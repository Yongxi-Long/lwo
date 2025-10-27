#' Calculate Shifted Multinomial Probabilities
#'
#' @description
#' Computes category probabilities given log odds shift(s) to a baseline probability vector
#' under either proportional odds or non-proportional odds
#'
#' @param lnORs Numeric vector of log odds ratio shifts to apply to the cumulative logits
#'   of the baseline distribution \code{p0}. If \code{common = TRUE}, only the first value
#'   is used and applied uniformly across all cumulative logits.
#' @param p0 p0 Numeric vector of baseline category probabilities (must sum to 1).
#' @param common Logical (default = \code{TRUE}). If \code{TRUE}, applies a common
#'   log-odds shift to all cumulative logits; if \code{FALSE}, allows category-specific shifts.
#'
#' @return Numeric vector of shifted category probabilities (same length as \code{p0})
#' @examples
#' # proportional odds
#' p1 <- calculate_multinomial_probs(lnORs=0.5,p0=c(0.1,0.2,0.3,0.4))
#' # non-proportional odds
#' p2 <- calculate_multinomial_probs(lnORs=c(-0.5,0.1,0.5),p0=c(0.1,0.2,0.3,0.4),common=FALSE)
#' @export
#'
calculate_multinomial_probs <- function(lnORs,p0,common=TRUE) {
  # ---- Input checks ----

  # Check baseline probability vector
  if (!is.numeric(p0) || any(p0 < 0) || abs(sum(p0) - 1) > 1e-6) {
    stop("`p0` must be a numeric probability vector that sums to 1.")
  }

  k <- length(p0)
  if (k < 2) {
    stop("`p0` must have at least two categories.")
  }

  # Check lnORs vector consistency
  if (common) {
    if (length(lnORs) < 1) {
      stop("When `common = TRUE`, `lnORs` must contain at least one value.")
    }
  } else {
    if (length(lnORs) != (k - 1)) {
      stop(paste0("When `common = FALSE`, `lnORs` must have length = length(p0) - 1 (",
                  k - 1, " in this case)."))
    }
  }
  # number of ordinal categories
  k <- length(p0)
  if(common) lnORs = rep(lnORs[1],k-1)
  # cumulative logit of p0
  logit0 <- -qlogis(cumsum(p0))[-k]
  # cumulative logit of p1 is shifted by supplied log ORs
  logit1 <- logit0 + lnORs
  # get cumulative probabilities for p1
  cump1 <- plogis(-logit1)
  # get multinomial probabilities for all categories
  p1 <- c(cump1,1)-c(0,cump1)
  return(p1)
}
