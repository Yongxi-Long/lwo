#' Convert individual-level data to pair-level data
#'
#' @description
#' Converts a dataset of individual-level observations into a pair-level dataset suitable for
#' win/loss/tie comparisons between subjects. Each row in the output corresponds to a unique subject pair
#' observed at a common visit.
#'
#' @param data A data frame of individual-level data to be converted
#' @param id.var A character string specifying the cluster id. For example, patient id for repeated measurements on the same patient
#'  Data are assumed to be sorted so that observations on each cluster appear as continuous rows in data. If data is not sorted this way,
#'  the function will not identify the clusters correctly.
#' @param visit.var A character string specifying the visit within each cluster. This is needed because subjects forming a pair must come from the same visit.
#' @param outcome.var Optional, a character string specifying the outcome variable. Pair-level outcome values are either a `win` (1), a `loss` (0), or a `tie` (0.5) from pairwise comparisons.
#' @param covariates A vector of character strings specifying the covariates. Pair-level covariate values are the difference between individual-level covariate values within the pairs.
#' @param time.vars Optional, needed if time-varying effects are modeled.
#' It should be a scalar or vector of variable names that represent time in the formula. This is IMPORTANT because
#' time variable will not be converted into pair-level.
#' @param larger larger Logical. If `TRUE` (default), larger values of the outcome indicate better outcomes.
#'
#' @return
#' A data frame containing the following columns
#' \describe{
#'   \item{L, R}{Identifiers of the left and right subjects forming each pair.}
#'   \item{pair.id}{Unique identifier for each subject pair.}
#'   \item{covariates}{Pair-level covariate differences (R - L).}
#'   \item{time.vars}{Time variables carried over from the left subject.}
#'   \item{outcome.var}{Pair-level win/loss/tie indicator (if \code{outcome.var} was supplied).}
#'   \item{visit.var}{Visit or time variable corresponding to each observation.}
#' }
#'
#' @examples
#' data("dat_SID")
#' dat_SID_long <- dat_SID$long_format
#' data_pair <- make_pairs(data=dat_SID_long,
#' id.var = "patient_ID",
#' visit.var = "week",
#' time.vars = "week",
#' outcome.var = "GBS_DS",
#' covariates = c("age","treat_ITT","pre_diarrhea","GBS_DS_baseline"),
#' larger = FALSE)
#'
#' @export
make_pairs <- function(data,
                       id.var,
                       visit.var,
                       outcome.var = NULL,   # allow missing
                       covariates,
                       time.vars = NULL,
                       larger = TRUE) {
  # ---- Input checks ----
  if (!is.data.frame(data)) stop("`data` must be a data frame.")
  if (!all(c(id.var, visit.var) %in% names(data)))
    stop("`id.var` and `visit.var` must be column names in `data`.")
  if (!all(covariates %in% names(data)))
    stop("All `covariates` must be present in `data`.")
  if (!is.null(outcome.var) && !(outcome.var %in% names(data)))
    stop("`outcome.var` not found in `data`.")
  if (!is.null(time.vars) && !all(time.vars %in% names(data)))
    stop("All `time.vars` must be present in `data`.")

  # Build a subject by visit table
  tab <- table(data[[id.var]], data[[visit.var]])

  # All possible subject pairs
  all.possible.pairs <- data.frame(t(combn(as.numeric(rownames(tab)), 2)))
  colnames(all.possible.pairs) <- c("L", "R")

  # Keep only pairs that ever meet in the same visit
  meet <- apply(all.possible.pairs, 1, function(i) {
    any(tab[rownames(tab) == as.character(i["L"]), ] +
          tab[rownames(tab) == as.character(i["R"]), ] == 2)
  })
  all.pairs <- all.possible.pairs[meet, ,drop=FALSE]
  all.pairs$pair.id <- seq_len(nrow(all.pairs))

  # Process visit by visit
  visits <- sort(unique(data[[visit.var]]))
  temp <- lapply(visits, function(v) {
    data.this.visit <- data[data[[visit.var]] == v, ,drop=FALSE]

    # Which pairs exist in this visit?
    pairs.this.visit <- intersect(
      which(all.pairs$L %in% data.this.visit[[id.var]]),
      which(all.pairs$R %in% data.this.visit[[id.var]])
    )
    dat.pairs.this.visit <- all.pairs[pairs.this.visit, ,drop=FALSE]

    # Variables to join
    join_vars <- c(id.var, covariates, time.vars)
    if (!is.null(outcome.var)) {
      join_vars <- c(join_vars, outcome.var)
    }

    # Join left subject
    dat.L <- dplyr::left_join(
      dat.pairs.this.visit,
      dplyr::select(data.this.visit, dplyr::all_of(join_vars)),
      by = c("L" = id.var)
    )
    # Join right subject
    dat.R <- dplyr::left_join(
      dat.pairs.this.visit,
      dplyr::select(data.this.visit, dplyr::all_of(join_vars)),
      by = c("R" = id.var)
    )

    # Compute differences (R - L) for covariates
    dat.RminusL <- dat.R
    for (cov in covariates) {
      dat.RminusL[[cov]] <- dat.R[[cov]] - dat.L[[cov]]
    }

    # Compute pair outcome (only if outcome.var provided)
    if (!is.null(outcome.var)) {
      diff_outcome <- dat.R[[outcome.var]] - dat.L[[outcome.var]]
      dat.RminusL[[outcome.var]] <- ifelse(diff_outcome > 0, 1,
                                           ifelse(diff_outcome < 0, 0, 0.5))
      if (!larger) {
        dat.RminusL[[outcome.var]] <- 1 - dat.RminusL[[outcome.var]]
      }
    }

    # Keep time.vars (take from L, same as R by construction)
    if (!is.null(time.vars)) {
      for (tv in time.vars) {
        dat.RminusL[[tv]] <- dat.L[[tv]]
      }
    }

    # Add visit
    dat.RminusL[[visit.var]] <- v

    # Add pair id
    dat.RminusL$pair.id <- dat.pairs.this.visit$pair.id

    return(dat.RminusL)
  })

  # Combine visits
  data.pair <- dplyr::bind_rows(temp)
  data.pair <- data.pair[order(data.pair$pair.id, data.pair[[visit.var]]), ]

  return(data.pair)
}

