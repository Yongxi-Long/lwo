#' Longitudinal Win Odds
#'
#' Fits a longitudinal probabilistic index model and estimate longitudinal win odds
#' for ordinal repeated measurements.
#'
#' @param formula A model formula specifying the outcome and the covariates.
#' @param data A data frame containing the variables in the model.
#' @param family Family of the link function, only logit link gives the win odds.
#' @param id.var A character string specifying the cluster id. For example, patient id for repeated measurements on the same patient
#'  Data are assumed to be sorted so that observations on each cluster appear as continuous rows in data. If data is not sorted this way,
#'  the function will not identify the clusters correctly.
#' @param visit.var A character string specifying the visit within each cluster. This is needed because subjects forming a pair must come from the same visit.
#' @param time.vars Optional, needed if time-varying effects are modeled.
#' It should be a scalar or vector of variable names that represent time in the formula. This is IMPORTANT because
#' time variable will not be converted into pair-level.
#' @param larger Logical. If `TRUE` (default), larger values of the outcome indicate better outcomes.
#' @param corstr Working correlation structure; passed to `geepack::geese.fit`.
#'   Options include `"independence"`, `"exchangeable"`, `"ar1"`, `"unstructured"`,
#'   `"userdefined"`, and `"fixed"`.
#' @param zcor A design matrix for correlation parameters.
#' @param std.err Type of standard error to compute.
#'   `"san.se"` for standard sandwich estimator
#'   `"san.se.modified"` for the modified version, see Mancle & DeRouen, 2001, Biometrics
#' @param weights Optional observation-level weights.
#' @param contrasts A list giving contrasts for some or all of the factors appearing in the model formula. The elements of the list should have the same name as the variable and should be either a contrast matrix (specifically, any full-rank matrix with as many rows as there are levels in the factor), or else a function to compute such a matrix given the number of levels.
#'
#' @importFrom stats binomial glm.fit is.empty.model model.frame model.matrix model.response model.weights plogis qlogis terms
#'
#' @details
#' The function internally converts individual-level
#' data to pair-level and performs a logistic regression on the binary pair-level
#' outcome and covariates, then uses a sandwich estimator to obtain the correct
#' variance-covariance matrix accounting for pairwise correlation. Setting std.err="san.se"
#' gives the sandwich estimate proposed in the paper. Setting std.err = "san.se.modified" gives
#' a modified sandwich estimate based on the modification of Mancle & DeRouen, 2001, Biometrics.
#'
#' @return
#' An object of class \code{"lwo"}, inheriting from
#' \code{"geeglm"}, \code{"gee"}, \code{"glm"}, and \code{"lm"}.
#' Components include:
#' \itemize{
#'   \item \code{coefficients} – estimated regression coefficients
#'   \item \code{varcov} – robust covariance matrix estimate
#'   \item \code{fitted.values} – model-fitted probabilities
#'   \item \code{model.matrix} – design matrix at the pair level
#'   \item \code{data.pair} – transformed pair-level data
#' }
#'
#' @examples
#' data("dat_SID")
#' dat_SID_long <- dat_SID$long_format
#' mod.lwo.spline <- lwo(
#'   GBS_DS ~ treat_ITT + age + pre_diarrhea + GBS_DS_baseline+
#'    treat_ITT:splines::ns(week,knots = 2),
#'  data = dat_SID_long,
#'  id.var = "patient_ID",
#'  visit.var = "week",
#'  time.vars = "week",
#'  corstr = "ar1",
#'  larger = FALSE
#' )
#'@export
lwo <- function(
    formula,
    data = parent.frame(),
    family = binomial(),
    id.var,
    visit.var,
    time.vars = NULL,
    larger=TRUE,
    corstr = "independence",
    zcor = NULL,
    std.err="san.se",
    weights,
    contrasts = NULL)
{
  if (!requireNamespace("geepack", quietly = TRUE)) {
    stop("Package 'geepack' is required but not installed.")
  }
  if (missing(data))
    data <- environment(formula) #If data is not explicitly provided, R uses the environment in which the formula was created to find the variables (e.g., the global environment or a specific data frame).
  ##---- Prepare model matrix and response vector
  # returns a full call to the function. e.g., glm(formula = y~x data=df, family = binomial())
  # expand.dot: expands any ... arguments to include all their contents explicitly. This is useful for saving the full call for reproducing results or diagnostics.
  call <- match.call(expand.dots = TRUE)
  mf.call <- match.call(expand.dots = FALSE)
  # find position of the following terms in mf
  m <- match(c("formula",
               "data",
               "weights"), names(mf.call), 0L) # returns 0 if not the term is not inputted
  # reduce mf to only contain the above terms
  mf.call <- mf.call[c(1L, m)] #c(1L, m) ensures the function name (first element) and the matched arguments are retained.
  # Other arguments, such as family or additional parameters, are excluded temporarily.
  mf.call$drop.unused.levels <- TRUE # If factors in the data frame have levels not present in the subset used for modeling, this avoids errors or inconsistencies during the fitting process.
  # The first element of mf (which was previously the function name, e.g., lwo) is replaced with stats::model.frame.
  # stats::model.frame(formula = y ~ x, data = df, ...)
  mf.call[[1L]] <- quote(stats::model.frame)

  # mf is the model frame
  # eval() executes the stats::model.frame() call in the appropriate context (parent.frame()), producing the model frame (mf).
  # The model frame includes the response variable, predictors, and any weights, offsets, or subsets specified in the arguments.
  mf <- eval(mf.call, parent.frame())
  # stats::model.frame() attaches a "terms" attribute to the resulting model frame. This "terms" object contains detailed information about the model formula, such as:
  # The response and predictor variables.
  # Their roles (response vs. predictors).
  # Any interactions or transformations specified in the formula.
  # mt is used later in the fitting process to handle predictors correctly.
  mt <- attr(mf, "terms")
  # return(list(mf,mt))

  # ###############################################################################
  # Create new data frame from input dataframe: convert outcome and covariate to pair level
  # ###############################################################################
  # we need input id name and time name, because time variable will not be converted to pair level
  relevant_vars <- extract_variable_names(formula)
  outcome.var <- relevant_vars[1]
  covariates <-  relevant_vars[2:length(relevant_vars)]

  data.pair <- make_pairs(data,
                          id.var,
                          visit.var,
                          outcome.var,
                          covariates,
                          time.vars,
                          larger)
  ##---- Modify model framework mf
  # since we have converted data frame from individual level to pair level
  # we have to modify the model framework accordingly
  mf.pair <- model.frame(formula, data.pair, drop.unused.levels = TRUE)
  mt.pair <- terms(mf.pair)
  # return(list(mf.pair,data.pair))
  #---- Construct X and Y
  Y  <- model.response(mf.pair, "numeric")
  N <- NROW(Y) # number of pairs
  X  <- if (!is.empty.model(mt.pair))
  {
    model.matrix(mt.pair, mf.pair, contrasts)
  }else
  {matrix(, N, 0)}
  # remove intercept
  X <- X[,colnames(X)!="(Intercept)",drop=FALSE]

  ##---- Clustering variable
  id.pair <- data.pair[,"pair.id"]
  if (is.null(id.pair)) stop("pair id variable not found.")

  ##---- Get family, for longitudinal win odds, only logit link is valid
  if (is.character(family))
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family))
    family <- family()
  if (is.null(family$family)) {
    print(family)
    stop("'family not recognized")
  }
  if (family$family!="binomial")
  {
    warning("Longitudinal win odds is only valid with a logit link!")
  }


  ##---- Working correlation structures
  if (corstr=="fixed" && is.null(zcor)){
    stop("When corstr is 'fixed' then 'zcor' must be given\n")
  }
  CORSTRS <- c("independence", "exchangeable", "ar1", "unstructured",
               "userdefined", "fixed")
  # do partial match if users put ind instead of independence
  corstrv <- pmatch(corstr, CORSTRS, -1)
  corstr  <- CORSTRS[corstrv]

  ##---- Weights for the observations (not used)
  weights <- model.weights(mf)
  if (is.null(weights)) weights <- rep(1, N)

  ## Check that factors in model do not have unused levels in data
  ## (otherwise R crashes).
  vars <- all.vars(formula)
  stopIt <- FALSE
  for(ii in seq_along(vars)){
    vv <- vars[ii]
    if(!is.na(match(vv, names(mf))) && is.factor(mf[,vv])){
      if (length(unique(mf[,vv])) != length(levels(mf[,vv]))){
        cat("Factors not allowed to have unused levels...\n")
        cat(" Levels of factor", vv,":", paste(levels(mf[, vv]), sep=' '),"\n")
        cat(" Used levels of factor", vv,":", paste(unique(mf[, vv]), sep=' '),"\n")
        stopIt <- TRUE
      }
    }
  }
  if (stopIt)
    stop("Can not continue...\n")


  ##---- Create dummy glm object in order for the model to fit in generic functions that have glm method
  glmFit <- suppressWarnings(glm.fit(x = X, y = Y, weights = rep(1,N), start = NULL, etastart = NULL,
                                     mustart = NULL, offset = rep(0,N), family = binomial(),
                                     control = list(), intercept = TRUE, singular.ok = TRUE))
  class(glmFit) <- "glm"
  glmFit$terms <- mt.pair
  glmFit$model <- mf.pair
  # modelmat <- model.matrix(glmFit)[-1] # no intercept
  # qqrr     <- qr(modelmat)
  # if (qqrr$rank < ncol(modelmat)){
  #   print(head(modelmat))
  #   stop("Model matrix is rank deficient; geeglm can not proceed\n")
  # }

  ##---- Use geese.fit to estimate parameters, the variances are wrong which will be
  # fixed later
  ans <- suppressWarnings(geepack::geese.fit(x=X, y=Y, id=id.pair,
                                             family=family,
                                             corstr=corstr,zcor = zcor,
                                             offset=rep(0,N), soffset=rep(0,N),
                                             weights=weights, waves = NULL, control = geepack::geese.control(),
                                             corp = NULL, b = NULL,
                                             alpha = NULL, gm = NULL,  mean.link = NULL, variance = NULL,
                                             cor.link = "identity", sca.link = "identity", link.same = TRUE,
                                             scale.fix = TRUE, scale.value = 1))
  ans <- c(ans, list(call = call, formula = formula))
  class(ans)  <- "geese"
  ans$model.matrix <- X
  ans$Y <- Y
  ans$data.pair <- data.pair
  # return(ans)
  ################################################################################
  # update to geese, modify the standard error vbeta
  ################################################################################
  # get the cluster-specifc score equations
  # a vector of cluster size for each pair
  clusz <- ans$clusz
  # model matrix
  # chop to cluster-level
  X_list <- chop_matrix(X,clusz)
  # outcomes
  # chop to cluster level
  Y_list <- chop_matrix(Y,clusz)
  # estimated parameters
  betas_hat <- c(ans$beta)
  rho_hat <- ans$alpha
  # fitted values
  mu_hat <- plogis(c(X%*%betas_hat))
  # chop fitted values
  mu_hat_list <- chop_matrix(mu_hat,clusz)
  # pair information
  pair_df <- dplyr::select(data.pair,c("L","R","pair.id"))
  IDs_subjectL <- unique(pair_df[,c("L","pair.id")])[,"L"]
  IDs_subjectR <- unique(pair_df[,c("R","pair.id")])[,"R"]

  # The derivative matrix
  D_list <- R_list <- A_list <- V_inv_list <- U_list <- Bread_list <- vector("list",length = length(clusz))
  for (i in 1:length(clusz))
  {
    csize <- clusz[i]
    X_i <- X_list[[i]]
    Y_i <- Y_list[[i]]
    mu_hat_i <- mu_hat_list[[i]]
    A_i <- diag(mu_hat_i*(1-mu_hat_i))
    A_list[[i]] <- A_i
    D_i <- A_i%*%X_i
    D_list[[i]] <- D_i
    # working correlation matrix
    R_i <- gen_corMatrix(csize,rho=rho_hat,corstr = corstr)
    R_list[[i]] <- R_i
    V_inv_i <- solve(sqrt(A_i)%*%R_i%*%sqrt(A_i))
    V_inv_list[[i]] <- V_inv_i
    U_i <- t(D_i)%*%V_inv_i%*%(Y_i - mu_hat_i)
    U_list[[i]] <- U_i |> as.vector()

    # bread part of the sandwich estimator
    Bread_list[[i]] <- t(D_i)%*%V_inv_i%*%D_i |>
      as.matrix()# p by p matrix
  }

  # the sum of score equations (should be zero)
  # rowSums(do.call(cbind,U_list))
  # matrix of score functions
  U_mat <- do.call(rbind,U_list)

  # the Bread part
  Bread <- Reduce("+", Bread_list)
  # Bread_array <- simplify2array(Bread_list)
  # Bread <- apply(Bread_array, c(1,2), sum)

  if(std.err == "san.se.modified")
  {
    M_list <- vector("list",length = length(clusz))
    for (i in 1:length(clusz))
    {
      csize <- clusz[i]
      D_i <- D_list[[i]]
      V_inv_i <- V_inv_list[[i]]
      Y_i <- Y_list[[i]]
      mu_hat_i <- mu_hat_list[[i]]
      # the H matrix from Mancle and DeRounen's modification
      H_ii <- D_i%*%solve(Bread)%*%t(D_i)%*%V_inv_i
      M_i <- t(D_i)%*%V_inv_i%*%solve((diag(csize)-H_ii))%*%(Y_i - mu_hat_i) |>
        as.vector()
      M_list[[i]] <- M_i
    }
    M_mat <- do.call(rbind,M_list)
    ## The Meat part
    mat <- M_mat
  } else
  {
    mat <- U_mat
  }
  # correlated pairs like (1,2) and (1,3), shared subject is on the same side
  shared.factor=1
  # correlated pairs like (1,2) and (2,4), shared subject is on either side
  switched.factor=1
  # self correlation, (1,2) and (1,2)
  self.factor=1
  # Compute column sums across rows of score matrix-like object for pairs with the same subject
  Usum1.tmp <- rowsum(mat,IDs_subjectL,reorder=FALSE) # compress mat with the same L subject, the 1st row = colSums(mat[IDs_subjectL==1,])
  Usum2.tmp <- rowsum(mat,IDs_subjectR,reorder=FALSE)
  Usum1  <- matrix(nrow = nrow(mat), ncol = ncol(mat),0)
  Usum2  <- matrix(nrow = nrow(mat), ncol = ncol(mat),0)
  Usum1[unique(IDs_subjectL),] <- Usum1.tmp
  Usum2[unique(IDs_subjectR),] <- Usum2.tmp
  UtUshared <- crossprod(Usum1) + crossprod(Usum2)
  UtUswitched <- crossprod(Usum1,Usum2) + crossprod(Usum2,Usum1)
  UDiag<-crossprod(mat) #Is counted twice as shared.factor, but needs to be counted as self.factor
  Meat <- shared.factor*UtUshared  +
    (switched.factor)*UtUswitched +
    (self.factor-2*shared.factor)*(UDiag)
  # sandwich estimate
  Bread_inv <- solve(Bread)
  varcov <- Bread_inv%*%Meat%*%Bread_inv
  colnames(varcov) <- rownames(varcov) <- names(betas_hat)
  # replace the wrong var-cov matrix in ans
  ans$vbeta <- varcov

  ##---- Prepare output
  out <- glmFit # build on the skeleton of glm
  toDelete    <- c("R", "deviance", "aic", "null.deviance", "iter",
                   "df.null", "converged", "boundary")
  out[match(toDelete,names(out))] <- NULL
  out$method <- "geese.fit"
  out$geese <- ans
  out$weights <- ans$weights
  out$coefficients <- ans$beta

  out$fitted.values <- family(out)$linkinv(out$linear.predictors)
  out$modelInfo <- ans$model
  out$data.pair <- ans$data.pair
  out$call <- ans$call
  out$corstr    <- ans$model$corstr
  out$cor.link  <- ans$model$cor.link
  out$control   <- ans$control
  out$std.err   <- std.err
  out$varcov <- varcov
  out$model.matrix <- ans$model.matrix
  class(out)    <- c("lwo","geeglm", "gee", "glm", "lm")
  return(out)
}


#' Calculate Variance-Covariance Matrix for a Fitted Model Object
#'
#' @description
#' This method returns the estimated variance-covariance matrix of model coefficients
#' from an object of class \code{"lwo"}.
#'
#' @param object An object of class \code{"lwo"}, typically the result of a call to \code{\link{lwo}}.
#' @param ... Additional arguments (currently ignored).
#'
#' @return
#' A symmetric matrix containing the estimated covariances between parameter estimates
#' from the fitted longitudinal win odds model.
#'
#' @details
#' The variance-covariance matrix is extracted from the fitted model object.
#' If \code{object$varcov} is missing or not a square matrix, an informative error is raised.
#' Column and row names are set to match the coefficient names.
#'
#' @examples
#' data("dat_SID")
#' dat_SID_long <- dat_SID$long_format
#' mod.lwo.spline <- lwo(
#'   GBS_DS ~ treat_ITT + age + pre_diarrhea + GBS_DS_baseline+
#'    treat_ITT:splines::ns(week,knots = 2),
#'  data = dat_SID_long,
#'  id.var = "patient_ID",
#'  visit.var = "week",
#'  time.vars = "week",
#'  corstr = "ar1",
#'  larger = FALSE
#' )
#' vcov(mod.lwo.spline)
#'
#' @export
vcov.lwo <- function(object, ...){
  if (!inherits(object, "lwo")) {
    stop("Input must be an object of class 'lwo'.")
  }

  if (is.null(object$varcov)) {
    stop("No variance-covariance matrix found in 'object$varcov'.")
  }

  out <- object$varcov

  if (!is.matrix(out)) {
    stop("'object$varcov' must be a matrix.")
  }
  coefs <- tryCatch(stats::coef(object), error = function(e) NULL)
  if (!is.null(coefs)) {
    pn <- names(coefs)
    n_coef <- length(pn)

    if (nrow(out) != n_coef || ncol(out) != n_coef) {
      warning("Dimension mismatch: resetting dimension names only if compatible.")
      if (length(pn) == nrow(out)) {
        dimnames(out) <- list(pn, pn)
      }
    } else {
      dimnames(out) <- list(pn, pn)
    }
  }
  cat("\nVariance is re-estimated to account for correlation between pseudo pairs\n")
  return(out)
}

#' Summarizing Longitudinal Win Odds Model Fits
#'
#' @description
#' The \code{summary} method for objects of class \code{"lwo"}.
#' Produces coefficient tables, correlation structure summaries, and model fit details
#' from a fitted longitudinal win odds model.
#'
#'
#' @param object object An object of class `lwo`, usually a resultof a call to \code{lwo}
#' @param ... Additional arguments (currently ignored).
#'
#' @return
#' An object of class \code{"summary.lwo"}, containing:
#' \describe{
#'   \item{call}{The function call that produced the \code{lwo} object.}
#'   \item{coefficients}{A data frame of estimated coefficients, standard errors, Wald statistics, and p-values.}
#'   \item{corstr}{A character string describing the assumed working correlation structure.}
#'   \item{corr}{A data frame summarizing estimated correlation parameters and their standard errors (if any).}
#' }
#'
#' @examples
#' data("dat_SID")
#' dat_SID_long <- dat_SID$long_format
#' mod.lwo.spline <- lwo(
#'   GBS_DS ~ treat_ITT + age + pre_diarrhea + GBS_DS_baseline+
#'    treat_ITT:splines::ns(week,knots = 2),
#'  data = dat_SID_long,
#'  id.var = "patient_ID",
#'  visit.var = "week",
#'  time.vars = "week",
#'  corstr = "ar1",
#'  larger = FALSE
#' )
#' summary(mod.lwo.spline)
#' @export
summary.lwo <- function(object, ...) {
  # ---- (1) Basic validation ----
  if (is.null(object$terms)) {
    stop("Invalid 'lwo' object: no 'terms' component found.")
  }
  if (is.null(object$geese) || !is.list(object$geese)) {
    stop("Invalid 'lwo' object: missing 'geese' component (expected from geepack).")
  }

  z <- object
  ans <- list(call = z$call)

  # ---- (2) Coefficient table ----
  beta <- z$geese$beta
  vbeta <- z$geese$vbeta

  if (is.null(beta) || is.null(vbeta)) {
    warning("Model appears to have no coefficients or variance estimates.")
    coefmat <- data.frame()
  } else {
    se <- sqrt(diag(vbeta))
    wald <- (beta / se)^2
    pval <- stats::pchisq(wald, df = 1, lower.tail = FALSE)
    coefmat <- data.frame(
      Estimate = beta,
      `Std.err` = se,
      Wald = wald,
      `Pr(>|W|)` = pval,
      check.names = FALSE
    )
    rownames(coefmat) <- names(beta)
  }

  ans$coefficients <- coefmat

  # ---- (3) Correlation structure ----
  ans$corstr <- if (!is.null(z$corstr)) z$corstr else NA_character_

  if (!is.null(z$geese$alpha)) {
    valpha <- z$geese$valpha
    corrmat <- data.frame(
      Estimate = z$geese$alpha,
      `Std. Err` = if (!is.null(valpha)) sqrt(diag(valpha)) else NA_real_,
      check.names = FALSE
    )
  } else {
    corrmat <- data.frame()
  }

  ans$corr <- corrmat

  # ---- (4) Model fit details ----
  if (!is.null(z$geese$criterion)) {
    ans$criterion <- z$geese$criterion
  }

  # ---- (5) Assign class ----
  class(ans) <- "summary.lwo"
  return(ans)
}

#' @export
print.summary.lwo <- function(x,
                              digits = max(3, getOption("digits") - 3),
                              quote = FALSE,
                              prefix = "",
                              ...) {
  # ---- (0) Input validation ----
  if (!inherits(x, "summary.lwo")) {
    stop("Input must be an object of class 'summary.lwo'.")
  }

  # ---- (1) Model call ----
  cat("\nCall:\n")
  if (!is.null(x$call)) {
    print(x$call)
  } else {
    cat("Call not available.\n")
  }

  # ---- (2) Coefficient table ----
  if (!is.null(x$coefficients)) {
    cat("\nCoefficients:\n")
    stats::printCoefmat(
      x$coefficients,
      digits = digits,
      signif.stars = getOption("show.signif.stars", TRUE)
    )
  } else {
    cat("\n(No coefficients available)\n")
  }

  # ---- (3) Correlation structure ----
  if (!is.null(x$corstr)) {
    cat("\nTemporal Correlation Structure:", x$corstr, "\n")
  } else {
    cat("\nTemporal Correlation Structure: (not specified)\n")
  }

  # ---- (4) Correlation parameter estimates ----
  if (!is.null(x$corr) &&
      !is.na(pmatch(x$corstr, "independence", nomatch = NA)) &&
      pmatch(x$corstr, "independence", nomatch = 0) == 0) {
    cat("\nEstimated Correlation Parameters:\n")
    print(x$corr, digits = digits)
  }

  invisible(x)
}

#' Predict Method for \code{lwo} Objects
#'
#' @description
#' Computes predicted linear predictors or response-scale win-probability estimates
#' from a fitted \code{lwo} model. Optionally returns confidence intervals.
#'
#' @details
#' The \code{predict.lwo()} method takes subject-level \code{newdata} and converts it
#' into pseudo-pair observations using the same pair construction used during model
#' fitting. Predictions are made on the resulting pair-level design matrix using the
#' estimated regression coefficients from the fitted model.
#'
#' If \code{type = "link"} (default), predictions are returned on the linear predictor
#' scale. If \code{type = "response"}, predictions are mapped through the inverse logit
#' and returned as estimated conditional win probabilities.
#'
#' If \code{conf.int = TRUE}, point estimates and corresponding Wald confidence intervals
#' are returned. Confidence intervals are computed using the estimated variance-covariance
#' matrix of the regression coefficients stored in the fitted model object.
#'
#' @param object A fitted \code{lwo} model object.
#' @param newdata Data frame of subject-level data at the visits for which predictions
#'   are desired. Must include all covariates appearing in the fitted model formula.
#' @param id.var Character string giving the subject identifier variable in \code{newdata}.
#' @param visit.var Character string giving the visit/time variable in \code{newdata}.
#' @param time.vars Optional character vector of time-varying covariate names. Defaults to \code{NULL}.
#' @param type Character string specifying the scale of prediction:
#'   \describe{
#'     \item{\code{"link"}}{Return predictions on the linear predictor scale (default).}
#'     \item{\code{"response"}}{Return predictions on the response scale (win probability).}
#'   }
#' @param conf.int Logical; if \code{TRUE}, return Wald confidence intervals. Default is \code{FALSE}.
#' @param alpha Significance level for confidence intervals. Default is \code{0.05}.
#' @param na.action Function indicating how to handle missing values. Default is \code{na.pass}.
#'
#' @importFrom stats binomial glm.fit is.empty.model model.frame model.matrix model.response model.weights plogis qlogis terms delete.response na.pass qnorm
#'
#' @return
#' If \code{conf.int = FALSE}, returns a numeric vector of predicted values.
#' If \code{conf.int = TRUE}, returns a matrix with columns:
#' \code{"lower.CI"}, \code{"estimate"}, and \code{"upper.CI"}.
#'
#' @examples
#' data("dat_SID")
#' dat_SID_long <- dat_SID$long_format
#' mod.lwo.spline <- lwo(
#'   GBS_DS ~ treat_ITT + age + pre_diarrhea + GBS_DS_baseline+
#'    treat_ITT:splines::ns(week,knots = 2),
#'  data = dat_SID_long,
#'  id.var = "patient_ID",
#'  visit.var = "week",
#'  time.vars = "week",
#'  corstr = "ar1",
#'  larger = FALSE
#' )
#' newdata <- dat_SID_long |>
#' dplyr::filter(patient_ID %in% c(1,4))
#' predict(object = mod.lwo.spline,
#'            newdata = newdata,
#'            id.var = "patient_ID",
#'            visit.var = "week",
#'            time.vars = "week",
#'            type = "response",
#'            conf.int = TRUE)
#'
#' @export
predict.lwo <- function(object,
                        newdata = NULL,
                        id.var,
                        visit.var,
                        time.vars = NULL,
                        type = c("link", "response"),
                        conf.int = FALSE, alpha = 0.05,
                        na.action = na.pass) {

  type <- match.arg(type)

  if (is.null(newdata))
    stop("newdata must be supplied for prediction.", call. = FALSE)

  # Original model formula and terms object (stores spline bases & contrasts)
  formula <- formula(object)
  terms.orig <- delete.response(terms(object))

  # Extract variable names
  relevant_vars <- extract_variable_names(formula)
  outcome.var <- relevant_vars[1]
  covariates <- relevant_vars[-1]

  # ---- Convert newdata to *pair-level*, but without generating outcomes ----
  data.pair <- make_pairs(
    data       = newdata,
    id.var     = id.var,
    visit.var  = visit.var,
    outcome.var = NULL,
    covariates = covariates,
    time.vars  = time.vars
  )

  # Add placeholder outcome variable (required by model.frame, not used)
  data.pair[[outcome.var]] <- 0

  # Construct model frame + model matrix using *stored* model structure
  mf.pair <- model.frame(terms.orig, data.pair, na.action = na.action)
  mod_mat <- model.matrix(terms.orig, mf.pair)

  # Remove intercept safely
  if ("(Intercept)" %in% colnames(mod_mat)) {
    mod_mat <- mod_mat[, colnames(mod_mat) != "(Intercept)", drop = FALSE]
  }

  # ---- Ensure columns match model coefficients ----
  coef <- object$coefficients

  missing_cols <- setdiff(names(coef), colnames(mod_mat))
  if (length(missing_cols) > 0) {
    stop(
      "The following required model terms are missing in newdata after processing:\n",
      paste(missing_cols, collapse = ", "),
      "\nThis typically occurs when spline bases or factor levels do not match.",
      call. = FALSE
    )
  }

  # Reorder mod_mat columns to match coefficient order
  mod_mat <- mod_mat[, names(coef), drop = FALSE]

  # ---- Compute linear predictor ----
  lp <- as.vector(mod_mat %*% coef)

  if (!conf.int) {
    if (type == "link") return(lp)
    else return(plogis(lp))
  }

  # ---- Compute standard errors and CIs ----
  var_beta <- object$geese$vbeta
  se_lp <- sqrt(rowSums((mod_mat %*% var_beta) * mod_mat))

  z <- qnorm(1 - alpha / 2)
  lp_lwr <- lp - z * se_lp
  lp_upr <- lp + z * se_lp

  if (type == "link") {
    out <- cbind(lower.CI = lp_lwr, estimate = lp, upper.CI = lp_upr)
  } else {
    out <- cbind(lower.CI = plogis(lp_lwr),
                 estimate = plogis(lp),
                 upper.CI = plogis(lp_upr))
  }

  class(out) <- "predict.lwo"
  return(out)
}

## helper functions not exported
# Function to extract variable names from a formula
extract_variable_names <- function(formula) {
  # Get all terms from the formula
  terms_object <- stats::terms(formula)

  # Extract the variables as they appear in the formula
  all_vars <- all.vars(formula)

  # Remove functions and transformations
  # Extract the variables directly involved in functions
  raw_vars <- attr(terms_object, "variables")[-1]  # Remove the response
  raw_vars <- sapply(raw_vars, function(x) {
    if (is.call(x)) all.vars(x) else as.character(x)
  })

  # Flatten and deduplicate the list of variables
  unique(unlist(raw_vars))
}

chop_matrix <- function(X, clusz) {
  sublist <- list()
  start <- 1
  for (size in clusz) {
    end <- start + size - 1
    if (is.matrix(X)) {
      sublist[[length(sublist) + 1]] <- X[start:end, , drop = FALSE]
    } else {
      sublist[[length(sublist) + 1]] <- X[start:end]
    }
    start <- end + 1
  }
  return(sublist)
}
