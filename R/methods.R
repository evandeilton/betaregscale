# ============================================================================ #
# S3 methods for objects of class "betaregscale"
#
# Output style follows the betareg package convention:
#   - coef(), vcov() accept  model = c("full", "mean", "precision")
#   - summary() produces separate tables for mean and precision
#   - print() shows the call + compact coefficient vectors
#   - Wald z-tests use pnorm (not pt)
# ============================================================================ #

# -- Class validation helper ------------------------------------------------ #

#' Validate a betaregscale object
#' @param x Object to validate.
#' @param call. Logical; passed to \code{stop()}.
#' @keywords internal
.check_class <- function(x, call. = FALSE) {
  if (!inherits(x, "brs")) {
    stop(
      "Expected an object of class 'brs', got '",
      paste(class(x), collapse = "', '"), "'.",
      call. = call.
    )
  }
}


# -- Display-name helpers --------------------------------------------------- #
#  These are used only in print/summary to produce clean output.
#  Internal names (names(est), rownames(hessian), etc.) are NOT changed,
#  so user code that indexes by name continues to work.

#' Strip (phi)_ prefix for summary display
#' @keywords internal
#' @noRd
.pretty_phi_names <- function(nms) {
  sub("^\\(phi\\)_", "", nms)
}

#' Convert internal Cholesky RE names to readable display names
#'
#' Mapping:
#'   (re_chol_logsd)_X|g  ->  logSD.X|g
#'   (re_chol)_X:Y|g      ->  cov.X:Y|g
#' @keywords internal
#' @noRd
.pretty_re_names <- function(nms) {
  nms <- sub("^\\(re_chol_logsd\\)_", "logSD.", nms)
  nms <- sub("^\\(re_chol\\)_", "cov.", nms)
  nms
}


# -- Extract coefficients --------------------------------------------------- #

#' Extract model coefficients
#'
#' @param object A fitted \code{"betaregscale"} object.
#' @param model  Character: which component to return.
#'   \code{"full"} (default) returns all parameters,
#'   \code{"mean"} returns only the mean-model coefficients,
#'   \code{"precision"} returns only the precision coefficients.
#' @param ... Ignored.
#'
#' @return Named numeric vector of estimated parameters.
#'
#' @seealso \code{\link{brs}}, \code{\link{brs_est}}, \code{\link{vcov.brs}}
#'
#' @examples
#' \donttest{
#' dat <- data.frame(
#'   y = c(
#'     0, 5, 20, 50, 75, 90, 100, 30, 60, 45,
#'     10, 40, 55, 70, 85, 25, 35, 65, 80, 15
#'   ),
#'   x1 = rep(c(1, 2), 10)
#' )
#' prep <- brs_prep(dat, ncuts = 100)
#' fit <- brs(y ~ x1, data = prep)
#' coef(fit)
#' coef(fit, model = "mean")
#' coef(fit, model = "precision")
#' }
#'
#' @method coef brs
#' @importFrom stats coef
#' @export
coef.brs <- function(object,
                     model = c("full", "mean", "precision"),
                     ...) {
  .check_class(object)
  model <- match.arg(model)
  switch(model,
    full      = object$par,
    mean      = object$coefficients$mean,
    precision = object$coefficients$precision
  )
}


# -- Variance-covariance matrix --------------------------------------------- #

#' Variance-covariance matrix of estimated coefficients
#'
#' @param object A fitted \code{"betaregscale"} object.
#' @param model  Character: which component (\code{"full"},
#'   \code{"mean"}, or \code{"precision"}).
#' @param ... Ignored.
#'
#' @return A square numeric matrix.
#'
#' @seealso \code{\link{brs}}, \code{\link{coef.brs}}, \code{\link{confint.brs}}
#'
#' @examples
#' \donttest{
#' dat <- data.frame(
#'   y = c(
#'     0, 5, 20, 50, 75, 90, 100, 30, 60, 45,
#'     10, 40, 55, 70, 85, 25, 35, 65, 80, 15
#'   ),
#'   x1 = rep(c(1, 2), 10)
#' )
#' prep <- brs_prep(dat, ncuts = 100)
#' fit <- brs(y ~ x1, data = prep)
#' vcov(fit)
#' vcov(fit, model = "mean")
#' }
#'
#' @method vcov brs
#' @importFrom stats vcov
#' @export
vcov.brs <- function(object,
                     model = c("full", "mean", "precision"),
                     ...) {
  .check_class(object)
  model <- match.arg(model)

  V <- tryCatch(
    solve(-object$hessian),
    error = function(e) {
      if (requireNamespace("MASS", quietly = TRUE)) {
        warning(
          "Hessian is computationally singular; returning a generalized inverse. ",
          "Standard errors may be unreliable.",
          call. = FALSE
        )
        MASS::ginv(-object$hessian)
      } else {
        # QUAL-L05: inform user why SEs will be NA
        warning(
          "Hessian is computationally singular and package 'MASS' is not available ",
          "for a generalized inverse. Install 'MASS' for approximate SEs. ",
          "Returning NA variance matrix.",
          call. = FALSE
        )
        matrix(NA_real_, nrow(object$hessian), ncol(object$hessian))
      }
    }
  )
  rownames(V) <- colnames(V) <- names(object$par)

  switch(model,
    full = V,
    mean = {
      idx <- seq_len(object$p)
      V[idx, idx, drop = FALSE]
    },
    precision = {
      idx <- object$p + seq_len(object$q)
      V[idx, idx, drop = FALSE]
    }
  )
}


# -- Log-likelihood --------------------------------------------------------- #

#' Extract log-likelihood
#'
#' @param object A fitted \code{"betaregscale"} object.
#' @param ... Ignored.
#'
#' @return An object of class \code{"logLik"} with attributes
#'   \code{df} (number of estimated parameters) and \code{nobs}
#'   (number of observations).
#'
#' @seealso \code{\link{brs}}, \code{\link{AIC.brs}}, \code{\link{BIC.brs}},
#'   \code{\link{brs_gof}}
#'
#' @examples
#' \donttest{
#' dat <- data.frame(
#'   y = c(
#'     0, 5, 20, 50, 75, 90, 100, 30, 60, 45,
#'     10, 40, 55, 70, 85, 25, 35, 65, 80, 15
#'   ),
#'   x1 = rep(c(1, 2), 10)
#' )
#' prep <- brs_prep(dat, ncuts = 100)
#' fit <- brs(y ~ x1, data = prep)
#' logLik(fit)
#' }
#'
#' @method logLik brs
#' @importFrom stats logLik
#' @export
logLik.brs <- function(object, ...) {
  .check_class(object)
  val <- object$value
  attr(val, "df") <- object$npar
  attr(val, "nobs") <- object$nobs
  class(val) <- "logLik"
  val
}


# -- AIC -------------------------------------------------------------------- #

#' Akaike information criterion
#'
#' @param object A fitted \code{"betaregscale"} object.
#' @param ... Ignored.
#' @param k    Penalty per parameter (default 2).
#'
#' @return Scalar AIC value.
#'
#' @seealso \code{\link{brs}}, \code{\link{logLik.brs}}, \code{\link{BIC.brs}},
#'   \code{\link{brs_gof}}
#'
#' @examples
#' \donttest{
#' dat <- data.frame(
#'   y = c(
#'     0, 5, 20, 50, 75, 90, 100, 30, 60, 45,
#'     10, 40, 55, 70, 85, 25, 35, 65, 80, 15
#'   ),
#'   x1 = rep(c(1, 2), 10)
#' )
#' prep <- brs_prep(dat, ncuts = 100)
#' fit <- brs(y ~ x1, data = prep)
#' AIC(fit)
#' }
#'
#' @method AIC brs
#' @importFrom stats AIC
#' @export
AIC.brs <- function(object, ..., k = 2) {
  .check_class(object)
  k * object$npar - 2 * object$value
}


# -- BIC -------------------------------------------------------------------- #

#' Bayesian information criterion
#'
#' @param object A fitted \code{"betaregscale"} object.
#' @param ... Ignored.
#'
#' @return Scalar BIC value.
#'
#' @seealso \code{\link{brs}}, \code{\link{logLik.brs}}, \code{\link{AIC.brs}},
#'   \code{\link{brs_gof}}
#'
#' @examples
#' \donttest{
#' dat <- data.frame(
#'   y = c(
#'     0, 5, 20, 50, 75, 90, 100, 30, 60, 45,
#'     10, 40, 55, 70, 85, 25, 35, 65, 80, 15
#'   ),
#'   x1 = rep(c(1, 2), 10)
#' )
#' prep <- brs_prep(dat, ncuts = 100)
#' fit <- brs(y ~ x1, data = prep)
#' BIC(fit)
#' }
#'
#' @method BIC brs
#' @importFrom stats BIC
#' @export
BIC.brs <- function(object, ...) {
  .check_class(object)
  log(object$nobs) * object$npar - 2 * object$value
}


# -- nobs ------------------------------------------------------------------- #

#' Number of observations
#'
#' @param object A fitted \code{"betaregscale"} object.
#' @param ... Ignored.
#'
#' @return Integer: number of observations.
#'
#' @seealso \code{\link{brs}}, \code{\link{fitted.brs}}, \code{\link{brs_gof}}
#'
#' @examples
#' \donttest{
#' dat <- data.frame(
#'   y = c(
#'     0, 5, 20, 50, 75, 90, 100, 30, 60, 45,
#'     10, 40, 55, 70, 85, 25, 35, 65, 80, 15
#'   ),
#'   x1 = rep(c(1, 2), 10)
#' )
#' prep <- brs_prep(dat, ncuts = 100)
#' fit <- brs(y ~ x1, data = prep)
#' nobs(fit)
#' }
#'
#' @method nobs brs
#' @importFrom stats nobs
#' @export
nobs.brs <- function(object, ...) {
  .check_class(object)
  object$nobs
}


# -- formula ---------------------------------------------------------------- #

#' Extract model formula
#'
#' @param x A fitted \code{"betaregscale"} object.
#' @param ... Ignored.
#'
#' @return The formula used to fit the model.
#'
#' @seealso \code{\link{brs}}, \code{\link{model.matrix.brs}},
#'   \code{\link{coef.brs}}
#'
#' @examples
#' \donttest{
#' dat <- data.frame(
#'   y = c(
#'     0, 5, 20, 50, 75, 90, 100, 30, 60, 45,
#'     10, 40, 55, 70, 85, 25, 35, 65, 80, 15
#'   ),
#'   x1 = rep(c(1, 2), 10)
#' )
#' prep <- brs_prep(dat, ncuts = 100)
#' fit <- brs(y ~ x1, data = prep)
#' formula(fit)
#' }
#'
#' @method formula brs
#' @importFrom stats formula
#' @export
formula.brs <- function(x, ...) {
  .check_class(x)
  x$formula
}


# -- model.matrix ----------------------------------------------------------- #

#' Extract design matrix
#'
#' @param object A fitted \code{"betaregscale"} object.
#' @param model  Character: \code{"mean"} (default) or
#'   \code{"precision"}.
#' @param ... Ignored.
#'
#' @return The design matrix for the specified submodel.
#'
#' @seealso \code{\link{brs}}, \code{\link{formula.brs}},
#'   \code{\link{coef.brs}}
#'
#' @examples
#' \donttest{
#' dat <- data.frame(
#'   y = c(
#'     0, 5, 20, 50, 75, 90, 100, 30, 60, 45,
#'     10, 40, 55, 70, 85, 25, 35, 65, 80, 15
#'   ),
#'   x1 = rep(c(1, 2), 10)
#' )
#' prep <- brs_prep(dat, ncuts = 100)
#' fit <- brs(y ~ x1, data = prep)
#' head(model.matrix(fit))
#' head(model.matrix(fit, model = "precision"))
#' }
#'
#' @method model.matrix brs
#' @importFrom stats model.matrix
#' @export
model.matrix.brs <- function(object,
                             model = c("mean", "precision"),
                             ...) {
  .check_class(object)
  model <- match.arg(model)
  # In the future, if we rebuild the matrix from data, we could pass ... to model.matrix
  # Currently object contains the matrix, so ... is truly not used here.
  # But we leave it available for consistency.
  switch(model,
    mean = object$model_matrices$X,
    precision = {
      if (!is.null(object$model_matrices$Z)) {
        object$model_matrices$Z
      } else {
        matrix(1,
          nrow = object$nobs, ncol = 1,
          dimnames = list(NULL, "(Intercept)")
        )
      }
    }
  )
}


# -- Summary ---------------------------------------------------------------- #

#' Summarize a fitted model (betareg style)
#'
#' @param object A fitted \code{"betaregscale"} object.
#' @param ...    Ignored.
#'
#' @return A list of class \code{"summary.betaregscale"}.
#'
#' @seealso \code{\link{brs}}, \code{\link{print.summary.brs}},
#'   \code{\link{brs_est}}, \code{\link{brs_gof}}
#'
#' @examples
#' \donttest{
#' dat <- data.frame(
#'   y = c(
#'     0, 5, 20, 50, 75, 90, 100, 30, 60, 45,
#'     10, 40, 55, 70, 85, 25, 35, 65, 80, 15
#'   ),
#'   x1 = rep(c(1, 2), 10)
#' )
#' prep <- brs_prep(dat, ncuts = 100)
#' fit <- brs(y ~ x1, data = prep)
#' s <- summary(fit)
#' s$coefficients$mean
#' }
#'
#' @method summary brs
#' @importFrom stats pnorm
#' @export
summary.brs <- function(object, ...) {
  .check_class(object)

  V <- vcov(object, model = "full")

  # Mean coefficients table
  cf_mu <- object$coefficients$mean
  se_mu <- sqrt(pmax(diag(V)[seq_len(object$p)], 0))
  z_mu <- cf_mu / se_mu
  p_mu <- 2 * stats::pnorm(-abs(z_mu))
  tab_mu <- cbind(
    Estimate = cf_mu,
    `Std. Error` = se_mu,
    `z value` = z_mu,
    `Pr(>|z|)` = p_mu
  )

  # Precision coefficients table
  cf_phi <- object$coefficients$precision
  idx_phi <- object$p + seq_len(object$q)
  se_phi <- sqrt(pmax(diag(V)[idx_phi], 0))
  z_phi <- cf_phi / se_phi
  p_phi <- 2 * stats::pnorm(-abs(z_phi))
  tab_phi <- cbind(
    Estimate = cf_phi,
    `Std. Error` = se_phi,
    `z value` = z_phi,
    `Pr(>|z|)` = p_phi
  )

  # Default residuals (RQR)
  rqr <- tryCatch(
    residuals(object, type = "rqr"),
    error = function(e) object$residuals
  )

  # Censoring summary
  delta <- object$delta
  cens_counts <- c(
    exact    = sum(delta == 0L),
    left     = sum(delta == 1L),
    right    = sum(delta == 2L),
    interval = sum(delta == 3L)
  )

  n_exact   <- cens_counts[["exact"]]
  n_total   <- sum(cens_counts)
  pct_exact <- if (n_total > 0L) n_exact / n_total else 1.0

  out <- list(
    call         = object$call,
    coefficients = list(mean = tab_mu, precision = tab_phi),
    residuals    = rqr,
    loglik       = object$value,
    AIC          = AIC(object),
    BIC          = BIC(object),
    df           = object$npar,
    nobs         = object$nobs,
    pseudo.r2    = object$pseudo.r.squared,
    pct_exact    = pct_exact,
    link         = object$link,
    link_phi     = object$link_phi,
    convergence  = object$convergence,
    iterations   = object$iterations,
    method       = object$optim_method,
    censoring    = cens_counts,
    repar        = object$repar
  )
  class(out) <- "summary.brs"
  out
}


#' Print a model summary (betareg style)
#'
#' @param x A \code{"summary.betaregscale"} object.
#' @param digits Number of digits.
#' @param ... Passed to \code{printCoefmat}.
#'
#' @return Invisibly returns the input object \code{x}. The function is called
#'   for its side effect of printing a comprehensive summary to the console,
#'   including the model call, quantile residuals, coefficient tables for mean
#'   and precision submodels with significance stars, goodness-of-fit statistics
#'   (log-likelihood, pseudo R-squared), optimization details, and censoring
#'   information.
#'
#' @seealso \code{\link{summary.brs}}, \code{\link{brs}},
#'   \code{\link{print.brs}}
#'
#' @examples
#' \donttest{
#' dat <- data.frame(
#'   y = c(
#'     0, 5, 20, 50, 75, 90, 100, 30, 60, 45,
#'     10, 40, 55, 70, 85, 25, 35, 65, 80, 15
#'   ),
#'   x1 = rep(c(1, 2), 10)
#' )
#' prep <- brs_prep(dat, ncuts = 100)
#' fit <- brs(y ~ x1, data = prep)
#' print(summary(fit))
#' }
#'
#' @method print summary.brs
#' @importFrom stats quantile printCoefmat
#' @export
print.summary.brs <- function(x,
                              digits = max(3, getOption("digits") - 3),
                              ...) {
  cat("\nCall:\n")
  print(x$call)
  cat("\n")

  # Quantile residuals summary
  rq <- quantile(x$residuals,
    probs = c(0, 0.25, 0.5, 0.75, 1),
    na.rm = TRUE
  )
  names(rq) <- c("Min", "1Q", "Median", "3Q", "Max")
  cat("Quantile residuals:\n")
  print(round(rq, digits))
  cat("\n")

  # Mean model
  cat("Coefficients (mean model with", x$link, "link):\n")
  stats::printCoefmat(x$coefficients$mean,
    digits = digits,
    P.values = TRUE, has.Pvalue = TRUE,
    signif.stars = TRUE,
    ...
  )
  cat("\n")

  # Precision model
  cat(paste0("Phi coefficients (precision model with ", x$link_phi, " link):\n"))
  tab_phi_display <- x$coefficients$precision
  rownames(tab_phi_display) <- .pretty_phi_names(rownames(tab_phi_display))
  stats::printCoefmat(tab_phi_display,
    digits = digits,
    P.values = TRUE, has.Pvalue = TRUE,
    signif.stars = TRUE,
    ...
  )

  cat("---\n")

  # Goodness-of-fit
  cat(
    "Log-likelihood:", formatC(x$loglik, format = "f", digits = 4),
    "on", x$df, "Df | AIC:", formatC(x$AIC, format = "f", digits = 4),
    "| BIC:", formatC(x$BIC, format = "f", digits = 4), "\n"
  )
  r2_note <- if (!is.null(x$pct_exact) && x$pct_exact < 0.5)
    " (midpoint approx.; interpret with caution for heavily censored data)" else ""
  cat("Pseudo R-squared:", formatC(x$pseudo.r2, format = "f", digits = 4),
      r2_note, "\n")
  cat(
    "Number of iterations:",
    if (!is.null(x$iterations)) x$iterations["function"] else "NA",
    paste0("(", x$method, ")"), "\n"
  )

  # Censoring info
  cc <- x$censoring
  parts <- character(0)
  if (cc["interval"] > 0) parts <- c(parts, paste(cc["interval"], "interval"))
  if (cc["left"] > 0) parts <- c(parts, paste(cc["left"], "left"))
  if (cc["right"] > 0) parts <- c(parts, paste(cc["right"], "right"))
  if (cc["exact"] > 0) parts <- c(parts, paste(cc["exact"], "exact"))
  if (length(parts) > 0) {
    cat("Censoring:", paste(parts, collapse = " | "), "\n")
  }

  cat("\n")
  invisible(x)
}


# -- Print ------------------------------------------------------------------ #

#' Print a fitted model (brief betareg style)
#'
#' @param x      A fitted \code{"betaregscale"} object.
#' @param digits Number of significant digits.
#' @param ... Included for consistency with generic methods. Currently
#'   passed to internal methods where applicable.
#'
#' @return Invisibly returns the input object \code{x}. The function is called
#'   for its side effect of printing a formatted summary of the fitted model
#'   to the console, including the model call, mean coefficients (with link
#'   function), and precision coefficients (with link function).
#'
#' @seealso \code{\link{summary.brs}}, \code{\link{print.summary.brs}},
#'   \code{\link{brs}}
#'
#' @examples
#' \donttest{
#' dat <- data.frame(
#'   y = c(
#'     0, 5, 20, 50, 75, 90, 100, 30, 60, 45,
#'     10, 40, 55, 70, 85, 25, 35, 65, 80, 15
#'   ),
#'   x1 = rep(c(1, 2), 10)
#' )
#' prep <- brs_prep(dat, ncuts = 100)
#' fit <- brs(y ~ x1, data = prep)
#' print(fit)
#' }
#'
#' @method print brs
#' @export
print.brs <- function(x,
                      digits = max(3, getOption("digits") - 3),
                      ...) {
  cat("\nCall:\n")
  print(x$call)
  cat("\n")

  cat("Coefficients (mean model with", x$link, "link):\n")
  print(round(x$coefficients$mean, digits))
  cat("\n")

  cat("Phi coefficients (precision model with", x$link_phi, "link):\n")
  prec_display <- x$coefficients$precision
  names(prec_display) <- .pretty_phi_names(names(prec_display))
  print(round(prec_display, digits))
  cat("\n")

  invisible(x)
}


# -- Fitted values ---------------------------------------------------------- #

#' Extract fitted values
#'
#' @param object A fitted \code{"betaregscale"} object.
#' @param type   Character: \code{"mu"} (default) or \code{"phi"}.
#' @param ...    Currently ignored.
#'
#' @return Numeric vector of fitted values.
#'
#' @seealso \code{\link{brs}}, \code{\link{residuals.brs}},
#'   \code{\link{predict.brs}}
#'
#' @examples
#' \donttest{
#' dat <- data.frame(
#'   y = c(
#'     0, 5, 20, 50, 75, 90, 100, 30, 60, 45,
#'     10, 40, 55, 70, 85, 25, 35, 65, 80, 15
#'   ),
#'   x1 = rep(c(1, 2), 10)
#' )
#' prep <- brs_prep(dat, ncuts = 100)
#' fit <- brs(y ~ x1, data = prep)
#' head(fitted(fit))
#' head(fitted(fit, type = "phi"))
#' }
#'
#' @method fitted brs
#' @importFrom stats fitted
#' @export
fitted.brs <- function(object, type = c("mu", "phi"), ...) {
  .check_class(object)
  type <- match.arg(type)
  if (type == "mu") {
    return(object$hatmu)
  }
  n <- length(object$hatmu)
  if (length(object$hatphi) == 1L) {
    rep(as.numeric(object$hatphi), n)
  } else {
    as.numeric(object$hatphi)
  }
}


# -- Residuals -------------------------------------------------------------- #

#' Extract residuals
#'
#' @param object A fitted \code{"betaregscale"} object.
#' @param type   Residual type. One of \code{"response"} (default),
#'   \code{"pearson"}, \code{"deviance"}, \code{"rqr"} (randomized
#'   quantile), \code{"weighted"}, or \code{"sweighted"}.
#' @param ...    Currently ignored.
#'
#' @return Numeric vector of residuals.
#'
#' @details
#' For Pearson residuals the variance formula depends on the
#' reparameterization stored in \code{object$repar}:
#' \describe{
#'   \item{repar = 1 (precision)}{V = mu(1 - mu) / (1 + phi)}
#'   \item{repar = 2 (mean-variance)}{V = mu(1 - mu) * phi}
#' }
#' The weighted and sweighted residuals use the digamma/trigamma
#' formulation from the precision parameterization (repar = 1),
#' so internal conversion is applied when \code{repar != 1}.
#'
#' @seealso \code{\link{brs}}, \code{\link{fitted.brs}}, \code{\link{plot.brs}}
#'
#' @examples
#' \donttest{
#' dat <- data.frame(
#'   y = c(
#'     0, 5, 20, 50, 75, 90, 100, 30, 60, 45,
#'     10, 40, 55, 70, 85, 25, 35, 65, 80, 15
#'   ),
#'   x1 = rep(c(1, 2), 10)
#' )
#' prep <- brs_prep(dat, ncuts = 100)
#' fit <- brs(y ~ x1, data = prep)
#' head(residuals(fit))
#' head(residuals(fit, type = "pearson"))
#' }
#'
#' @method residuals brs
#' @importFrom stats residuals qnorm pbeta dbeta qlogis
#' @export
residuals.brs <- function(object,
                          type = c(
                            "response", "pearson",
                            "deviance", "rqr",
                            "weighted", "sweighted"
                          ),
                          ...) {
  .check_class(object)
  type <- match.arg(type)

  y <- object$Y[, "yt"]
  mu <- object$hatmu
  phi <- object$hatphi
  repar <- object$repar

  if (type == "response") {
    return(object$residuals)
  }

  # Helper: get shape parameters (a, b) from (mu, phi) per repar
  get_shapes <- function(mu, phi, repar) {
    rp <- brs_repar(mu, phi, repar = repar)
    list(a = rp$shape1, b = rp$shape2)
  }

  # Helper: convert to precision scale (repar=1) phi
  to_precision <- function(mu, phi, repar) {
    if (repar == 1L) {
      return(phi)
    }
    if (repar == 2L) {
      return((1 - phi) / phi)
    }
    # repar == 0: shapes are (mu, phi) directly, precision = a + b = mu + phi
    mu + phi
  }

  switch(type,
    pearson = {
      # Variance depends on reparameterization
      if (repar == 1L) {
        # phi is precision: V[Y] = mu(1-mu)/(1+phi)
        v <- mu * (1 - mu) / (1 + phi)
      } else if (repar == 2L) {
        # phi is dispersion in (0,1): V[Y] = mu(1-mu)*phi
        v <- mu * (1 - mu) * phi
      } else {
        # repar = 0: direct shape parameters, compute variance
        sh <- get_shapes(mu, phi, repar)
        s <- sh$a + sh$b
        v <- (sh$a * sh$b) / (s^2 * (s + 1))
      }
      (y - mu) / sqrt(v)
    },
    deviance = {
      sh <- get_shapes(mu, phi, repar)
      # BUG-C02: deviance uses the SATURATED log-likelihood where mu_sat = y.
      # The fitted contribution evaluates density at y with fitted shapes (a,b).
      # The saturated contribution uses shapes derived from y itself.
      y_safe <- pmin(pmax(y, 1e-7), 1 - 1e-7)
      sh_sat  <- brs_repar(y_safe, phi, repar = repar)
      ll_sat  <- stats::dbeta(y_safe, sh_sat$shape1, sh_sat$shape2, log = TRUE)
      ll_fit  <- stats::dbeta(y_safe, sh$a, sh$b, log = TRUE)
      sign(y - mu) * sqrt(abs(2 * (ll_sat - ll_fit)))
    },
    rqr = {
      sh <- get_shapes(mu, phi, repar)
      left <- object$Y[, "left"]
      right <- object$Y[, "right"]
      delta <- as.integer(object$Y[, "delta"])

      f_left <- stats::pbeta(left, sh$a, sh$b)
      f_right <- stats::pbeta(right, sh$a, sh$b)
      f_y <- stats::pbeta(y, sh$a, sh$b)

      # Randomized PIT to respect censoring intervals.
      lo <- ifelse(delta == 0L, f_y,
        ifelse(delta == 1L, 0, f_left)
      )
      hi <- ifelse(delta == 0L, f_y,
        ifelse(delta == 2L, 1, f_right)
      )
      hi <- pmax(hi, lo)

      u <- stats::runif(length(lo), min = lo, max = hi)
      u <- pmin(pmax(u, 1e-10), 1 - 1e-10)
      stats::qnorm(u)
    },
    weighted = {
      prec <- to_precision(mu, phi, repar)
      ystar <- stats::qlogis(y)
      mustar <- digamma(mu * prec) - digamma((1 - mu) * prec)
      v <- trigamma(mu * prec) + trigamma((1 - mu) * prec)
      (ystar - mustar) / sqrt(prec * v)
    },
    sweighted = {
      prec <- to_precision(mu, phi, repar)
      ystar <- stats::qlogis(y)
      mustar <- digamma(mu * prec) - digamma((1 - mu) * prec)
      v <- trigamma(mu * prec) + trigamma((1 - mu) * prec)
      (ystar - mustar) / sqrt(v)
    }
  )
}


# -- Confidence intervals --------------------------------------------------- #

#' Wald confidence intervals
#'
#' @description
#' Computes Wald confidence intervals for model parameters using the
#' normal approximation.
#'
#' @param object A fitted \code{"betaregscale"} object.
#' @param parm   Character or integer: which parameters. If missing,
#'   all parameters are returned.
#' @param level  Confidence level (default 0.95).
#' @param model  Character: \code{"full"}, \code{"mean"}, or
#'   \code{"precision"}.
#' @param ...    Currently ignored.
#'
#' @return Matrix with columns for lower and upper confidence bounds.
#'
#' @seealso \code{\link{brs}}, \code{\link{coef.brs}}, \code{\link{vcov.brs}},
#'   \code{\link{brs_est}}
#'
#' @examples
#' \donttest{
#' dat <- data.frame(
#'   y = c(
#'     0, 5, 20, 50, 75, 90, 100, 30, 60, 45,
#'     10, 40, 55, 70, 85, 25, 35, 65, 80, 15
#'   ),
#'   x1 = rep(c(1, 2), 10)
#' )
#' prep <- brs_prep(dat, ncuts = 100)
#' fit <- brs(y ~ x1, data = prep)
#' confint(fit)
#' confint(fit, model = "mean")
#' }
#'
#' @method confint brs
#' @importFrom stats confint qnorm
#' @export
confint.brs <- function(object, parm, level = 0.95,
                        model = c("full", "mean", "precision"),
                        ...) {
  .check_class(object)
  model <- match.arg(model)

  cf <- coef(object, model = model)
  se <- sqrt(pmax(diag(vcov(object, model = model)), 0))
  z <- stats::qnorm(1 - (1 - level) / 2)

  ci <- cbind(cf - z * se, cf + z * se)
  colnames(ci) <- paste0(
    format(100 * c((1 - level) / 2, 1 - (1 - level) / 2), digits = 3),
    " %"
  )

  if (!missing(parm)) {
    ci <- ci[parm, , drop = FALSE]
  }

  ci
}


# -- Predict ---------------------------------------------------------------- #

#' Predict from a fitted model
#'
#' @param object  A fitted \code{"betaregscale"} object.
#' @param newdata Optional data frame for prediction.
#' @param type    Prediction type: \code{"response"} (default),
#'   \code{"link"}, \code{"precision"}, \code{"variance"}, or
#'   \code{"quantile"}.
#' @param at      Numeric vector of probabilities for quantile
#'   predictions (default 0.5).
#' @param ...     Currently ignored.
#'
#' @return Numeric vector or matrix.
#'
#' @seealso \code{\link{brs}}, \code{\link{fitted.brs}},
#'   \code{\link{brs_predict_scoreprob}}
#'
#' @examples
#' \donttest{
#' dat <- data.frame(
#'   y = c(
#'     0, 5, 20, 50, 75, 90, 100, 30, 60, 45,
#'     10, 40, 55, 70, 85, 25, 35, 65, 80, 15
#'   ),
#'   x1 = rep(c(1, 2), 10)
#' )
#' prep <- brs_prep(dat, ncuts = 100)
#' fit <- brs(y ~ x1, data = prep)
#' head(predict(fit))
#' head(predict(fit, type = "precision"))
#' newdat <- data.frame(x1 = c(1, 2))
#' predict(fit, newdata = newdat)
#' }
#'
#' @method predict brs
#' @importFrom stats predict qbeta model.matrix terms model.frame
#' @export
predict.brs <- function(object, newdata = NULL,
                        type = c(
                          "response", "link",
                          "precision", "variance",
                          "quantile"
                        ),
                        at = 0.5, ...) {
  .check_class(object)
  type <- match.arg(type)

  if (is.null(newdata)) {
    mu <- object$hatmu
    phi <- if (length(object$hatphi) == 1L) {
      rep(as.numeric(object$hatphi), length(mu))
    } else {
      as.numeric(object$hatphi)
    }
    eta_mu <- apply_link(pmin(pmax(mu, 1e-7), 1 - 1e-7), object$link)
  } else {
    # Build X from newdata
    mt_mu <- stats::delete.response(object$terms$mean)
    mf <- stats::model.frame(mt_mu, data = newdata, ...)
    X <- stats::model.matrix(mt_mu, mf)
    eta_mu <- as.numeric(X %*% object$coefficients$mean)
    mu <- apply_inv_link(eta_mu, object$link)

    # BUG-H02: detect variable-dispersion by presence of Z matrix with
    # non-intercept columns, not by q > 1 (which misclassifies y ~ x | 1).
    has_var_phi <- !is.null(object$model_matrices$Z) &&
                   !is.null(object$terms$precision) &&
                   length(attr(object$terms$precision, "term.labels")) > 0L
    # Build Z from newdata (variable dispersion)
    if (has_var_phi) {
      mt_phi <- object$terms$precision
      mf_z <- stats::model.frame(mt_phi, data = newdata, ...)
      Z <- stats::model.matrix(mt_phi, mf_z)
      eta_phi <- as.numeric(Z %*% object$coefficients$precision)
      phi <- apply_inv_link(eta_phi, object$link_phi)
    } else {
      phi_scalar <- apply_inv_link(
        as.numeric(object$coefficients$precision),
        object$link_phi
      )
      phi <- rep(phi_scalar, nrow(X))
    }
  }

  switch(type,
    response = mu,
    link = eta_mu,
    precision = phi,
    variance = {
      repar <- object$repar
      if (repar == 1L) {
        mu * (1 - mu) / (1 + phi)
      } else if (repar == 2L) {
        mu * (1 - mu) * phi
      } else {
        sh <- brs_repar(mu, phi, repar = repar)
        s <- sh$shape1 + sh$shape2
        (sh$shape1 * sh$shape2) / (s^2 * (s + 1))
      }
    },
    quantile = {
      rp <- brs_repar(mu, phi, repar = object$repar)
      rval <- sapply(at, function(p) {
        stats::qbeta(p, rp$shape1, rp$shape2)
      })
      if (length(at) > 1L) {
        if (NCOL(rval) == 1L) {
          rval <- matrix(rval,
            ncol = length(at),
            dimnames = list(NULL, paste0("q_", at))
          )
        } else {
          colnames(rval) <- paste0("q_", at)
        }
      } else {
        rval <- drop(rval)
      }
      rval
    }
  )
}


# -- Convenience extractors ------------------------------------------------ #

#' Goodness-of-fit measures
#'
#' @param object A fitted \code{"brs"} or \code{"brsmm"} object.
#'
#' @return Data frame with logLik, AIC, BIC, and pseudo-R-squared.
#'
#' @seealso \code{\link{brs}}, \code{\link{brs_est}}, \code{\link{brs_hessian}}
#'
#' @examples
#' \donttest{
#' dat <- data.frame(
#'   y = c(
#'     0, 5, 20, 50, 75, 90, 100, 30, 60, 45,
#'     10, 40, 55, 70, 85, 25, 35, 65, 80, 15
#'   ),
#'   x1 = rep(c(1, 2), 10)
#' )
#' prep <- brs_prep(dat, ncuts = 100)
#' fit <- brs(y ~ x1, data = prep)
#' brs_gof(fit)
#' }
#' @references
#' Lopes, J. E. (2023). \emph{Modelos de regressao beta para dados de escala}.
#' Master's dissertation, Universidade Federal do Parana, Curitiba.
#' URI: \url{https://hdl.handle.net/1884/86624}.
#'
#' Hawker, G. A., Mian, S., Kendzerska, T., and French, M. (2011).
#' Measures of adult pain: Visual Analog Scale for Pain (VAS Pain),
#' Numeric Rating Scale for Pain (NRS Pain), McGill Pain Questionnaire (MPQ),
#' Short-Form McGill Pain Questionnaire (SF-MPQ), Chronic Pain Grade Scale
#' (CPGS), Short Form-36 Bodily Pain Scale (SF-36 BPS), and Measure of
#' Intermittent and Constant Osteoarthritis Pain (ICOAP).
#' Arthritis Care and Research, 63(S11), S240-S252.
#' \doi{10.1002/acr.20543}
#'
#' Hjermstad, M. J., Fayers, P. M., Haugen, D. F., et al. (2011).
#' Studies comparing Numerical Rating Scales, Verbal Rating Scales, and
#' Visual Analogue Scales for assessment of pain intensity in adults:
#' a systematic literature review.
#' Journal of Pain and Symptom Management, 41(6), 1073-1093.
#' \doi{10.1016/j.jpainsymman.2010.08.016}
#' @rdname brs_gof
#' @export
brs_gof <- function(object) {
  if (!inherits(object, c("brs", "brsmm"))) {
    stop("Expected a 'brs' or 'brsmm' object.", call. = FALSE)
  }
  data.frame(
    logLik    = as.numeric(logLik(object)),
    AIC       = AIC(object),
    BIC       = BIC(object),
    pseudo_r2 = object$pseudo.r.squared
  )
}

#' Coefficient estimates with inference
#'
#' @param object A fitted \code{"betaregscale"} object.
#' @param alpha  Significance level (default 0.05).
#'
#' @return Data frame of estimates, standard errors, z-values, and
#'   p-values.
#'
#' @seealso \code{\link{brs}}, \code{\link{brs_gof}}, \code{\link{brs_hessian}},
#'   \code{\link{summary.brs}}
#'
#' @examples
#' \donttest{
#' dat <- data.frame(
#'   y = c(
#'     0, 5, 20, 50, 75, 90, 100, 30, 60, 45,
#'     10, 40, 55, 70, 85, 25, 35, 65, 80, 15
#'   ),
#'   x1 = rep(c(1, 2), 10)
#' )
#' prep <- brs_prep(dat, ncuts = 100)
#' fit <- brs(y ~ x1, data = prep)
#' brs_est(fit)
#' }
#'
#' @references
#' Lopes, J. E. (2023). \emph{Modelos de regressao beta para dados de escala}.
#' Master's dissertation, Universidade Federal do Parana, Curitiba.
#' URI: \url{https://hdl.handle.net/1884/86624}.
#'
#' Ferrari, S. L. P., and Cribari-Neto, F. (2004).
#' Beta regression for modelling rates and proportions.
#' \emph{Journal of Applied Statistics}, \bold{31}(7), 799--815.
#' \doi{10.1080/0266476042000214501}
#'
#' Hawker, G. A., Mian, S., Kendzerska, T., and French, M. (2011).
#' Measures of adult pain: Visual Analog Scale for Pain (VAS Pain),
#' Numeric Rating Scale for Pain (NRS Pain), McGill Pain Questionnaire (MPQ),
#' Short-Form McGill Pain Questionnaire (SF-MPQ), Chronic Pain Grade Scale
#' (CPGS), Short Form-36 Bodily Pain Scale (SF-36 BPS), and Measure of
#' Intermittent and Constant Osteoarthritis Pain (ICOAP).
#' Arthritis Care and Research, 63(S11), S240-S252.
#' \doi{10.1002/acr.20543}
#'
#' Hjermstad, M. J., Fayers, P. M., Haugen, D. F., et al. (2011).
#' Studies comparing Numerical Rating Scales, Verbal Rating Scales, and
#' Visual Analogue Scales for assessment of pain intensity in adults:
#' a systematic literature review.
#' Journal of Pain and Symptom Management, 41(6), 1073-1093.
#' \doi{10.1016/j.jpainsymman.2010.08.016}
#'
#' @importFrom stats pnorm
#' @rdname brs_est
#' @export
brs_est <- function(object, alpha = 0.05) {
  if (!inherits(object, c("brs", "brsmm"))) {
    stop("Expected a 'brs' or 'brsmm' object.", call. = FALSE)
  }
  V <- vcov(object)
  se <- sqrt(pmax(diag(V), 0))
  z <- object$par / se
  p <- 2 * stats::pnorm(-abs(z))
  z_alpha <- stats::qnorm(1 - alpha / 2)

  data.frame(
    variable  = names(object$par),
    estimate  = unname(object$par),
    se        = unname(se),
    z_value   = unname(z),
    p_value   = unname(p),
    ci_lower  = unname(object$par - z_alpha * se),
    ci_upper  = unname(object$par + z_alpha * se),
    row.names = NULL
  )
}

#' Internal coefficient table (deprecated, use brs_est() or summary())
#'
#' @description
#' Deprecated convenience wrapper. Use \code{\link{brs_est}} for coefficient
#' estimates or \code{\link{summary.brs}} for a full model summary.
#'
#' @param fit   A fitted \code{"brs"} object.
#' @param alpha Significance level.
#'
#' @return A list with components \code{est} (from \code{\link{brs_est}})
#'   and \code{gof} (from \code{\link{brs_gof}}).
#'
#' @seealso \code{\link{brs_est}}, \code{\link{brs_gof}}, \code{\link{summary.brs}}
#'
#' @examples
#' \donttest{
#' dat <- data.frame(
#'   y = c(
#'     0, 5, 20, 50, 75, 90, 100, 30, 60, 45,
#'     10, 40, 55, 70, 85, 25, 35, 65, 80, 15
#'   ),
#'   x1 = rep(c(1, 2), 10)
#' )
#' prep <- brs_prep(dat, ncuts = 100)
#' fit <- brs(y ~ x1, data = prep)
#' suppressWarnings(brs_coef(fit))  # deprecated; use brs_est()
#' }
#'
#' @references
#' Lopes, J. E. (2023). \emph{Modelos de regressao beta para dados de escala}.
#' Master's dissertation, Universidade Federal do Parana, Curitiba.
#' URI: \url{https://hdl.handle.net/1884/86624}.
#'
#' Ferrari, S. L. P., and Cribari-Neto, F. (2004).
#' Beta regression for modelling rates and proportions.
#' \emph{Journal of Applied Statistics}, \bold{31}(7), 799--815.
#' \doi{10.1080/0266476042000214501}
#'
#' Hawker, G. A., Mian, S., Kendzerska, T., and French, M. (2011).
#' Measures of adult pain: Visual Analog Scale for Pain (VAS Pain),
#' Numeric Rating Scale for Pain (NRS Pain), McGill Pain Questionnaire (MPQ),
#' Short-Form McGill Pain Questionnaire (SF-MPQ), Chronic Pain Grade Scale
#' (CPGS), Short Form-36 Bodily Pain Scale (SF-36 BPS), and Measure of
#' Intermittent and Constant Osteoarthritis Pain (ICOAP).
#' Arthritis Care and Research, 63(S11), S240-S252.
#' \doi{10.1002/acr.20543}
#'
#' Hjermstad, M. J., Fayers, P. M., Haugen, D. F., et al. (2011).
#' Studies comparing Numerical Rating Scales, Verbal Rating Scales, and
#' Visual Analogue Scales for assessment of pain intensity in adults:
#' a systematic literature review.
#' Journal of Pain and Symptom Management, 41(6), 1073-1093.
#' \doi{10.1016/j.jpainsymman.2010.08.016}
#' @export
brs_coef <- function(fit, alpha = 0.05) {
  .Deprecated("brs_est")
  .check_class(fit)
  list(est = brs_est(fit, alpha = alpha), gof = brs_gof(fit))
}

#' Extract the Hessian matrix
#'
#' @param object A fitted \code{"brs"} or \code{"brsmm"} object.
#'
#' @return Numeric Hessian matrix.
#'
#' @seealso \code{\link{brs}}, \code{\link{vcov.brs}}, \code{\link{brs_est}}
#'
#' @examples
#' \donttest{
#' dat <- data.frame(
#'   y = c(
#'     0, 5, 20, 50, 75, 90, 100, 30, 60, 45,
#'     10, 40, 55, 70, 85, 25, 35, 65, 80, 15
#'   ),
#'   x1 = rep(c(1, 2), 10)
#' )
#' prep <- brs_prep(dat, ncuts = 100)
#' fit <- brs(y ~ x1, data = prep)
#' brs_hessian(fit)
#' }
#'
#' @references
#' Lopes, J. E. (2023). \emph{Modelos de regressao beta para dados de escala}.
#' Master's dissertation, Universidade Federal do Parana, Curitiba.
#' URI: \url{https://hdl.handle.net/1884/86624}.
#'
#' Hawker, G. A., Mian, S., Kendzerska, T., and French, M. (2011).
#' Measures of adult pain: Visual Analog Scale for Pain (VAS Pain),
#' Numeric Rating Scale for Pain (NRS Pain), McGill Pain Questionnaire (MPQ),
#' Short-Form McGill Pain Questionnaire (SF-MPQ), Chronic Pain Grade Scale
#' (CPGS), Short Form-36 Bodily Pain Scale (SF-36 BPS), and Measure of
#' Intermittent and Constant Osteoarthritis Pain (ICOAP).
#' Arthritis Care and Research, 63(S11), S240-S252.
#' \doi{10.1002/acr.20543}
#'
#' Hjermstad, M. J., Fayers, P. M., Haugen, D. F., et al. (2011).
#' Studies comparing Numerical Rating Scales, Verbal Rating Scales, and
#' Visual Analogue Scales for assessment of pain intensity in adults:
#' a systematic literature review.
#' Journal of Pain and Symptom Management, 41(6), 1073-1093.
#' \doi{10.1016/j.jpainsymman.2010.08.016}
#'
#' @rdname brs_hessian
#' @export
brs_hessian <- function(object) {
  .check_class(object)
  object$hessian
}
