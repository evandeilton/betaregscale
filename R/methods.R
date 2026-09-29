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

#' Validate a brs object
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
#' @param object A fitted \code{"brs"} object.
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
#' @param object A fitted \code{"brs"} object.
#' @param model  Character: which component (\code{"full"},
#'   \code{"mean"}, or \code{"precision"}).
#' @param ... Ignored.
#'
#' @details
#' \eqn{(-H)^{-1}} with \eqn{H} the Hessian of the log-likelihood at the
#' estimate. No generalised inverse is used: a singular Hessian gives an
#' \code{NA} matrix, and negative or non-finite variances become \code{NA}
#' (row and column); both cases warn (see \code{fit$diagnostics}).
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

  # No generalised inverse: singular or indefinite Hessians give NA, with a warning
  V <- .brs_vcov(object$hessian, names(object$par))

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
#' @param object A fitted \code{"brs"} object.
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
#' @param object A fitted \code{"brs"} object.
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
#' @param object A fitted \code{"brs"} object.
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
#' @param object A fitted \code{"brs"} object.
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
#' @param x A fitted \code{"brs"} object.
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
#' @param object A fitted \code{"brs"} object.
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

#' Summarize a fitted beta interval model
#'
#' @description
#' Wald tables for the mean and precision (or dispersion) coefficients,
#' information criteria, a pseudo \eqn{R^2}, the censoring counts and a
#' summary of the randomized quantile residuals.
#'
#' @details
#' For each coefficient \eqn{\hat\theta_j}, the standard error is
#' \eqn{SE_j = \sqrt{[(-H)^{-1}]_{jj}}} from \code{\link{vcov.brs}}, the Wald
#' statistic is \eqn{z_j = \hat\theta_j / SE_j} and the two-sided p-value is
#' \eqn{2\Phi(-|z_j|)} (Lopes, 2023, "Inferencia"). The test is on the link
#' scale (\eqn{H_0: \theta_j = 0}). When \eqn{-H} is singular, or its inverse
#' has negative variances, the affected standard errors, statistics and
#' p-values are \code{NA} (see 'Fit diagnostics' in \code{\link{brs}}).
#'
#' \eqn{\mathrm{AIC} = -2\ell + 2k} and \eqn{\mathrm{BIC} = -2\ell + k\log n},
#' with \eqn{\ell} the maximised log-likelihood, \eqn{k} the number of
#' coefficients and \eqn{n} the number of observations. The pseudo
#' \eqn{R^2} is the squared correlation between the fitted linear predictor
#' of the mean and \eqn{g_1(y_i)} at the cell centres \code{yt} (Ferrari and
#' Cribari-Neto, 2004); under \code{repar = 0} both sides are on the logit
#' scale. It uses the cell centres, so it is rough when most observations are
#' censored (the print says so).
#'
#' The randomized quantile residuals (\code{\link{residuals.brs}},
#' \code{type = "rqr"}) are drawn without changing the caller's RNG state.
#'
#' @param object A fitted \code{"brs"} object.
#' @param ... Currently ignored.
#'
#' @return A list of class \code{"summary.brs"} with \code{coefficients}
#'   (tables \code{mean} and \code{precision} with columns \code{Estimate},
#'   \code{Std. Error}, \code{z value}, \code{Pr(>|z|)}), \code{residuals}
#'   (RQR), \code{loglik}, \code{AIC}, \code{BIC}, \code{df}, \code{nobs},
#'   \code{pseudo.r2}, \code{censoring} (counts by type), \code{link},
#'   \code{link_phi}, \code{repar}, \code{convergence} and \code{iterations}.
#'
#' @seealso \code{\link{brs}}, \code{\link{confint.brs}},
#'   \code{\link{anova.brs}}, \code{\link{brs_gof}}
#'
#' @references
#' Lopes, J. E. (2023). \emph{Modelos de regressao beta para dados de escala}.
#' Master's dissertation, Universidade Federal do Parana, Curitiba.
#' URI: https://hdl.handle.net/1884/86624.
#'
#' Ferrari, S. L. P., and Cribari-Neto, F. (2004).
#' Beta regression for modelling rates and proportions.
#' \emph{Journal of Applied Statistics}, \bold{31}(7), 799--815.
#' \doi{10.1080/0266476042000214501}
#'
#' @examples
#' set.seed(2023)
#' d <- data.frame(time = factor(rep(c("6h", "12h", "24h"), each = 60),
#'                               levels = c("6h", "12h", "24h")))
#' shp <- brs_repar(mu = plogis(-1.3 + c(0, 0.75, 0.3)[d$time]), phi = 0.3)
#' d$y <- round(10 * rbeta(nrow(d), shp$shape1, shp$shape2))
#' fit <- brs(y ~ time, data = d, ncuts = 10)
#' s <- summary(fit)
#' s
#' s$coefficients$mean
#' c(AIC = s$AIC, BIC = s$BIC, pseudo_R2 = s$pseudo.r2)
#'
#' @method summary brs
#' @importFrom stats pnorm
#' @export
summary.brs <- function(object, ...) {
  .check_class(object)

  V <- vcov(object, model = "full")

  # Mean coefficients table
  cf_mu <- object$coefficients$mean
  se_mu <- sqrt(diag(V)[seq_len(object$p)])
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
  se_phi <- sqrt(diag(V)[idx_phi])
  z_phi <- cf_phi / se_phi
  p_phi <- 2 * stats::pnorm(-abs(z_phi))
  tab_phi <- cbind(
    Estimate = cf_phi,
    `Std. Error` = se_phi,
    `z value` = z_phi,
    `Pr(>|z|)` = p_phi
  )

  # Default residuals (RQR); their random draws leave the user's RNG state intact
  rqr <- .brs_keep_seed(tryCatch(
    residuals(object, type = "rqr"),
    error = function(e) object$residuals
  ))

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
#' @param x A \code{"summary.brs"} object.
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
#' @param x      A fitted \code{"brs"} object.
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
#' @param object A fitted \code{"brs"} object.
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
  n <- length(object$hatmu)
  phi <- rep_len(as.numeric(object$hatphi), n)
  if (type == "mu") {
    # E[Y] = a / (a + b): hatmu itself under repar 1/2, p / (p + q) under 0
    return(.brs_mean(object$hatmu, phi, object$repar))
  }
  phi
}


# -- Residuals -------------------------------------------------------------- #

#' Residuals of a fitted beta interval model
#'
#' @description
#' Residuals of a \code{"brs"} fit. Randomized quantile residuals
#' (\code{type = "rqr"}) use the censoring of each observation and are the
#' recommended ones; the other types are evaluated at one point of the cell
#' (see Details).
#'
#' @details
#' Let \eqn{(a_i, b_i)} be the fitted shapes, \eqn{\hat\mu_i = a_i/(a_i + b_i)}
#' the fitted mean and \eqn{V_i} the fitted variance (\code{\link{brs_repar}}).
#' All types except \code{"rqr"} use \eqn{y_i =} \code{yt}, the centre of the
#' cell: \eqn{s/K} under \code{interval = "mid"}, so that the border scores
#' sit at \eqn{10^{-5}} and \eqn{1 - 10^{-5}}, and \eqn{(s + 0.5)/(K + 1)}
#' under \code{"right"}/\code{"left"}; exact values are used as they are. This
#' is the midpoint convention of Lopes (2023, "Analise de residuos"); it makes
#' these residuals unreliable at the borders of the scale.
#' \describe{
#'   \item{\code{"response"}}{\eqn{y_i - \hat\mu_i}.}
#'   \item{\code{"pearson"}}{\eqn{(y_i - \hat\mu_i)/\sqrt{V_i}}. Lopes (2023)
#'     writes \eqn{V_i} as \eqn{\mu(1 - \mu)/(1 + \phi)}
#'     (parameterisation 1); it is the same variance under every
#'     \code{repar}.}
#'   \item{\code{"deviance"}}{\eqn{\mathrm{sign}(y_i - \hat\mu_i)
#'     \sqrt{|2\{\ell_i(y_i) - \ell_i(\hat\mu_i)\}|}}, where \eqn{\ell_i(m)} is
#'     the beta log-density at \eqn{y_i} with mean \eqn{m} and precision
#'     \eqn{a_i + b_i} (Ferrari and Cribari-Neto, 2004). The saturated mean
#'     is taken as \eqn{y_i}, as in \pkg{betareg}. For a small precision
#'     (U- or J-shaped densities) \eqn{\ell_i(y_i)} can be below
#'     \eqn{\ell_i(\hat\mu_i)}; the absolute value is then used and only the
#'     sign carries information.}
#'   \item{\code{"rqr"}}{\eqn{\Phi^{-1}(u_i)}, with \eqn{u_i} uniform on
#'     \eqn{(F(l_i), F(u_i))} for \eqn{\delta_i = 3}, on \eqn{(0, F(u_i))} for
#'     \eqn{\delta_i = 1} and on \eqn{(F(l_i), 1)} for \eqn{\delta_i = 2}, and
#'     \eqn{u_i = F(y_i)} for exact values (Dunn and Smyth, 1996); \eqn{u_i}
#'     is kept in \eqn{[10^{-10}, 1 - 10^{-10}]}. They are standard normal
#'     under the model whatever the censoring. They are random: set a seed
#'     to reproduce them; \code{summary()} draws them without changing the
#'     caller's RNG state.}
#'   \item{\code{"weighted"}, \code{"sweighted"}}{\eqn{(y_i^* - \mu_i^*)/
#'     \sqrt{(a_i + b_i) v_i}} and \eqn{(y_i^* - \mu_i^*)/\sqrt{v_i}}, with
#'     \eqn{y_i^* = \mathrm{logit}(y_i)}, \eqn{\mu_i^* = \psi(a_i) - \psi(b_i)}
#'     and \eqn{v_i = \psi'(a_i) + \psi'(b_i)} (Espinheira, Ferrari and
#'     Cribari-Neto, 2008).}
#' }
#' Lopes (2023) also recommends the adjusted quantile residuals of Pereira
#' (2019), which are not implemented.
#'
#' @param object A fitted \code{"brs"} object.
#' @param type Residual type: \code{"response"} (default), \code{"pearson"},
#'   \code{"deviance"}, \code{"rqr"}, \code{"weighted"} or
#'   \code{"sweighted"}.
#' @param ... Currently ignored.
#'
#' @return Numeric vector of residuals, one per observation.
#'
#' @seealso \code{\link{brs}}, \code{\link{fitted.brs}}, \code{\link{plot.brs}}
#'
#' @references
#' Lopes, J. E. (2023). \emph{Modelos de regressao beta para dados de escala}.
#' Master's dissertation, Universidade Federal do Parana, Curitiba.
#' URI: https://hdl.handle.net/1884/86624.
#'
#' Dunn, P. K., and Smyth, G. K. (1996). Randomized quantile residuals.
#' \emph{Journal of Computational and Graphical Statistics}, \bold{5}(3),
#' 236--244.
#'
#' Espinheira, P. L., Ferrari, S. L. P., and Cribari-Neto, F. (2008). On beta
#' regression residuals. \emph{Journal of Applied Statistics}, \bold{35}(4),
#' 407--419.
#'
#' Ferrari, S. L. P., and Cribari-Neto, F. (2004).
#' Beta regression for modelling rates and proportions.
#' \emph{Journal of Applied Statistics}, \bold{31}(7), 799--815.
#' \doi{10.1080/0266476042000214501}
#'
#' Pereira, G. H. A. (2019). On quantile residuals in beta regression.
#' \emph{Communications in Statistics - Simulation and Computation},
#' \bold{48}(1), 302--316.
#'
#' @examples
#' # Synthetic NRS-11 scores: 3 post-operative times. Simulated, not real data.
#' set.seed(2023)
#' nrs <- data.frame(time = factor(rep(c("6h", "12h", "24h"), each = 80),
#'                                 levels = c("6h", "12h", "24h")))
#' shp <- brs_repar(mu = plogis(-1.3 + c(0, 0.75, 0.3)[nrs$time]), phi = 0.3,
#'                  repar = 2)
#' nrs$y <- round(10 * rbeta(nrow(nrs), shp$shape1, shp$shape2))
#' fit <- brs(y ~ time, data = nrs, ncuts = 10)
#'
#' # Randomized quantile residuals: approximately N(0, 1), borders included
#' set.seed(1)
#' r_q <- residuals(fit, type = "rqr")
#' qqnorm(r_q); qqline(r_q)
#'
#' # Midpoint-based residuals are extreme at the border scores 0 and 10
#' r_p <- residuals(fit, type = "pearson")
#' tapply(r_p, cut(nrs$y, c(-1, 0, 9, 10), labels = c("0", "1-9", "10")), mean)
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
  # `mu` is the FIRST parameter (shape p under repar 0); `ey` is E[Y].
  mu <- object$hatmu
  phi <- rep_len(as.numeric(object$hatphi), length(mu))
  repar <- object$repar
  ey <- .brs_mean(mu, phi, repar)

  if (type == "response") {
    return(object$residuals)
  }

  # Helper: get shape parameters (a, b) from (mu, phi) per repar
  get_shapes <- function(mu, phi, repar) {
    rp <- brs_repar(mu, phi, repar = repar)
    list(a = rp$shape1, b = rp$shape2)
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
      (y - ey) / sqrt(v)
    },
    deviance = {
      sh <- get_shapes(mu, phi, repar)
      # BUG-C02: deviance uses the SATURATED log-likelihood where mu_sat = y.
      # The fitted contribution evaluates density at y with fitted shapes (a,b).
      # Saturated model: mean y, same precision a + b (= brs_repar(y, phi) under
      # repar 1/2; defined under repar 0 too).
      y_safe <- pmin(pmax(y, 1e-7), 1 - 1e-7)
      prec <- sh$a + sh$b
      ll_sat  <- stats::dbeta(y_safe, y_safe * prec, (1 - y_safe) * prec, log = TRUE)
      ll_fit  <- stats::dbeta(y_safe, sh$a, sh$b, log = TRUE)
      sign(y - ey) * sqrt(abs(2 * (ll_sat - ll_fit)))
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
    weighted = ,
    sweighted = {
      # Espinheira et al. (2008) from the shapes: a = mu*prec, b = (1-mu)*prec
      sh <- get_shapes(mu, phi, repar)
      prec <- sh$a + sh$b
      ystar <- stats::qlogis(y)
      mustar <- digamma(sh$a) - digamma(sh$b)
      v <- trigamma(sh$a) + trigamma(sh$b)
      if (type == "weighted") {
        (ystar - mustar) / sqrt(prec * v)
      } else {
        (ystar - mustar) / sqrt(v)
      }
    }
  )
}


# -- Confidence intervals --------------------------------------------------- #

#' Wald confidence intervals
#'
#' @description
#' Wald intervals \eqn{\hat\theta_j \pm z_{1 - \alpha/2} SE_j} on the link
#' scale, with \eqn{SE_j} from \code{\link{vcov.brs}} (Lopes, 2023,
#' "Inferencia").
#'
#' @details
#' Intervals for a mean or precision on the response scale follow by the
#' inverse link of the limits (monotone links). A limit is \code{NA} when the
#' variance is not estimable (see 'Fit diagnostics' in \code{\link{brs}}).
#' For small samples or parameters near the border of the scale,
#' \code{\link{brs_bootstrap}} gives intervals that do not rely on the normal
#' approximation.
#'
#' @param object A fitted \code{"brs"} object.
#' @param parm Character or integer: which parameters. If missing, all
#'   parameters of \code{model} are returned.
#' @param level Confidence level (default 0.95).
#' @param model Character: \code{"full"}, \code{"mean"} or
#'   \code{"precision"}.
#' @param ... Currently ignored.
#'
#' @return Matrix with the lower and upper limits.
#'
#' @seealso \code{\link{brs}}, \code{\link{vcov.brs}},
#'   \code{\link{brs_bootstrap}}, \code{\link{brs_est}}
#'
#' @references
#' Lopes, J. E. (2023). \emph{Modelos de regressao beta para dados de escala}.
#' Master's dissertation, Universidade Federal do Parana, Curitiba.
#' URI: https://hdl.handle.net/1884/86624.
#'
#' @examples
#' set.seed(2023)
#' d <- data.frame(x = runif(150))
#' s <- brs_sim(~ x, data = d, beta = c(-0.5, 1), phi = qlogis(0.3), ncuts = 10)
#' fit <- brs(y ~ x, data = s)
#' confint(fit)
#' # Mean at x = 0 on (0, 1): inverse logit of the intercept limits
#' plogis(confint(fit, parm = "(Intercept)"))
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
  se <- sqrt(diag(vcov(object, model = model)))
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

# Expected score sum_s s P(S = s) from the first parameter and phi, using the
# cells of the fit's interval (shared by predict.brs and predict.brsmm).
.brs_expected_score <- function(mu, phi, object) {
  K <- as.integer(object$ncuts)
  P <- .brs_score_prob_matrix(
    mu = mu, phi = phi, repar = object$repar, ncuts = K, lim = object$lim,
    scores = 0:K, interval = .brs_interval_of(object)
  )
  as.numeric(P %*% (0:K))
}

#' Predict from a fitted model
#'
#' @param object  A fitted \code{"brs"} object.
#' @param newdata Optional data frame for prediction.
#' @param type    Prediction type: \code{"response"} (default; the mean
#'   \eqn{E[Y] = a / (a + b)}), \code{"link"} (linear predictor of the
#'   first parameter), \code{"precision"} (second parameter on its own
#'   scale), \code{"variance"}, \code{"quantile"}, \code{"score"} or
#'   \code{"expected_score"}. \code{"score"} is the latent continuous score
#'   of the fit's \code{interval} at \eqn{y^* = E[Y]}: \eqn{K y^*} on
#'   \eqn{(0, K)} for \code{"mid"}, \eqn{(K + 1) y^*} on \eqn{(0, K + 1)} for
#'   \code{"right"} and \eqn{(K + 1) y^* - 1} on \eqn{(-1, K)} for
#'   \code{"left"} (it can be negative); under \code{"right"}/\code{"left"} it
#'   is about 0.5 above/below the expected recorded score, because a recorded
#'   score is the lower/upper end of its latent interval.
#'   \code{"expected_score"} is the expected recorded score
#'   \eqn{\sum_s s\, P(S = s)} (\code{\link{brs_predict_scoreprob}}), the
#'   same for \code{"right"} and \code{"left"}, and a proper expectation only
#'   when the cells partition \eqn{[0, 1]}.
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
                          "quantile", "score", "expected_score"
                        ),
                        at = 0.5, ...) {
  .check_class(object)
  type <- match.arg(type)

  # mu is the FIRST parameter (shape p under repar 0); E[Y] via .brs_mean().
  if (is.null(newdata)) {
    mu <- object$hatmu
    phi <- rep_len(as.numeric(object$hatphi), length(mu))
    eta_mu <- as.numeric(object$model_matrices$X %*% object$coefficients$mean)
  } else {
    # xlev = levels of the fit: a new factor level errors instead of misaligning X
    mt_mu <- stats::delete.response(object$terms$mean)
    mf <- stats::model.frame(mt_mu, data = newdata,
                             xlev = object$xlevels$mean, ...)
    X <- stats::model.matrix(mt_mu, mf)
    eta_mu <- as.numeric(X %*% object$coefficients$mean)
    mu <- .clamp_mu_by_repar(apply_inv_link(eta_mu, object$link), object$repar)

    # BUG-H02: detect variable-dispersion by presence of Z matrix with
    # non-intercept columns, not by q > 1 (which misclassifies y ~ x | 1).
    has_var_phi <- !is.null(object$model_matrices$Z) &&
                   !is.null(object$terms$precision) &&
                   length(attr(object$terms$precision, "term.labels")) > 0L
    # Build Z from newdata (variable dispersion)
    if (has_var_phi) {
      mt_phi <- object$terms$precision
      mf_z <- stats::model.frame(mt_phi, data = newdata,
                                 xlev = object$xlevels$precision, ...)
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
    # Same clamp as the compiled likelihood and as hatphi.
    phi <- .clamp_phi_by_repar(phi, object$repar)
  }

  switch(type,
    response = .brs_mean(mu, phi, object$repar),
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
    # Latent score of the mean on the original scale (direction-specific)
    score = .brs_latent_score(.brs_mean(mu, phi, object$repar), object$ncuts,
                              .brs_interval_of(object)),
    # sum_s s P(S = s) from the score cells of the fit's interval
    expected_score = .brs_expected_score(mu, phi, object),
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
#' URI: https://hdl.handle.net/1884/86624.
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
#' @param object A fitted \code{"brs"} object.
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
#' URI: https://hdl.handle.net/1884/86624.
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
  se <- sqrt(diag(V))
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
#' URI: https://hdl.handle.net/1884/86624.
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
#' URI: https://hdl.handle.net/1884/86624.
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
  if (!inherits(object, c("brs", "brsmm"))) {
    stop("Expected a 'brs' or 'brsmm' object.", call. = FALSE)
  }
  object$hessian
}
