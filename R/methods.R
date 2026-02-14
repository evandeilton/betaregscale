# ============================================================================ #
# S3 methods for objects of class "betaregscale"
#
# Output style follows the betareg package convention:
#   - coef(), vcov() accept  model = c("full", "mean", "precision")
#   - summary() produces separate tables for mean and precision
#   - print() shows the call + compact coefficient vectors
#   - Wald z-tests use pnorm (not pt) per Lopes (2024, Eq. 2.34--2.35)
# ============================================================================ #

# -- Class validation helper ------------------------------------------------ #

#' Validate a betaregscale object
#' @param x Object to validate.
#' @param call. Logical; passed to \code{stop()}.
#' @keywords internal
.check_class <- function(x, call. = FALSE) {
  if (!inherits(x, "betaregscale")) {
    stop(
      "Expected an object of class 'betaregscale', got '",
      paste(class(x), collapse = "', '"), "'.",
      call. = call.
    )
  }
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
#' @method coef betaregscale
#' @importFrom stats coef
#' @export
coef.betaregscale <- function(object,
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
#' @method vcov betaregscale
#' @importFrom stats vcov
#' @export
vcov.betaregscale <- function(object,
                              model = c("full", "mean", "precision"),
                              ...) {
  .check_class(object)
  model <- match.arg(model)

  V <- tryCatch(
    solve(-object$hessian),
    error = function(e) {
      warning(
        "Hessian is computationally singular; ",
        "returning a generalized inverse.",
        call. = FALSE
      )
      if (requireNamespace("MASS", quietly = TRUE)) {
        MASS::ginv(-object$hessian)
      } else {
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
#' @method logLik betaregscale
#' @importFrom stats logLik
#' @export
logLik.betaregscale <- function(object, ...) {
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
#' @method AIC betaregscale
#' @importFrom stats AIC
#' @export
AIC.betaregscale <- function(object, ..., k = 2) {
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
#' @method BIC betaregscale
#' @importFrom stats BIC
#' @export
BIC.betaregscale <- function(object, ...) {
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
#' @method nobs betaregscale
#' @importFrom stats nobs
#' @export
nobs.betaregscale <- function(object, ...) {
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
#' @method formula betaregscale
#' @importFrom stats formula
#' @export
formula.betaregscale <- function(x, ...) {
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
#' @method model.matrix betaregscale
#' @importFrom stats model.matrix
#' @export
model.matrix.betaregscale <- function(object,
                                      model = c("mean", "precision"),
                                      ...) {
  .check_class(object)
  model <- match.arg(model)
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
#' @method summary betaregscale
#' @importFrom stats pnorm
#' @export
summary.betaregscale <- function(object, ...) {
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

  out <- list(
    call         = object$call,
    coefficients = list(mean = tab_mu, precision = tab_phi),
    residuals    = rqr,
    loglik       = object$value,
    df           = object$npar,
    nobs         = object$nobs,
    pseudo.r2    = object$pseudo.r.squared,
    link         = object$link,
    link_phi     = object$link_phi,
    convergence  = object$convergence,
    iterations   = object$iterations,
    method       = object$optim_method,
    censoring    = cens_counts,
    repar        = object$repar
  )
  class(out) <- "summary.betaregscale"
  out
}


#' Print a model summary (betareg style)
#'
#' @param x A \code{"summary.betaregscale"} object.
#' @param digits Number of digits.
#' @param ... Ignored.
#'
#' @method print summary.betaregscale
#' @importFrom stats quantile printCoefmat
#' @export
print.summary.betaregscale <- function(x,
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
    signif.stars = TRUE
  )
  cat("\n")

  # Precision model
  phi_label <- if (nrow(x$coefficients$precision) > 1L) {
    paste0("Phi coefficients (precision model with ", x$link_phi, " link):\n")
  } else {
    paste0("Phi coefficients (precision model with ", x$link_phi, " link):\n")
  }
  cat(phi_label)
  stats::printCoefmat(x$coefficients$precision,
    digits = digits,
    P.values = TRUE, has.Pvalue = TRUE,
    signif.stars = TRUE
  )

  cat("---\n")

  # Goodness-of-fit
  cat(
    "Log-likelihood:", formatC(x$loglik, format = "f", digits = 4),
    "on", x$df, "Df\n"
  )
  cat("Pseudo R-squared:", formatC(x$pseudo.r2, format = "f", digits = 4), "\n")
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
#' @param ...    Ignored.
#'
#' @method print betaregscale
#' @export
print.betaregscale <- function(x,
                               digits = max(3, getOption("digits") - 3),
                               ...) {
  cat("\nCall:\n")
  print(x$call)
  cat("\n")

  cat("Coefficients (mean model with", x$link, "link):\n")
  print(round(x$coefficients$mean, digits))
  cat("\n")

  cat("Phi coefficients (precision model with", x$link_phi, "link):\n")
  print(round(x$coefficients$precision, digits))
  cat("\n")

  invisible(x)
}


# -- Fitted values ---------------------------------------------------------- #

#' Extract fitted values
#'
#' @param object A fitted \code{"betaregscale"} object.
#' @param type   Character: \code{"mu"} (default) or \code{"phi"}.
#' @param ...    Ignored.
#'
#' @return Numeric vector of fitted values.
#'
#' @method fitted betaregscale
#' @importFrom stats fitted
#' @export
fitted.betaregscale <- function(object, type = c("mu", "phi"), ...) {
  .check_class(object)
  type <- match.arg(type)
  if (type == "mu") object$hatmu else object$hatphi
}


# -- Residuals -------------------------------------------------------------- #

#' Extract residuals
#'
#' @param object A fitted \code{"betaregscale"} object.
#' @param type   Residual type. One of \code{"response"} (default),
#'   \code{"pearson"}, \code{"deviance"}, \code{"rqr"} (randomized
#'   quantile), \code{"weighted"}, or \code{"sweighted"}.
#' @param ...    Ignored.
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
#' @method residuals betaregscale
#' @importFrom stats residuals qnorm pbeta dbeta qlogis
#' @export
residuals.betaregscale <- function(object,
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
    rp <- beta_reparam(mu, phi, repar = repar)
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
      # Individual log-density at observed vs fitted
      ll_obs <- stats::dbeta(y, sh$a, sh$b, log = TRUE)
      ll_fit <- stats::dbeta(mu, sh$a, sh$b, log = TRUE)
      sign(y - mu) * sqrt(2 * pmax(ll_obs - ll_fit, 0))
    },
    rqr = {
      sh <- get_shapes(mu, phi, repar)
      u <- stats::pbeta(y, sh$a, sh$b)
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
#' normal approximation (Lopes, 2024, Eq. 2.30--2.31).
#'
#' @param object A fitted \code{"betaregscale"} object.
#' @param parm   Character or integer: which parameters. If missing,
#'   all parameters are returned.
#' @param level  Confidence level (default 0.95).
#' @param model  Character: \code{"full"}, \code{"mean"}, or
#'   \code{"precision"}.
#' @param ...    Ignored.
#'
#' @return Matrix with columns for lower and upper confidence bounds.
#'
#' @method confint betaregscale
#' @importFrom stats confint qnorm
#' @export
confint.betaregscale <- function(object, parm, level = 0.95,
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
#' @param ...     Ignored.
#'
#' @return Numeric vector or matrix.
#'
#' @method predict betaregscale
#' @importFrom stats predict qbeta model.matrix make.link terms model.frame
#' @export
predict.betaregscale <- function(object, newdata = NULL,
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
    phi <- object$hatphi
    eta_mu <- stats::make.link(object$link)$linkfun(mu)
  } else {
    # Build X from newdata
    mt_mu <- stats::delete.response(object$terms$mean)
    mf <- stats::model.frame(mt_mu, data = newdata)
    X <- stats::model.matrix(mt_mu, mf)
    eta_mu <- as.numeric(X %*% object$coefficients$mean)
    mu <- stats::make.link(object$link)$linkinv(eta_mu)

    # Build Z from newdata (variable dispersion)
    if (!is.null(object$model_matrices$Z) && object$q > 1L) {
      mt_phi <- object$terms$precision
      mf_z <- stats::model.frame(mt_phi, data = newdata)
      Z <- stats::model.matrix(mt_phi, mf_z)
      eta_phi <- as.numeric(Z %*% object$coefficients$precision)
      phi <- stats::make.link(object$link_phi)$linkinv(eta_phi)
    } else {
      phi <- stats::make.link(object$link_phi)$linkinv(
        as.numeric(object$coefficients$precision)
      )
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
        sh <- beta_reparam(mu, phi, repar = repar)
        s <- sh$shape1 + sh$shape2
        (sh$shape1 * sh$shape2) / (s^2 * (s + 1))
      }
    },
    quantile = {
      rp <- beta_reparam(mu, phi, repar = object$repar)
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
#' @param object A fitted \code{"betaregscale"} object.
#'
#' @return Data frame with logLik, AIC, BIC, and pseudo-R-squared.
#'
#' @examples
#' \dontrun{
#' gof(fit)
#' }
#' @export
gof <- function(object) {
  .check_class(object)
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
#' @examples
#' \donttest{
#' sim <- betaregscale_simulate(
#'   formula = ~x1, data = data.frame(x1 = rnorm(50)),
#'   beta = c(0, 0.5), phi = 0.1, ncuts = 10, repar = 2
#' )
#' fit <- betaregscale(y ~ x1, data = sim, repar = 2)
#' est(fit)
#' }
#'
#' @importFrom stats pnorm
#' @export
est <- function(object, alpha = 0.05) {
  .check_class(object)
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

#' Internal coefficient table (deprecated, use est() or summary())
#'
#' @param fit   A fitted \code{"betaregscale"} object.
#' @param alpha Significance level.
#' @return A list with \code{est} and \code{gof}.
#' @keywords internal
#' @export
betaregscale_coef <- function(fit, alpha = 0.05) {
  .check_class(fit)
  list(est = est(fit, alpha = alpha), gof = gof(fit))
}

#' Extract the Hessian matrix
#'
#' @param object A fitted \code{"betaregscale"} object.
#'
#' @return Numeric Hessian matrix.
#'
#' @examples
#' \donttest{
#' sim <- betaregscale_simulate(
#'   formula = ~x1, data = data.frame(x1 = rnorm(50)),
#'   beta = c(0, 0.5), phi = 0.1, ncuts = 10, repar = 2
#' )
#' fit <- betaregscale(y ~ x1, data = sim, repar = 2)
#' hessian_matrix(fit)
#' }
#'
#' @export
hessian_matrix <- function(object) {
  .check_class(object)
  object$hessian
}
