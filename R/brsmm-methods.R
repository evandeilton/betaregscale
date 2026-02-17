# ============================================================================ #
# S3 methods for brsmm objects
# ============================================================================ #

#' Validate a brsmm object
#' @param x Object to validate.
#' @param call. Logical; passed to \code{stop()}.
#' @keywords internal
.check_class_mm <- function(x, call. = FALSE) {
  if (!inherits(x, "brsmm")) {
    stop(
      "Expected an object of class 'brsmm', got '",
      paste(class(x), collapse = "', '"), "'.",
      call. = call.
    )
  }
}


#' Extract coefficients from a brsmm fit
#'
#' @param object A fitted \code{"brsmm"} object.
#' @param model Character: \code{"full"} (default), \code{"mean"},
#'   \code{"precision"}, or \code{"random"}.
#' @param ... Currently ignored.
#'
#' @return Named numeric vector.
#'
#' @method coef brsmm
#' @importFrom stats coef
#' @export
coef.brsmm <- function(object,
                       model = c("full", "mean", "precision", "random"),
                       ...) {
  .check_class_mm(object)
  model <- match.arg(model)
  switch(model,
    full = object$par,
    mean = object$coefficients$mean,
    precision = object$coefficients$precision,
    random = object$coefficients$random
  )
}


#' Variance-covariance matrix for brsmm coefficients
#'
#' @param object A fitted \code{"brsmm"} object.
#' @param model Character: \code{"full"}, \code{"mean"},
#'   \code{"precision"}, or \code{"random"}.
#' @param ... Currently ignored.
#'
#' @return Numeric matrix.
#'
#' @method vcov brsmm
#' @importFrom stats vcov
#' @export
vcov.brsmm <- function(object,
                       model = c("full", "mean", "precision", "random"),
                       ...) {
  .check_class_mm(object)
  model <- match.arg(model)

  V <- tryCatch(
    solve(-object$hessian),
    error = function(e) {
      warning(
        "Hessian is computationally singular; returning generalized inverse.",
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

  p <- object$p
  q <- object$q
  k_re <- object$k_re
  idx_mean <- seq_len(p)
  idx_precision <- p + seq_len(q)
  idx_random <- p + q + seq_len(k_re)

  switch(model,
    full = V,
    mean = V[idx_mean, idx_mean, drop = FALSE],
    precision = V[idx_precision, idx_precision, drop = FALSE],
    random = V[idx_random, idx_random, drop = FALSE]
  )
}


#' Log-likelihood for brsmm models
#'
#' @param object A fitted \code{"brsmm"} object.
#' @param ... Currently ignored.
#'
#' @return Object of class \code{"logLik"}.
#'
#' @method logLik brsmm
#' @importFrom stats logLik
#' @export
logLik.brsmm <- function(object, ...) {
  .check_class_mm(object)
  val <- object$value
  attr(val, "df") <- object$npar
  attr(val, "nobs") <- object$nobs
  class(val) <- "logLik"
  val
}


#' AIC for brsmm models
#'
#' @param object A fitted \code{"brsmm"} object.
#' @param ... Currently ignored.
#' @param k Numeric penalty per parameter.
#'
#' @return Numeric scalar.
#'
#' @method AIC brsmm
#' @importFrom stats AIC
#' @export
AIC.brsmm <- function(object, ..., k = 2) {
  .check_class_mm(object)
  -2 * object$value + k * object$npar
}


#' BIC for brsmm models
#'
#' @param object A fitted \code{"brsmm"} object.
#' @param ... Currently ignored.
#'
#' @return Numeric scalar.
#'
#' @method BIC brsmm
#' @importFrom stats BIC
#' @export
BIC.brsmm <- function(object, ...) {
  .check_class_mm(object)
  -2 * object$value + log(object$nobs) * object$npar
}


#' Number of observations in a brsmm fit
#'
#' @param object A fitted \code{"brsmm"} object.
#' @param ... Currently ignored.
#'
#' @return Integer.
#'
#' @method nobs brsmm
#' @importFrom stats nobs
#' @export
nobs.brsmm <- function(object, ...) {
  .check_class_mm(object)
  object$nobs
}


#' Fitted values from a brsmm model
#'
#' @param object A fitted \code{"brsmm"} object.
#' @param type Character: \code{"mu"} (default) or \code{"phi"}.
#' @param ... Currently ignored.
#'
#' @return Numeric vector.
#'
#' @method fitted brsmm
#' @importFrom stats fitted
#' @export
fitted.brsmm <- function(object, type = c("mu", "phi"), ...) {
  .check_class_mm(object)
  type <- match.arg(type)
  if (identical(type, "mu")) {
    return(object$fitted_mu)
  }
  object$fitted_phi
}


#' Predict from a brsmm model
#'
#' @param object A fitted \code{"brsmm"} object.
#' @param newdata Optional data frame.
#' @param type Character: \code{"response"} (default), \code{"link"},
#'   \code{"precision"}, or \code{"variance"}.
#' @param ... Currently ignored.
#'
#' @return Numeric vector.
#'
#' @method predict brsmm
#' @importFrom stats predict model.frame model.matrix delete.response
#' @export
predict.brsmm <- function(object,
                          newdata = NULL,
                          type = c("response", "link", "precision", "variance"),
                          ...) {
  .check_class_mm(object)
  type <- match.arg(type)

  p <- object$p
  q <- object$q
  beta <- object$par[seq_len(p)]
  gamma <- object$par[p + seq_len(q)]

  if (is.null(newdata)) {
    if (is.matrix(object$random$mode_b)) {
      b_obs <- object$random$mode_b[object$group_index, , drop = FALSE]
      eta_mu <- as.numeric(object$model_matrices$X %*% beta) +
        rowSums(object$model_matrices$Xr * b_obs)
    } else {
      b_obs <- object$random$mode_b[object$group_index]
      eta_mu <- as.numeric(object$model_matrices$X %*% beta) +
        object$model_matrices$Xr[, 1L] * b_obs
    }
    eta_phi <- as.numeric(object$model_matrices$Z %*% gamma)
  } else {
    if (!is.data.frame(newdata)) {
      stop("'newdata' must be a data.frame.", call. = FALSE)
    }

    tm_mu <- stats::delete.response(object$terms$mean)
    tm_mu <- stats::delete.response(object$terms$mean)
    mf_mu <- stats::model.frame(tm_mu, data = newdata, ...)
    Xn <- stats::model.matrix(tm_mu, mf_mu)

    mf_phi <- stats::model.frame(object$terms$precision, data = newdata, ...)
    Zn <- stats::model.matrix(object$terms$precision, mf_phi)

    mf_r <- stats::model.frame(object$random$re_terms, data = newdata, ...)
    Xrn <- stats::model.matrix(object$random$re_terms, mf_r)
    if (nrow(Xrn) != nrow(Xn)) {
      stop(
        "Rows used by fixed and random design matrices in 'newdata' do not match.",
        call. = FALSE
      )
    }

    if (is.matrix(object$random$mode_b)) {
      bnew <- matrix(0, nrow = nrow(Xrn), ncol = ncol(Xrn))
      if (object$random$group %in% names(newdata)) {
        gnew <- as.character(newdata[[object$random$group]])
        idx <- match(gnew, rownames(object$random$mode_b))
        ok <- !is.na(idx)
        if (any(ok)) {
          bnew[ok, ] <- object$random$mode_b[idx[ok], , drop = FALSE]
        }
      }
      eta_mu <- as.numeric(Xn %*% beta + rowSums(Xrn * bnew))
    } else {
      bnew <- rep(0, nrow(Xrn))
      if (object$random$group %in% names(newdata)) {
        gnew <- as.character(newdata[[object$random$group]])
        map <- object$random$mode_b
        bnew <- as.numeric(map[gnew])
        bnew[is.na(bnew)] <- 0
      }
      eta_mu <- as.numeric(Xn %*% beta + Xrn[, 1L] * bnew)
    }
    eta_phi <- as.numeric(Zn %*% gamma)
  }

  mu <- apply_inv_link(eta_mu, object$link)
  phi <- apply_inv_link(eta_phi, object$link_phi)

  switch(type,
    response = mu,
    link = eta_mu,
    precision = phi,
    variance = {
      shp <- brs_repar(mu = mu, phi = phi, repar = object$repar)
      s <- shp$shape1 + shp$shape2
      (shp$shape1 * shp$shape2) / (s^2 * (s + 1))
    }
  )
}


#' Residuals from a brsmm model
#'
#' @param object A fitted \code{"brsmm"} object.
#' @param type Character: \code{"response"} (default) or \code{"pearson"}.
#' @param ... Currently ignored.
#'
#' @return Numeric vector.
#'
#' @method residuals brsmm
#' @importFrom stats residuals
#' @export
residuals.brsmm <- function(object, type = c("response", "pearson"), ...) {
  .check_class_mm(object)
  type <- match.arg(type)

  y <- as.numeric(object$Y[, "yt"])
  mu <- as.numeric(object$fitted_mu)
  r <- y - mu
  if (identical(type, "response")) {
    return(r)
  }

  v <- predict(object, type = "variance")
  v <- pmax(v, 1e-12)
  r / sqrt(v)
}


#' Summarize a fitted brsmm model
#'
#' @param object A fitted \code{"brsmm"} object.
#' @param ... Currently ignored.
#'
#' @return Object of class \code{"summary.brsmm"}.
#'
#' @method summary brsmm
#' @importFrom stats pnorm
#' @export
summary.brsmm <- function(object, ...) {
  .check_class_mm(object)
  V <- vcov(object, model = "full")
  se <- sqrt(diag(V))
  est <- object$par
  z <- est / se
  p <- 2 * (1 - stats::pnorm(abs(z)))

  tab <- cbind(
    Estimate = est,
    `Std. Error` = se,
    `z value` = z,
    `Pr(>|z|)` = p
  )
  out <- list(
    call = object$call,
    coefficients = tab,
    logLik = object$value,
    AIC = AIC(object),
    BIC = BIC(object),
    nobs = object$nobs,
    ngroups = object$ngroups,
    integration = object$int_method
  )
  class(out) <- "summary.brsmm"
  out
}


#' Print summary for brsmm models
#'
#' @param x A \code{"summary.brsmm"} object.
#' @param digits Number of digits.
#' @param ... Passed to \code{printCoefmat}.
#'
#' @return Invisibly returns \code{x}.
#'
#' @method print summary.brsmm
#' @importFrom stats printCoefmat
#' @export
print.summary.brsmm <- function(x,
                                digits = max(3, getOption("digits") - 3),
                                ...) {
  cat("\nCall:\n")
  print(x$call)
  method_name <- switch(x$integration,
    "laplace" = "Laplace",
    "aghq" = "AGHQ",
    "qmc" = "QMC",
    x$integration
  )
  cat("\nMixed beta interval model (", method_name, ")\n", sep = "")
  cat("Observations:", x$nobs, " | Groups:", x$ngroups, "\n")
  cat("logLik =", formatC(x$logLik, digits = digits, format = "f"),
    " | AIC =", formatC(x$AIC, digits = digits, format = "f"),
    " | BIC =", formatC(x$BIC, digits = digits, format = "f"), "\n\n",
    sep = ""
  )
  stats::printCoefmat(x$coefficients, digits = digits, P.values = TRUE, has.Pvalue = TRUE, ...)
  invisible(x)
}


#' Print a fitted brsmm model
#'
#' @param x A fitted \code{"brsmm"} object.
#' @param digits Number of digits.
#' @param ... Included for consistency with generic methods.
#'
#' @return Invisibly returns \code{x}.
#'
#' @method print brsmm
#' @export
print.brsmm <- function(x,
                        digits = max(3, getOption("digits") - 3),
                        ...) {
  .check_class_mm(x)
  cat("\nCall:\n")
  print(x$call)
  method_name <- switch(x$int_method,
    "laplace" = "Laplace",
    "aghq" = "AGHQ",
    "qmc" = "QMC",
    x$int_method
  )
  cat("\nMixed beta interval model (", method_name, ")\n", sep = "")
  cat("Observations:", x$nobs, " | Groups:", x$ngroups, "\n")
  cat("Log-likelihood:", formatC(x$value, digits = digits, format = "f"), "\n")
  re_sd <- x$random$sd_b
  if (length(re_sd) == 1L) {
    cat("Random SD:", formatC(re_sd, digits = digits, format = "f"), "\n")
  } else {
    cat("Random SDs:\n")
    print(round(re_sd, digits))
  }
  cat("Convergence code:", x$convergence, "\n")
  invisible(x)
}


#' Extract random effects from a brsmm model
#'
#' @param object A fitted \code{"brsmm"} object.
#' @param ... Currently ignored.
#'
#' @return Named numeric vector of group-specific random-intercept modes.
#' @export
ranef.brsmm <- function(object, ...) {
  .check_class_mm(object)
  object$random$mode_b
}
