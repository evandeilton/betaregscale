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
#' @seealso \code{\link{brsmm}}, \code{\link{vcov.brsmm}},
#'   \code{\link{confint.brsmm}}
#'
#' @examples
#' \donttest{
#' dat <- data.frame(
#'   y = c(
#'     0, 5, 20, 50, 75, 90, 100, 30, 60, 45,
#'     10, 40, 55, 70, 85, 25, 35, 65, 80, 15
#'   ),
#'   x1 = rep(c(1, 2), 10),
#'   id = factor(rep(1:4, each = 5))
#' )
#' prep <- brs_prep(dat, ncuts = 100)
#' fit <- brsmm(y ~ x1, random = ~ 1 | id, data = prep)
#' coef(fit)
#' coef(fit, model = "mean")
#' coef(fit, model = "random")
#' }
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
#' @seealso \code{\link{brsmm}}, \code{\link{coef.brsmm}},
#'   \code{\link{confint.brsmm}}
#'
#' @examples
#' \donttest{
#' dat <- data.frame(
#'   y = c(
#'     0, 5, 20, 50, 75, 90, 100, 30, 60, 45,
#'     10, 40, 55, 70, 85, 25, 35, 65, 80, 15
#'   ),
#'   x1 = rep(c(1, 2), 10),
#'   id = factor(rep(1:4, each = 5))
#' )
#' prep <- brs_prep(dat, ncuts = 100)
#' fit <- brsmm(y ~ x1, random = ~ 1 | id, data = prep)
#' vcov(fit, model = "mean")
#' }
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


#' Extract model formula
#'
#' @param x A fitted \code{"brsmm"} object.
#' @param ... Ignored.
#'
#' @return The formula used to fit the model.
#'
#' @seealso \code{\link{brsmm}}, \code{\link{model.matrix.brsmm}}
#'
#' @examples
#' \donttest{
#' dat <- data.frame(
#'   y = c(
#'     0, 5, 20, 50, 75, 90, 100, 30, 60, 45,
#'     10, 40, 55, 70, 85, 25, 35, 65, 80, 15
#'   ),
#'   x1 = rep(c(1, 2), 10),
#'   id = factor(rep(1:4, each = 5))
#' )
#' prep <- brs_prep(dat, ncuts = 100)
#' fit <- brsmm(y ~ x1, random = ~ 1 | id, data = prep)
#' formula(fit)
#' }
#'
#' @method formula brsmm
#' @importFrom stats formula
#' @export
formula.brsmm <- function(x, ...) {
  .check_class_mm(x)
  x$formula
}


#' Extract design matrix
#'
#' @param object A fitted \code{"brsmm"} object.
#' @param model  Character: \code{"mean"} (default), \code{"precision"}, or \code{"random"}.
#' @param ... Ignored.
#'
#' @return The design matrix for the specified submodel.
#'
#' @seealso \code{\link{brsmm}}, \code{\link{formula.brsmm}}
#'
#' @examples
#' \donttest{
#' dat <- data.frame(
#'   y = c(
#'     0, 5, 20, 50, 75, 90, 100, 30, 60, 45,
#'     10, 40, 55, 70, 85, 25, 35, 65, 80, 15
#'   ),
#'   x1 = rep(c(1, 2), 10),
#'   id = factor(rep(1:4, each = 5))
#' )
#' prep <- brs_prep(dat, ncuts = 100)
#' fit <- brsmm(y ~ x1, random = ~ 1 | id, data = prep)
#' head(model.matrix(fit))
#' head(model.matrix(fit, model = "random"))
#' }
#'
#' @method model.matrix brsmm
#' @importFrom stats model.matrix
#' @export
model.matrix.brsmm <- function(object,
                               model = c("mean", "precision", "random"),
                               ...) {
  .check_class_mm(object)
  model <- match.arg(model)
  switch(model,
    mean = object$model_matrices$X,
    precision = object$model_matrices$Z,
    random = object$model_matrices$Xr
  )
}


#' Wald confidence intervals for brsmm models
#'
#' @param object A fitted \code{"brsmm"} object.
#' @param parm   Character or integer: which parameters.
#' @param level  Confidence level (default 0.95).
#' @param model  Character: \code{"full"}, \code{"mean"}, \code{"precision"}, or \code{"random"}.
#' @param ...    Currently ignored.
#'
#' @return Matrix with columns for lower and upper confidence bounds.
#'
#' @seealso \code{\link{brsmm}}, \code{\link{coef.brsmm}},
#'   \code{\link{vcov.brsmm}}
#'
#' @examples
#' \donttest{
#' dat <- data.frame(
#'   y = c(
#'     0, 5, 20, 50, 75, 90, 100, 30, 60, 45,
#'     10, 40, 55, 70, 85, 25, 35, 65, 80, 15
#'   ),
#'   x1 = rep(c(1, 2), 10),
#'   id = factor(rep(1:4, each = 5))
#' )
#' prep <- brs_prep(dat, ncuts = 100)
#' fit <- brsmm(y ~ x1, random = ~ 1 | id, data = prep)
#' confint(fit, model = "mean")
#' }
#'
#' @method confint brsmm
#' @importFrom stats confint qnorm
#' @export
confint.brsmm <- function(object, parm, level = 0.95,
                          model = c("full", "mean", "precision", "random"),
                          ...) {
  .check_class_mm(object)
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


#' Log-likelihood for brsmm models
#'
#' @param object A fitted \code{"brsmm"} object.
#' @param ... Currently ignored.
#'
#' @return Object of class \code{"logLik"}.
#'
#' @seealso \code{\link{brsmm}}, \code{\link{AIC.brsmm}},
#'   \code{\link{BIC.brsmm}}, \code{\link{brs_gof}}
#'
#' @examples
#' \donttest{
#' dat <- data.frame(
#'   y = c(
#'     0, 5, 20, 50, 75, 90, 100, 30, 60, 45,
#'     10, 40, 55, 70, 85, 25, 35, 65, 80, 15
#'   ),
#'   x1 = rep(c(1, 2), 10),
#'   id = factor(rep(1:4, each = 5))
#' )
#' prep <- brs_prep(dat, ncuts = 100)
#' fit <- brsmm(y ~ x1, random = ~ 1 | id, data = prep)
#' logLik(fit)
#' }
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
#' @seealso \code{\link{brsmm}}, \code{\link{logLik.brsmm}},
#'   \code{\link{BIC.brsmm}}, \code{\link{brs_gof}}
#'
#' @examples
#' \donttest{
#' dat <- data.frame(
#'   y = c(
#'     0, 5, 20, 50, 75, 90, 100, 30, 60, 45,
#'     10, 40, 55, 70, 85, 25, 35, 65, 80, 15
#'   ),
#'   x1 = rep(c(1, 2), 10),
#'   id = factor(rep(1:4, each = 5))
#' )
#' prep <- brs_prep(dat, ncuts = 100)
#' fit <- brsmm(y ~ x1, random = ~ 1 | id, data = prep)
#' AIC(fit)
#' }
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
#' @seealso \code{\link{brsmm}}, \code{\link{logLik.brsmm}},
#'   \code{\link{AIC.brsmm}}, \code{\link{brs_gof}}
#'
#' @examples
#' \donttest{
#' dat <- data.frame(
#'   y = c(
#'     0, 5, 20, 50, 75, 90, 100, 30, 60, 45,
#'     10, 40, 55, 70, 85, 25, 35, 65, 80, 15
#'   ),
#'   x1 = rep(c(1, 2), 10),
#'   id = factor(rep(1:4, each = 5))
#' )
#' prep <- brs_prep(dat, ncuts = 100)
#' fit <- brsmm(y ~ x1, random = ~ 1 | id, data = prep)
#' BIC(fit)
#' }
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
#' @seealso \code{\link{brsmm}}, \code{\link{fitted.brsmm}}
#'
#' @examples
#' \donttest{
#' dat <- data.frame(
#'   y = c(
#'     0, 5, 20, 50, 75, 90, 100, 30, 60, 45,
#'     10, 40, 55, 70, 85, 25, 35, 65, 80, 15
#'   ),
#'   x1 = rep(c(1, 2), 10),
#'   id = factor(rep(1:4, each = 5))
#' )
#' prep <- brs_prep(dat, ncuts = 100)
#' fit <- brsmm(y ~ x1, random = ~ 1 | id, data = prep)
#' nobs(fit)
#' }
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
#' @seealso \code{\link{brsmm}}, \code{\link{residuals.brsmm}},
#'   \code{\link{predict.brsmm}}
#'
#' @examples
#' \donttest{
#' dat <- data.frame(
#'   y = c(
#'     0, 5, 20, 50, 75, 90, 100, 30, 60, 45,
#'     10, 40, 55, 70, 85, 25, 35, 65, 80, 15
#'   ),
#'   x1 = rep(c(1, 2), 10),
#'   id = factor(rep(1:4, each = 5))
#' )
#' prep <- brs_prep(dat, ncuts = 100)
#' fit <- brsmm(y ~ x1, random = ~ 1 | id, data = prep)
#' head(fitted(fit))
#' head(fitted(fit, type = "phi"))
#' }
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
#'   \code{"precision"}, \code{"variance"}, or \code{"quantile"}.
#' @param at Numeric vector of probabilities for quantile
#'   predictions (default 0.5).
#' @param ... Currently ignored.
#'
#' @return Numeric vector.
#'
#' @seealso \code{\link{brsmm}}, \code{\link{fitted.brsmm}},
#'   \code{\link{brs_predict_scoreprob}}
#'
#' @examples
#' \donttest{
#' dat <- data.frame(
#'   y = c(
#'     0, 5, 20, 50, 75, 90, 100, 30, 60, 45,
#'     10, 40, 55, 70, 85, 25, 35, 65, 80, 15
#'   ),
#'   x1 = rep(c(1, 2), 10),
#'   id = factor(rep(1:4, each = 5))
#' )
#' prep <- brs_prep(dat, ncuts = 100)
#' fit <- brsmm(y ~ x1, random = ~ 1 | id, data = prep)
#' head(predict(fit))
#' head(predict(fit, type = "precision"))
#' }
#'
#' @method predict brsmm
#' @importFrom stats predict model.frame model.matrix delete.response qbeta
#' @export
predict.brsmm <- function(object,
                          newdata = NULL,
                          type = c("response", "link", "precision", "variance", "quantile"),
                          at = 0.5,
                          ...) {
  .check_class_mm(object)
  type <- match.arg(type)

  p <- object$p
  q <- object$q
  beta <- object$par[seq_len(p)]
  gamma <- object$par[p + seq_len(q)]

  if (is.null(newdata)) {
    eta_fixed <- as.numeric(object$model_matrices$X %*% beta)
    if (is.matrix(object$random$mode_b)) {
      b_obs <- object$random$mode_b[object$group_index, , drop = FALSE]
      eta_mu <- eta_fixed + rowSums(object$model_matrices$Xr * b_obs)
    } else {
      b_obs <- object$random$mode_b[object$group_index]
      eta_mu <- eta_fixed + object$model_matrices$Xr[, 1L] * b_obs
    }
    eta_phi <- as.numeric(object$model_matrices$Z %*% gamma)
  } else {
    if (!is.data.frame(newdata)) {
      stop("'newdata' must be a data.frame.", call. = FALSE)
    }

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

    eta_fixed <- as.numeric(Xn %*% beta)
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
      eta_mu <- eta_fixed + rowSums(Xrn * bnew)
    } else {
      bnew <- rep(0, nrow(Xrn))
      if (object$random$group %in% names(newdata)) {
        gnew <- as.character(newdata[[object$random$group]])
        map <- object$random$mode_b
        bnew <- as.numeric(map[gnew])
        bnew[is.na(bnew)] <- 0
      }
      eta_mu <- eta_fixed + Xrn[, 1L] * bnew
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


#' Residuals from a brsmm model
#'
#' @param object A fitted \code{"brsmm"} object.
#' @param type Character: \code{"response"} (default), \code{"pearson"},
#'   \code{"deviance"}, \code{"rqr"}, \code{"weighted"}, or \code{"sweighted"}.
#' @param ... Currently ignored.
#'
#' @return Numeric vector.
#'
#' @seealso \code{\link{brsmm}}, \code{\link{fitted.brsmm}},
#'   \code{\link{plot.brsmm}}
#'
#' @examples
#' \donttest{
#' dat <- data.frame(
#'   y = c(
#'     0, 5, 20, 50, 75, 90, 100, 30, 60, 45,
#'     10, 40, 55, 70, 85, 25, 35, 65, 80, 15
#'   ),
#'   x1 = rep(c(1, 2), 10),
#'   id = factor(rep(1:4, each = 5))
#' )
#' prep <- brs_prep(dat, ncuts = 100)
#' fit <- brsmm(y ~ x1, random = ~ 1 | id, data = prep)
#' head(residuals(fit))
#' head(residuals(fit, type = "pearson"))
#' }
#'
#' @method residuals brsmm
#' @importFrom stats residuals qnorm pbeta dbeta qlogis runif
#' @export
residuals.brsmm <- function(object, type = c(
                              "response", "pearson",
                              "deviance", "rqr",
                              "weighted", "sweighted"
                            ), ...) {
  .check_class_mm(object)
  type <- match.arg(type)

  y <- as.numeric(object$Y[, "yt"])
  mu <- as.numeric(object$fitted_mu)
  r <- y - mu
  if (type == "response") {
    return(r)
  }

  phi <- as.numeric(object$fitted_phi)
  repar <- object$repar

  get_shapes <- function(mu, phi, repar) {
    rp <- brs_repar(mu, phi, repar = repar)
    list(a = rp$shape1, b = rp$shape2)
  }

  to_precision <- function(mu, phi, repar) {
    if (repar == 1L) {
      return(phi)
    }
    if (repar == 2L) {
      return((1 - phi) / phi)
    }
    mu + phi
  }

  switch(type,
    pearson = {
      v <- predict(object, type = "variance")
      v <- pmax(v, 1e-12)
      r / sqrt(v)
    },
    deviance = {
      sh <- get_shapes(mu, phi, repar)
      ll_obs <- stats::dbeta(y, sh$a, sh$b, log = TRUE)
      ll_fit <- stats::dbeta(mu, sh$a, sh$b, log = TRUE)
      sign(y - mu) * sqrt(2 * pmax(ll_obs - ll_fit, 0))
    },
    rqr = {
      sh <- get_shapes(mu, phi, repar)
      left <- as.numeric(object$Y[, "left"])
      right <- as.numeric(object$Y[, "right"])
      delta <- object$delta

      f_left <- stats::pbeta(left, sh$a, sh$b)
      f_right <- stats::pbeta(right, sh$a, sh$b)
      f_y <- stats::pbeta(y, sh$a, sh$b)

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
      prec <- to_precision(mu, phi, repar)
      ystar <- stats::qlogis(y)
      mustar <- digamma(mu * prec) - digamma((1 - mu) * prec)
      v <- trigamma(mu * prec) + trigamma((1 - mu) * prec)
      if (type == "weighted") {
        (ystar - mustar) / sqrt(prec * v)
      } else {
        (ystar - mustar) / sqrt(v)
      }
    }
  )
}


#' Summarize a fitted brsmm model
#'
#' @param object A fitted \code{"brsmm"} object.
#' @param ... Currently ignored.
#'
#' @return Object of class \code{"summary.brsmm"}.
#'
#' @seealso \code{\link{brsmm}}, \code{\link{print.summary.brsmm}},
#'   \code{\link{brs_gof}}, \code{\link{brsmm_re_study}}
#'
#' @examples
#' \donttest{
#' dat <- data.frame(
#'   y = c(
#'     0, 5, 20, 50, 75, 90, 100, 30, 60, 45,
#'     10, 40, 55, 70, 85, 25, 35, 65, 80, 15
#'   ),
#'   x1 = rep(c(1, 2), 10),
#'   id = factor(rep(1:4, each = 5))
#' )
#' prep <- brs_prep(dat, ncuts = 100)
#' fit <- brsmm(y ~ x1, random = ~ 1 | id, data = prep)
#' s <- summary(fit)
#' s$coefficients$mean
#' }
#'
#' @method summary brsmm
#' @importFrom stats pnorm residuals
#' @export
summary.brsmm <- function(object, ...) {
  .check_class_mm(object)

  V <- vcov(object, model = "full")

  # Setup arrays
  est <- object$par
  se <- sqrt(pmax(diag(V), 0))
  z <- est / se
  p <- 2 * stats::pnorm(-abs(z))

  tab <- cbind(
    Estimate = est,
    `Std. Error` = se,
    `z value` = z,
    `Pr(>|z|)` = p
  )

  idx_beta <- seq_len(object$p)
  idx_gamma <- object$p + seq_len(object$q)
  idx_re <- object$p + object$q + seq_len(object$k_re)

  # Check residuals (Prioritizing randomized quantile residuals for censored data)
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
    call = object$call,
    coefficients = list(
      mean      = tab[idx_beta, , drop = FALSE],
      precision = tab[idx_gamma, , drop = FALSE],
      random    = tab[idx_re, , drop = FALSE]
    ),
    residuals = rqr,
    loglik = object$value,
    AIC = AIC(object),
    BIC = BIC(object),
    df = object$npar,
    nobs = object$nobs,
    pseudo.r2 = object$pseudo.r.squared,
    link = object$link,
    link_phi = object$link_phi,
    convergence = object$convergence,
    iterations = object$iterations,
    method = object$method,
    censoring = cens_counts,
    repar = object$repar,
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
#' @seealso \code{\link{summary.brsmm}}, \code{\link{brsmm}},
#'   \code{\link{print.brsmm}}
#'
#' @examples
#' \donttest{
#' dat <- data.frame(
#'   y = c(
#'     0, 5, 20, 50, 75, 90, 100, 30, 60, 45,
#'     10, 40, 55, 70, 85, 25, 35, 65, 80, 15
#'   ),
#'   x1 = rep(c(1, 2), 10),
#'   id = factor(rep(1:4, each = 5))
#' )
#' prep <- brs_prep(dat, ncuts = 100)
#' fit <- brsmm(y ~ x1, random = ~ 1 | id, data = prep)
#' print(summary(fit))
#' }
#'
#' @method print summary.brsmm
#' @importFrom stats printCoefmat quantile
#' @export
print.summary.brsmm <- function(x,
                                digits = max(3, getOption("digits") - 3),
                                ...) {
  cat("\nCall:\n")
  print(x$call)
  cat("\n")

  method_name <- switch(x$integration,
    "laplace" = "Laplace",
    "aghq" = "AGHQ",
    "qmc" = "QMC",
    x$integration
  )

  # Quantile residuals summary
  rq <- stats::quantile(x$residuals,
    probs = c(0, 0.25, 0.5, 0.75, 1),
    na.rm = TRUE
  )
  names(rq) <- c("Min", "1Q", "Median", "3Q", "Max")
  cat("Randomized Quantile Residuals:\n")
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
  phi_label <- paste0("Phi coefficients (precision model with ", x$link_phi, " link):\n")
  cat(phi_label)
  stats::printCoefmat(x$coefficients$precision,
    digits = digits,
    P.values = TRUE, has.Pvalue = TRUE,
    signif.stars = TRUE,
    ...
  )
  cat("\n")

  # Random effects
  cat("Random-effects parameters (Cholesky scale):\n")
  stats::printCoefmat(x$coefficients$random,
    digits = digits,
    P.values = TRUE, has.Pvalue = TRUE,
    signif.stars = TRUE,
    ...
  )

  cat("---\n")
  cat("Mixed beta interval model (", method_name, ")\n", sep = "")
  cat("Observations:", x$nobs, " | Groups:", x$ngroups, "\n")

  # Goodness-of-fit
  cat(
    "Log-likelihood:", formatC(x$loglik, format = "f", digits = 4),
    "on", x$df, "Df | AIC:", formatC(x$AIC, format = "f", digits = 4),
    "| BIC:", formatC(x$BIC, format = "f", digits = 4), "\n"
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


#' Print a fitted brsmm model
#'
#' @param x A fitted \code{"brsmm"} object.
#' @param digits Number of digits.
#' @param ... Included for consistency with generic methods.
#'
#' @return Invisibly returns \code{x}.
#'
#' @seealso \code{\link{summary.brsmm}}, \code{\link{print.summary.brsmm}},
#'   \code{\link{brsmm}}
#'
#' @examples
#' \donttest{
#' dat <- data.frame(
#'   y = c(
#'     0, 5, 20, 50, 75, 90, 100, 30, 60, 45,
#'     10, 40, 55, 70, 85, 25, 35, 65, 80, 15
#'   ),
#'   x1 = rep(c(1, 2), 10),
#'   id = factor(rep(1:4, each = 5))
#' )
#' prep <- brs_prep(dat, ncuts = 100)
#' fit <- brsmm(y ~ x1, random = ~ 1 | id, data = prep)
#' print(fit)
#' }
#'
#' @method print brsmm
#' @export
print.brsmm <- function(x,
                        digits = max(3, getOption("digits") - 3),
                        ...) {
  .check_class_mm(x)
  cat("\nCall:\n")
  print(x$call)
  cat("\n")

  method_name <- switch(x$int_method,
    "laplace" = "Laplace",
    "aghq" = "AGHQ",
    "qmc" = "QMC",
    x$int_method
  )

  cat("Coefficients (mean model with", x$link, "link):\n")
  print(round(x$coefficients$mean, digits))
  cat("\n")

  cat("Phi coefficients (precision model with", x$link_phi, "link):\n")
  print(round(x$coefficients$precision, digits))
  cat("\n")

  cat("Random-effects parameters:\n")
  print(round(x$coefficients$random, digits))
  cat("\n")

  re_sd <- x$random$sd_b
  if (length(re_sd) == 1L) {
    cat("Random SD:", formatC(re_sd, digits = digits, format = "f"), "\n")
  } else {
    cat("Random SDs:\n")
    print(round(re_sd, digits))
  }

  cat("---\n")
  cat("Mixed beta interval model (", method_name, ")\n", sep = "")
  cat("Observations:", x$nobs, " | Groups:", x$ngroups, "\n")
  cat("Log-likelihood:", formatC(x$value, digits = digits, format = "f"), "\n")
  cat("Convergence code:", x$convergence, "\n")
  invisible(x)
}


#' Extract random effects
#'
#' @description Generic function for extracting random effects.
#' @param object A fitted model object.
#' @param ... Additional arguments passed to methods.
#'
#' @return Method-specific; for \code{"brsmm"} objects, a matrix or named
#'   numeric vector of group-specific random-effect modes.
#'
#' @seealso \code{\link{ranef.brsmm}}, \code{\link{brsmm_re_study}}
#'
#' @examples
#' \donttest{
#' dat <- data.frame(
#'   y = c(
#'     0, 5, 20, 50, 75, 90, 100, 30, 60, 45,
#'     10, 40, 55, 70, 85, 25, 35, 65, 80, 15
#'   ),
#'   x1 = rep(c(1, 2), 10),
#'   id = factor(rep(1:4, each = 5))
#' )
#' prep <- brs_prep(dat, ncuts = 100)
#' fit <- brsmm(y ~ x1, random = ~ 1 | id, data = prep)
#' ranef(fit)
#' }
#'
#' @export
ranef <- function(object, ...) UseMethod("ranef")


#' Extract random effects from a brsmm model
#'
#' @param object A fitted \code{"brsmm"} object.
#' @param ... Currently ignored.
#'
#' @return A matrix or named numeric vector of group-specific random-effect
#'   posterior modes.
#'
#' @method ranef brsmm
#'
#' @seealso \code{\link{brsmm}}, \code{\link{brsmm_re_study}},
#'   \code{\link{ranef}}
#'
#' @examples
#' \donttest{
#' dat <- data.frame(
#'   y = c(
#'     0, 5, 20, 50, 75, 90, 100, 30, 60, 45,
#'     10, 40, 55, 70, 85, 25, 35, 65, 80, 15
#'   ),
#'   x1 = rep(c(1, 2), 10),
#'   id = factor(rep(1:4, each = 5))
#' )
#' prep <- brs_prep(dat, ncuts = 100)
#' fit <- brsmm(y ~ x1, random = ~ 1 | id, data = prep)
#' ranef(fit)
#' }
#'
#' @export
ranef.brsmm <- function(object, ...) {
  .check_class_mm(object)
  object$random$mode_b
}
