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

  V <- .brs_vcov(object$hessian, names(object$par))

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
#' @description
#' \eqn{-2\ell + \log(n)\,k} with \eqn{n} = the number of observations
#' (\code{nobs(object)}), as \code{lme4} does, and \eqn{k} the number of
#' parameters (fixed effects, precision and packed random-effect
#' parameters). There is no single sample size for a mixed model: counting
#' groups instead (\eqn{n = } \code{object$ngroups}) penalises more and is a
#' common alternative, \code{-2 * logLik(object) + log(object$ngroups) * k}.
#' Compare BIC values only between fits with the same convention.
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
    # E[Y] = a / (a + b): fitted_mu itself under repar 1/2, p/(p+q) under 0
    return(.brs_mean(object$fitted_mu, object$fitted_phi, object$repar))
  }
  object$fitted_phi
}


#' Predict from a brsmm model
#'
#' @param object A fitted \code{"brsmm"} object.
#' @param newdata Optional data frame.
#' @param type Character: \code{"response"} (default), \code{"link"},
#'   \code{"precision"}, \code{"variance"}, \code{"quantile"},
#'   \code{"score"} or \code{"expected_score"}. \code{"score"} is the latent
#'   score of the fit's \code{interval} at the conditional mean (support
#'   \eqn{(0, K)}, \eqn{(0, K + 1)} or \eqn{(-1, K)}; about 0.5 above/below
#'   the expected recorded score under \code{"right"}/\code{"left"});
#'   \code{"expected_score"} is the expected recorded score
#'   \eqn{\sum_s s\, P(S = s)}. Details: \code{\link{predict.brs}}.
#' @param at Numeric vector of probabilities for quantile
#'   predictions (default 0.5).
#' @param ... Currently ignored.
#'
#' @return Numeric vector, except when \code{type = "quantile"} and
#'   \code{at} has length greater than 1, in which case a numeric matrix
#'   with one column per requested quantile (named \code{q_<value>}, e.g.
#'   \code{"q_0.5"}) and one row per observation.
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
                          type = c("response", "link", "precision", "variance",
                                   "quantile", "score", "expected_score"),
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

    # xlev: factor levels of the fit, so a new level is a clear error
    tm_mu <- stats::delete.response(object$terms$mean)
    mf_mu <- stats::model.frame(tm_mu, data = newdata,
                                xlev = object$xlevels$mean, ...)
    Xn <- stats::model.matrix(tm_mu, mf_mu)

    mf_phi <- stats::model.frame(object$terms$precision, data = newdata,
                                 xlev = object$xlevels$precision, ...)
    Zn <- stats::model.matrix(object$terms$precision, mf_phi)

    mf_r <- stats::model.frame(object$random$re_terms, data = newdata,
                               xlev = object$xlevels$random, ...)
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

  # mu is the FIRST parameter (shape p under repar 0); clamps mirror the C++.
  mu <- .clamp_mu_by_repar(apply_inv_link(eta_mu, object$link), object$repar)
  phi <- .clamp_phi_by_repar(apply_inv_link(eta_phi, object$link_phi), object$repar)

  switch(type,
    response = .brs_mean(mu, phi, object$repar),
    link = eta_mu,
    precision = phi,
    variance = {
      shp <- brs_repar(mu = mu, phi = phi, repar = object$repar)
      s <- shp$shape1 + shp$shape2
      (shp$shape1 * shp$shape2) / (s^2 * (s + 1))
    },
    # Latent score of the mean and sum_s s P(S = s), as in predict.brs
    score = .brs_latent_score(.brs_mean(mu, phi, object$repar), object$ncuts,
                              .brs_interval_of(object)),
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
  # `mu` is the FIRST parameter (shape p under repar 0); `ey` is E[Y].
  mu <- as.numeric(object$fitted_mu)
  phi <- as.numeric(object$fitted_phi)
  repar <- object$repar
  ey <- .brs_mean(mu, phi, repar)
  r <- y - ey
  if (type == "response") {
    return(r)
  }

  get_shapes <- function(mu, phi, repar) {
    rp <- brs_repar(mu, phi, repar = repar)
    list(a = rp$shape1, b = rp$shape2)
  }

  switch(type,
    pearson = {
      v <- predict(object, type = "variance")
      v <- pmax(v, 1e-12)
      r / sqrt(v)
    },
    deviance = {
      # Ferrari & Cribari-Neto, identical to residuals.brs: saturated model with
      # mean y and precision a + b (the old form evaluated the density at E[Y]).
      sh <- get_shapes(mu, phi, repar)
      y_safe <- pmin(pmax(y, 1e-7), 1 - 1e-7)
      prec <- sh$a + sh$b
      ll_sat <- stats::dbeta(y_safe, y_safe * prec, (1 - y_safe) * prec, log = TRUE)
      ll_fit <- stats::dbeta(y_safe, sh$a, sh$b, log = TRUE)
      sign(y - ey) * sqrt(abs(2 * (ll_sat - ll_fit)))
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


#' Summarize a fitted brsmm model
#'
#' @description
#' Wald tests for the fixed effects. The random effects are reported as
#' standard deviations and correlations (\code{varcorr}) with Wald intervals
#' built on a transformed scale and mapped back: \eqn{\exp} of the interval
#' for \eqn{\log SD}, \eqn{\tanh} of the interval for
#' \eqn{\mathrm{atanh}(\rho)} (delta method from the packed Cholesky
#' parameters). No test or p-value is given for them: a z-test of
#' \eqn{\log SD} tests \eqn{SD = 1}, and \eqn{SD = 0} lies on the boundary;
#' use \code{\link{anova.brsmm}} (chi-bar-square mixture) against the model
#' without the term. The randomized quantile residuals are drawn without
#' changing the caller's RNG state.
#'
#' @param object A fitted \code{"brsmm"} object.
#' @param level Confidence level of the \code{varcorr} intervals.
#' @param ... Currently ignored.
#'
#' @return Object of class \code{"summary.brsmm"}; \code{coefficients$random}
#'   holds the packed Cholesky parameters (estimate and standard error only)
#'   and \code{varcorr} the SD/correlation table.
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
summary.brsmm <- function(object, level = 0.95, ...) {
  .check_class_mm(object)

  V <- vcov(object, model = "full")

  # Setup arrays
  est <- object$par
  se <- sqrt(diag(V))
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
  # RQR draws leave the user's RNG state intact
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

  out <- list(
    call = object$call,
    coefficients = list(
      mean      = tab[idx_beta, , drop = FALSE],
      precision = tab[idx_gamma, , drop = FALSE],
      # Cholesky scale: no z/p (a z-test of log SD tests SD = 1)
      random    = tab[idx_re, c("Estimate", "Std. Error"), drop = FALSE]
    ),
    varcorr = .brsmm_varcorr(object, V[idx_re, idx_re, drop = FALSE], level),
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
  cat(paste0("Phi coefficients (precision model with ", x$link_phi, " link):\n"))
  tab_phi_display <- x$coefficients$precision
  rownames(tab_phi_display) <- .pretty_phi_names(rownames(tab_phi_display))
  stats::printCoefmat(tab_phi_display,
    digits = digits,
    P.values = TRUE, has.Pvalue = TRUE,
    signif.stars = TRUE,
    ...
  )
  cat("\n")

  # Random effects: SD / Corr with transformed Wald intervals, no tests
  vc <- x$varcorr
  cat(sprintf(paste0("Random effects (SD and Corr; %s%% Wald CI on the log / ",
                     "atanh scale; no tests, see anova()):\n"),
              format(100 * attr(vc, "level"))))
  vc_tab <- as.matrix(vc[, c("estimate", "lower", "upper")])
  dimnames(vc_tab) <- list(vc$term, c("Estimate", "Lower", "Upper"))
  print(round(vc_tab, digits))

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
  prec_display <- x$coefficients$precision
  names(prec_display) <- .pretty_phi_names(names(prec_display))
  print(round(prec_display, digits))
  cat("\n")

  cat("Random-effects parameters:\n")
  re_display <- x$coefficients$random
  names(re_display) <- .pretty_re_names(names(re_display))
  print(round(re_display, digits))
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


#' SD and correlation of the random effects with transformed Wald intervals
#'
#' @description
#' From the packed Cholesky parameters \eqn{\theta} (log-diagonal,
#' off-diagonal as is): \eqn{D = LL^\top}, \eqn{SD_r = \sqrt{D_{rr}}},
#' \eqn{\rho_{rs} = D_{rs}/(SD_r SD_s)}. Intervals are Wald intervals for
#' \eqn{\log SD} and \eqn{\mathrm{atanh}\,\rho} (analytic Jacobian, delta
#' method) mapped back by \eqn{\exp} and \eqn{\tanh}.
#'
#' @param object A \code{"brsmm"} fit.
#' @param V_re Covariance matrix of the packed Cholesky parameters.
#' @param level Confidence level.
#' @return A data frame (term, type, estimate, lower, upper, se_transformed).
#' @keywords internal
#' @noRd
.brsmm_varcorr <- function(object, V_re, level = 0.95) {
  theta <- as.numeric(object$coefficients$random)
  q <- object$q_re
  nm <- object$random$terms
  tr <- .brsmm_varcorr_transform(theta, q)
  J <- .brsmm_varcorr_jacobian(theta, q)
  # Only the parameters each quantity depends on: an NA elsewhere does not spread
  se <- vapply(seq_len(nrow(J)), function(r) {
    nz <- which(J[r, ] != 0)
    sqrt(drop(J[r, nz] %*% V_re[nz, nz, drop = FALSE] %*% J[r, nz]))
  }, numeric(1))
  z <- stats::qnorm(1 - (1 - level) / 2)
  n_sd <- q
  is_sd <- seq_along(tr) <= n_sd
  pairs <- if (q > 1L) which(lower.tri(diag(q)), arr.ind = TRUE) else NULL
  terms <- c(
    paste0("SD ", nm),
    if (q > 1L) paste0("Corr ", nm[pairs[, 1L]], ",", nm[pairs[, 2L]])
  )
  back <- function(v) ifelse(is_sd, exp(v), tanh(v))
  out <- data.frame(
    term = terms,
    type = ifelse(is_sd, "sd", "corr"),
    estimate = back(tr),
    lower = back(tr - z * se),
    upper = back(tr + z * se),
    se_transformed = se,
    row.names = NULL,
    stringsAsFactors = FALSE
  )
  attr(out, "level") <- level
  out
}

# (log SD_1..q, atanh rho_rs for r > s) from the packed Cholesky parameters.
.brsmm_varcorr_transform <- function(theta, q) {
  L <- .brsmm_unpack_chol(theta, q)
  D <- L %*% t(L)
  sdv <- sqrt(diag(D))
  out <- log(sdv)
  if (q > 1L) {
    pr <- which(lower.tri(D), arr.ind = TRUE)
    out <- c(out, atanh(D[pr] / (sdv[pr[, 1L]] * sdv[pr[, 2L]])))
  }
  out
}

# Analytic Jacobian of .brsmm_varcorr_transform(): dD = dL L' + L dL'.
.brsmm_varcorr_jacobian <- function(theta, q) {
  L <- .brsmm_unpack_chol(theta, q)
  D <- L %*% t(L)
  dd <- diag(D)
  pr <- if (q > 1L) which(lower.tri(D), arr.ind = TRUE) else matrix(0L, 0L, 2L)
  rho <- if (q > 1L) D[pr] / sqrt(dd[pr[, 1L]] * dd[pr[, 2L]]) else numeric(0)
  idx <- which(lower.tri(D, diag = TRUE), arr.ind = TRUE)
  idx <- idx[order(idx[, 2L], idx[, 1L]), , drop = FALSE]  # column-wise packing
  J <- matrix(0, q + nrow(pr), length(theta))
  for (k in seq_along(theta)) {
    i <- idx[k, 1L]
    j <- idx[k, 2L]
    dL <- matrix(0, q, q)
    # Diagonal entries are exp(theta): d/dtheta = L_ii
    dL[i, j] <- if (i == j) L[i, i] else 1
    dD <- dL %*% t(L) + L %*% t(dL)
    g_sd <- 0.5 * diag(dD) / dd
    g_rho <- if (q > 1L) {
      (dD[pr] / sqrt(dd[pr[, 1L]] * dd[pr[, 2L]]) -
        rho / 2 * (diag(dD)[pr[, 1L]] / dd[pr[, 1L]] + diag(dD)[pr[, 2L]] / dd[pr[, 2L]])) /
        (1 - rho^2)
    } else {
      numeric(0)
    }
    J[, k] <- c(g_sd, g_rho)
  }
  J
}

# Lower-triangular L from the column-wise packed Cholesky parameters.
.brsmm_unpack_chol <- function(theta, q) {
  L <- matrix(0, q, q)
  k <- 1L
  for (j in seq_len(q)) {
    for (i in j:q) {
      L[i, j] <- if (i == j) exp(theta[k]) else theta[k]
      k <- k + 1L
    }
  }
  L
}
