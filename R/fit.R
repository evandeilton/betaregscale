# ============================================================================ #
# Model-fitting functions
# ============================================================================ #

.validate_brs_common_args <- function(data, ncuts, lim, repar,
                                      link = NULL, link_phi = NULL,
                                      interval = NULL) {
  if (!is.data.frame(data)) {
    stop("`data` must be a data.frame.", call. = FALSE)
  }
  # The endpoints were built with brs_prep()'s ncuts/lim/interval: the
  # attribute wins (scoreprob/autoplot need the same grid).
  lim_from_attr <- isTRUE(attr(data, "is_prepared", exact = TRUE)) &&
    !is.null(attr(data, "lim", exact = TRUE))
  # Explicit interval: partial matching, as match.arg() in brs_check()/brs_prep()
  if (!is.null(interval)) {
    interval <- tryCatch(match.arg(interval, .brs_intervals), error = function(e) {
      stop("`interval` must be one of \"mid\", \"right\", \"left\".", call. = FALSE)
    })
  }
  # Prepared data without an interval attribute (made before `interval`
  # existed, or by hand) hold mid cells: that is the stored value
  if (isTRUE(attr(data, "is_prepared", exact = TRUE)) &&
    is.null(attr(data, "interval", exact = TRUE))) {
    attr(data, "interval") <- "mid"
  }
  ncuts <- .brs_prep_attr_arg(ncuts, data, "ncuts", 100L)
  lim <- .brs_prep_attr_arg(lim, data, "lim", 0.5)
  interval <- .brs_prep_attr_arg(interval, data, "interval", "mid")
  if (!is.character(interval) || length(interval) != 1L ||
    !(interval %in% .brs_intervals)) {
    stop("`interval` must be one of \"mid\", \"right\", \"left\".", call. = FALSE)
  }
  ncuts <- as.integer(ncuts)
  if (length(ncuts) != 1L || !is.finite(ncuts) || ncuts < 2L) {
    stop("`ncuts` must be an integer >= 2.", call. = FALSE)
  }
  # brs_prep() already warned about a stored lim; warn only for explicit values
  .brs_lim_check(lim, interval, warn = !lim_from_attr)
  repar <- as.integer(repar)
  if (length(repar) != 1L || is.na(repar) || !(repar %in% 0:2)) {
    stop("`repar` must be one of 0, 1, or 2.", call. = FALSE)
  }
  links <- .resolve_links(link, link_phi, repar)
  list(
    ncuts = ncuts, lim = as.numeric(lim), repar = repar,
    link = links$link, link_phi = links$link_phi, interval = interval
  )
}

# Value of `ncuts`/`lim`/`interval`: the brs_prep() attribute when present
# (warning if an explicit different value was passed), else the explicit
# value, else default.
.brs_prep_attr_arg <- function(value, data, name, default) {
  # Not prepared data: explicit value, else default.
  stored <- attr(data, name, exact = TRUE)
  if (!isTRUE(attr(data, "is_prepared", exact = TRUE)) || is.null(stored)) {
    return(if (is.null(value)) default else value)
  }
  if (is.null(value)) {
    return(stored)
  }
  # Numeric attributes compare by value, character ones (interval) exactly
  same <- if (is.numeric(stored)) {
    isTRUE(all.equal(as.numeric(value), as.numeric(stored)))
  } else {
    identical(as.character(value), as.character(stored))
  }
  if (!same) {
    warning(
      "`", name, " = ", format(value), "` differs from the value used by ",
      "brs_prep() (", format(stored), "); the prepared endpoints were built ",
      "with ", format(stored), ", which is used instead.",
      call. = FALSE
    )
  }
  stored
}

# Pseudo R2 = cor(eta, g(y))^2 (betareg style). Under repar 0 eta is for the
# shape p, so logit(E[Y]) vs logit(y) is used instead.
.brs_pseudo_r2 <- function(eta, ey, y, link, repar) {
  y_safe <- pmin(pmax(y, 1e-7), 1 - 1e-7)
  if (as.integer(repar) == 0L) {
    ey_safe <- pmin(pmax(ey, 1e-7), 1 - 1e-7)
    u <- stats::qlogis(ey_safe)
  } else {
    u <- as.numeric(eta)
  }
  v <- apply_link(y_safe, if (as.integer(repar) == 0L) "logit" else link)
  # An intercept-only model (or a constant response) has no correlation to
  # report: NA rather than cor()'s zero-sd warning.
  if (length(u) < 2L || !all(is.finite(u)) || !all(is.finite(v)) ||
    stats::sd(u) == 0 || stats::sd(v) == 0) {
    return(NA_real_)
  }
  r2 <- stats::cor(u, v)^2
  if (!is.finite(r2)) NA_real_ else as.numeric(r2)
}

# The sqrt link maps eta <= 0 to 0 (flat): estimates on that plateau are not
# a proper maximum of the likelihood in that parameter.
.warn_sqrt_plateau <- function(eta, link, arg) {
  if (identical(link, "sqrt")) {
    n_flat <- sum(eta <= 0, na.rm = TRUE)
    if (n_flat > 0L) {
      warning(
        "With `", arg, " = \"sqrt\"` the fitted linear predictor is <= 0 for ",
        n_flat, " observation(s): the inverse link is flat there (parameter ",
        "clamped at 0) and the estimates sit on a plateau. Consider `", arg,
        " = \"log\"`.",
        call. = FALSE
      )
    }
  }
  invisible(NULL)
}

#' Fit a fixed-dispersion beta interval regression model
#'
#' @description
#' Estimates the parameters of a beta regression model with a single
#' (scalar) dispersion parameter using maximum likelihood.  The
#' log-likelihood and its gradient are evaluated by the compiled C++
#' backend supporting the complete likelihood with mixed censoring
#' types.
#'
#' @param formula Two-sided formula \code{y ~ x1 + x2 + ...}.
#' @param data   Data frame.
#' @param link   Link for the first parameter (the mean under
#'   \code{repar = 1, 2}; the shape \eqn{p} under \code{repar = 0}).
#'   \code{NULL} (default) selects the link implied by \code{repar}:
#'   \code{"logit"} for \code{repar = 1, 2}, \code{"log"} for
#'   \code{repar = 0}. See the 'Reparameterizations and links' section of
#'   \code{\link{brs}} for the admissible values.
#' @param link_phi Link for the second parameter. \code{NULL} (default)
#'   selects \code{"logit"} for \code{repar = 2} (dispersion on
#'   \eqn{(0, 1)}) and \code{"log"} for \code{repar = 0, 1} (positive
#'   shape/precision).
#' @param ncuts  Number of scale categories. \code{NULL} (default) uses the
#'   value stored by \code{\link{brs_prep}} in \code{attr(data, "ncuts")},
#'   or 100 when \code{data} was not prepared. A value that differs from
#'   the stored one is ignored with a warning (the endpoints were built
#'   with the stored value).
#' @param lim    Half-width of the score cell in \eqn{(0, 0.5]}
#'   (\code{interval = "mid"} only). \code{NULL} (default) uses
#'   \code{attr(data, "lim")} from \code{\link{brs_prep}}, or 0.5; same
#'   rule as \code{ncuts}. Values below 0.5 warn (partial coarsening).
#' @param hessian_method Character: \code{"cpp"} (default),
#'   \code{"numDeriv"} or \code{"optim"}. \code{"cpp"} uses the compiled
#'   chain-rule Hessian (per-observation second derivatives in the linear
#'   predictors); \code{"numDeriv"} differentiates the log-likelihood with
#'   \code{\link[numDeriv]{hessian}}; \code{"optim"} keeps the optimizer's
#'   own approximation.
#' @param repar  Reparameterization scheme (default 2); see
#'   \code{\link{brs_repar}}.
#' @param method Optimization method: \code{"BFGS"} (default) or
#'   \code{"L-BFGS-B"}.
#' @param start Optional numeric vector of starting values (mean
#'   coefficients, then dispersion coefficients). \code{NULL} (default) uses
#'   \code{compute_start()}. Refits (bootstrap, jackknife) pass the parent
#'   estimate here.
#' @param control Control list for \code{\link[stats]{optim}}; its entries
#'   are merged into the default \code{list(maxit = 5000L)}.
#' @param interval Direction of the uncertainty interval, \code{"mid"},
#'   \code{"right"} or \code{"left"} (see \code{\link{brs_check}}).
#'   \code{NULL} (default) uses \code{attr(data, "interval")} from
#'   \code{\link{brs_prep}}, or \code{"mid"}; same rule as \code{ncuts}.
#'
#' @return An object of class \code{"brs"}.
#'
#' @examples
#' \donttest{
#' dat <- data.frame(
#'   y = c(
#'     0, 5, 20, 50, 75, 90, 100, 30, 60, 45,
#'     10, 40, 55, 70, 85, 25, 35, 65, 80, 15
#'   ),
#'   x1 = rep(c(1, 2), 10),
#'   x2 = rep(c(0, 0, 1, 1), 5)
#' )
#' prep <- brs_prep(dat, ncuts = 100)
#' fit <- brs_fit_fixed(y ~ x1 + x2, data = prep)
#' print(fit)
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
#' @importFrom stats optim cor model.frame model.matrix model.response terms
#' @importFrom numDeriv hessian
#' @keywords internal
#' @export
brs_fit_fixed <- function(formula, data,
                          link = NULL,
                          link_phi = NULL,
                          ncuts = NULL,
                          lim = NULL,
                          hessian_method = c("cpp", "numDeriv", "optim"),
                          repar = 2L,
                          method = c("BFGS", "L-BFGS-B"),
                          interval = NULL,
                          start = NULL,
                          control = list()) {
  cl <- match.call()
  method <- match.arg(method)
  hessian_method <- match.arg(hessian_method)
  validated <- .validate_brs_common_args(data, ncuts, lim, repar, link, link_phi,
                                        interval)
  ncuts <- validated$ncuts
  lim <- validated$lim
  repar <- validated$repar
  link <- validated$link
  link_phi <- validated$link_phi
  interval <- validated$interval

  # Build matrices
  mf <- stats::model.frame(formula, data = data)
  mtX <- stats::terms(formula, data = data, rhs = 1L)
  Y <- .extract_response(mf, data, ncuts = ncuts, lim = lim, interval = interval)
  X <- stats::model.matrix(mtX, mf)
  n <- nrow(X)
  p <- ncol(X)
  # Aliased columns stop here instead of a non-finite value inside optim
  .brs_check_design(X, "mean")

  # Extract delta from brs_check output
  delta <- as.integer(Y[, "delta"])

  # Starting values: user / parent estimate, else compute_start()
  ini <- if (is.null(start)) {
    compute_start(
      formula = formula, data = data, link = link,
      link_phi = link_phi, ncuts = ncuts,
      lim = lim, repar = repar, interval = interval
    )
  } else {
    .brs_check_start(start, p + 1L)
  }

  # Pre-compute link codes for C++
  lc_mu <- link_to_code(link)
  lc_phi <- link_to_code(link_phi)

  # Objective: -loglik (we minimize)
  fn_obj <- function(par) {
    -.brs_loglik_fixed_cpp(
      par, X, Y[, "left"], Y[, "right"], Y[, "yt"],
      delta, lc_mu, lc_phi, repar
    )
  }

  # Gradient: -grad
  gr_obj <- function(par) {
    -.brs_grad_fixed_cpp(
      par, X, Y[, "left"], Y[, "right"], Y[, "yt"],
      delta, lc_mu, lc_phi, repar
    )
  }

  # Optimize
  opt <- stats::optim(
    par     = ini,
    fn      = fn_obj,
    gr      = gr_obj,
    method  = method,
    hessian = (hessian_method == "optim"),
    control = .brs_merge_control(list(maxit = 5000L), control)
  )

  # BUG-H04: warn if optimizer did not converge
  if (opt$convergence != 0L) {
    warning(
      "Optimizer did not converge (code ", opt$convergence, ")",
      if (!is.null(opt$message)) paste0(": ", opt$message) else ".",
      "\nResults may be unreliable. Try increasing 'maxit' or changing 'method'.",
      call. = FALSE
    )
  }

  # Hessian (on the log-likelihood scale)
  if (hessian_method == "cpp") {
    opt$hessian <- .brs_hessian_fixed_cpp(
      opt$par, X, Y[, "left"], Y[, "right"], Y[, "yt"],
      delta, lc_mu, lc_phi, repar
    )
  } else if (hessian_method == "numDeriv") {
    fn_ll <- function(par) {
      .brs_loglik_fixed_cpp(
        par, X, Y[, "left"], Y[, "right"], Y[, "yt"],
        delta, lc_mu, lc_phi, repar
      )
    }
    opt$hessian <- numDeriv::hessian(fn_ll, opt$par)
  } else {
    opt$hessian <- -opt$hessian
  }

  # hatmu is the FIRST parameter (shape p under repar 0); E[Y] via .brs_mean().
  est <- opt$par
  eta_mu <- X %*% est[1:p]
  mu_raw <- apply_inv_link(eta_mu, link)
  phi_raw <- apply_inv_link(est[p + 1L], link_phi)
  # Same clamps as the compiled likelihood (src/brs_common.h).
  hatmu <- .clamp_mu_by_repar(mu_raw, repar)
  hatphi <- .clamp_phi_by_repar(phi_raw, repar)
  # Gradient, Hessian and clamp checks at the estimate (warn + stored)
  diagnostics <- .brs_fit_diagnostics(-gr_obj(est), opt$hessian, mu_raw, phi_raw,
                                      repar, opt$convergence,
                                      badly_scaled = .brs_badly_scaled(X))
  y_mid <- Y[, "yt"]
  ey <- .brs_mean(as.numeric(hatmu), hatphi, repar)
  resid <- as.numeric(y_mid - ey)

  pseudo_r2 <- .brs_pseudo_r2(eta_mu, ey, y_mid, link, repar)
  .warn_sqrt_plateau(eta_mu, link, "link")
  .warn_sqrt_plateau(est[p + 1L], link_phi, "link_phi")

  # --- betareg-style parameter naming ---
  # Mean coefficients: use column names of X
  mean_names <- colnames(X)
  # Precision coefficient: single scalar for fixed model
  phi_names <- "(phi)"

  par_names <- c(mean_names, phi_names)
  names(est) <- par_names

  # Named coefficient lists (betareg style)
  coefficients <- list(
    mean      = est[seq_len(p)],
    precision = est[p + 1L]
  )
  names(coefficients$precision) <- phi_names

  # Name the hessian
  rownames(opt$hessian) <- colnames(opt$hessian) <- par_names

  # Build result object
  result <- list(
    call             = cl,
    par              = est,
    coefficients     = coefficients,
    # Exactly the log-likelihood at the returned estimate
    value            = -fn_obj(opt$par),
    hessian          = opt$hessian,
    convergence      = opt$convergence,
    message          = opt$message,
    iterations       = opt$counts,
    hatmu            = as.numeric(hatmu),
    hatphi           = as.numeric(hatphi),
    residuals        = resid,
    pseudo.r.squared = as.numeric(pseudo_r2),
    link             = link,
    link_phi         = link_phi,
    formula          = formula,
    formula_x        = formula,
    formula_z        = ~1,
    terms            = list(mean = mtX, full = mtX),
    xlevels          = list(mean = stats::.getXlevels(mtX, mf)),
    model_matrices   = list(X = X),
    Y                = Y,
    delta            = delta,
    data             = data,
    nobs             = n,
    npar             = length(est),
    p                = p,
    q                = 1L,
    repar            = repar,
    ncuts            = ncuts,
    lim              = lim,
    interval         = interval,
    method           = method,
    optim_method     = method,
    diagnostics      = diagnostics
  )

  class(result) <- "brs"
  invisible(result)
}


#' Fit a variable-dispersion beta interval regression model
#'
#' @description
#' Estimates the parameters of a beta regression model with
#' observation-specific dispersion governed by a second linear
#' predictor.  Both submodels are estimated jointly via maximum
#' likelihood, using the complete likelihood with mixed censoring.
#'
#' @param formula A \code{\link[Formula]{Formula}}-style formula with
#'   two parts: \code{y ~ x1 + x2 | z1 + z2}.
#' @param data   Data frame.
#' @param link   Link for the first parameter (the mean under
#'   \code{repar = 1, 2}; the shape \eqn{p} under \code{repar = 0}).
#'   \code{NULL} (default) selects the link implied by \code{repar}:
#'   \code{"logit"} for \code{repar = 1, 2}, \code{"log"} for
#'   \code{repar = 0}. See the 'Reparameterizations and links' section of
#'   \code{\link{brs}} for the admissible values.
#' @param link_phi Link for the second parameter. \code{NULL} (default)
#'   selects \code{"logit"} for \code{repar = 2} (dispersion on
#'   \eqn{(0, 1)}) and \code{"log"} for \code{repar = 0, 1} (positive
#'   shape/precision).
#' @param hessian_method Character: \code{"cpp"} (default),
#'   \code{"numDeriv"} or \code{"optim"} (see \code{\link{brs_fit_fixed}}).
#' @param ncuts  Number of scale categories. \code{NULL} (default) uses the
#'   value stored by \code{\link{brs_prep}} in \code{attr(data, "ncuts")},
#'   or 100 when \code{data} was not prepared. A value that differs from
#'   the stored one is ignored with a warning (the endpoints were built
#'   with the stored value).
#' @param lim    Half-width of the score cell in \eqn{(0, 0.5]}
#'   (\code{interval = "mid"} only). \code{NULL} (default) uses
#'   \code{attr(data, "lim")} from \code{\link{brs_prep}}, or 0.5; same
#'   rule as \code{ncuts}. Values below 0.5 warn (partial coarsening).
#' @param repar  Reparameterization scheme (default 2); see
#'   \code{\link{brs_repar}}.
#' @param method Optimization method (default \code{"BFGS"}).
#' @param start Optional numeric vector of starting values (mean
#'   coefficients, then dispersion coefficients). \code{NULL} (default) uses
#'   \code{compute_start()}. Refits (bootstrap, jackknife) pass the parent
#'   estimate here.
#' @param control Control list for \code{\link[stats]{optim}}; its entries
#'   are merged into the default \code{list(maxit = 5000L)}.
#' @param interval Direction of the uncertainty interval, \code{"mid"},
#'   \code{"right"} or \code{"left"} (see \code{\link{brs_check}}).
#'   \code{NULL} (default) uses \code{attr(data, "interval")} from
#'   \code{\link{brs_prep}}, or \code{"mid"}; same rule as \code{ncuts}.
#'
#' @return An object of class \code{"brs"}.
#'
#' @examples
#' \donttest{
#' dat <- data.frame(
#'   y = c(
#'     0, 5, 20, 50, 75, 90, 100, 30, 60, 45,
#'     10, 40, 55, 70, 85, 25, 35, 65, 80, 15
#'   ),
#'   x1 = rep(c(1, 2), 10),
#'   x2 = rep(c(0, 0, 1, 1), 5)
#' )
#' prep <- brs_prep(dat, ncuts = 100)
#' fit <- brs_fit_var(y ~ x1 | x2, data = prep)
#' print(fit)
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
#' @importFrom Formula as.Formula Formula
#' @importFrom stats optim cor delete.response
#' @importFrom numDeriv hessian
#' @keywords internal
#' @export
brs_fit_var <- function(formula, data,
                        link = NULL,
                        link_phi = NULL,
                        hessian_method = c("cpp", "numDeriv", "optim"),
                        ncuts = NULL,
                        lim = NULL,
                        repar = 2L,
                        method = c("BFGS", "L-BFGS-B"),
                        interval = NULL,
                        start = NULL,
                        control = list()) {
  cl <- match.call()
  method <- match.arg(method)
  hessian_method <- match.arg(hessian_method)
  validated <- .validate_brs_common_args(data, ncuts, lim, repar, link, link_phi,
                                        interval)
  ncuts <- validated$ncuts
  lim <- validated$lim
  repar <- validated$repar
  link <- validated$link
  link_phi <- validated$link_phi
  interval <- validated$interval

  # Parse multi-part formula
  formula_orig <- formula
  formula <- Formula::as.Formula(formula)
  if (length(formula)[2L] < 2L) {
    formula <- Formula::as.Formula(formula(formula), ~1)
  } else if (length(formula)[2L] > 2L) {
    formula <- Formula::Formula(formula(formula, rhs = 1:2))
  }

  mf <- stats::model.frame(formula, data = data)
  mtX <- stats::terms(formula, data = data, rhs = 1L)
  mtZ <- stats::delete.response(
    stats::terms(formula, data = data, rhs = 2L)
  )
  Y <- .extract_response(mf, data, ncuts = ncuts, lim = lim, interval = interval)
  X <- stats::model.matrix(mtX, mf)
  Z <- stats::model.matrix(mtZ, mf)
  n <- nrow(X)
  p <- ncol(X)
  q <- ncol(Z)
  # Aliased columns stop here instead of a non-finite value inside optim
  .brs_check_design(X, "mean")
  .brs_check_design(Z, "precision")

  # Extract delta from brs_check output
  delta <- as.integer(Y[, "delta"])

  # Starting values: user / parent estimate, else compute_start()
  ini <- if (is.null(start)) {
    compute_start(
      formula = formula, data = data, link = link,
      link_phi = link_phi, ncuts = ncuts,
      lim = lim, repar = repar, interval = interval
    )
  } else {
    .brs_check_start(start, p + q)
  }

  # Link codes
  lc_mu <- link_to_code(link)
  lc_phi <- link_to_code(link_phi)

  # Objective
  fn_obj <- function(par) {
    -.brs_loglik_variable_cpp(
      par, X, Z, Y[, "left"], Y[, "right"], Y[, "yt"],
      delta, lc_mu, lc_phi, repar
    )
  }

  gr_obj <- function(par) {
    -.brs_grad_variable_cpp(
      par, X, Z, Y[, "left"], Y[, "right"], Y[, "yt"],
      delta, lc_mu, lc_phi, repar
    )
  }

  opt <- stats::optim(
    par     = ini,
    fn      = fn_obj,
    gr      = gr_obj,
    method  = method,
    hessian = (hessian_method == "optim"),
    control = .brs_merge_control(list(maxit = 5000L), control)
  )

  # BUG-H04: warn if optimizer did not converge
  if (opt$convergence != 0L) {
    warning(
      "Optimizer did not converge (code ", opt$convergence, ")",
      if (!is.null(opt$message)) paste0(": ", opt$message) else ".",
      "\nResults may be unreliable. Try increasing 'maxit' or changing 'method'.",
      call. = FALSE
    )
  }

  # Hessian
  if (hessian_method == "cpp") {
    opt$hessian <- .brs_hessian_variable_cpp(
      opt$par, X, Z, Y[, "left"], Y[, "right"], Y[, "yt"],
      delta, lc_mu, lc_phi, repar
    )
  } else if (hessian_method == "numDeriv") {
    fn_ll <- function(par) {
      .brs_loglik_variable_cpp(
        par, X, Z, Y[, "left"], Y[, "right"], Y[, "yt"],
        delta, lc_mu, lc_phi, repar
      )
    }
    opt$hessian <- numDeriv::hessian(fn_ll, opt$par)
  } else {
    opt$hessian <- -opt$hessian
  }

  # Fitted values
  est <- opt$par
  idx_beta <- seq_len(p)
  idx_zeta <- p + seq_len(q)

  # `hatmu` is the FIRST parameter (see brs_fit_fixed); means via .brs_mean().
  eta_mu <- X %*% est[idx_beta]
  eta_phi <- Z %*% est[idx_zeta]
  mu_raw <- apply_inv_link(eta_mu, link)
  phi_raw <- apply_inv_link(eta_phi, link_phi)
  # Same clamps as the compiled likelihood (src/brs_common.h).
  hatmu <- .clamp_mu_by_repar(mu_raw, repar)
  hatphi <- .clamp_phi_by_repar(phi_raw, repar)
  # Gradient, Hessian and clamp checks at the estimate (warn + stored)
  diagnostics <- .brs_fit_diagnostics(-gr_obj(est), opt$hessian, mu_raw, phi_raw,
                                      repar, opt$convergence,
                                      badly_scaled = .brs_badly_scaled(X, Z))
  y_mid <- Y[, "yt"]
  ey <- .brs_mean(as.numeric(hatmu), hatphi, repar)
  resid <- as.numeric(y_mid - ey)

  pseudo_r2 <- .brs_pseudo_r2(eta_mu, ey, y_mid, link, repar)
  .warn_sqrt_plateau(eta_mu, link, "link")
  .warn_sqrt_plateau(eta_phi, link_phi, "link_phi")

  # --- betareg-style parameter naming ---
  # Mean coefficients: use column names of X
  mean_names <- colnames(X)
  # Precision coefficients: prefix with "(phi)_"
  phi_names <- paste0("(phi)_", colnames(Z))

  par_names <- c(mean_names, phi_names)
  names(est) <- par_names

  # Named coefficient lists (betareg style)
  coefficients <- list(
    mean      = est[idx_beta],
    precision = est[idx_zeta]
  )
  names(coefficients$mean) <- mean_names
  names(coefficients$precision) <- phi_names

  # Name the hessian
  rownames(opt$hessian) <- colnames(opt$hessian) <- par_names

  # Store formula components
  formula_x <- Formula::as.Formula(formula(formula, rhs = 1L))
  formula_z <- Formula::as.Formula(
    stats::delete.response(stats::terms(formula, data = data, rhs = 2L))
  )

  result <- list(
    call             = cl,
    par              = est,
    coefficients     = coefficients,
    # Exactly the log-likelihood at the returned estimate
    value            = -fn_obj(opt$par),
    hessian          = opt$hessian,
    convergence      = opt$convergence,
    message          = opt$message,
    iterations       = opt$counts,
    hatmu            = as.numeric(hatmu),
    hatphi           = as.numeric(hatphi),
    residuals        = resid,
    pseudo.r.squared = as.numeric(pseudo_r2),
    link             = link,
    link_phi         = link_phi,
    formula          = formula,
    formula_x        = formula_x,
    formula_z        = formula_z,
    terms            = list(mean = mtX, precision = mtZ, full = mtX),
    xlevels          = list(
      mean = stats::.getXlevels(mtX, mf),
      precision = stats::.getXlevels(mtZ, mf)
    ),
    model_matrices   = list(X = X, Z = Z),
    Y                = Y,
    delta            = delta,
    data             = data,
    nobs             = n,
    npar             = length(est),
    p                = p,
    q                = q,
    repar            = repar,
    ncuts            = ncuts,
    lim              = lim,
    interval         = interval,
    method           = method,
    optim_method     = method,
    diagnostics      = diagnostics
  )

  class(result) <- "brs"
  invisible(result)
}


#' Fit a beta interval regression model
#'
#' @description
#' Unified interface that dispatches to \code{\link{brs_fit_fixed}}
#' (fixed dispersion) or \code{\link{brs_fit_var}} (variable
#' dispersion) based on the formula structure.
#'
#' @details
#' If the formula contains a \code{|} separator
#' (e.g., \code{y ~ x1 + x2 | z1}), the variable-dispersion model is
#' fitted; otherwise, a fixed-dispersion model is used.
#'
#' @section Reparameterizations and links:
#' The three schemes of \code{\link{brs_repar}} model different parameters,
#' so the admissible links differ: parameters on \eqn{(0, 1)} use a
#' \code{(0, 1)}-link, parameters on \eqn{(0, \infty)} use \code{"log"} or
#' \code{"sqrt"}. \code{link = NULL} and \code{link_phi = NULL} (the
#' defaults) select the first entry of each cell; any other combination is
#' rejected with an error.
#' \tabular{lll}{
#'   \code{repar} \tab \code{link} (first parameter) \tab
#'     \code{link_phi} (second parameter) \cr
#'   0 (shapes \eqn{p, q}) \tab \code{log}, \code{sqrt} \tab
#'     \code{log}, \code{sqrt} \cr
#'   1 (mean, precision) \tab \code{logit}, \code{probit}, \code{cauchit},
#'     \code{cloglog} \tab \code{log}, \code{sqrt} \cr
#'   2 (mean, dispersion) \tab \code{logit}, \code{probit}, \code{cauchit},
#'     \code{cloglog} \tab \code{logit}, \code{probit}, \code{cauchit},
#'     \code{cloglog}
#' }
#' \code{"identity"}, \code{"inverse"} and \code{"1/mu^2"} are not accepted
#' for positive parameters (their inverse does not map the real line onto
#' \eqn{(0, \infty)}). With \code{"sqrt"} the inverse link is flat for
#' \eqn{\eta \le 0}; a warning is issued after the fit when a fitted linear
#' predictor lies on that plateau.
#'
#' Under \code{repar = 0} the fitted object stores the shape \eqn{p} in
#' \code{hatmu} (and \code{predict(type = "link")} is its linear
#' predictor), while \code{fitted()}, \code{predict(type = "response")},
#' residuals and marginal effects use the mean \eqn{E[Y] = p / (p + q)}.
#'
#' @section Interval direction:
#' \code{interval} selects how a score \eqn{s} is coarsened into a cell of
#' \eqn{(0, 1)}: \code{"mid"} \eqn{[s - \mathrm{lim}, s + \mathrm{lim}] / K}
#' (default), \code{"right"} and \code{"left"} \eqn{[s, s + 1] / (K + 1)}
#' (the dissertation's \eqn{r} and \eqn{l} directions; equal cells, a
#' package normalisation). \code{"right"} and \code{"left"} give the same
#' likelihood and coefficients and differ only in the latent score read
#' back by \code{predict(type = "score")} (one unit). The modes are
#' different coarsening models, so their log-likelihoods are not
#' comparable; see \code{\link{brs_check}} and \code{\link{anova.brs}}.
#'
#' @section Fit diagnostics:
#' Before optimisation, an exactly rank-deficient mean or precision model
#' matrix is an error naming the aliased columns, and a nearly collinear one
#' (condition number of the unit-column matrix above \eqn{10^4}) a warning.
#' After it, \code{fit$diagnostics} stores the largest absolute gradient
#' (\code{grad_norm}), the log-likelihood gain \eqn{\frac12 g^\top (-H)^{-1} g}
#' of the Newton step left to the optimum (\code{grad_gain}) and that step in
#' SE units (\code{grad_step}), the extreme eigenvalues of \eqn{-H} (\code{min_eig},
#' \code{max_eig}, \code{min_eig_scaled} in correlation form,
#' \code{hessian_nd}) and the number of observations whose mean, second
#' parameter or beta shapes sit on the clamps of the likelihood
#' (\code{n_clamped}, \code{clamped}). One-line warnings flag a gradient
#' not \eqn{\approx 0} (gain above 0.01), a Hessian that is not negative
#' definite, and clamped observations (possible non-identifiability, e.g.
#' every observation left-censored). \code{\link{vcov.brs}} then returns
#' \code{NA} for variances it cannot estimate instead of a generalised
#' inverse.
#'
#' @inheritParams brs_fit_var
#'
#' @return An object of class \code{"brs"}; \code{diagnostics} holds the
#'   post-fit checks described above.
#'
#' @examples
#' \donttest{
#' dat <- data.frame(
#'   y = c(
#'     0, 5, 20, 50, 75, 90, 100, 30, 60, 45,
#'     10, 40, 55, 70, 85, 25, 35, 65, 80, 15
#'   ),
#'   x1 = rep(c(1, 2), 10),
#'   x2 = rep(c(0, 0, 1, 1), 5)
#' )
#' prep <- brs_prep(dat, ncuts = 100)
#' # Fixed dispersion
#' fit1 <- brs(y ~ x1, data = prep)
#' print(fit1)
#' # Variable dispersion
#' fit2 <- brs(y ~ x1 | x2, data = prep)
#' print(fit2)
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
#' @importFrom Formula as.Formula Formula
#' @export
brs <- function(formula, data,
                link = NULL,
                link_phi = NULL,
                ncuts = NULL,
                lim = NULL,
                repar = 2L,
                method = c("BFGS", "L-BFGS-B"),
                hessian_method = c("cpp", "numDeriv", "optim"),
                interval = NULL,
                start = NULL,
                control = list()) {
  cl <- match.call()
  formula_parsed <- Formula::as.Formula(formula)

  if (length(formula_parsed)[2L] < 2L) {
    fit <- brs_fit_fixed(
      formula = formula, data = data,
      link = link, link_phi = link_phi,
      ncuts = ncuts, lim = lim,
      hessian_method = hessian_method,
      repar = repar, method = method,
      interval = interval, start = start, control = control
    )
  } else {
    fit <- brs_fit_var(
      formula = formula, data = data,
      link = link, link_phi = link_phi,
      hessian_method = hessian_method,
      ncuts = ncuts, lim = lim,
      repar = repar, method = method,
      interval = interval, start = start, control = control
    )
  }

  # Override the call with the unified interface call
  fit$call <- cl
  fit
}

# Validate a user / parent starting vector (length and finiteness).
.brs_check_start <- function(start, npar) {
  start <- as.numeric(start)
  if (length(start) != npar || any(!is.finite(start))) {
    stop("'start' must be a finite numeric vector of length ", npar, ".", call. = FALSE)
  }
  start
}
