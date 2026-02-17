# ============================================================================ #
# Model-fitting functions
# ============================================================================ #

.validate_brs_common_args <- function(data, ncuts, lim, repar) {
  if (!is.data.frame(data)) {
    stop("`data` must be a data.frame.", call. = FALSE)
  }
  ncuts <- as.integer(ncuts)
  if (!is.finite(ncuts) || ncuts < 2L) {
    stop("`ncuts` must be an integer >= 2.", call. = FALSE)
  }
  if (!is.numeric(lim) || length(lim) != 1L || !is.finite(lim) || lim <= 0) {
    stop("`lim` must be a positive finite scalar.", call. = FALSE)
  }
  repar <- as.integer(repar)
  if (!(repar %in% 0:2)) {
    stop("`repar` must be one of 0, 1, or 2.", call. = FALSE)
  }
  list(ncuts = ncuts, lim = as.numeric(lim), repar = repar)
}

#' Fit a fixed-dispersion beta interval regression model
#'
#' @description
#' Estimates the parameters of a beta regression model with a single
#' (scalar) dispersion parameter using maximum likelihood.  The
#' log-likelihood and its gradient are evaluated by the compiled C++
#' backend supporting the complete likelihood with mixed censoring
#' types (Lopes, 2024, Eq. 2.24).
#'
#' @param formula Two-sided formula \code{y ~ x1 + x2 + ...}.
#' @param data   Data frame.
#' @param link   Mean link function (default \code{"logit"}).
#' @param link_phi Dispersion link function (default \code{"logit"}).
#' @param ncuts  Number of scale categories (default 100).
#' @param lim    Uncertainty half-width (default 0.5).
#' @param hessian_method Character: \code{"numDeriv"} (default) or
#'   \code{"optim"}.  With \code{"numDeriv"} the Hessian is computed
#'   after convergence using \code{\link[numDeriv]{hessian}}, which is
#'   typically more accurate than the built-in optim Hessian.
#' @param repar  Reparameterization scheme (default 2).
#' @param method Optimization method: \code{"BFGS"} (default) or
#'   \code{"L-BFGS-B"}.
#'
#' @return An object of class \code{"brs"}.
#'
#' @examples
#' set.seed(42)
#' n <- 100
#' dat <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
#' sim <- brs_sim(
#'   formula = ~ x1 + x2, data = dat,
#'   beta = c(0.2, -0.5, 0.3), phi = 1 / 5
#' )
#' fit <- brs_fit_fixed(
#'   formula = y ~ x1 + x2, data = sim,
#'   link = "logit", link_phi = "logit"
#' )
#' print(fit)
#'
#' @references
#' Hawker, G. A., Mian, S., Kendzerska, T., and French, M. (2011).
#' Measures of adult pain: Visual Analog Scale for Pain (VAS Pain),
#' Numeric Rating Scale for Pain (NRS Pain), McGill Pain Questionnaire (MPQ),
#' Short-Form McGill Pain Questionnaire (SF-MPQ), Chronic Pain Grade Scale
#' (CPGS), Short Form-36 Bodily Pain Scale (SF-36 BPS), and Measure of
#' Intermittent and Constant Osteoarthritis Pain (ICOAP).
#' Arthritis Care and Research, 63(S11), S240-S252.
#' doi:10.1002/acr.20543.
#'
#' Hjermstad, M. J., Fayers, P. M., Haugen, D. F., et al. (2011).
#' Studies comparing Numerical Rating Scales, Verbal Rating Scales, and
#' Visual Analogue Scales for assessment of pain intensity in adults:
#' a systematic literature review.
#' Journal of Pain and Symptom Management, 41(6), 1073-1093.
#' doi:10.1016/j.jpainsymman.2010.08.016.
#'
#' @importFrom stats optim cor model.frame model.matrix model.response terms
#' @importFrom stats make.link
#' @importFrom numDeriv hessian
#' @keywords internal
#' @export
brs_fit_fixed <- function(formula, data,
                          link = "logit",
                          link_phi = "logit",
                          ncuts = 100L,
                          lim = 0.5,
                          hessian_method = c("numDeriv", "optim"),
                          repar = 2L,
                          method = c("BFGS", "L-BFGS-B")) {
  cl <- match.call()
  method <- match.arg(method)
  hessian_method <- match.arg(hessian_method)
  link <- match.arg(link, .mu_links)
  link_phi <- match.arg(link_phi, .phi_links)
  validated <- .validate_brs_common_args(data, ncuts, lim, repar)
  ncuts <- validated$ncuts
  lim <- validated$lim
  repar <- validated$repar

  # Build matrices
  mf <- stats::model.frame(formula, data = data)
  mtX <- stats::terms(formula, data = data, rhs = 1L)
  Y <- .extract_response(mf, data, ncuts = ncuts, lim = lim)
  X <- stats::model.matrix(mtX, mf)
  n <- nrow(X)
  p <- ncol(X)

  # Extract delta from brs_check output
  delta <- as.integer(Y[, "delta"])

  # Starting values
  ini <- compute_start(
    formula = formula, data = data, link = link,
    link_phi = link_phi, ncuts = ncuts,
    lim = lim
  )

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
    control = list(maxit = 5000L)
  )

  # Hessian (on the log-likelihood scale)
  if (hessian_method == "numDeriv") {
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

  # Fitted values
  est <- opt$par
  hatmu <- apply_inv_link(X %*% est[1:p], link)
  hatphi <- apply_inv_link(est[p + 1L], link_phi)
  y_mid <- Y[, "yt"]
  resid <- as.numeric(y_mid - hatmu)

  pseudo_r2 <- stats::cor(
    X %*% est[1:p],
    stats::make.link(link)$linkfun(y_mid)
  )^2

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
    value            = -opt$value,
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
    method           = method,
    optim_method     = method
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
#' likelihood, using the complete likelihood with mixed censoring
#' (Lopes, 2024, Eq. 2.24).
#'
#' @param formula A \code{\link[Formula]{Formula}}-style formula with
#'   two parts: \code{y ~ x1 + x2 | z1 + z2}.
#' @param data   Data frame.
#' @param link   Mean link function (default \code{"logit"}).
#' @param link_phi Dispersion link function (default \code{"logit"}).
#' @param hessian_method Character: \code{"numDeriv"} or
#'   \code{"optim"}.
#' @param ncuts  Number of scale categories (default 100).
#' @param lim    Uncertainty half-width (default 0.5).
#' @param repar  Reparameterization scheme (default 2).
#' @param method Optimization method (default \code{"BFGS"}).
#'
#' @return An object of class \code{"brs"}.
#'
#' @examples
#' set.seed(42)
#' n <- 100
#' dat <- data.frame(
#'   x1 = rnorm(n), x2 = rnorm(n),
#'   z1 = runif(n)
#' )
#' sim <- brs_sim(
#'   formula = ~ x1 + x2 | z1,
#'   data = dat,
#'   beta = c(0.2, -0.5, 0.3),
#'   zeta = c(1, 1.2)
#' )
#' fit <- brs_fit_var(
#'   formula = y ~ x1 + x2 | z1, data = sim,
#'   link = "logit", link_phi = "logit"
#' )
#' print(fit)
#'
#' @references
#' Hawker, G. A., Mian, S., Kendzerska, T., and French, M. (2011).
#' Measures of adult pain: Visual Analog Scale for Pain (VAS Pain),
#' Numeric Rating Scale for Pain (NRS Pain), McGill Pain Questionnaire (MPQ),
#' Short-Form McGill Pain Questionnaire (SF-MPQ), Chronic Pain Grade Scale
#' (CPGS), Short Form-36 Bodily Pain Scale (SF-36 BPS), and Measure of
#' Intermittent and Constant Osteoarthritis Pain (ICOAP).
#' Arthritis Care and Research, 63(S11), S240-S252.
#' doi:10.1002/acr.20543.
#'
#' Hjermstad, M. J., Fayers, P. M., Haugen, D. F., et al. (2011).
#' Studies comparing Numerical Rating Scales, Verbal Rating Scales, and
#' Visual Analogue Scales for assessment of pain intensity in adults:
#' a systematic literature review.
#' Journal of Pain and Symptom Management, 41(6), 1073-1093.
#' doi:10.1016/j.jpainsymman.2010.08.016.
#'
#' @importFrom Formula as.Formula Formula
#' @importFrom stats optim cor make.link delete.response
#' @importFrom numDeriv hessian
#' @keywords internal
#' @export
brs_fit_var <- function(formula, data,
                        link = "logit",
                        link_phi = "logit",
                        hessian_method = c("numDeriv", "optim"),
                        ncuts = 100L,
                        lim = 0.5,
                        repar = 2L,
                        method = c("BFGS", "L-BFGS-B")) {
  cl <- match.call()
  method <- match.arg(method)
  hessian_method <- match.arg(hessian_method)
  link <- match.arg(link, .mu_links)
  link_phi <- match.arg(link_phi, .phi_links)
  validated <- .validate_brs_common_args(data, ncuts, lim, repar)
  ncuts <- validated$ncuts
  lim <- validated$lim
  repar <- validated$repar

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
  Y <- .extract_response(mf, data, ncuts = ncuts, lim = lim)
  X <- stats::model.matrix(mtX, mf)
  Z <- stats::model.matrix(mtZ, mf)
  n <- nrow(X)
  p <- ncol(X)
  q <- ncol(Z)

  # Extract delta from brs_check output
  delta <- as.integer(Y[, "delta"])

  # Starting values
  ini <- compute_start(
    formula = formula, data = data, link = link,
    link_phi = link_phi, ncuts = ncuts,
    lim = lim
  )

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
    control = list(maxit = 5000L)
  )

  # Hessian
  if (hessian_method == "numDeriv") {
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

  hatmu <- apply_inv_link(X %*% est[idx_beta], link)
  hatphi <- apply_inv_link(Z %*% est[idx_zeta], link_phi)
  y_mid <- Y[, "yt"]
  resid <- as.numeric(y_mid - hatmu)

  pseudo_r2 <- stats::cor(
    X %*% est[idx_beta],
    stats::make.link(link)$linkfun(y_mid)
  )^2

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
    value            = -opt$value,
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
    method           = method,
    optim_method     = method
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
#' @inheritParams brs_fit_var
#'
#' @return An object of class \code{"brs"}.
#'
#' @examples
#' set.seed(42)
#' n <- 100
#' dat <- data.frame(
#'   x1 = rnorm(n), x2 = rnorm(n),
#'   z1 = runif(n)
#' )
#' sim <- brs_sim(
#'   formula = ~ x1 + x2 | z1,
#'   data = dat,
#'   beta = c(0.2, -0.5, 0.3),
#'   zeta = c(1, 1.2)
#' )
#'
#' # Fixed dispersion
#' fit1 <- brs(y ~ x1 + x2, data = sim)
#' print(fit1)
#'
#' # Variable dispersion
#' fit2 <- brs(y ~ x1 + x2 | z1, data = sim)
#' print(fit2)
#'
#' @references
#' Hawker, G. A., Mian, S., Kendzerska, T., and French, M. (2011).
#' Measures of adult pain: Visual Analog Scale for Pain (VAS Pain),
#' Numeric Rating Scale for Pain (NRS Pain), McGill Pain Questionnaire (MPQ),
#' Short-Form McGill Pain Questionnaire (SF-MPQ), Chronic Pain Grade Scale
#' (CPGS), Short Form-36 Bodily Pain Scale (SF-36 BPS), and Measure of
#' Intermittent and Constant Osteoarthritis Pain (ICOAP).
#' Arthritis Care and Research, 63(S11), S240-S252.
#' doi:10.1002/acr.20543.
#'
#' Hjermstad, M. J., Fayers, P. M., Haugen, D. F., et al. (2011).
#' Studies comparing Numerical Rating Scales, Verbal Rating Scales, and
#' Visual Analogue Scales for assessment of pain intensity in adults:
#' a systematic literature review.
#' Journal of Pain and Symptom Management, 41(6), 1073-1093.
#' doi:10.1016/j.jpainsymman.2010.08.016.
#'
#' @importFrom Formula as.Formula Formula
#' @export
brs <- function(formula, data,
                link = "logit",
                link_phi = "logit",
                ncuts = 100L,
                lim = 0.5,
                repar = 2L,
                method = c("BFGS", "L-BFGS-B"),
                hessian_method = c("numDeriv", "optim")) {
  cl <- match.call()
  formula_parsed <- Formula::as.Formula(formula)

  if (length(formula_parsed)[2L] < 2L) {
    fit <- brs_fit_fixed(
      formula = formula, data = data,
      link = link, link_phi = link_phi,
      ncuts = ncuts, lim = lim,
      hessian_method = hessian_method,
      repar = repar, method = method
    )
  } else {
    fit <- brs_fit_var(
      formula = formula, data = data,
      link = link, link_phi = link_phi,
      hessian_method = hessian_method,
      ncuts = ncuts, lim = lim,
      repar = repar, method = method
    )
  }

  # Override the call with the unified interface call
  fit$call <- cl
  fit
}
