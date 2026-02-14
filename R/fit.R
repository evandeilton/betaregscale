# ============================================================================ #
# Model-fitting functions
# ============================================================================ #

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
#' @param type   \strong{Deprecated.}
#'   Interval type (default \code{"m"}).
#'   Use \code{\link{bs_prepare}} to control interval geometry instead.
#' @param lim    Uncertainty half-width (default 0.5).
#' @param hessian_method Character: \code{"numDeriv"} (default) or
#'   \code{"optim"}.  With \code{"numDeriv"} the Hessian is computed
#'   after convergence using \code{\link[numDeriv]{hessian}}, which is
#'   typically more accurate than the built-in optim Hessian.
#' @param repar  Reparameterization scheme (default 2).
#' @param method Optimization method: \code{"BFGS"} (default) or
#'   \code{"L-BFGS-B"}.
#'
#' @return An object of class \code{"betaregscale"}.
#'
#' @examples
#' set.seed(42)
#' n <- 100
#' dat <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
#' sim <- betaregscale_simulate(
#'   formula = ~ x1 + x2, data = dat,
#'   beta = c(0.2, -0.5, 0.3), phi = 1 / 5
#' )
#' fit <- betaregscale_fit(
#'   formula = y ~ x1 + x2, data = sim,
#'   link = "logit", link_phi = "logit"
#' )
#' print(fit)
#'
#' @importFrom stats optim cor model.frame model.matrix model.response terms
#' @importFrom stats make.link
#' @importFrom numDeriv hessian
#' @export
betaregscale_fit <- function(formula, data,
                             link = "logit",
                             link_phi = "logit",
                             ncuts = 100L,
                             type = "m",
                             lim = 0.5,
                             hessian_method = c("numDeriv", "optim"),
                             repar = 2L,
                             method = c("BFGS", "L-BFGS-B")) {
  cl <- match.call()
  if (!missing(type)) {
    .Deprecated(msg = paste0(
      "The 'type' argument of betaregscale_fit() is deprecated ",
      "and will be removed in a future version. ",
      "Use bs_prepare() to control interval geometry."
    ))
  }
  method <- match.arg(method)
  hessian_method <- match.arg(hessian_method)
  link <- match.arg(link, .mu_links)
  link_phi <- match.arg(link_phi, .phi_links)
  repar <- as.integer(repar)

  # Build matrices
  mf <- stats::model.frame(formula, data = data)
  mtX <- stats::terms(formula, data = data, rhs = 1L)
  Y <- .extract_response(mf, data, ncuts = ncuts, type = type, lim = lim)
  X <- stats::model.matrix(mtX, mf)
  n <- nrow(X)
  p <- ncol(X)

  # Extract delta from check_response output
  delta <- as.integer(Y[, "delta"])

  # Starting values
  ini <- compute_start(
    formula = formula, data = data, link = link,
    link_phi = link_phi, ncuts = ncuts, type = type,
    lim = lim
  )

  # Pre-compute link codes for C++
  lc_mu <- link_to_code(link)
  lc_phi <- link_to_code(link_phi)

  # Objective: -loglik (we minimize)
  fn_obj <- function(par) {
    -.betaregscale_loglik_fixed_cpp(
      par, X, Y[, "left"], Y[, "right"], Y[, "yt"],
      delta, lc_mu, lc_phi, repar
    )
  }

  # Gradient: -grad
  gr_obj <- function(par) {
    -.betaregscale_grad_fixed_cpp(
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
      .betaregscale_loglik_fixed_cpp(
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
    type             = type,
    lim              = lim,
    method           = method,
    optim_method     = method
  )

  class(result) <- "betaregscale"
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
#' @param type   \strong{Deprecated.}
#'   Interval type (default \code{"m"}).
#'   Use \code{\link{bs_prepare}} to control interval geometry instead.
#' @param lim    Uncertainty half-width (default 0.5).
#' @param repar  Reparameterization scheme (default 2).
#' @param method Optimization method (default \code{"BFGS"}).
#'
#' @return An object of class \code{"betaregscale"}.
#'
#' @examples
#' set.seed(42)
#' n <- 100
#' dat <- data.frame(
#'   x1 = rnorm(n), x2 = rnorm(n),
#'   z1 = runif(n)
#' )
#' sim <- betaregscale_simulate_z(
#'   formula_x = ~ x1 + x2, formula_z = ~z1,
#'   data = dat,
#'   beta = c(0.2, -0.5, 0.3),
#'   zeta = c(1, 1.2)
#' )
#' fit <- betaregscale_fit_z(
#'   formula = y ~ x1 + x2 | z1, data = sim,
#'   link = "logit", link_phi = "logit"
#' )
#' print(fit)
#'
#' @importFrom Formula as.Formula Formula
#' @importFrom stats optim cor make.link delete.response
#' @importFrom numDeriv hessian
#' @export
betaregscale_fit_z <- function(formula, data,
                               link = "logit",
                               link_phi = "logit",
                               hessian_method = c("numDeriv", "optim"),
                               ncuts = 100L,
                               type = "m",
                               lim = 0.5,
                               repar = 2L,
                               method = c("BFGS", "L-BFGS-B")) {
  cl <- match.call()
  if (!missing(type)) {
    .Deprecated(msg = paste0(
      "The 'type' argument of betaregscale_fit_z() is deprecated ",
      "and will be removed in a future version. ",
      "Use bs_prepare() to control interval geometry."
    ))
  }
  method <- match.arg(method)
  hessian_method <- match.arg(hessian_method)
  link <- match.arg(link, .mu_links)
  link_phi <- match.arg(link_phi, .phi_links)
  repar <- as.integer(repar)

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
  Y <- .extract_response(mf, data, ncuts = ncuts, type = type, lim = lim)
  X <- stats::model.matrix(mtX, mf)
  Z <- stats::model.matrix(mtZ, mf)
  n <- nrow(X)
  p <- ncol(X)
  q <- ncol(Z)

  # Extract delta from check_response output
  delta <- as.integer(Y[, "delta"])

  # Starting values
  ini <- compute_start(
    formula = formula, data = data, link = link,
    link_phi = link_phi, ncuts = ncuts, type = type,
    lim = lim
  )

  # Link codes
  lc_mu <- link_to_code(link)
  lc_phi <- link_to_code(link_phi)

  # Objective
  fn_obj <- function(par) {
    -.betaregscale_loglik_variable_cpp(
      par, X, Z, Y[, "left"], Y[, "right"], Y[, "yt"],
      delta, lc_mu, lc_phi, repar
    )
  }

  gr_obj <- function(par) {
    -.betaregscale_grad_variable_cpp(
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
      .betaregscale_loglik_variable_cpp(
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
    type             = type,
    lim              = lim,
    method           = method,
    optim_method     = method
  )

  class(result) <- "betaregscale"
  invisible(result)
}


#' Fit a beta interval regression model
#'
#' @description
#' Unified interface that dispatches to \code{\link{betaregscale_fit}}
#' (fixed dispersion) or \code{\link{betaregscale_fit_z}} (variable
#' dispersion) based on the formula structure.
#'
#' @details
#' If the formula contains a \code{|} separator
#' (e.g., \code{y ~ x1 + x2 | z1}), the variable-dispersion model is
#' fitted; otherwise, a fixed-dispersion model is used.
#'
#' @inheritParams betaregscale_fit_z
#'
#' @return An object of class \code{"betaregscale"}.
#'
#' @examples
#' set.seed(42)
#' n <- 100
#' dat <- data.frame(
#'   x1 = rnorm(n), x2 = rnorm(n),
#'   z1 = runif(n)
#' )
#' sim <- betaregscale_simulate_z(
#'   formula_x = ~ x1 + x2, formula_z = ~z1,
#'   data = dat,
#'   beta = c(0.2, -0.5, 0.3),
#'   zeta = c(1, 1.2)
#' )
#'
#' # Fixed dispersion
#' fit1 <- betaregscale(y ~ x1 + x2, data = sim)
#' print(fit1)
#'
#' # Variable dispersion
#' fit2 <- betaregscale(y ~ x1 + x2 | z1, data = sim)
#' print(fit2)
#'
#' @importFrom Formula as.Formula Formula
#' @export
betaregscale <- function(formula, data,
                         link = "logit",
                         link_phi = "logit",
                         ncuts = 100L,
                         type = "m",
                         lim = 0.5,
                         repar = 2L,
                         method = c("BFGS", "L-BFGS-B"),
                         hessian_method = c("numDeriv", "optim")) {
  cl <- match.call()
  if (!missing(type)) {
    .Deprecated(msg = paste0(
      "The 'type' argument of betaregscale() is deprecated ",
      "and will be removed in a future version. ",
      "Use bs_prepare() to control interval geometry."
    ))
  }
  formula_parsed <- Formula::as.Formula(formula)

  if (length(formula_parsed)[2L] < 2L) {
    fit <- betaregscale_fit(
      formula = formula, data = data,
      link = link, link_phi = link_phi,
      ncuts = ncuts, type = type, lim = lim,
      hessian_method = hessian_method,
      repar = repar, method = method
    )
  } else {
    fit <- betaregscale_fit_z(
      formula = formula, data = data,
      link = link, link_phi = link_phi,
      hessian_method = hessian_method,
      ncuts = ncuts, type = type, lim = lim,
      repar = repar, method = method
    )
  }

  # Override the call with the unified interface call
  fit$call <- cl
  fit
}
