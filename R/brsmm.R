# ============================================================================ #
# Mixed-effects beta interval regression
# ============================================================================ #

#' Fit a mixed-effects beta interval regression model
#'
#' @description
#' Fits a beta interval-censored mixed model with Gaussian random
#' intercepts/slopes using marginal maximum likelihood. The implementation supports
#' random-effects formulas such as \code{~ 1 | group} and \code{~ 1 + x | group},
#' and offers three integration methods for the
#' random effects: Laplace approximation, Adaptive Gauss-Hermite Quadrature
#' (AGHQ), and Quasi-Monte Carlo (QMC).
#'
#' @details
#' The conditional contribution for each observation follows the same mixed
#' censoring likelihood used by \code{\link{brs}}:
#'
#' \enumerate{
#'   \item \eqn{\delta=0}: exact contribution via beta density,
#'   \item \eqn{\delta=1}: left-censored contribution via beta CDF,
#'   \item \eqn{\delta=2}: right-censored contribution via survival CDF,
#'   \item \eqn{\delta=3}: interval contribution via CDF difference.
#' }
#'
#' For group \eqn{i}, the random-effects vector
#' \eqn{\mathbf{b}_i \sim N(\mathbf{0}, D)} is integrated out numerically.
#'
#' \itemize{
#'   \item \code{"laplace"}: Uses a second-order Laplace approximation at the
#'     conditional mode. Fast and generally accurate for \eqn{n_i} large.
#'   \item \code{"aghq"}: Adaptive Gauss-Hermite Quadrature. Uses \code{n_points}
#'     quadrature nodes centered and scaled by the conditional mode and curvature.
#'     More accurate than Laplace, especially for small \eqn{n_i}.
#'   \item \code{"qmc"}: Quasi-Monte Carlo integration using a Halton sequence.
#'     Uses \code{qmc_points} evaluation points. Suitable for high-dimensional
#'     integration (future proofing) or checking robustness.
#' }
#'
#' @param formula Model formula. Supports one- or two-part formulas:
#'   \code{y ~ x1 + x2} or \code{y ~ x1 + x2 | z1 + z2}.
#' @param random Random-effects specification of the form
#'   \code{~ terms | group}, e.g. \code{~ 1 | id} or \code{~ 1 + x | id}.
#' @param data Data frame.
#' @param link Mean link function.
#' @param link_phi Precision link function.
#' @param repar Beta reparameterization code (0, 1, 2).
#' @param ncuts Number of categories on the original scale.
#' @param lim Half-width used to construct interval endpoints.
#' @param int_method Integration method: \code{"laplace"} (default),
#'   \code{"aghq"}, or \code{"qmc"}.
#' @param n_points Number of quadrature points for \code{int_method="aghq"}.
#'   Ignored for other methods. Default is 11.
#' @param qmc_points Number of QMC points for \code{int_method="qmc"}.
#'   Default is 1024.
#' @param start Optional numeric vector of starting values
#'   (\code{beta}, \code{gamma}, and packed lower-Cholesky random parameters).
#' @param method Optimizer passed to \code{\link[stats]{optim}}.
#' @param hessian_method \code{"numDeriv"} (default) or \code{"optim"}.
#' @param control Control list for \code{\link[stats]{optim}}.
#'
#' @return An object of class \code{"brsmm"}.
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
#' fit_mm <- brsmm(y ~ x1, random = ~ 1 | id, data = prep)
#' fit_mm
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
#' @importFrom Formula as.Formula Formula
#' @importFrom stats model.frame terms delete.response model.matrix make.link
#' @importFrom stats cor optim
#' @importFrom numDeriv hessian
#' @export
brsmm <- function(formula,
                  random = ~ 1 | id,
                  data,
                  link = "logit",
                  link_phi = "logit",
                  repar = 2L,
                  ncuts = 100L,
                  lim = 0.5,
                  int_method = c("laplace", "aghq", "qmc"),
                  n_points = 11L,
                  qmc_points = 1024L,
                  start = NULL,
                  method = c("BFGS", "L-BFGS-B"),
                  hessian_method = c("numDeriv", "optim"),
                  control = list(maxit = 2000L)) {
  cl <- match.call()
  method <- match.arg(method)
  hessian_method <- match.arg(hessian_method)
  link <- match.arg(link, .mu_links)
  link_phi <- match.arg(link_phi, .phi_links)
  int_method <- match.arg(int_method)
  repar <- as.integer(repar)
  n_points <- as.integer(n_points)
  qmc_points <- as.integer(qmc_points)

  if (!is.data.frame(data)) {
    stop("'data' must be a data.frame.", call. = FALSE)
  }
  # int_method check removed to support aghq and qmc
  if (!is.finite(n_points) || n_points < 1L) {
    stop("'n_points' must be >= 1.", call. = FALSE)
  }
  if (!is.finite(qmc_points) || qmc_points < 16L) {
    stop("'qmc_points' must be >= 16.", call. = FALSE)
  }

  random_spec <- .brsmm_parse_random(random)
  group_var <- random_spec$group_var

  formula_parsed <- Formula::as.Formula(formula)
  if (length(formula_parsed)[2L] < 2L) {
    formula_parsed <- Formula::as.Formula(formula(formula_parsed), ~1)
  } else if (length(formula_parsed)[2L] > 2L) {
    formula_parsed <- Formula::Formula(formula(formula_parsed, rhs = 1:2))
  }

  mf <- stats::model.frame(formula_parsed, data = data)
  mtX <- stats::terms(formula_parsed, data = data, rhs = 1L)
  mtZ <- stats::delete.response(stats::terms(formula_parsed, data = data, rhs = 2L))

  X <- stats::model.matrix(mtX, mf)
  Z <- stats::model.matrix(mtZ, mf)
  Y <- .extract_response(mf, data, ncuts = ncuts, lim = lim)
  delta <- as.integer(Y[, "delta"])

  rows_idx <- .brsmm_row_index(mf = mf, data = data)
  data_sub <- data[rows_idx, , drop = FALSE]

  group <- .brsmm_extract_group(mf = mf, data = data, group_var = group_var)
  group <- factor(group)
  if (nlevels(group) < 2L) {
    stop("Random intercept requires at least 2 groups.", call. = FALSE)
  }
  group_index <- as.integer(group)

  p <- ncol(X)
  q <- ncol(Z)
  Xr <- stats::model.matrix(random_spec$re_terms, data_sub)
  q_re <- ncol(Xr)
  k_re <- q_re * (q_re + 1L) / 2L
  n <- nrow(X)
  g <- nlevels(group)

  if (q_re > 1L && int_method != "laplace") {
    stop(
      "For random-effects dimension > 1, only int_method = 'laplace' is currently supported.",
      call. = FALSE
    )
  }

  if (is.null(start)) {
    start_fix <- compute_start(
      formula = formula_parsed,
      data = data,
      link = link,
      link_phi = link_phi,
      ncuts = ncuts,
      lim = lim
    )
    if (length(start_fix) != (p + q)) {
      stop(
        "Internal error: starting vector from compute_start() has unexpected length.",
        call. = FALSE
      )
    }
    theta_re_start <- numeric(k_re)
    k <- 1L
    for (j in seq_len(q_re)) {
      for (i in j:q_re) {
        theta_re_start[k] <- if (i == j) log(0.3) else 0
        k <- k + 1L
      }
    }
    start <- c(start_fix, theta_re_start)
  } else {
    start <- as.numeric(start)
    if (length(start) != (p + q + k_re)) {
      stop(
        "'start' must have length p + q + q_re * (q_re + 1) / 2.",
        call. = FALSE
      )
    }
  }

  lc_mu <- link_to_code(link)
  lc_phi <- link_to_code(link_phi)

  # Map string method to integer code
  method_code <- match(int_method, c("laplace", "aghq", "qmc")) - 1L

  # Determine number of points
  n_pts <- if (int_method == "qmc") qmc_points else n_points

  fn_ll <- function(par) {
    # Call the new Eigen-based backend
    brsmm_loglik_eigen(
      param = as.numeric(par),
      X = X,
      Z = Z,
      Xr = Xr,
      y_left = as.numeric(Y[, "left"]),
      y_right = as.numeric(Y[, "right"]),
      yt = as.numeric(Y[, "yt"]),
      delta = delta,
      group = group_index,
      link_mu = lc_mu,
      link_phi = lc_phi,
      repar = repar,
      method = method_code,
      n_points = n_pts
    )
  }

  fn_obj <- function(par) -fn_ll(par)

  opt <- stats::optim(
    par = start,
    fn = fn_obj,
    method = method,
    hessian = (hessian_method == "optim"),
    control = control
  )

  if (hessian_method == "numDeriv") {
    hess <- numDeriv::hessian(fn_ll, opt$par)
  } else {
    hess <- -opt$hessian
  }

  est <- as.numeric(opt$par)
  idx_beta <- seq_len(p)
  idx_gamma <- p + seq_len(q)
  idx_re <- p + q + seq_len(k_re)

  beta_hat <- est[idx_beta]
  gamma_hat <- est[idx_gamma]
  theta_re_hat <- est[idx_re]

  L <- matrix(0, nrow = q_re, ncol = q_re)
  k <- 1L
  for (j in seq_len(q_re)) {
    for (i in j:q_re) {
      L[i, j] <- if (i == j) exp(theta_re_hat[k]) else theta_re_hat[k]
      k <- k + 1L
    }
  }
  D <- L %*% t(L)
  sd_b_terms <- sqrt(diag(D))

  gm <- brsmm_group_modes_eigen(
    param = est,
    X = X,
    Z = Z,
    Xr = Xr,
    y_left = as.numeric(Y[, "left"]),
    y_right = as.numeric(Y[, "right"]),
    yt = as.numeric(Y[, "yt"]),
    delta = delta,
    group = group_index,
    link_mu = lc_mu,
    link_phi = lc_phi,
    repar = repar
  )
  mode_b <- as.matrix(gm)
  if (ncol(mode_b) != q_re) {
    stop("Internal error while computing group modes.", call. = FALSE)
  }
  eta_phi <- as.numeric(Z %*% gamma_hat)
  y_mid <- as.numeric(Y[, "yt"])

  pseudo_r2 <- suppressWarnings(stats::cor(
    as.numeric(X %*% beta_hat),
    stats::make.link(link)$linkfun(y_mid)
  )^2)
  if (!is.finite(pseudo_r2)) {
    pseudo_r2 <- NA_real_
  }

  mean_names <- colnames(X)
  phi_names <- paste0("(phi)_", colnames(Z))
  re_colnames <- colnames(Xr)
  re_colnames[is.na(re_colnames) | re_colnames == ""] <- paste0("re", seq_len(q_re))
  re_param_names <- character(k_re)
  k <- 1L
  for (j in seq_len(q_re)) {
    for (i in j:q_re) {
      re_param_names[k] <- if (i == j) {
        paste0("(re_chol_logsd)_", re_colnames[i], "|", group_var)
      } else {
        paste0("(re_chol)_", re_colnames[i], ":", re_colnames[j], "|", group_var)
      }
      k <- k + 1L
    }
  }
  par_names <- c(mean_names, phi_names, re_param_names)
  names(est) <- par_names
  rownames(hess) <- colnames(hess) <- par_names

  coefficients <- list(
    mean = est[idx_beta],
    precision = est[idx_gamma],
    random = stats::setNames(est[idx_re], re_param_names)
  )
  levels_group <- levels(group)
  rownames(mode_b) <- levels_group
  colnames(mode_b) <- re_colnames
  names(sd_b_terms) <- re_colnames

  if (q_re == 1L) {
    mode_store <- as.numeric(mode_b[, 1L])
    names(mode_store) <- levels_group
    b_obs <- mode_store[group_index]
    eta_mu <- as.numeric(X %*% beta_hat + Xr[, 1L] * b_obs)
    sigma_b_hat <- as.numeric(sd_b_terms[1L])
  } else {
    mode_store <- mode_b
    b_obs <- mode_b[group_index, , drop = FALSE]
    eta_mu <- as.numeric(X %*% beta_hat + rowSums(Xr * b_obs))
    sigma_b_hat <- NA_real_
  }

  hatmu <- apply_inv_link(eta_mu, link)
  hatphi <- apply_inv_link(eta_phi, link_phi)

  out <- list(
    call = cl,
    par = est,
    coefficients = coefficients,
    value = -opt$value,
    hessian = hess,
    convergence = opt$convergence,
    message = opt$message,
    iterations = opt$counts,
    fitted_mu = as.numeric(hatmu),
    fitted_phi = as.numeric(hatphi),
    residuals = as.numeric(y_mid - hatmu),
    pseudo.r.squared = pseudo_r2,
    random = list(
      group = group_var,
      levels = levels_group,
      terms = re_colnames,
      re_terms = random_spec$re_terms,
      mode_b = mode_store,
      sd_b = sd_b_terms,
      D = D,
      L = L,
      sigma_b = sigma_b_hat
    ),
    link = link,
    link_phi = link_phi,
    formula = formula_parsed,
    random_formula = random,
    terms = list(mean = mtX, precision = mtZ, full = mtX),
    model_matrices = list(X = X, Z = Z, Xr = Xr),
    Y = Y,
    delta = delta,
    group = group,
    group_index = group_index,
    data = data,
    nobs = n,
    ngroups = g,
    npar = length(est),
    p = p,
    q = q,
    q_re = q_re,
    k_re = k_re,
    repar = repar,
    ncuts = ncuts,
    lim = lim,
    method = method,
    hessian_method = hessian_method,
    int_method = int_method,
    n_points = n_points,
    qmc_points = qmc_points,
    diagnostics = list(
      integration = list(
        method = int_method,
        n_groups = g
      )
    )
  )

  class(out) <- "brsmm"
  out
}

#' Parse random-effect specification for brsmm
#' @keywords internal
#' @noRd
.brsmm_parse_random <- function(random) {
  if (!inherits(random, "formula")) {
    stop("'random' must be a formula like ~ 1 | id or ~ 1 + x | id.", call. = FALSE)
  }
  if (length(random) != 2L) {
    stop("'random' must be one-sided, e.g. ~ 1 | id.", call. = FALSE)
  }

  rhs <- random[[2L]]
  if (!is.call(rhs) || !identical(rhs[[1L]], as.name("|")) || length(rhs) != 3L) {
    stop("'random' must have the form ~ terms | group.", call. = FALSE)
  }

  re_part <- rhs[[2L]]
  group_vars <- all.vars(rhs[[3L]])
  if (length(group_vars) != 1L) {
    stop("'random' must define exactly one grouping variable.", call. = FALSE)
  }
  re_formula <- stats::as.formula(paste("~", deparse(re_part)))
  re_terms <- stats::terms(re_formula)
  list(
    group_var = group_vars[[1L]],
    re_formula = re_formula,
    re_terms = re_terms
  )
}

#' Row index from model.frame to data
#' @keywords internal
#' @noRd
.brsmm_row_index <- function(mf, data) {
  rows_num <- suppressWarnings(as.integer(rownames(mf)))
  if (all(!is.na(rows_num))) {
    return(rows_num)
  }
  rows <- match(rownames(mf), rownames(data))
  if (anyNA(rows)) {
    stop("Could not map model.frame rows back to 'data'.", call. = FALSE)
  }
  rows
}

#' Extract grouping variable aligned with model.frame rows
#' @keywords internal
#' @noRd
.brsmm_extract_group <- function(mf, data, group_var) {
  if (!(group_var %in% names(data))) {
    stop("Grouping variable '", group_var, "' not found in data.", call. = FALSE)
  }
  rows <- .brsmm_row_index(mf, data)
  grp <- data[[group_var]][rows]
  if (anyNA(grp)) {
    stop("Grouping variable contains missing values after subsetting.", call. = FALSE)
  }
  grp
}
