# ============================================================================ #
# Mixed-effects beta interval regression
# ============================================================================ #

#' Fit a mixed-effects beta interval regression model
#'
#' @description
#' Fits a beta interval-censored mixed model with Gaussian random
#' intercepts using marginal maximum likelihood. The implementation supports
#' \code{random = ~ 1 | group} and offers three integration methods for the
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
#' For group \eqn{i}, the random intercept \eqn{b_i \sim N(0, \sigma_b^2)} is
#' integrated out numerically.
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
#' @param random Random-effects specification. The supported format is
#'   \code{~ 1 | group}.
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
#'   (\code{beta}, \code{gamma}, \code{log_sigma_b}).
#' @param method Optimizer passed to \code{\link[stats]{optim}}.
#' @param hessian_method \code{"numDeriv"} (default) or \code{"optim"}.
#' @param control Control list for \code{\link[stats]{optim}}.
#' @param seed Optional seed used by integration methods that depend on
#'   randomized points (reserved for future use).
#'
#' @return An object of class \code{"brsmm"}.
#'
#' @examples
#' \donttest{
#' set.seed(123)
#' g <- 15
#' ni <- 8
#' id <- factor(rep(seq_len(g), each = ni))
#' n <- length(id)
#' x1 <- rnorm(n)
#' b <- rnorm(g, sd = 0.5)
#' eta_mu <- 0.2 + 0.6 * x1 + b[as.integer(id)]
#' mu <- plogis(eta_mu)
#' phi <- plogis(-0.2 + 0.2 * x1)
#' shp <- brs_repar(mu = mu, phi = phi, repar = 2)
#' y <- round(stats::rbeta(n, shp$shape1, shp$shape2) * 100)
#' d <- data.frame(y = y, x1 = x1, id = id)
#'
#' fit_mm <- brsmm(y ~ x1, random = ~ 1 | id, data = d, repar = 2)
#' fit_mm
#' }
#'
#' @references
#' Lopes, J. E. (2024). \emph{Beta Regression for Interval-Censored
#' Scale-Derived Outcomes}. MSc Dissertation, PPGMNE/UFPR.
#'
#' Ferrari, S. and Cribari-Neto, F. (2004). Beta regression for modelling
#' rates and proportions. \emph{Journal of Applied Statistics},
#' \bold{31}(7), 799--815.
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
                  control = list(maxit = 2000L),
                  seed = NULL) {
  cl <- match.call()
  method <- match.arg(method)
  hessian_method <- match.arg(hessian_method)
  link <- match.arg(link, .mu_links)
  link_phi <- match.arg(link_phi, .phi_links)
  int_method <- match.arg(int_method)
  repar <- as.integer(repar)
  n_points <- as.integer(n_points)
  qmc_points <- as.integer(qmc_points)

  if (!is.null(seed)) {
    set.seed(seed)
  }
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

  group_var <- .brsmm_parse_random(random)

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

  group <- .brsmm_extract_group(mf = mf, data = data, group_var = group_var)
  group <- factor(group)
  if (nlevels(group) < 2L) {
    stop("Random intercept requires at least 2 groups.", call. = FALSE)
  }
  group_index <- as.integer(group)

  p <- ncol(X)
  q <- ncol(Z)
  n <- nrow(X)
  g <- nlevels(group)

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
    start <- c(start_fix, log(0.3))
  } else {
    start <- as.numeric(start)
    if (length(start) != (p + q + 1L)) {
      stop(
        "'start' must have length p + q + 1 (beta, gamma, log_sigma_b).",
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
  idx_sd <- p + q + 1L

  beta_hat <- est[idx_beta]
  gamma_hat <- est[idx_gamma]
  sigma_b_hat <- exp(est[idx_sd])

  gm <- .brsmm_group_modes_cpp(
    param = est,
    X = X,
    Z = Z,
    y_left = as.numeric(Y[, "left"]),
    y_right = as.numeric(Y[, "right"]),
    yt = as.numeric(Y[, "yt"]),
    delta = delta,
    group = group_index,
    link_mu_code = lc_mu,
    link_phi_code = lc_phi,
    repar = repar
  )
  mode_b <- as.numeric(gm[, 1L])
  sd_b <- as.numeric(gm[, 2L])

  b_obs <- mode_b[group_index]
  eta_mu <- as.numeric(X %*% beta_hat + b_obs)
  eta_phi <- as.numeric(Z %*% gamma_hat)

  hatmu <- apply_inv_link(eta_mu, link)
  hatphi <- apply_inv_link(eta_phi, link_phi)
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
  sd_name <- paste0("(sd)_", group_var)
  par_names <- c(mean_names, phi_names, sd_name)
  names(est) <- par_names
  rownames(hess) <- colnames(hess) <- par_names

  coefficients <- list(
    mean = est[idx_beta],
    precision = est[idx_gamma],
    random = stats::setNames(est[idx_sd], sd_name)
  )
  levels_group <- levels(group)
  names(mode_b) <- levels_group
  names(sd_b) <- levels_group

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
      mode_b = mode_b,
      sd_b = sd_b,
      sigma_b = sigma_b_hat
    ),
    link = link,
    link_phi = link_phi,
    formula = formula_parsed,
    random_formula = random,
    terms = list(mean = mtX, precision = mtZ, full = mtX),
    model_matrices = list(X = X, Z = Z),
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
    stop("'random' must be a formula like ~ 1 | id.", call. = FALSE)
  }
  if (length(random) != 2L) {
    stop("'random' must be one-sided, e.g. ~ 1 | id.", call. = FALSE)
  }

  rhs <- random[[2L]]
  if (!is.call(rhs) || !identical(rhs[[1L]], as.name("|")) || length(rhs) != 3L) {
    stop("'random' must have the form ~ 1 | group.", call. = FALSE)
  }

  re_part <- rhs[[2L]]
  if (!(identical(re_part, 1) || identical(re_part, 1L) || identical(as.character(re_part), "1"))) {
    stop("Only random intercept is currently supported: use ~ 1 | group.", call. = FALSE)
  }

  group_vars <- all.vars(rhs[[3L]])
  if (length(group_vars) != 1L) {
    stop("'random' must define exactly one grouping variable.", call. = FALSE)
  }
  group_vars[[1L]]
}

#' Extract grouping variable aligned with model.frame rows
#' @keywords internal
#' @noRd
.brsmm_extract_group <- function(mf, data, group_var) {
  if (!(group_var %in% names(data))) {
    stop("Grouping variable '", group_var, "' not found in data.", call. = FALSE)
  }
  rows_num <- suppressWarnings(as.integer(rownames(mf)))
  if (all(!is.na(rows_num))) {
    grp <- data[[group_var]][rows_num]
  } else {
    rows <- match(rownames(mf), rownames(data))
    grp <- data[[group_var]][rows]
  }
  if (anyNA(grp)) {
    stop("Grouping variable contains missing values after subsetting.", call. = FALSE)
  }
  grp
}
