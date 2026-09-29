# ============================================================================ #
# Mixed-effects beta interval regression
# ============================================================================ #

#' Fit a mixed-effects beta interval regression model
#'
#' @description
#' Beta interval regression (\code{\link{brs}}) with Gaussian random effects
#' in the linear predictor of the first parameter (the mean, or the shape
#' \eqn{p} under \code{repar = 0}), fitted by marginal maximum likelihood.
#' \code{random = ~ 1 | id} gives a random intercept per group,
#' \code{~ 1 + x | id} a random intercept and slope with a free correlation.
#'
#' @details
#' For group \eqn{g} with observations \eqn{i}, the model is
#' \deqn{g_1(\mu_{gi}) = x_{gi}^\top \beta + x_{r,gi}^\top b_g, \qquad
#'   g_2(\phi_{gi}) = z_{gi}^\top \gamma, \qquad b_g \sim N(0, D),}
#' with the censored contributions of \code{\link{brs}}. The marginal
#' log-likelihood is \eqn{\sum_g \log \int \exp\{h_g(b)\}\, db}, where
#' \eqn{h_g(b) = \sum_i \ell_{gi}(b) + \log \varphi(b; 0, D)}. The integral
#' is computed around the mode \eqn{\hat b_g} of \eqn{h_g}, with
#' \eqn{H_g = -\partial^2 h_g / \partial b\, \partial b^\top} there:
#' \describe{
#'   \item{\code{"laplace"}}{\eqn{h_g(\hat b_g) + \frac{q}{2}\log(2\pi) -
#'     \frac12 \log |H_g|}; fast, accurate when groups are not tiny.}
#'   \item{\code{"aghq"}}{adaptive Gauss-Hermite quadrature: a product grid of
#'     \code{n_points} nodes per dimension at
#'     \eqn{\hat b_g + \sqrt{2}\, H_g^{-1/2} z}, with the symmetric square
#'     root \eqn{H_g^{-1/2}} (\code{n_points}\eqn{^q} nodes, at most
#'     500000).}
#'   \item{\code{"qmc"}}{importance sampling from \eqn{N(\hat b_g, H_g^{-1})}
#'     (nodes \eqn{\hat b_g + H_g^{-1/2} z}) on \code{qmc_points} Halton
#'     points. It is deterministic and, with two or more random effects,
#'     underestimates the log-likelihood (about \eqn{-0.05} at 1024 points
#'     in the package's two-effect checks); prefer \code{"aghq"} for up to
#'     three random effects.}
#' }
#' The inner mode is found by a Levenberg--Marquardt Newton method,
#' warm-started from the modes of the previous evaluation (the cache is
#' cleared at the start of each fit). A group whose curvature is not positive
#' definite at its mode adds the penalty value \eqn{-10^6} instead of a
#' silently regularised term, and \code{brsmm()} warns.
#'
#' \eqn{D = LL^\top} is parameterised by its lower Cholesky factor: the
#' log of each diagonal entry and the off-diagonal entries as they are
#' (\code{(re_chol_logsd)_} and \code{(re_chol)_} in \code{coef()}).
#' \code{\link{summary.brsmm}} reports the standard deviations and
#' correlations with intervals; \code{\link{brsmm_re_study}} the intraclass
#' correlation.
#'
#' @section Estimation:
#' \code{\link[stats]{optim}} maximises the marginal log-likelihood with its
#' compiled gradient: the derivative of the chosen approximation by the chain
#' rule on the linear predictor and the implicit-function theorem at the group
#' modes (including the movement of the quadrature nodes), with
#' per-observation derivatives by central differences. The Hessian for
#' \code{vcov()} (\code{hessian_method = "cpp"}, the default) is a Richardson
#' central difference of that gradient. Starting values: \code{start} when
#' given; otherwise those of \code{\link{brs}} for the fixed effects and
#' \eqn{\log} of the between-group SD of the cell centres (at least 0.1) for
#' the random effects.
#'
#' @section Diagnostics:
#' The checks of \code{\link{brs}} ('Fit diagnostics') apply, with the
#' compiled gradient (a central difference with step \eqn{10^{-3}} if it is
#' not finite). In addition:
#' \describe{
#'   \item{\dQuote{Variance component at the boundary}}{A random-effect term
#'     has log SD below \eqn{-6}, or raises the log-likelihood by less than
#'     \eqn{10^{-3}} over the same fit with that SD at zero. Its variance is
#'     essentially zero; the standard error of its log SD is meaningless. Test
#'     the term with \code{\link{anova.brsmm}} (chi-bar-square mixture) and
#'     drop it if not needed. When the gain is negative the message adds that
#'     SD \eqn{\approx 0} has a higher log-likelihood: \code{optim} stopped
#'     short of the maximum.}
#'   \item{\dQuote{group(s) have no positive-definite random-effect mode}}{Those
#'     groups contribute the penalty value; the fit is not reliable. Simplify
#'     the random-effects structure or check the data of those groups.}
#' }
#' \code{fit$diagnostics} stores \code{re_boundary}, \code{re_gain} and
#' \code{inner} (number of such groups, largest gradient norm at the modes).
#' Rank-deficient fixed-effect or random-effect design matrices are an error.
#'
#' @param formula Model formula: \code{y ~ x1 + x2} or
#'   \code{y ~ x1 + x2 | z1 + z2} (see \code{\link{brs}}).
#' @param random Random-effects formula \code{~ terms | group}, e.g.
#'   \code{~ 1 | id} or \code{~ 1 + x | id}.
#' @param data Data frame (raw scores, or the output of
#'   \code{\link{brs_prep}}).
#' @param link,link_phi Links for the first and second parameter;
#'   \code{NULL} (default) selects those implied by \code{repar} (see
#'   \code{\link{brs}}).
#' @param repar Parameterisation (0, 1 or 2); see \code{\link{brs_repar}}.
#' @param ncuts Integer \eqn{K}: the maximum score (scale
#'   \eqn{0, \ldots, K}). \code{NULL} (default) uses the value stored by
#'   \code{\link{brs_prep}}, or 100; an explicit different value is ignored
#'   with a warning.
#' @param lim Half-width of the cell in \eqn{(0, 0.5]} (\code{interval =
#'   "mid"} only); \code{NULL} uses the stored value, or 0.5.
#' @param interval \code{"mid"}, \code{"right"} or \code{"left"} (see
#'   \code{\link{brs_check}}); \code{NULL} uses the stored value, or
#'   \code{"mid"}.
#' @param int_method \code{"laplace"} (default), \code{"aghq"} or
#'   \code{"qmc"}; see Details.
#' @param n_points Nodes per dimension for \code{"aghq"} (default 11).
#' @param qmc_points Halton points for \code{"qmc"} (default 1024).
#' @param start Optional starting vector: fixed effects, precision
#'   coefficients, then the packed Cholesky parameters.
#' @param method \code{"BFGS"} (default) or \code{"L-BFGS-B"}.
#' @param hessian_method \code{"cpp"} (default; Richardson differences of
#'   the compiled gradient), \code{"numDeriv"} or \code{"optim"}.
#' @param control Control list for \code{\link[stats]{optim}}, merged into
#'   the default \code{list(maxit = 2000L)}.
#'
#' @return An object of class \code{"brsmm"} with the components of a
#'   \code{"brs"} fit (\code{par}, \code{coefficients} with a \code{random}
#'   part, \code{value}, \code{hessian}, \code{diagnostics}, ...) and
#'   \code{random} (group variable, levels, conditional modes \code{mode_b},
#'   \code{D}, \code{L} and the SDs \code{sd_b}), \code{ngroups},
#'   \code{int_method}.
#'
#' @seealso \code{\link{brs}}, \code{\link{summary.brsmm}},
#'   \code{\link{anova.brsmm}}, \code{\link{brsmm_re_study}},
#'   \code{\link{ranef.brsmm}}, \code{\link{predict.brsmm}}
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
#' Pinheiro, J. C., and Bates, D. M. (1995). Approximations to the
#' log-likelihood function in the nonlinear mixed-effects model.
#' \emph{Journal of Computational and Graphical Statistics}, \bold{4}(1),
#' 12--35. \doi{10.1080/10618600.1995.10474663}
#'
#' @examples
#' # Synthetic NRS-11 scores (0-10) of 40 patients at 6h, 12h and 24h;
#' # intercepts and time slopes vary by patient. Simulated, not real data.
#' set.seed(21)
#' nrs <- expand.grid(id = 1:40, time = c("6h", "12h", "24h"))
#' nrs$tc <- c(-1, 0, 1)[nrs$time]                 # centred time for the slope
#' b0 <- rnorm(40, sd = 0.8)
#' b1 <- rnorm(40, sd = 0.5)
#' eta <- -1.3 + c(0, 0.75, 0.3)[nrs$time] + b0[nrs$id] + b1[nrs$id] * nrs$tc
#' shp <- brs_repar(mu = plogis(eta), phi = 0.2, repar = 2)
#' nrs$y <- round(10 * rbeta(nrow(nrs), shp$shape1, shp$shape2))
#'
#' # Random intercept per patient (Laplace): SD with its interval, no z-test
#' m1 <- brsmm(y ~ time, random = ~ 1 | id, data = nrs, ncuts = 10)
#' summary(m1)
#'
#' # Random intercept and slope: SDs and their correlation
#' m2 <- brsmm(y ~ time, random = ~ 1 + tc | id, data = nrs, ncuts = 10)
#' summary(m2)$varcorr
#' head(ranef(m2))
#'
#' # Is the slope needed? One added random term: the p-value uses the
#' # mixture 1/2 chi2(1) + 1/2 chi2(2) (see the printed heading)
#' anova(m1, m2)
#'
#' # A known patient uses its conditional mode; a new one gets b = 0
#' nd <- data.frame(time = "12h", tc = 0, id = c(1, 999))
#' predict(m2, newdata = nd)
#'
#' @importFrom Formula as.Formula Formula
#' @importFrom stats model.frame terms delete.response model.matrix
#' @importFrom stats cor optim
#' @importFrom numDeriv hessian
#' @export
brsmm <- function(formula,
                  random = ~ 1 | id,
                  data,
                  link = NULL,
                  link_phi = NULL,
                  repar = 2L,
                  ncuts = NULL,
                  lim = NULL,
                  int_method = c("laplace", "aghq", "qmc"),
                  n_points = 11L,
                  qmc_points = 1024L,
                  start = NULL,
                  method = c("BFGS", "L-BFGS-B"),
                  hessian_method = c("cpp", "numDeriv", "optim"),
                  control = list(maxit = 2000L),
                  interval = NULL) {
  cl <- match.call()
  .brsmm_reset_cache()   # warm starts never carry over from earlier fits
  method <- match.arg(method)
  hessian_method <- match.arg(hessian_method)
  int_method <- match.arg(int_method)
  n_points <- as.integer(n_points)
  qmc_points <- as.integer(qmc_points)

  # Same checks as brs(); an unknown repar used to reach the C++ `default:` branch.
  validated <- .validate_brs_common_args(data, ncuts, lim, repar, link, link_phi,
                                        interval)
  ncuts <- validated$ncuts
  lim <- validated$lim
  repar <- validated$repar
  link <- validated$link
  link_phi <- validated$link_phi
  interval <- validated$interval
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
  Y <- .extract_response(mf, data, ncuts = ncuts, lim = lim, interval = interval)
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
  if (nrow(Xr) != nrow(X)) {
    # model.matrix() drops rows with NA in the random-effect variables while
    # X/Y/group keep them; the compiled code assumes equal lengths.
    stop(
      "Random-effect variables in 'random' contain missing values (",
      nrow(X) - nrow(Xr), " row(s)). Remove or impute them before fitting.",
      call. = FALSE
    )
  }
  q_re <- ncol(Xr)
  k_re <- q_re * (q_re + 1L) / 2L
  # Aliased columns stop before optim (a flat direction gave non-finite values)
  .brs_check_design(X, "mean")
  .brs_check_design(Z, "precision")
  .brs_check_design(Xr, "random-effects")
  n <- nrow(X)
  g <- nlevels(group)


  if (is.null(start)) {
    start_fix <- compute_start(
      formula = formula_parsed,
      data = data,
      link = link,
      link_phi = link_phi,
      ncuts = ncuts,
      lim = lim,
      repar = repar,
      interval = interval
    )
    if (length(start_fix) != (p + q)) {
      stop(
        "Internal error: starting vector from compute_start() has unexpected length.",
        call. = FALSE
      )
    }
    # QUAL-L04: estimate sigma_b from between-group variance of y_mid
    y_mid_start <- as.numeric(Y[, "yt"])
    group_means <- tapply(y_mid_start, group_index, mean)
    sigma_b_init <- max(stats::sd(as.numeric(group_means)), 0.1, na.rm = TRUE)

    theta_re_start <- numeric(k_re)
    k <- 1L
    for (j in seq_len(q_re)) {
      for (i in j:q_re) {
        theta_re_start[k] <- if (i == j) log(sigma_b_init) else 0
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
    .brsmm_loglik_eigen(
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

  # Gradient of the chosen approximation: chain rule + implicit-function theorem,
  # per-observation finite differences in eta (see ?brsmm, hessian_method).
  gr_ll <- function(par) {
    .brsmm_grad_cpp(
      param = as.numeric(par), X = X, Z = Z, Xr = Xr,
      y_left = as.numeric(Y[, "left"]), y_right = as.numeric(Y[, "right"]),
      yt = as.numeric(Y[, "yt"]), delta = delta, group = group_index,
      link_mu = lc_mu, link_phi = lc_phi, repar = repar,
      method = method_code, n_points = n_pts
    )
  }
  gr_obj <- function(par) -gr_ll(par)

  opt <- stats::optim(
    par = start,
    fn = fn_obj,
    gr = gr_obj,
    method = method,
    hessian = (hessian_method == "optim"),
    # User entries override; the default maxit stays unless given
    control = .brs_merge_control(list(maxit = 2000L), control)
  )

  # BUG-H04: warn if optimizer did not converge
  if (opt$convergence != 0L) {
    warning(
      "Optimizer did not converge (code ", opt$convergence, ")",
      if (!is.null(opt$message)) paste0(": ", opt$message) else ".",
      "\nResults may be unreliable. Try increasing 'control$maxit' or changing 'method'.",
      call. = FALSE
    )
  }

  if (hessian_method == "cpp") {
    hess <- .brsmm_hessian_cpp(
      param = opt$par, X = X, Z = Z, Xr = Xr,
      y_left = as.numeric(Y[, "left"]), y_right = as.numeric(Y[, "right"]),
      yt = as.numeric(Y[, "yt"]), delta = delta, group = group_index,
      link_mu = lc_mu, link_phi = lc_phi, repar = repar,
      method = method_code, n_points = n_pts
    )
  } else if (hessian_method == "numDeriv") {
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

  gm <- .brsmm_group_modes_eigen(
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
  # Inner-mode diagnostics (reported by .brs_fit_diagnostics): a group without a
  # positive-definite mode is penalised in the likelihood.
  diag_re <- .brsmm_mode_diag_cpp(
    param = est, X = X, Z = Z, Xr = Xr,
    y_left = as.numeric(Y[, "left"]), y_right = as.numeric(Y[, "right"]),
    yt = as.numeric(Y[, "yt"]), delta = delta, group = group_index,
    link_mu = lc_mu, link_phi = lc_phi, repar = repar, warm = TRUE
  )
  eta_phi <- as.numeric(Z %*% gamma_hat)
  y_mid <- as.numeric(Y[, "yt"])

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

  # fitted_mu is the FIRST parameter (shape p under repar 0); E[Y] via .brs_mean().
  # Both clamps mirror the compiled likelihood (src/brs_common.h).
  mu_raw <- apply_inv_link(eta_mu, link)
  phi_raw <- apply_inv_link(eta_phi, link_phi)
  hatmu <- .clamp_mu_by_repar(mu_raw, repar)
  hatphi <- .clamp_phi_by_repar(phi_raw, repar)
  # Compiled gradient of the marginal likelihood; central differences only if
  # it is not finite. With the fit's Hessian it gives the Newton-gain criterion.
  g_hat <- tryCatch(gr_ll(est), error = function(e) NULL)
  if (is.null(g_hat) || any(!is.finite(g_hat))) g_hat <- .brs_num_grad(fn_ll, est)
  fit_diag <- .brs_fit_diagnostics(
    g_hat, hess, mu_raw, phi_raw, repar, opt$convergence,
    re_logsd = c(log(sd_b_terms), theta_re_hat[.brsmm_chol_diag_index(q_re)]),
    re_gain = .brsmm_re_gain(fn_ll, -opt$value, est, idx_re, q_re),
    badly_scaled = .brs_badly_scaled(X, Z, Xr),
    inner = list(bad_groups = sum(diag_re$ok == 0), max_grad = max(diag_re$grad_inf))
  )
  ey <- .brs_mean(as.numeric(hatmu), hatphi, repar)

  pseudo_r2 <- suppressWarnings(
    .brs_pseudo_r2(X %*% beta_hat, ey, y_mid, link, repar)
  )
  .warn_sqrt_plateau(eta_mu, link, "link")
  .warn_sqrt_plateau(eta_phi, link_phi, "link_phi")

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
    residuals = as.numeric(y_mid - ey),
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
    xlevels = list(
      mean = stats::.getXlevels(mtX, mf),
      precision = stats::.getXlevels(mtZ, mf),
      random = stats::.getXlevels(
        random_spec$re_terms,
        stats::model.frame(random_spec$re_terms, data_sub)
      )
    ),
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
    interval = interval,
    method = method,
    hessian_method = hessian_method,
    int_method = int_method,
    n_points = n_points,
    qmc_points = qmc_points,
    diagnostics = c(fit_diag, list(
      integration = list(
        method = int_method,
        n_groups = g
      )
    ))
  )

  class(out) <- "brsmm"
  out
}

#' Parse random-effect specification for brsmm
#' @keywords internal
#' @noRd
# Named entries of `user` replace those of `default` (base-R modifyList()).
.brs_merge_control <- function(default, user) {
  user <- as.list(user)
  if (length(user) && is.null(names(user))) {
    stop("'control' must be a named list.", call. = FALSE)
  }
  default[names(user)] <- user
  default
}

# Log-likelihood gain of each random-effect term over the same parameters with
# that term's Cholesky row at ~0 (log-diagonal -10: smaller SDs lose Laplace accuracy).
.brsmm_re_gain <- function(fn_ll, ll_hat, est, idx_re, q_re) {
  # ll_hat = -opt$value: the log-likelihood the optimizer reached (no re-evaluation)
  theta <- est[idx_re]
  vapply(seq_len(q_re), function(i) {
    th <- theta
    k <- 1L
    for (j in seq_len(q_re)) {
      for (r in j:q_re) {
        if (r == i) th[k] <- if (r == j) -10 else 0
        k <- k + 1L
      }
    }
    par0 <- est
    par0[idx_re] <- th
    ll_hat - fn_ll(par0)
  }, numeric(1))
}

# Positions of the log-diagonal Cholesky entries in the column-wise theta_re.
.brsmm_chol_diag_index <- function(q_re) {
  cumsum(c(1L, if (q_re > 1L) rev(seq_len(q_re - 1L)) + 1L))[seq_len(q_re)]
}

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
  # Row names are labels, not positions (see .extract_response); a permuted
  # or subsetted `data` silently misaligned group/Xr with X/Y otherwise.
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
