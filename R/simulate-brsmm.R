# ============================================================================ #
# Simulation of clustered scores and simulate() methods
# ============================================================================ #

#' Simulate clustered scores from a mixed beta interval model
#'
#' @description
#' Simulates scores on \eqn{0, 1, \ldots, K} (\eqn{K =} \code{ncuts}) from the
#' model of \code{\link{brsmm}}: Gaussian random effects in the linear
#' predictor of the first parameter, a latent beta value per row, and the
#' coarsening of \code{\link{brs_sim}}. The output can be passed to
#' \code{\link{brsmm}} as it is.
#'
#' @details
#' For group \eqn{g} of the grouping variable of \code{random},
#' \eqn{b_g \sim N(0, D)} with \eqn{D = \mathrm{diag}(s)\, R\,
#' \mathrm{diag}(s)}, \eqn{s =} \code{re_sd} and \eqn{R =} \code{re_cor}
#' (the standard deviations and correlations reported by
#' \code{\link{summary.brsmm}}). Row \eqn{i} of group \eqn{g} has
#' \deqn{g_1(\mu_{gi}) = x_{gi}^\top \beta + x_{r,gi}^\top b_g, \qquad
#'   g_2(\phi_{gi}) = z_{gi}^\top \zeta,}
#' a latent \eqn{Y_{gi} \sim \mathrm{Beta}} with the shapes of
#' \code{\link{brs_repar}}, and the score of \code{\link{brs_sim}}:
#' \eqn{\mathrm{round}(K Y)} under \code{"mid"}, \eqn{\lfloor (K + 1) Y
#' \rfloor} under \code{"right"}/\code{"left"}, with the boundary cutpoints
#' \code{c(c_lo, c_hi)} when given. Score 0 is left-censored, \eqn{K}
#' right-censored and the rest interval-censored.
#'
#' The random effects are drawn first (a \eqn{G \times q} matrix of standard
#' normals, \eqn{G} groups in the order of \code{levels(factor(group))},
#' times \eqn{L^\top}, with \eqn{D = LL^\top} the lower Cholesky factor),
#' then the \eqn{n} beta values; \code{set.seed()} before the call makes the
#' output reproducible, and \code{\link{simulate.brsmm}} uses the same
#' sequence.
#'
#' @param formula Model formula with one (mean) or two parts
#'   (mean \code{|} dispersion); a left-hand side is ignored.
#' @param random Random-effects formula \code{~ terms | group}, as in
#'   \code{\link{brsmm}}.
#' @param data Data frame with the covariates and the grouping variable.
#' @param beta Fixed-effect coefficients of the first parameter, on its link
#'   scale.
#' @param phi Second parameter on its link scale (one-part formulas).
#' @param zeta Coefficients of the second parameter, on its link scale
#'   (two-part formulas).
#' @param re_sd Standard deviations of the random effects, one per column
#'   of the random-effects design (e.g. intercept, then slope); positive.
#' @param re_cor Correlation matrix of the random effects; \code{NULL}
#'   (default) for independent effects. With two random effects a single
#'   correlation is accepted.
#' @param link,link_phi Links; \code{NULL} (default) selects those implied by
#'   \code{repar} (see \code{\link{brs}}).
#' @param ncuts Integer \eqn{K}: the maximum score (default 100).
#' @param lim Half-width of the cell under \code{"mid"} (default 0.5).
#' @param repar Parameterization (0, 1 or 2); see \code{\link{brs_repar}}.
#' @param interval \code{"mid"} (default), \code{"right"} or \code{"left"}.
#' @param cutpoints \code{NULL} (default) or \code{c(c_lo, c_hi)}, boundary
#'   cutpoints on \eqn{(0, 1)}, as in \code{\link{brs_sim}}.
#'
#' @return \code{data} with the columns \code{left}, \code{right}, \code{yt},
#'   \code{y} (score) and \code{delta} in front (existing columns with these
#'   names are replaced), and the attributes \code{"is_prepared"},
#'   \code{"ncuts"}, \code{"lim"} and \code{"interval"} that
#'   \code{\link{brsmm}} reuses. Attribute \code{"ranef"}: the simulated
#'   \eqn{b_g}, one row per group level. Attribute \code{"truth"}: a list with
#'   \code{beta}, \code{phi} or \code{zeta}, \code{re_sd}, \code{re_cor},
#'   \code{D}, \code{L}, \code{cutpoints} and \code{par}, the true parameters
#'   in the order and on the scale of \code{coef()} of a \code{brsmm} fit.
#'
#' @seealso \code{\link{brs_sim}}, \code{\link{brsmm}},
#'   \code{\link{simulate.brsmm}}
#'
#' @examples
#' # 100 patients with 5 visits; random intercept and time slope
#' d <- expand.grid(visit = 1:5, id = 1:100)
#' d$tc <- (d$visit - 3) / 2
#' set.seed(7)
#' s <- brsmm_sim(~ tc, random = ~ 1 + tc | id, data = d,
#'                beta = c(-0.5, 0.4), phi = qlogis(0.2),
#'                re_sd = c(0.8, 0.3), re_cor = 0.3, ncuts = 10)
#' head(s)
#' head(attr(s, "ranef"))
#' attr(s, "truth")$par
#'
#' # Random intercept only, boundary cutpoints (heaping at 0 and 10)
#' set.seed(8)
#' s1 <- brsmm_sim(~ tc, random = ~ 1 | id, data = d, beta = c(-0.5, 0.4),
#'                 phi = qlogis(0.2), re_sd = 0.8, ncuts = 10,
#'                 cutpoints = c(0.1, 0.9))
#' table(s1$y)
#' \donttest{
#' m <- brsmm(y ~ tc, random = ~ 1 | id, data = s1, cutpoints = "boundary")
#' round(cbind(true = attr(s1, "truth")$par, est = coef(m)), 2)
#' }
#'
#' @export
brsmm_sim <- function(formula,
                      random = ~ 1 | id,
                      data,
                      beta,
                      phi = 1 / 5,
                      zeta = NULL,
                      re_sd,
                      re_cor = NULL,
                      link = NULL,
                      link_phi = NULL,
                      ncuts = 100L,
                      lim = 0.5,
                      repar = 2L,
                      interval = c("mid", "right", "left"),
                      cutpoints = NULL) {
  if (!is.data.frame(data)) {
    stop("'data' must be a data.frame.", call. = FALSE)
  }
  interval <- match.arg(interval)
  .brs_lim_check(lim, interval, warn = TRUE, generator = TRUE)
  ncuts <- as.integer(ncuts)
  if (length(ncuts) != 1L || is.na(ncuts) || ncuts < 1L) {
    stop("'ncuts' must be a positive integer.", call. = FALSE)
  }
  if (!is.null(cutpoints)) {
    # Same rules as brs_sim(): score data on cells that partition (0, 1)
    if (abs(lim - 0.5) > sqrt(.Machine$double.eps)) {
      stop("`cutpoints` needs lim = 0.5.", call. = FALSE)
    }
    cutpoints <- .brs_cut_check_values(cutpoints, ncuts, interval)
  }
  repar <- as.integer(repar)
  links <- .resolve_links(link, link_phi, repar)
  link <- links$link
  link_phi <- links$link_phi

  spec <- .brsmm_parse_random(random)
  gv <- spec$group_var
  if (!(gv %in% names(data))) {
    stop("Grouping variable '", gv, "' not found in data.", call. = FALSE)
  }
  if (anyNA(data[[gv]])) {
    stop("Grouping variable '", gv, "' contains missing values.", call. = FALSE)
  }
  group <- factor(data[[gv]])
  gidx <- as.integer(group)

  design <- .sim_design_matrices(formula, data)
  X <- design$X
  Z <- design$Z
  n <- nrow(data)
  if (nrow(X) != n || (!is.null(Z) && nrow(Z) != n)) {
    stop("Covariates of 'formula' contain missing values.", call. = FALSE)
  }
  mf_re <- stats::model.frame(spec$re_terms, data, na.action = stats::na.pass)
  Xr <- stats::model.matrix(spec$re_terms, mf_re)
  if (nrow(Xr) != n || anyNA(Xr)) {
    stop("Random-effect variables in 'random' contain missing values.", call. = FALSE)
  }
  q_re <- ncol(Xr)
  re_names <- colnames(Xr)

  if (!is.numeric(beta) || length(beta) != ncol(X) || anyNA(beta)) {
    stop("'beta' length (", length(beta), ") must equal ncol(X) (", ncol(X), ").",
         call. = FALSE)
  }
  if (is.null(Z)) {
    if (!is.null(zeta)) {
      stop("'zeta' is only valid for formulas with '|'.", call. = FALSE)
    }
    if (length(phi) != 1L) {
      stop("'phi' must be a scalar for one-part formulas.", call. = FALSE)
    }
    eta_phi <- rep(phi, n)
    disp <- c("(phi)_(Intercept)" = phi)
  } else {
    if (is.null(zeta) || length(zeta) != ncol(Z)) {
      stop("'zeta' must have length ncol(Z) (", ncol(Z), ").", call. = FALSE)
    }
    eta_phi <- drop(Z %*% zeta)
    disp <- stats::setNames(zeta, paste0("(phi)_", colnames(Z)))
  }
  re <- .brsmm_sim_cov(re_sd, re_cor, q_re, re_names)

  eta_fix <- drop(X %*% beta)
  # Validate the parameters on their ranges (as brs_sim() does) at b = 0
  brs_repar(mu = apply_inv_link(eta_fix, link),
            phi = apply_inv_link(eta_phi, link_phi), repar = repar)

  sim <- .brsmm_sim_draw(eta_fix, eta_phi, Xr, gidx, nlevels(group), re$L,
                         link, link_phi, repar)
  out_y <- .build_simulated_response(sim$y, delta = NULL, ncuts = ncuts, lim = lim,
                                     interval = interval, cuts = cutpoints)
  d_out <- as.integer(out_y[, "delta"])
  if (n > 0L && (all(d_out == 1L) || all(d_out == 2L))) {
    warning("Every simulated observation is censored on the same side (delta = ",
            d_out[1L], "): brsmm() has no finite MLE for these data.", call. = FALSE)
  }

  resp <- colnames(out_y)
  keep <- data[, setdiff(names(data), resp), drop = FALSE]
  result <- data.frame(out_y, keep, check.names = FALSE)
  rownames(result) <- rownames(data)
  attr(result, "is_prepared") <- TRUE
  attr(result, "ncuts") <- ncuts
  attr(result, "lim") <- lim
  attr(result, "interval") <- interval

  b <- sim$b
  dimnames(b) <- list(levels(group), re_names)
  # True parameters in the order and scale of coef() of a brsmm fit
  par <- c(stats::setNames(beta, colnames(X)), disp,
           .brsmm_pack_chol(re$L, re_names, gv))
  if (!is.null(cutpoints)) {
    par <- c(par, stats::setNames(
      .brs_cut_par(cutpoints, ncuts, .brs_cut_anchors(ncuts, interval)),
      c("(cut_lo)", "(cut_hi)")))
  }
  attr(result, "ranef") <- b
  attr(result, "truth") <- list(
    beta = stats::setNames(beta, colnames(X)),
    phi = if (is.null(Z)) phi else NULL,
    zeta = if (is.null(Z)) NULL else stats::setNames(zeta, colnames(Z)),
    re_sd = stats::setNames(sqrt(diag(re$D)), re_names),
    re_cor = re$R, D = re$D, L = re$L,
    cutpoints = cutpoints, par = par
  )
  result
}

# Covariance of the random effects from SDs and correlations: validated D, its
# correlation matrix R and lower Cholesky factor L.
.brsmm_sim_cov <- function(re_sd, re_cor, q, nm) {
  if (missing(re_sd) || !is.numeric(re_sd) || length(re_sd) != q ||
      any(!is.finite(re_sd)) || any(re_sd <= 0)) {
    stop("'re_sd' must be ", q, " positive finite number(s), one per random ",
         "effect (", paste(nm, collapse = ", "), ").", call. = FALSE)
  }
  R <- if (is.null(re_cor)) {
    diag(q)
  } else if (q == 2L && is.numeric(re_cor) && length(re_cor) == 1L) {
    matrix(c(1, re_cor, re_cor, 1), 2L)
  } else {
    as.matrix(re_cor)
  }
  ok <- is.numeric(R) && identical(dim(R), c(q, q)) && all(is.finite(R)) &&
    isSymmetric(unname(R), tol = 1e-10) && all(abs(diag(R) - 1) < 1e-10) &&
    all(abs(R) <= 1)
  if (!ok) {
    stop("'re_cor' must be NULL, a ", q, " x ", q, " correlation matrix",
         if (q == 2L) " or one correlation in [-1, 1]" else "", ".", call. = FALSE)
  }
  D <- diag(re_sd, q) %*% R %*% diag(re_sd, q)
  U <- tryCatch(chol(D), error = function(e) NULL)
  # brsmm() parameterizes D by a Cholesky factor with positive diagonal
  if (is.null(U) || any(diag(U) < 1e-8 * max(re_sd))) {
    stop("The random-effect covariance from 're_sd' and 're_cor' is not ",
         "positive definite.", call. = FALSE)
  }
  dimnames(R) <- dimnames(D) <- list(nm, nm)
  L <- t(U)
  dimnames(L) <- list(nm, nm)
  list(D = D, R = R, L = L)
}

# Packed Cholesky parameters (log diagonal, off-diagonal as is), named as in
# coef() of a brsmm fit.
.brsmm_pack_chol <- function(L, re_names, group_var) {
  q <- nrow(L)
  val <- numeric(0)
  nm <- character(0)
  for (j in seq_len(q)) {
    for (i in j:q) {
      if (i == j) {
        val <- c(val, log(L[i, i]))
        nm <- c(nm, paste0("(re_chol_logsd)_", re_names[i], "|", group_var))
      } else {
        val <- c(val, L[i, j])
        nm <- c(nm, paste0("(re_chol)_", re_names[i], ":", re_names[j], "|", group_var))
      }
    }
  }
  stats::setNames(val, nm)
}

# One draw of the random effects (G x q standard normals times t(L)) and of the
# latent beta values; `b` given skips the first step (conditional simulation).
.brsmm_sim_draw <- function(eta_fix, eta_phi, Xr, gidx, G, L, link, link_phi,
                            repar, b = NULL) {
  q <- ncol(Xr)
  if (is.null(b)) {
    b <- matrix(stats::rnorm(G * q), G, q) %*% t(L)
  }
  eta <- eta_fix + rowSums(Xr * b[gidx, , drop = FALSE])
  sh <- .brs_shapes_cpp(eta, eta_phi, link, link_phi, repar)
  y <- stats::rbeta(length(eta), sh$a, sh$b)
  list(b = b, y = y)
}

# Scores of latent values with the cells of a fit (estimated cutpoints for
# boundary fits); rows with exact values (delta = 0) keep the latent value.
.brs_sim_coarsen <- function(y, object) {
  K <- as.integer(object$ncuts)
  s <- .brs_score_from_unit_cut(y, K, .brs_interval_of(object), .brs_cuts_of(object))
  exact <- as.integer(object$Y[, "delta"]) == 0L
  s[exact] <- y[exact]
  s
}

# RNG state handling of stats::simulate(): the "seed" attribute value and, with
# a seed, the caller's state to restore (`old`).
.brs_sim_seed <- function(seed) {
  if (!exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    stats::runif(1)
  }
  if (is.null(seed)) {
    return(list(rng = get(".Random.seed", envir = globalenv()), old = NULL))
  }
  old <- get(".Random.seed", envir = globalenv())
  set.seed(seed)
  list(rng = structure(seed, kind = as.list(RNGkind())), old = old)
}


#' Simulate responses from a fitted beta interval model
#'
#' @description
#' Draws new responses from a \code{\link{brs}} or \code{\link{brsmm}} fit at
#' the estimated parameters, coarsened with the cells of the fit, as
#' \code{\link{brs_sim}} and \code{\link{brsmm_sim}} do.
#'
#' @details
#' Each simulation draws one latent beta value per row with the fitted shapes
#' and returns its score on \eqn{0, \ldots, K} (cells of the fit's
#' \code{interval} map, or its estimated boundary cutpoints); rows fitted as
#' exact values (\eqn{\delta = 0}) get the latent value itself.
#'
#' For a \code{brsmm} fit, \code{re = "marginal"} (default) draws new random
#' effects \eqn{b_g \sim N(0, \hat D)} for each simulation (new groups with
#' the same design: the marginal distribution of the data), in the order of
#' \code{\link{brsmm_sim}}, so \code{simulate(fit, seed = s)} and
#' \code{set.seed(s); brsmm_sim(...)} with the estimates as true values give
#' the same scores. \code{re = "conditional"} keeps the random effects at their
#' conditional modes (\code{\link{ranef.brsmm}}) and draws only the latent
#' values: the same groups again.
#'
#' @param object A fitted \code{"brs"} or \code{"brsmm"} object.
#' @param nsim Number of simulated response vectors (default 1).
#' @param seed \code{NULL} or an integer seed, as in
#'   \code{\link[stats]{simulate}}: with a seed the caller's RNG state is
#'   restored afterwards.
#' @param re For \code{brsmm} fits: \code{"marginal"} (default; new random
#'   effects in each simulation) or \code{"conditional"} (random effects at
#'   the conditional modes).
#' @param ... Currently ignored.
#'
#' @return A data frame with \code{nsim} columns \code{sim_1}, \ldots, one row
#'   per observation of the fit, and the attribute \code{"seed"} (as
#'   \code{\link[stats]{simulate}}). For \code{brsmm} fits with
#'   \code{re = "marginal"}, the attribute \code{"ranef"} holds the random
#'   effects of each simulation (a list of matrices, one row per group).
#'
#' @seealso \code{\link{brs_sim}}, \code{\link{brsmm_sim}}
#'
#' @examples
#' set.seed(2023)
#' d <- data.frame(x = runif(150))
#' s <- brs_sim(~ x, data = d, beta = c(-0.5, 1), phi = qlogis(0.3), ncuts = 10)
#' fit <- brs(y ~ x, data = s)
#' sims <- simulate(fit, nsim = 200, seed = 1)
#' # Observed and simulated frequencies of the scores
#' rbind(observed = tabulate(s$y + 1, 11) / nrow(s),
#'       simulated = tabulate(unlist(sims) + 1, 11) / (200 * nrow(s)))
#'
#' @name simulate.brs
#' @importFrom stats simulate
#' @method simulate brs
#' @export
simulate.brs <- function(object, nsim = 1, seed = NULL, ...) {
  .check_class(object)
  nsim <- .brs_check_nsim(nsim)
  st <- .brs_sim_seed(seed)
  if (!is.null(st$old)) {
    on.exit(assign(".Random.seed", st$old, envir = globalenv()), add = TRUE)
  }
  p <- object$p
  par <- unname(object$par)
  mm <- object$model_matrices
  eta_mu <- drop(mm$X %*% par[seq_len(p)])
  eta_phi <- if (is.null(mm$Z)) {
    rep(par[p + 1L], length(eta_mu))
  } else {
    drop(mm$Z %*% par[p + seq_len(ncol(mm$Z))])
  }
  # Shapes as the likelihood builds them
  sh <- .brs_shapes_cpp(eta_mu, eta_phi, object$link, object$link_phi, object$repar)
  n <- length(eta_mu)
  out <- lapply(seq_len(nsim), function(r) {
    .brs_sim_coarsen(stats::rbeta(n, sh$a, sh$b), object)
  })
  .brs_sim_frame(out, object, st$rng)
}

#' @rdname simulate.brs
#' @method simulate brsmm
#' @export
simulate.brsmm <- function(object, nsim = 1, seed = NULL,
                           re = c("marginal", "conditional"), ...) {
  .check_class_mm(object)
  re <- match.arg(re)
  nsim <- .brs_check_nsim(nsim)
  st <- .brs_sim_seed(seed)
  if (!is.null(st$old)) {
    on.exit(assign(".Random.seed", st$old, envir = globalenv()), add = TRUE)
  }
  p <- object$p
  par <- unname(object$par)
  mm <- object$model_matrices
  eta_fix <- drop(mm$X %*% par[seq_len(p)])
  eta_phi <- drop(mm$Z %*% par[p + seq_len(ncol(mm$Z))])
  G <- length(object$random$levels)
  L <- object$random$L
  b_mode <- matrix(object$random$mode_b, nrow = G)
  draws <- lapply(seq_len(nsim), function(r) {
    .brsmm_sim_draw(eta_fix, eta_phi, mm$Xr, object$group_index, G, L,
                    object$link, object$link_phi, object$repar,
                    b = if (re == "conditional") b_mode else NULL)
  })
  out <- .brs_sim_frame(lapply(draws, function(d) .brs_sim_coarsen(d$y, object)),
                        object, st$rng)
  if (re == "marginal") {
    attr(out, "ranef") <- lapply(draws, function(d) {
      dimnames(d$b) <- list(object$random$levels, object$random$terms)
      d$b
    })
  }
  out
}

# nsim as a positive integer.
.brs_check_nsim <- function(nsim) {
  nsim <- as.integer(nsim)
  if (length(nsim) != 1L || is.na(nsim) || nsim < 1L) {
    stop("'nsim' must be a positive integer.", call. = FALSE)
  }
  nsim
}

# Data frame sim_1..sim_nsim with the fit's row names and the "seed" attribute.
.brs_sim_frame <- function(sims, object, rng) {
  out <- as.data.frame(stats::setNames(sims, paste0("sim_", seq_along(sims))))
  rn <- rownames(object$Y)
  if (!is.null(rn) && !anyDuplicated(rn)) rownames(out) <- rn
  attr(out, "seed") <- rng
  out
}
