# ============================================================================ #
# Parametric bootstrap for beta interval regression
#
# Provides bootstrap-based confidence intervals for model parameters,
# complementing the asymptotic (Wald) intervals from the Hessian.
# ============================================================================ #

#' Parametric bootstrap confidence intervals for brs models
#'
#' @description
#' Computes bootstrap-based confidence intervals for the parameters of a
#' fitted \code{"brs"} model by repeatedly simulating data from the fitted
#' model and re-estimating parameters. Only \code{"brs"} (fixed or
#' variable-dispersion) objects are supported; \code{"brsmm"} is not supported.
#'
#' @details
#' For each replicate, data are simulated via \code{\link{brs_sim}} using
#' the estimated coefficients (on the link scale) and the original
#' design. The model is then re-fitted with \code{\link{brs}}. Replicates
#' that fail to converge are discarded; if the number of successful replicates
#' is too low, a warning is issued. Intervals are the empirical quantiles of
#' the bootstrap distribution of each parameter.
#'
#' @param object A fitted \code{"brs"} object (fixed or variable dispersion).
#' @param R Integer: number of bootstrap replicates (default 199).
#' @param level Numeric: confidence level (default 0.95).
#' @param seed Optional integer: random seed for reproducibility.
#' @param ci_type Character: type of confidence interval. One of
#'   \code{"percentile"} (default), \code{"basic"}, \code{"normal"},
#'   or \code{"bca"}.
#' @param max_tries Optional integer: maximum number of bootstrap attempts
#'   to obtain converged replicates. If \code{NULL}, uses \code{max(3 * R, 50)}.
#' @param keep_draws Logical: if \code{TRUE}, stores successful bootstrap
#'   parameter draws in attribute \code{"boot_draws"}.
#'
#' @return A data frame with columns \code{parameter}, \code{estimate}
#'   (original point estimate), \code{se_boot} (bootstrap standard error),
#'   \code{ci_lower}, \code{ci_upper}, \code{mcse_lower}, \code{mcse_upper},
#'   \code{wald_lower}, \code{wald_upper}, and \code{level}. The attribute
#'   \code{"n_success"} gives the number of replicates that converged.
#'   Additional attributes include \code{"R"}, \code{"n_attempted"},
#'   \code{"ci_type"}, and optionally \code{"boot_draws"}.
#'
#' @examples
#' set.seed(42)
#' n <- 80
#' dat <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
#' sim <- brs_sim(
#'   formula = ~ x1 + x2, data = dat,
#'   beta = c(0.2, -0.5, 0.3), phi = 1 / 5, ncuts = 100
#' )
#' fit <- brs(y ~ x1 + x2, data = sim)
#' \donttest{
#' boot <- brs_bootstrap(fit, R = 99, level = 0.95, seed = 1)
#' print(boot)
#' }
#'
#' @references
#' Lopes, J. E. (2024). \emph{Beta Regression for Interval-Censored
#' Scale-Derived Outcomes}. MSc Dissertation, PPGMNE/UFPR.
#'
#' @seealso \code{\link{confint.brs}} for Wald intervals;
#'   \code{\link{brs_sim}} for simulation; \code{\link{brs}} for fitting.
#'
#' @rdname brs_bootstrap
#' @export
brs_bootstrap <- function(object,
                         R = 199L,
                         level = 0.95,
                         seed = NULL,
                         ci_type = c("percentile", "basic", "normal", "bca"),
                         max_tries = NULL,
                         keep_draws = FALSE) {
  if (!inherits(object, "brs")) {
    stop("'object' must be a fitted 'brs' object.", call. = FALSE)
  }
  if (inherits(object, "brsmm")) {
    stop("'brs_bootstrap' does not support 'brsmm' objects.", call. = FALSE)
  }

  R <- as.integer(R)
  if (length(R) != 1L || is.na(R) || R < 10L) {
    stop("'R' must be at least 10.", call. = FALSE)
  }
  if (length(level) != 1L || is.na(level) || level <= 0 || level >= 1) {
    stop("'level' must be in (0, 1).", call. = FALSE)
  }
  ci_type <- match.arg(ci_type)
  keep_draws <- isTRUE(keep_draws)
  if (is.null(max_tries)) {
    max_tries <- max(3L * R, 50L)
  }
  max_tries <- as.integer(max_tries)
  if (length(max_tries) != 1L || is.na(max_tries) || max_tries < R) {
    stop("'max_tries' must be a single integer >= R.", call. = FALSE)
  }

  if (!is.null(seed)) {
    has_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (has_seed) {
      old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    }
    on.exit({
      if (has_seed) {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    }, add = TRUE)
    set.seed(as.integer(seed))
  }

  p <- object$p
  q <- object$q
  par_orig <- object$par
  formula <- object$formula
  data <- object$data
  link <- object$link
  link_phi <- object$link_phi
  ncuts <- object$ncuts
  lim <- object$lim
  repar <- object$repar

  alpha <- 1 - level
  probs <- c(alpha / 2, 1 - alpha / 2)

  # Build parameter vector for simulation: beta and phi (scalar) or zeta (vector)
  beta <- par_orig[seq_len(p)]
  if (q == 1L) {
    phi <- par_orig[p + 1L]
    zeta <- NULL
  } else {
    phi <- NULL
    zeta <- par_orig[p + seq_len(q)]
  }

  boot_par <- matrix(NA_real_, nrow = R, ncol = length(par_orig))
  n_ok <- 0L
  n_attempted <- 0L

  while (n_ok < R && n_attempted < max_tries) {
    n_attempted <- n_attempted + 1L
    sim_r <- tryCatch(
      brs_sim(
        formula = formula,
        data = data,
        beta = beta,
        phi = phi,
        zeta = zeta,
        link = link,
        link_phi = link_phi,
        ncuts = ncuts,
        lim = lim,
        repar = repar
      ),
      error = function(e) NULL
    )
    if (is.null(sim_r)) next

    fit_r <- tryCatch(
      brs(
        formula = formula,
        data = sim_r,
        link = link,
        link_phi = link_phi,
        ncuts = ncuts,
        lim = lim,
        repar = repar
      ),
      error = function(e) NULL
    )
    if (is.null(fit_r) || fit_r$convergence != 0L) next
    if (length(fit_r$par) != length(par_orig)) next
    if (any(!is.finite(fit_r$par))) next

    boot_par[n_ok + 1L, ] <- fit_r$par
    n_ok <- n_ok + 1L
  }

  min_success <- max(10L, ceiling(0.6 * R))
  if (n_ok < min_success) {
    stop(
      "Too few successful bootstrap replicates (", n_ok, "). ",
      "Need at least ", min_success, " successes. ",
      "Increase 'max_tries', simplify the model, or check convergence.",
      call. = FALSE
    )
  }

  boot_par <- boot_par[seq_len(n_ok), , drop = FALSE]
  par_names <- names(par_orig)

  se_boot <- sqrt(.colVars(boot_par))
  q_lo <- apply(boot_par, 2L, stats::quantile, probs = probs[1L], names = FALSE)
  q_hi <- apply(boot_par, 2L, stats::quantile, probs = probs[2L], names = FALSE)
  z <- stats::qnorm(probs[2L])
  mcse <- matrix(NA_real_, nrow = 2L, ncol = length(par_orig))
  ci <- switch(ci_type,
    percentile = {
      for (j in seq_len(ncol(boot_par))) {
        mcse[, j] <- .boot_mcse_limits(boot_par[, j], probs = probs)
      }
      rbind(q_lo, q_hi)
    },
    basic = {
      for (j in seq_len(ncol(boot_par))) {
        mcse_raw <- .boot_mcse_limits(boot_par[, j], probs = rev(1 - probs))
        mcse[, j] <- mcse_raw
      }
      rbind(2 * par_orig - q_hi, 2 * par_orig - q_lo)
    },
    normal = {
      rbind(par_orig - z * se_boot, par_orig + z * se_boot)
    },
    bca = {
      bca <- .boot_bca_ci(
        object = object,
        boot_par = boot_par,
        par_orig = par_orig,
        probs = probs
      )
      mcse <- bca$mcse
      bca$ci
    }
  )
  ci <- pmin(pmax(ci, -Inf), Inf)
  V_wald <- vcov(object, model = "full")
  se_wald <- sqrt(pmax(diag(V_wald), 0))
  wald_ci <- rbind(par_orig - z * se_wald, par_orig + z * se_wald)

  out <- data.frame(
    parameter = par_names,
    estimate  = unname(par_orig),
    se_boot   = unname(se_boot),
    ci_lower  = unname(ci[1L, ]),
    ci_upper  = unname(ci[2L, ]),
    mcse_lower = unname(mcse[1L, ]),
    mcse_upper = unname(mcse[2L, ]),
    wald_lower = unname(wald_ci[1L, ]),
    wald_upper = unname(wald_ci[2L, ]),
    level     = level,
    row.names = NULL
  )

  if (n_ok < R) {
    warning(
      "Only ", n_ok, " successful replicates were obtained in ",
      n_attempted, " attempts (target R = ", R, "). ",
      "Consider increasing 'max_tries' or checking model convergence.",
      call. = FALSE
    )
  }

  attr(out, "n_success") <- n_ok
  attr(out, "R") <- R
  attr(out, "n_attempted") <- n_attempted
  attr(out, "ci_type") <- ci_type
  if (keep_draws) {
    colnames(boot_par) <- names(par_orig)
    attr(out, "boot_draws") <- boot_par
  }
  class(out) <- c("brs_bootstrap", "data.frame")
  out
}


#' @describeIn brs_bootstrap Print method for bootstrap results
#' @param x Object returned by \code{brs_bootstrap}.
#' @param ... Ignored.
#' @export
print.brs_bootstrap <- function(x, ...) {
  cat("Bootstrap confidence intervals\n")
  cat(
    "  Level:", unique(x$level),
    "| CI:", attr(x, "ci_type"),
    "| Successful replicates:", attr(x, "n_success"), "/", attr(x, "R"),
    "| Attempts:", attr(x, "n_attempted"),
    "\n\n"
  )
  print(as.data.frame(x))
  invisible(x)
}


# Column variances (no external dependency)
.colVars <- function(x) {
  n <- nrow(x)
  if (n < 2L) return(rep(NA_real_, ncol(x)))
  cent <- x - rep(colMeans(x), each = n)
  colSums(cent^2) / (n - 1L)
}

# Monte Carlo error approximation for CI endpoints (quantile-based)
.boot_mcse_limits <- function(x, probs) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  n <- length(x)
  if (n < 30L) return(c(NA_real_, NA_real_))
  dens <- stats::density(x, na.rm = TRUE, n = 512)
  out <- rep(NA_real_, length(probs))
  for (k in seq_along(probs)) {
    p <- probs[k]
    q <- as.numeric(stats::quantile(x, probs = p, names = FALSE))
    f_q <- stats::approx(dens$x, dens$y, xout = q, rule = 2)$y
    if (!is.finite(f_q) || f_q <= 0) next
    out[k] <- sqrt((p * (1 - p)) / n) / f_q
  }
  out
}

# BCa intervals with jackknife acceleration
.boot_bca_ci <- function(object, boot_par, par_orig, probs) {
  B <- nrow(boot_par)
  p <- ncol(boot_par)
  n <- nrow(object$data)

  # Bias-correction
  z0 <- vapply(seq_len(p), function(j) {
    prop <- mean(boot_par[, j] < par_orig[j], na.rm = TRUE)
    prop <- min(max(prop, 1 / (2 * B)), 1 - 1 / (2 * B))
    stats::qnorm(prop)
  }, numeric(1))

  # Jackknife acceleration
  jack <- matrix(NA_real_, nrow = n, ncol = p)
  for (i in seq_len(n)) {
    d_i <- object$data[-i, , drop = FALSE]
    fit_i <- tryCatch(
      brs(
        formula = object$formula,
        data = d_i,
        link = object$link,
        link_phi = object$link_phi,
        ncuts = object$ncuts,
        lim = object$lim,
        repar = object$repar
      ),
      error = function(e) NULL
    )
    if (is.null(fit_i) || fit_i$convergence != 0L || length(fit_i$par) != p) next
    jack[i, ] <- fit_i$par
  }
  a <- rep(0, p)
  for (j in seq_len(p)) {
    jj <- jack[, j]
    jj <- jj[is.finite(jj)]
    if (length(jj) < max(20L, ceiling(0.7 * n))) next
    u <- mean(jj) - jj
    num <- sum(u^3)
    den <- 6 * (sum(u^2)^(3 / 2))
    if (is.finite(den) && den > 0) {
      a[j] <- num / den
    }
  }

  z_alpha <- stats::qnorm(probs)
  ci <- matrix(NA_real_, nrow = 2L, ncol = p)
  mcse <- matrix(NA_real_, nrow = 2L, ncol = p)
  for (j in seq_len(p)) {
    adj <- stats::pnorm(z0[j] + (z0[j] + z_alpha) / (1 - a[j] * (z0[j] + z_alpha)))
    adj <- pmin(pmax(adj, 1 / (B + 1)), B / (B + 1))
    ci[, j] <- as.numeric(stats::quantile(boot_par[, j], probs = adj, names = FALSE, na.rm = TRUE))
    mcse[, j] <- .boot_mcse_limits(boot_par[, j], probs = adj)
  }
  list(ci = ci, mcse = mcse)
}
