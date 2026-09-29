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
#' Each replicate draws a new response \eqn{y^*_i \sim \mathrm{Beta}(a_i, b_i)}
#' at the fitted shapes of the rows used by the fit and substitutes it into a
#' copy of \code{object$data}; covariates, factor levels, transformations and
#' the formula stay those of the original fit, which is then re-fitted with
#' \code{\link{brs}} under the same \code{ncuts}, \code{lim}, \code{interval},
#' links, \code{repar} and optimizer. The observation mechanism of each row is
#' reproduced: exact observations (\eqn{\delta = 0}) stay continuous; scores
#' are re-coarsened on the fit's grid (the mapping of \code{\link{brs_sim}}).
#' Rows whose bounds came from the analyst (\code{\link{brs_prep}} Modes 2--4,
#' bounds that are not the cell of the row's score) keep their thresholds as
#' fixed and non-informative (independent of \eqn{Y}), and \eqn{\delta} is
#' re-drawn by the cell of the partition they induce where \eqn{y^*} falls:
#' for an upper bound \eqn{c}, \eqn{y^* \le c} gives \eqn{\delta = 1} on
#' \eqn{[\epsilon, c]}, otherwise \eqn{\delta = 2} on \eqn{[c, 1 - \epsilon]}
#' (an interval \eqn{[l, u]} gives three cells). When the original design had
#' more thresholds than a row records (e.g. Mode 4 intervals cut from a finer
#' instrument) the bootstrap is conservative, and for Mode 2 rows with a
#' forced \eqn{\delta} (a threshold built from the score itself, which is
#' informative) it is only an approximation. The response must be a variable
#' (not an expression such as \code{I(y / 10)}).
#'
#' Each refit starts from the estimate of \code{object} and uses the compiled
#' Hessian (\code{hessian_method = "cpp"}), so replicates are cheap.
#'
#' Replicates that fail (refit error, non-convergence, non-finite estimates)
#' are discarded and counted: attributes \code{"n_failed"} and
#' \code{"fail_rate"}, also printed. Intervals are computed from the bootstrap
#' distribution of each parameter, with the method controlled by
#' \code{ci_type}: \code{"percentile"} (default) uses the raw empirical
#' quantiles; \code{"basic"} uses reflected empirical quantiles;
#' \code{"normal"} uses a normal approximation from the bootstrap standard
#' error (no quantiles); \code{"bca"} uses bias-corrected-and-accelerated
#' adjusted quantiles. With parametric resampling and a nonparametric
#' (leave-one-out) jackknife acceleration, \code{"bca"} is an approximation;
#' a warning says so once per session, and failed jackknife refits are
#' reported in \code{"n_jack_failed"}.
#'
#' @section Cost of \code{ci_type = "bca"}:
#' The bias-corrected and accelerated interval needs an acceleration constant,
#' which is obtained here by a leave-one-out jackknife. That requires \code{n}
#' additional model fits, one per observation, on top of the \code{R} bootstrap
#' replicates: the total is \code{R + n} fits rather than \code{R}. The cost is
#' therefore driven by the sample size, not by \code{R}, and grows quickly --
#' for \code{n = 1000} the jackknife alone dominates the run time by an order of
#' magnitude. The other three interval types need only the \code{R} replicates.
#' Prefer \code{"percentile"} or \code{"basic"} for exploratory work on large
#' samples, and reserve \code{"bca"} for a final result.
#'
#' @param object A fitted \code{"brs"} object (fixed or variable dispersion).
#' @param R Integer: number of bootstrap replicates (default 199).
#' @param level Numeric: confidence level (default 0.95).
#' @param ci_type Character: type of confidence interval. One of
#'   \code{"percentile"} (default), \code{"basic"}, \code{"normal"},
#'   or \code{"bca"}. See the section on the cost of \code{"bca"} below.
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
#'   \code{"n_failed"}, \code{"fail_rate"}, \code{"n_jack_failed"} (BCa
#'   only), \code{"ci_type"}, and optionally \code{"boot_draws"}.
#'
#' @examples
#' # Synthetic NRS-11 scores at 6h, 12h, 24h (time is a factor)
#' set.seed(3)
#' nrs <- data.frame(time = factor(rep(c("6h", "12h", "24h"), each = 40),
#'                                 levels = c("6h", "12h", "24h")))
#' shp <- brs_repar(mu = plogis(-1.3 + c(0, 0.75, 0.3)[nrs$time]), phi = 0.3)
#' nrs$y <- round(10 * rbeta(nrow(nrs), shp$shape1, shp$shape2))
#' fit <- brs(y ~ time, data = nrs, ncuts = 10)
#'
#' # Percentile intervals from 30 parametric replicates (use R >= 199 in practice)
#' set.seed(4)
#' bt <- brs_bootstrap(fit, R = 30)
#' bt
#' # Bootstrap next to Wald limits, and the bootstrap/Wald SE ratio
#' cols <- c("parameter", "ci_lower", "ci_upper", "wald_lower", "wald_upper")
#' as.data.frame(bt)[, cols]
#' round(bt$se_boot / sqrt(diag(vcov(fit))), 2)
#'
#' @seealso \code{\link{confint.brs}} for Wald intervals;
#'   \code{\link{brs_sim}} for simulation; \code{\link{brs}} for fitting.
#'
#' @rdname brs_bootstrap
#' @export
brs_bootstrap <- function(object,
                          R = 199L,
                          level = 0.95,
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

  par_orig <- object$par
  # Rows, response and observation mechanism of the fit (once, not per replicate)
  setup <- .brs_boot_setup(object)

  alpha <- 1 - level
  probs <- c(alpha / 2, 1 - alpha / 2)

  boot_par <- matrix(NA_real_, nrow = R, ncol = length(par_orig))
  n_ok <- 0L
  n_attempted <- 0L

  while (n_ok < R && n_attempted < max_tries) {
    n_attempted <- n_attempted + 1L
    # Only the response changes; covariates, factors and formula stay the fit's own
    data_r <- .brs_boot_data(setup, object)
    fit_r <- .brs_refit(object, data_r)
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
      # Lower basic limit = 2 theta - upper quantile (and vice versa): reversed MCSE
      for (j in seq_len(ncol(boot_par))) {
        mcse[, j] <- .boot_mcse_limits(boot_par[, j], probs = rev(probs))
      }
      rbind(2 * par_orig - q_hi, 2 * par_orig - q_lo)
    },
    normal = {
      rbind(par_orig - z * se_boot, par_orig + z * se_boot)
    },
    bca = {
      .brs_warn_once("bca_approx", paste0(
        "ci_type = \"bca\" is an approximation here: parametric resampling with a ",
        "nonparametric (leave-one-out) jackknife acceleration. Shown once per session."
      ))
      bca <- .boot_bca_ci(
        object = object,
        boot_par = boot_par,
        par_orig = par_orig,
        probs = probs,
        setup = setup
      )
      mcse <- bca$mcse
      n_jack_failed <- bca$n_failed
      bca$ci
    }
  )
  V_wald <- vcov(object, model = "full")
  se_wald <- sqrt(diag(V_wald))
  wald_ci <- rbind(par_orig - z * se_wald, par_orig + z * se_wald)

  out <- data.frame(
    parameter = par_names,
    estimate = unname(par_orig),
    se_boot = unname(se_boot),
    ci_lower = unname(ci[1L, ]),
    ci_upper = unname(ci[2L, ]),
    mcse_lower = unname(mcse[1L, ]),
    mcse_upper = unname(mcse[2L, ]),
    wald_lower = unname(wald_ci[1L, ]),
    wald_upper = unname(wald_ci[2L, ]),
    level = level,
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
  # Failed replicates: simulation, refit error, non-convergence or bad estimates
  attr(out, "n_failed") <- n_attempted - n_ok
  attr(out, "fail_rate") <- (n_attempted - n_ok) / n_attempted
  if (identical(ci_type, "bca")) attr(out, "n_jack_failed") <- n_jack_failed
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
    "\n"
  )
  nf <- attr(x, "n_failed")
  if (!is.null(nf)) {
    cat(sprintf("  Failed replicates: %d (%.1f%% of attempts)", nf,
                100 * attr(x, "fail_rate")))
    nj <- attr(x, "n_jack_failed")
    if (!is.null(nj)) cat(" | Failed jackknife refits:", nj)
    cat("\n")
  }
  cat("\n")
  print(as.data.frame(x))
  invisible(x)
}


# Column variances (no external dependency)
.colVars <- function(x) {
  n <- nrow(x)
  if (n < 2L) {
    return(rep(NA_real_, ncol(x)))
  }
  cent <- x - rep(colMeans(x), each = n)
  colSums(cent^2) / (n - 1L)
}

# Monte Carlo error approximation for CI endpoints (quantile-based)
.boot_mcse_limits <- function(x, probs) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  n <- length(x)
  if (n < 30L) {
    return(c(NA_real_, NA_real_))
  }
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

# BCa intervals with jackknife acceleration (over the rows used by the fit)
.boot_bca_ci <- function(object, boot_par, par_orig, probs, setup) {
  B <- nrow(boot_par)
  p <- ncol(boot_par)
  rows <- setup$rows
  n <- length(rows)

  # Bias-correction
  z0 <- vapply(seq_len(p), function(j) {
    prop <- mean(boot_par[, j] < par_orig[j], na.rm = TRUE)
    prop <- min(max(prop, 1 / (2 * B)), 1 - 1 / (2 * B))
    stats::qnorm(prop)
  }, numeric(1))

  # Jackknife acceleration: leave one fitted row out of the data (response included)
  jack <- matrix(NA_real_, nrow = n, ncol = p)
  for (i in seq_len(n)) {
    fit_i <- .brs_refit(object, setup$data[-rows[i], , drop = FALSE])
    if (is.null(fit_i) || fit_i$convergence != 0L || length(fit_i$par) != p) next
    jack[i, ] <- fit_i$par
  }
  n_failed <- sum(!stats::complete.cases(jack))
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
  list(ci = ci, mcse = mcse, n_failed = n_failed)
}


# -- Refit and response simulation ---------------------------------------- #

# THE refit of a brs model on new data (bootstrap replicates and jackknife):
# the fit's own settings; advisory warnings muffled; NULL on error.
.brs_refit <- function(object, data) {
  meth <- if (is.null(object$method)) "BFGS" else object$method
  # Warm start from the parent estimate; compiled (cpp) Hessian, the brs() default
  tryCatch(
    suppressMessages(.brs_quiet_advisory(brs(
      formula = object$formula,
      data = data,
      link = object$link,
      link_phi = object$link_phi,
      ncuts = object$ncuts,
      lim = object$lim,
      repar = object$repar,
      method = meth,
      hessian_method = "cpp",
      interval = .brs_interval_of(object),
      start = unname(object$par)
    ))),
    error = function(e) NULL
  )
}

#' Observation mechanism of a brs fit, for the parametric bootstrap
#'
#' @description
#' Rows of \code{object$data} used by the fit, the response column (added to
#' the data when it lives in the formula environment), fitted beta shapes, and
#' the mechanism of each row: \code{"exact"} (\eqn{\delta = 0}, stays
#' continuous), \code{"score"} (cell of a score on the fit's grid, re-coarsened
#' with \code{ncuts}/\code{lim}/\code{interval}) or \code{"analyst"} (brs_prep
#' thresholds that are not the cell of the row's score: thresholds kept, the
#' new value reported by its side of them).
#'
#' @param object A \code{"brs"} fit.
#' @return A list used by \code{.brs_boot_data()}.
#' @keywords internal
#' @noRd
.brs_boot_setup <- function(object) {
  data <- object$data
  f <- Formula::as.Formula(object$formula)
  lhs <- formula(f)[[2L]]
  if (!is.name(lhs)) {
    stop("brs_bootstrap() needs the response to be a variable (found '",
         deparse(lhs), "' on the left-hand side); create that column first.",
         call. = FALSE)
  }
  resp <- as.character(lhs)
  # Materialise a response kept outside `data` (formula environment)
  if (!resp %in% names(data)) {
    v <- eval(lhs, data, environment(f))
    if (length(v) != nrow(data)) {
      stop("The response '", resp, "' does not match the rows of the fitted data.",
           call. = FALSE)
    }
    data[[resp]] <- v
  }
  mf <- stats::model.frame(f, data = data)
  rows <- match(rownames(mf), rownames(data))
  if (anyNA(rows) || length(rows) != object$nobs) {
    stop("Unable to align the fitted rows with `object$data`.", call. = FALSE)
  }

  prepared <- isTRUE(attr(data, "is_prepared")) &&
    all(c("left", "right", "yt", "delta") %in% names(data))
  K <- as.integer(object$ncuts)
  interval <- .brs_interval_of(object)
  Y <- object$Y
  delta <- as.integer(object$delta)
  mech <- ifelse(delta == 0L, "exact", "score")
  if (prepared) {
    # A censored row is a score row only if it is exactly the cell of its score
    y <- as.numeric(Y[, "y"])
    is_score <- delta != 0L & !is.na(y) & abs(y - round(y)) < 1e-8 &
      y >= 0 & y <= K
    same <- rep(FALSE, length(y))
    if (any(is_score)) {
      auto <- .brs_check_core(round(y[is_score]), ncuts = K, lim = object$lim,
                              delta = NULL, interval = interval)
      same[is_score] <- auto[, "delta"] == delta[is_score] &
        abs(auto[, "left"] - Y[is_score, "left"]) < 1e-9 &
        abs(auto[, "right"] - Y[is_score, "right"]) < 1e-9
    }
    mech[delta != 0L & !same] <- "analyst"
  }

  # Fitted shapes, clamped as in the compiled likelihood
  n <- object$nobs
  sh <- brs_repar(object$hatmu, rep_len(object$hatphi, n), repar = object$repar)
  list(
    data = data, rows = rows, resp = resp, prepared = prepared, mech = mech,
    shape1 = pmin(pmax(sh$shape1, 1e-12), 1e8),
    shape2 = pmin(pmax(sh$shape2, 1e-12), 1e8),
    y_unit = as.numeric(Y[, "y"]) > 0 & as.numeric(Y[, "y"]) < 1,
    left = as.numeric(Y[, "left"]), right = as.numeric(Y[, "right"]),
    delta = delta, K = K, lim = object$lim, interval = interval
  )
}

# One bootstrap data set: a new response for the fitted rows, same covariates.
.brs_boot_data <- function(setup, object) {
  eps <- 1e-5
  n <- length(setup$rows)
  ys <- stats::rbeta(n, setup$shape1, setup$shape2)
  yc <- pmin(pmax(ys, eps), 1 - eps)
  K <- setup$K
  iv <- setup$interval
  ex <- setup$mech == "exact"
  sc <- setup$mech == "score"
  an <- setup$mech == "analyst"
  # Scores by the likelihood's own coarsening (Lote 4 mapping)
  s <- .brs_score_from_unit(yc, K, iv)
  data <- setup$data
  rows <- setup$rows

  if (!setup$prepared) {
    # Raw data: exact rows stay in (0, 1), score rows get a new score
    v <- data[[setup$resp]]
    v[rows[ex]] <- yc[ex]
    v[rows[sc]] <- s[sc]
    data[[setup$resp]] <- v
    return(data)
  }

  left <- setup$left
  right <- setup$right
  yt <- yc
  d <- setup$delta
  ycol <- ifelse(setup$y_unit, yc, .brs_latent_score(yc, K, iv))
  left[ex] <- yc[ex]
  right[ex] <- yc[ex]
  if (any(sc)) {
    cells <- .brs_check_core(s[sc], ncuts = K, lim = setup$lim, delta = NULL,
                             interval = iv)
    left[sc] <- cells[, "left"]
    right[sc] <- cells[, "right"]
    yt[sc] <- cells[, "yt"]
    d[sc] <- as.integer(cells[, "delta"])
    ycol[sc] <- s[sc]
  }
  if (any(an)) {
    # Fixed analyst thresholds: report the side (or the interval) holding y*
    l0 <- setup$left[an]
    u0 <- setup$right[an]
    d0 <- setup$delta[an]
    y0 <- yc[an]
    i1 <- d0 == 1L
    i2 <- d0 == 2L
    i3 <- d0 == 3L
    cut <- ifelse(i2, l0, u0)
    below <- (i1 & y0 <= cut) | (i2 & y0 < cut) | (i3 & y0 < l0)
    above <- (i1 & y0 > cut) | (i2 & y0 >= cut) | (i3 & y0 > u0)
    # below: [eps, lower threshold]; above: [upper threshold, 1 - eps]; else [l0, u0]
    new_d <- ifelse(below, 1L, ifelse(above, 2L, 3L))
    new_l <- ifelse(below, eps, ifelse(above, ifelse(i3, u0, cut), l0))
    new_r <- ifelse(below, ifelse(i3, l0, cut), ifelse(above, 1 - eps, u0))
    new_yt <- ifelse(new_d == 1L, new_r / 2,
      ifelse(new_d == 2L, (new_l + 1) / 2, (new_l + new_r) / 2))
    left[an] <- new_l
    right[an] <- new_r
    yt[an] <- new_yt
    d[an] <- new_d
    ycol[an] <- .brs_latent_score(new_yt, K, iv)
  }
  data[["left"]][rows] <- left
  data[["right"]][rows] <- right
  data[["yt"]][rows] <- yt
  data[["delta"]][rows] <- d
  data[[setup$resp]][rows] <- ycol
  data
}
