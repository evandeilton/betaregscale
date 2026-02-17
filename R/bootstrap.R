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
#'
#' @return A data frame with columns \code{parameter}, \code{estimate}
#'   (original point estimate), \code{se_boot} (bootstrap standard error),
#'   \code{ci_lower}, \code{ci_upper}, and \code{level}. The attribute
#'   \code{"n_success"} gives the number of replicates that converged.
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
                         seed = NULL) {
  if (!inherits(object, "brs")) {
    stop("'object' must be a fitted 'brs' object.", call. = FALSE)
  }
  if (inherits(object, "brsmm")) {
    stop("'brs_bootstrap' does not support 'brsmm' objects.", call. = FALSE)
  }

  R <- as.integer(R)
  if (R < 10L) {
    stop("'R' must be at least 10.", call. = FALSE)
  }
  if (level <= 0 || level >= 1) {
    stop("'level' must be in (0, 1).", call. = FALSE)
  }

  if (!is.null(seed)) {
    set.seed(seed)
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

  for (r in seq_len(R)) {
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

    boot_par[n_ok + 1L, ] <- fit_r$par
    n_ok <- n_ok + 1L
  }

  if (n_ok < 10L) {
    stop(
      "Too few successful bootstrap replicates (", n_ok, "). ",
      "Increase R or check model convergence.",
      call. = FALSE
    )
  }

  boot_par <- boot_par[seq_len(n_ok), , drop = FALSE]
  par_names <- names(par_orig)

  se_boot <- sqrt(.colVars(boot_par))
  ci <- apply(boot_par, 2L, stats::quantile, probs = probs, names = FALSE)

  out <- data.frame(
    parameter = par_names,
    estimate  = unname(par_orig),
    se_boot   = unname(se_boot),
    ci_lower  = unname(ci[1L, ]),
    ci_upper  = unname(ci[2L, ]),
    level     = level,
    row.names = NULL
  )

  if (n_ok < R) {
    warning(
      "Only ", n_ok, " of ", R, " replicates converged. ",
      "Consider increasing R or checking the model.",
      call. = FALSE
    )
  }

  attr(out, "n_success") <- n_ok
  attr(out, "R") <- R
  class(out) <- c("brs_bootstrap", "data.frame")
  out
}


#' @describeIn brs_bootstrap Print method for bootstrap results
#' @param x Object returned by \code{brs_bootstrap}.
#' @param ... Ignored.
#' @export
print.brs_bootstrap <- function(x, ...) {
  cat("Bootstrap confidence intervals\n")
  cat("  Level:", unique(x$level), "| Successful replicates:", attr(x, "n_success"), "/", attr(x, "R"), "\n\n")
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
