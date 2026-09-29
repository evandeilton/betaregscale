# ============================================================================ #
# Fit diagnostics shared by brs() and brsmm(): design checks before optim,
# gradient/Hessian/clamp checks after it, the variance matrix, and the
# advisory-warning plumbing used by refits (bootstrap, jackknife).
# ============================================================================ #

# Advisory warning with a class, so refits can muffle it (see .brs_quiet_advisory).
.brs_advisory <- function(msg, class = "brs_fit_diagnostic") {
  warning(structure(
    class = c(class, "brs_advisory", "warning", "condition"),
    list(message = msg, call = NULL)
  ))
}

# Evaluate `expr` without the package's advisory warnings (lim, mixed values,
# fit diagnostics); real errors and other warnings still propagate.
.brs_quiet_advisory <- function(expr) {
  withCallingHandlers(expr, brs_advisory = function(w) {
    invokeRestart("muffleWarning")
  })
}

# One warning per R session for a given id (documented approximations).
.brs_once <- new.env(parent = emptyenv())
.brs_warn_once <- function(id, msg) {
  if (isTRUE(.brs_once[[id]])) return(invisible(FALSE))
  assign(id, TRUE, envir = .brs_once)
  warning(msg, call. = FALSE)
  invisible(TRUE)
}

# Run `expr` and restore the RNG state afterwards (summary() must not consume it).
.brs_keep_seed <- function(expr) {
  env <- globalenv()
  had <- exists(".Random.seed", envir = env, inherits = FALSE)
  old <- if (had) get(".Random.seed", envir = env, inherits = FALSE) else NULL
  on.exit({
    if (had) {
      assign(".Random.seed", old, envir = env)
    } else if (exists(".Random.seed", envir = env, inherits = FALSE)) {
      rm(".Random.seed", envir = env)
    }
  })
  expr
}

#' Design-matrix checks before fitting
#'
#' @description
#' Stops on exact rank deficiency (pivoted QR, the \code{lm()} tolerance
#' \code{1e-7}) naming the aliased columns; warns on near collinearity when
#' the condition number of the column-equilibrated matrix (unit-length
#' columns, Belsley's scaling) exceeds \code{1e4}.
#'
#' @param M Numeric model matrix.
#' @param what Label used in the message (e.g. \code{"mean"}).
#' @return \code{M}, invisibly.
#' @keywords internal
#' @noRd
.brs_check_design <- function(M, what) {
  k <- ncol(M)
  if (is.null(k) || k == 0L || nrow(M) == 0L) return(invisible(M))
  if (any(!is.finite(M))) {
    stop("The ", what, " model matrix contains non-finite values ",
         "(check log(), division or missing codes in the covariates).", call. = FALSE)
  }
  nms <- colnames(M)
  if (is.null(nms)) nms <- paste0("V", seq_len(k))
  qx <- qr(M, tol = 1e-7)
  if (qx$rank < k) {
    alias <- nms[qx$pivot[seq.int(qx$rank + 1L, k)]]
    stop("The ", what, " model matrix is rank deficient (rank ", qx$rank, " < ", k,
         "): column(s) ", paste0("'", alias, "'", collapse = ", "),
         " are linear combinations of the others; remove them.", call. = FALSE)
  }
  if (k >= 2L) {
    # Unit-length columns: the condition number no longer depends on covariate units
    Ms <- sweep(M, 2L, sqrt(colSums(M^2)), "/")
    sv <- svd(Ms, nu = 0L)
    kappa <- sv$d[1L] / sv$d[k]
    if (!is.finite(kappa) || kappa > 1e4) {
      v <- abs(sv$v[, k])
      involved <- nms[v >= 0.1 * max(v)]
      .brs_advisory(paste0(
        "The ", what, " model matrix is nearly collinear (condition number ",
        format(signif(kappa, 3)), "; columns ",
        paste0("'", involved, "'", collapse = ", "),
        "): estimates and SEs may be unreliable."
      ))
    }
  }
  invisible(M)
}

#' Post-fit diagnostics (one helper for brs and brsmm)
#'
#' @description
#' Uses only what the fit carries: the gradient and the Hessian of the
#' log-likelihood at the estimate (any source) and the inverse-linked
#' parameters before the clamps of the compiled likelihood.
#' \itemize{
#'   \item gradient: log-likelihood gain of the Newton step left to the
#'     optimum, \eqn{\frac12 g^\top (-H)^{-1} g} (scale free, the same rule
#'     for brs and brsmm); warns above \eqn{10^{-2}} when optim reported
#'     convergence (otherwise optim already warned). The largest step in SE
#'     units is stored too;
#'   \item Hessian: \eqn{-H} must be positive definite, also in correlation
#'     form (smallest eigenvalue of \eqn{D^{-1/2}(-H)D^{-1/2}} above
#'     \code{sqrt(.Machine$double.eps)}, a scale-free singularity check);
#'   \item clamps: observations whose mean/shape p, dispersion/precision or
#'     beta shapes sit on the clamps of \code{src/brs_common.h};
#'   \item brsmm: a random-effect log SD below -6, or a log-likelihood gain
#'     below 1e-3 over the same fit with that term's SD at 0 (a flat ridge
#'     where optim stops before the boundary).
#' }
#'
#' @param grad Gradient of the log-likelihood at the estimate.
#' @param hessian Hessian of the log-likelihood at the estimate.
#' @param mu_raw,phi_raw Inverse-linked first/second parameter, unclamped.
#' @param repar Parameterisation.
#' @param convergence optim convergence code.
#' @param re_logsd Optional log SDs of the random effects (brsmm).
#' @param re_gain Optional log-likelihood gains of each random-effect term
#'   over its removal (brsmm).
#' @param badly_scaled Logical: the design has columns of very different
#'   scales (the gradient warning then suggests rescaling).
#' @param warn Emit the warnings.
#' @return A list stored as \code{fit$diagnostics}.
#' @keywords internal
#' @noRd
.brs_fit_diagnostics <- function(grad, hessian, mu_raw, phi_raw, repar,
                                 convergence = 0L, re_logsd = NULL, re_gain = NULL,
                                 badly_scaled = FALSE, warn = TRUE) {
  H <- -hessian
  H <- (H + t(H)) / 2
  finite <- all(is.finite(H)) && all(is.finite(grad))
  ev <- if (finite) eigen(H, symmetric = TRUE, only.values = TRUE)$values else NA_real_
  dg <- if (finite) diag(H) else NA_real_
  # Correlation form: eigenvalues free of the parameter (covariate) scales
  sev <- if (finite && all(dg > 0)) {
    min(eigen(H / sqrt(outer(dg, dg)), symmetric = TRUE, only.values = TRUE)$values)
  } else {
    -Inf
  }
  nd_ok <- finite && min(ev) > 0 && sev > sqrt(.Machine$double.eps)

  # Newton step left to the optimum: its log-likelihood gain (the criterion)
  # and its size in SE units; meaningful only when -H is positive definite
  step <- NA_real_
  gain <- NA_real_
  if (nd_ok) {
    V <- tryCatch(solve(H), error = function(e) NULL)
    if (!is.null(V) && all(diag(V) > 0)) {
      s_newton <- drop(V %*% grad)
      step <- max(abs(s_newton) / sqrt(diag(V)))
      gain <- 0.5 * sum(grad * s_newton)
    }
  }

  # Clamps of the compiled likelihood (mean/shape p, second parameter, shapes)
  n <- length(mu_raw)
  phi_raw <- rep_len(as.numeric(phi_raw), n)
  mu_c <- .clamp_mu_by_repar(mu_raw, repar)
  phi_c <- .clamp_phi_by_repar(phi_raw, repar)
  on_mu <- !is.finite(mu_raw) | mu_raw != mu_c
  on_phi <- !is.finite(phi_raw) | phi_raw != phi_c
  sh <- .brs_shapes_unclamped(mu_c, phi_c, repar)
  on_sh <- !is.finite(sh$a) | !is.finite(sh$b) |
    sh$a < 1e-12 | sh$b < 1e-12 | sh$a > 1e8 | sh$b > 1e8
  on_any <- on_mu | on_phi | on_sh

  # Boundary: tiny SD, or the term adds nothing to the likelihood
  re_bnd <- any(re_logsd < -6) || any(re_gain < 1e-3)

  out <- list(
    grad_norm = if (length(grad)) max(abs(grad)) else 0,
    grad_step = step,
    grad_gain = gain,
    min_eig = min(ev),
    max_eig = max(ev),
    min_eig_scaled = sev,
    hessian_nd = nd_ok,
    n_clamped = sum(on_any),
    clamped = c(mu = sum(on_mu), phi = sum(on_phi), shape = sum(on_sh)),
    re_boundary = re_bnd,
    re_gain = re_gain
  )

  if (warn) {
    # On the boundary a non-zero gradient is expected: that warning explains it
    if (identical(as.integer(convergence), 0L) && is.finite(gain) && gain > 1e-2 &&
        !re_bnd) {
      .brs_advisory(sprintf(paste0(
        "Gradient not ~0 at the estimate (a Newton step would gain %.3g in ",
        "log-likelihood): the optimizer stopped early; %s."), gain,
        if (badly_scaled) "rescale the covariates (their scales differ widely) and refit"
        else "refit with the other method or rescaled covariates"))
    }
    if (!nd_ok) {
      .brs_advisory(sprintf(paste0(
        "Hessian not negative definite (SEs unreliable): smallest eigenvalue ",
        "of -H = %.3g (%.3g in correlation form)."), min(ev), sev))
    }
    if (out$n_clamped > 0L) {
      .brs_advisory(sprintf(paste0(
        "%d of %d observations have mu, phi or the beta shapes on the clamp ",
        "boundary (possible non-identifiability)."), out$n_clamped, n))
    }
    if (re_bnd) {
      # A negative gain: sd ~ 0 beats the estimate, so optim stopped short of the MLE
      short <- any(re_gain < 0)
      .brs_advisory(paste0(
        "Variance component at the boundary (log sd < -6 or no likelihood ",
        "gain over sd = 0)",
        if (short) paste0(": sd ~ 0 has a higher log-likelihood than the estimate ",
                          "(the interior estimate is not the MLE; the fit stopped short)") else "",
        ". Use the LRT with the 1/2 chi2(Df - 1) + 1/2 chi2(Df) mixture (anova())."))
    }
  }
  out
}

# Columns whose root mean squares differ by more than 1e3 (e.g. x^3 of x ~ 1e5).
.brs_badly_scaled <- function(...) {
  rms <- unlist(lapply(list(...), function(M) {
    if (is.null(M)) NULL else sqrt(colMeans(M^2))
  }))
  rms <- rms[is.finite(rms) & rms > 0]
  length(rms) > 1L && max(rms) / min(rms) > 1e3
}

# Beta shapes from (clamped) mu and phi before their own clamp (src beta_shapes()).
.brs_shapes_unclamped <- function(mu, phi, repar) {
  switch(as.character(as.integer(repar)),
    "0" = list(a = mu, b = phi),
    "1" = list(a = mu * phi, b = (1 - mu) * phi),
    list(a = mu * (1 - phi) / phi, b = (1 - mu) * (1 - phi) / phi)
  )
}

# Central-difference gradient for brsmm (no compiled gradient). Step 1e-3: the
# marginal likelihood has an inner mode search, so tiny steps measure its tolerance.
.brs_num_grad <- function(fn, par, h_rel = 1e-3) {
  vapply(seq_along(par), function(j) {
    h <- h_rel * max(1, abs(par[j]))
    e <- replace(numeric(length(par)), j, h)
    (fn(par + e) - fn(par - e)) / (2 * h)
  }, numeric(1))
}

#' Variance matrix from the Hessian (vcov.brs and vcov.brsmm)
#'
#' @description
#' \code{solve(-H)}; no generalised inverse. A singular Hessian gives an
#' \code{NA} matrix; negative or non-finite variances become \code{NA}
#' (row and column). Both cases warn.
#'
#' @param hessian Hessian of the log-likelihood.
#' @param nms Parameter names.
#' @return Covariance matrix.
#' @keywords internal
#' @noRd
.brs_vcov <- function(hessian, nms) {
  V <- tryCatch(solve(-hessian), error = function(e) NULL)
  if (is.null(V)) {
    warning("Hessian is singular: the variance matrix is not available (NA). ",
            "See fit$diagnostics.", call. = FALSE)
    V <- matrix(NA_real_, nrow(hessian), ncol(hessian))
  } else {
    bad <- !is.finite(diag(V)) | diag(V) <= 0
    if (any(bad)) {
      warning(sum(bad), " negative or non-finite variance(s) set to NA (",
              paste(nms[bad], collapse = ", "),
              "): the Hessian is not negative definite.", call. = FALSE)
      V[bad, ] <- NA_real_
      V[, bad] <- NA_real_
    }
  }
  dimnames(V) <- list(nms, nms)
  V
}
