# ============================================================================ #
# Random-effects numeric study utilities for brsmm
# ============================================================================ #

#' Random-effects study for brsmm models
#'
#' @description
#' Provides a compact numeric study of random effects, including:
#' estimated covariance matrix, correlation matrix, per-term standard
#' deviations, empirical mean/SD of posterior modes, shrinkage ratio, and
#' a normality check by Shapiro-Wilk (when applicable).
#'
#' @details
#' \code{icc} is the intraclass correlation of \eqn{\mathrm{logit}(Y)} implied
#' by the fitted model: for two observations of the same group with the
#' covariates of observation \eqn{i},
#' \deqn{\mathrm{ICC}_i = \frac{\mathrm{Var}_b[\psi(a_i) - \psi(b_i)]}
#'   {\mathrm{Var}_b[\psi(a_i) - \psi(b_i)] + E_b[\psi_1(a_i) + \psi_1(b_i)]},}
#' where \eqn{a_i(b), b_i(b)} are the beta shapes with random part
#' \eqn{b \sim N(0, x_{r,i}^\top D x_{r,i})}, and \eqn{\psi},
#' \eqn{\psi_1} are the digamma and trigamma functions
#' (\eqn{E[\mathrm{logit}\,Y] = \psi(a) - \psi(b)},
#' \eqn{\mathrm{Var}[\mathrm{logit}\,Y] = \psi_1(a) + \psi_1(b)}). The
#' expectations over \eqn{b} use 40-point Gauss-Hermite quadrature and the
#' reported value is the mean of \eqn{\mathrm{ICC}_i} over the observations.
#' The level-1 variance is that of the beta, so the value depends on the
#' precision; with the logit link (only) and a large precision it approaches
#' \eqn{\sigma_b^2 / (\sigma_b^2 + \psi_1(a) + \psi_1(b))}. It replaces the
#' logistic-latent formula \eqn{\sigma_b^2 / (\sigma_b^2 + \pi^2/3)}, which
#' does not describe a beta response.
#'
#' The moments of \eqn{\mathrm{logit}(Y)} over \eqn{b} can be infinite: with a
#' probit link when \eqn{\sigma_b^2 \ge 1/2}, with a cloglog link for every
#' \eqn{\sigma_b > 0}, and numerically with any link when \eqn{\sigma_b} is
#' very large. The value would then be set by the clamp of the mean
#' (\eqn{10^{-5}}), so \code{icc} is \code{NA}, with a warning, whenever
#' \eqn{E_b[\psi_1(a) + \psi_1(b)]} changes by more than 10\% when that clamp
#' is tightened to \eqn{10^{-8}}.
#'
#' @param object A fitted \code{"brsmm"} object.
#' @param ... Currently ignored.
#'
#' @return A list with class \code{"brsmm_re_study"}.
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
#' @seealso \code{\link{brsmm}}, \code{\link{ranef.brsmm}}
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
#' fit <- brsmm(y ~ x1, random = ~ 1 | id, data = prep)
#' rs <- brsmm_re_study(fit)
#' print(rs)
#' rs$summary
#' }
#'
#' @importFrom stats cov2cor sd shapiro.test
#' @export
brsmm_re_study <- function(object, ...) {
  .check_class_mm(object)

  re <- object$random$mode_b
  if (is.matrix(re)) {
    B <- re
  } else {
    B <- matrix(as.numeric(re), ncol = 1L)
    rownames(B) <- names(re)
    cn <- object$random$terms
    if (is.null(cn) || length(cn) == 0L) cn <- "(Intercept)"
    colnames(B) <- cn[1L]
  }

  D <- object$random$D
  if (is.null(D)) {
    sd_single <- object$random$sd_b
    D <- matrix(as.numeric(sd_single)^2, nrow = 1L, ncol = 1L)
  }
  # Term names on D (the fit stores it without them; print showed re1, re2)
  dimnames(D) <- list(colnames(B), colnames(B))
  Corr <- stats::cov2cor(D)

  mode_mean <- colMeans(B)
  mode_sd <- apply(B, 2, stats::sd)
  mode_var <- pmax(mode_sd^2, 0)
  model_var <- pmax(diag(D), 1e-12)
  shrinkage_ratio <- pmin(pmax(mode_var / model_var, 0), 1)

  shapiro_p <- rep(NA_real_, ncol(B))
  if (nrow(B) >= 3L && nrow(B) <= 5000L) {
    shapiro_p <- apply(B, 2, function(x) {
      x_num <- as.numeric(x)
      # Safe-guard against singular fits (identical modes = 0 variance)
      if (stats::sd(x_num) > 1e-6) {
        stats::shapiro.test(x_num)$p.value
      } else {
        NA_real_
      }
    })
  }

  summary_df <- data.frame(
    term = colnames(B),
    sd_model = sqrt(model_var),
    mean_mode = as.numeric(mode_mean),
    sd_mode = as.numeric(mode_sd),
    shrinkage_ratio = as.numeric(shrinkage_ratio),
    shapiro_p = as.numeric(shapiro_p),
    row.names = NULL
  )

  # ICC of logit(Y) implied by the fitted beta model (not pi^2/3, a logistic value)
  icc <- .brsmm_icc(object, D)

  out <- list(
    summary = summary_df,
    D = D,
    Corr = Corr,
    icc = icc,
    n_groups = nrow(B),
    modes = B
  )
  class(out) <- "brsmm_re_study"
  out
}

#' Print a random-effects study
#'
#' @description
#' Prints a compact summary of the random-effects study returned by
#' \code{\link{brsmm_re_study}}, including per-term standard deviations,
#' shrinkage ratios, Shapiro-Wilk p-values, and the estimated covariance
#' and correlation matrices.
#'
#' @param x A \code{"brsmm_re_study"} object returned by
#'   \code{\link{brsmm_re_study}}.
#' @param digits Integer: number of significant digits for rounding
#'   (default \code{max(3, getOption("digits") - 3)}).
#' @param ... Currently ignored.
#'
#' @return Invisibly returns \code{x}. Called for its side-effect of
#'   printing the study to the console.
#'
#' @method print brsmm_re_study
#'
#' @seealso \code{\link{brsmm_re_study}}, \code{\link{brsmm}},
#'   \code{\link{ranef.brsmm}}
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
#' fit <- brsmm(y ~ x1, random = ~ 1 | id, data = prep)
#' rs <- brsmm_re_study(fit)
#' print(rs)
#' }
#'
#' @export
print.brsmm_re_study <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  cat("\nRandom-effects study\n")
  cat("Groups:", x$n_groups, "\n\n")

  # --- VarCorr block (lme4 style) ---
  cat("Random-effects (VarCorr):\n")
  D <- x$D
  Corr <- x$Corr
  nms <- colnames(D)
  q <- nrow(D)
  sd_vals <- sqrt(pmax(diag(D), 0))

  # Header
  corr_header <- if (q > 1L) paste0("  Corr") else ""
  cat(sprintf("  %-22s  %10s%s\n", "Name", "Std.Dev.", corr_header))
  for (i in seq_len(q)) {
    nm <- if (is.null(nms)) paste0("re", i) else nms[i]
    sdv <- formatC(sd_vals[i], format = "f", digits = digits)
    corr_part <- ""
    if (q > 1L && i > 1L) {
      cors <- formatC(Corr[i, seq_len(i - 1L)], format = "f", digits = digits)
      corr_part <- paste0("  ", paste(cors, collapse = "  "))
    }
    cat(sprintf("  %-22s  %10s%s\n", nm, sdv, corr_part))
  }

  # --- ICC ---
  cat(sprintf("\nICC (logit(Y) scale, beta level-1 variance): %.4f\n", x$icc))

  # --- Per-term summary ---
  cat("\nSummary by term (SD_model = model SD; shrinkage = Var(modes)/Var(model)):\n")
  sm <- x$summary
  is_num <- vapply(sm, is.numeric, logical(1L))
  sm[is_num] <- lapply(sm[is_num], round, digits = digits)
  print(sm, row.names = FALSE)

  invisible(x)
}


# ICC of logit(Y) implied by a brsmm fit (see ?brsmm_re_study): averaged over rows.
.brsmm_icc <- function(object, D) {
  mm <- object$model_matrices
  eta0 <- as.numeric(mm$X %*% object$coefficients$mean)
  phi <- .clamp_phi_by_repar(
    apply_inv_link(as.numeric(mm$Z %*% object$coefficients$precision), object$link_phi),
    object$repar
  )
  # Variance of the random part x_r' b of each row
  s2 <- rowSums((mm$Xr %*% D) * mm$Xr)
  .brs_icc_logit(eta0, phi, s2, object$link, object$repar)
}

# ICC_i = Var_b[m] / (Var_b[m] + E_b[v]), m = psi(a) - psi(b), v = psi1(a) + psi1(b),
# b ~ N(0, s2_i) by Gauss-Hermite; mean over i. NA (warning) when the clamp drives it.
.brs_icc_logit <- function(eta0, phi, s2, link, repar, n_gh = 40L) {
  at_pkg <- .brs_icc_parts(eta0, phi, s2, link, repar, n_gh, eps = 1e-5)
  at_tight <- .brs_icc_parts(eta0, phi, s2, link, repar, n_gh, eps = 1e-8)
  # Infinite logit(Y) moments (probit, cloglog tails): E_b[v] then follows the clamp
  rel <- abs(at_tight$ev - at_pkg$ev) / at_pkg$ev
  if (!is.finite(rel) || rel > 0.1) {
    warning("ICC not available: the variance of logit(Y) is driven by the clamp of ",
            "the mean (link '", link, "', E_b[psi1(a) + psi1(b)] changes by ",
            format(signif(100 * rel, 3)), "% between the 1e-5 and 1e-8 clamps).",
            call. = FALSE)
    return(NA_real_)
  }
  at_pkg$icc
}

# Mean ICC_i and mean E_b[v_i] with the mean clamped to [eps, 1 - eps]
# ([eps, 1e8] for the shape p under repar 0); beta shapes as in the likelihood.
.brs_icc_parts <- function(eta0, phi, s2, link, repar, n_gh, eps) {
  gh <- .brs_gh_normal(n_gh)
  n <- length(eta0)
  eta <- eta0 + outer(sqrt(pmax(rep_len(s2, n), 0)), gh$x)
  mu <- as.numeric(apply_inv_link(eta, link))
  top <- if (as.integer(repar) == 0L) 1e8 else 1 - eps
  mu[!is.finite(mu)] <- top
  mu <- pmin(pmax(mu, eps), top)
  # phi does not depend on b: one column per node, repeated
  sh <- .brs_shapes_unclamped(mu, rep(rep_len(phi, n), times = n_gh), repar)
  a <- matrix(pmin(pmax(sh$a, 1e-12), 1e8), nrow = n)
  b <- matrix(pmin(pmax(sh$b, 1e-12), 1e8), nrow = n)
  m <- digamma(a) - digamma(b)
  v <- trigamma(a) + trigamma(b)
  em <- drop(m %*% gh$w)
  vm <- pmax(drop((m - em)^2 %*% gh$w), 0)
  ev <- drop(v %*% gh$w)
  list(icc = mean(vm / (vm + ev)), ev = mean(ev))
}

# Gauss-Hermite nodes/weights for E[f(Z)], Z ~ N(0, 1) (Golub-Welsch).
.brs_gh_normal <- function(n) {
  J <- matrix(0, n, n)
  off <- sqrt(seq_len(n - 1L))
  J[cbind(seq_len(n - 1L), 2:n)] <- off
  J[cbind(2:n, seq_len(n - 1L))] <- off
  e <- eigen(J, symmetric = TRUE)
  list(x = e$values, w = e$vectors[1L, ]^2)
}
