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
#' @param object A fitted \code{"brsmm"} object.
#' @param ... Currently ignored.
#'
#' @return A list with class \code{"brsmm_re_study"}.
#'
#' @examples
#' \donttest{
#' set.seed(123)
#' g <- 10
#' ni <- 8
#' id <- factor(rep(seq_len(g), each = ni))
#' n <- length(id)
#' x1 <- rnorm(n)
#' b0 <- rnorm(g, sd = 0.4)
#' b1 <- rnorm(g, sd = 0.2)
#' mu <- plogis(0.2 + 0.5 * x1 + b0[id] + b1[id] * x1)
#' phi <- plogis(-0.2)
#' shp <- brs_repar(mu = mu, phi = rep(phi, n), repar = 2)
#' y <- round(stats::rbeta(n, shp$shape1, shp$shape2) * 100)
#' d <- data.frame(y = y, x1 = x1, id = id)
#' fit <- brsmm(y ~ x1, random = ~ 1 + x1 | id, data = d, repar = 2)
#'
#' rs <- brsmm_re_study(fit)
#' print(rs)
#' knitr::kable(rs$summary, digits = 4)
#' }
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
    colnames(D) <- rownames(D) <- colnames(B)
  }
  Corr <- stats::cov2cor(D)

  mode_mean <- colMeans(B)
  mode_sd <- apply(B, 2, stats::sd)
  mode_var <- pmax(mode_sd^2, 0)
  model_var <- pmax(diag(D), 1e-12)
  shrinkage_ratio <- pmin(pmax(mode_var / model_var, 0), 1)

  shapiro_p <- rep(NA_real_, ncol(B))
  if (nrow(B) >= 3L && nrow(B) <= 5000L) {
    shapiro_p <- apply(B, 2, function(x) {
      stats::shapiro.test(as.numeric(x))$p.value
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

  out <- list(
    summary = summary_df,
    D = D,
    Corr = Corr,
    n_groups = nrow(B),
    modes = B
  )
  class(out) <- "brsmm_re_study"
  out
}

#' Print method for random-effects study
#'
#' @param x A \code{"brsmm_re_study"} object.
#' @param digits Number of digits.
#' @param ... Currently ignored.
#' @return Invisibly returns \code{x}.
#' @export
print.brsmm_re_study <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  cat("\nRandom-effects study\n")
  cat("Groups:", x$n_groups, "\n\n")
  cat("Summary by term:\n")
  sm <- x$summary
  is_num <- vapply(sm, is.numeric, logical(1))
  sm[is_num] <- lapply(sm[is_num], round, digits = digits)
  print(sm, row.names = FALSE)
  cat("\nEstimated covariance matrix D:\n")
  print(round(x$D, digits))
  cat("\nEstimated correlation matrix:\n")
  print(round(x$Corr, digits))
  invisible(x)
}
