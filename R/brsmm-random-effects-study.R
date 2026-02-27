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

  # ICC on latent logistic scale: ICC = sigma_b^2 / (sigma_b^2 + pi^2/3)
  # Uses only the intercept RE variance (column 1) as the between-group component.
  sigma_b2_intercept <- as.numeric(model_var[1L])
  icc_latent <- sigma_b2_intercept / (sigma_b2_intercept + pi^2 / 3)

  out <- list(
    summary = summary_df,
    D = D,
    Corr = Corr,
    icc = icc_latent,
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
  cat(sprintf("\nICC (latent logistic scale): %.4f\n", x$icc))

  # --- Per-term summary ---
  cat("\nSummary by term (SD_model = model SD; shrinkage = Var(modes)/Var(model)):\n")
  sm <- x$summary
  is_num <- vapply(sm, is.numeric, logical(1L))
  sm[is_num] <- lapply(sm[is_num], round, digits = digits)
  print(sm, row.names = FALSE)

  invisible(x)
}
