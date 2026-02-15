# ============================================================================ #
# ggplot2 autoplot methods for brsmm
# ============================================================================ #

#' ggplot2 autoplot for brsmm models
#'
#' @description
#' Produces ggplot2 diagnostics tailored to mixed beta interval models.
#'
#' @param object A fitted \code{"brsmm"} object.
#' @param type Plot type:
#'   \code{"calibration"}, \code{"score_dist"}, \code{"ranef_qq"}, or
#'   \code{"residuals_by_group"}.
#' @param bins Number of bins used in calibration plots.
#' @param scores Optional integer vector of scores for \code{"score_dist"}.
#'   Defaults to all scores from \code{0} to \code{ncuts}.
#' @param residual_type Residual type passed to \code{\link{residuals.brsmm}}
#'   for \code{type = "residuals_by_group"}.
#' @param max_groups Maximum number of groups displayed in
#'   \code{"residuals_by_group"}.
#' @param ... Currently ignored.
#'
#' @return A \code{ggplot2} object.
#'
#' @examples
#' \donttest{
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   set.seed(123)
#'   g <- 10
#'   ni <- 8
#'   id <- factor(rep(seq_len(g), each = ni))
#'   n <- length(id)
#'   x1 <- rnorm(n)
#'   b <- rnorm(g, sd = 0.4)
#'
#'   mu <- plogis(0.1 + 0.5 * x1 + b[as.integer(id)])
#'   phi <- plogis(-0.2)
#'   shp <- brs_repar(mu = mu, phi = rep(phi, n), repar = 2)
#'   y <- round(stats::rbeta(n, shp$shape1, shp$shape2) * 100)
#'   d <- data.frame(y = y, x1 = x1, id = id)
#'
#'   fit_mm <- brsmm(y ~ x1, random = ~ 1 | id, data = d, repar = 2)
#'
#'   autoplot.brsmm(fit_mm, type = "calibration")
#'   autoplot.brsmm(fit_mm, type = "ranef_qq")
#' }
#' }
#'
#' @export
autoplot.brsmm <- function(object,
                           type = c(
                             "calibration",
                             "score_dist",
                             "ranef_qq",
                             "residuals_by_group"
                           ),
                           bins = 10L,
                           scores = NULL,
                           residual_type = c("response", "pearson"),
                           max_groups = 25L,
                           ...) {
  .check_class_mm(object)
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for autoplot.brsmm().", call. = FALSE)
  }

  type <- match.arg(type)
  residual_type <- match.arg(residual_type)
  bins <- as.integer(bins)
  max_groups <- as.integer(max_groups)

  if (!is.finite(bins) || bins < 3L) {
    stop("'bins' must be an integer >= 3.", call. = FALSE)
  }
  if (!is.finite(max_groups) || max_groups < 2L) {
    stop("'max_groups' must be an integer >= 2.", call. = FALSE)
  }

  switch(type,
    calibration = .brsmm_autoplot_calibration(object, bins = bins),
    score_dist = .brsmm_autoplot_score_dist(object, scores = scores),
    ranef_qq = .brsmm_autoplot_ranef_qq(object),
    residuals_by_group = .brsmm_autoplot_resid_group(
      object,
      residual_type = residual_type,
      max_groups = max_groups
    )
  )
}

#' @keywords internal
.brsmm_autoplot_calibration <- function(object, bins = 10L) {
  df <- data.frame(
    observed = as.numeric(object$Y[, "yt"]),
    predicted = as.numeric(object$fitted_mu)
  )
  probs <- seq(0, 1, length.out = bins + 1L)
  breaks <- unique(stats::quantile(df$predicted, probs = probs, na.rm = TRUE))
  if (length(breaks) < 3L) {
    breaks <- seq(min(df$predicted), max(df$predicted), length.out = bins + 1L)
  }
  df$bin <- cut(df$predicted, breaks = breaks, include.lowest = TRUE, ordered_result = TRUE)

  cal <- stats::aggregate(df[, c("predicted", "observed")], by = list(bin = df$bin), FUN = mean)
  cal$n <- as.integer(table(df$bin)[as.character(cal$bin)])

  ggplot2::ggplot(cal, ggplot2::aes(x = .data$predicted, y = .data$observed, size = .data$n)) +
    ggplot2::geom_point(color = "#1b9e77", alpha = 0.9) +
    ggplot2::geom_line(color = "#1b9e77", alpha = 0.6) +
    ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray35") +
    ggplot2::labs(
      title = "Calibration Plot (brsmm)",
      x = "Mean predicted response (bin average)",
      y = "Mean observed response (bin average)",
      size = "Bin n"
    ) +
    ggplot2::theme_minimal()
}

#' @keywords internal
.brsmm_autoplot_score_dist <- function(object, scores = NULL) {
  K <- as.integer(object$ncuts)
  if (is.null(scores)) {
    scores <- 0:K
  }
  scores <- sort(unique(as.integer(scores)))
  if (any(!is.finite(scores)) || any(scores < 0L) || any(scores > K)) {
    stop("'scores' must be integers in [0, ncuts].", call. = FALSE)
  }

  obs_scores <- .brs_observed_scores(object$Y[, "y"], K = K)
  obs_counts <- as.numeric(table(factor(obs_scores, levels = scores)))

  probs <- .brs_score_prob_matrix(
    mu = object$fitted_mu,
    phi = object$fitted_phi,
    repar = object$repar,
    ncuts = K,
    lim = object$lim,
    scores = scores
  )
  exp_counts <- colSums(probs)

  df <- rbind(
    data.frame(score = scores, count = obs_counts, source = "Observed"),
    data.frame(score = scores, count = exp_counts, source = "Expected")
  )

  ggplot2::ggplot(df, ggplot2::aes(x = .data$score, y = .data$count, fill = .data$source)) +
    ggplot2::geom_col(position = "dodge", alpha = 0.85) +
    ggplot2::scale_fill_manual(values = c(Observed = "#1b9e77", Expected = "#7570b3")) +
    ggplot2::labs(
      title = "Observed vs Expected Score Distribution (brsmm)",
      x = "Score",
      y = "Count",
      fill = ""
    ) +
    ggplot2::theme_minimal()
}

#' @keywords internal
.brsmm_autoplot_ranef_qq <- function(object) {
  re <- sort(as.numeric(object$random$mode_b))
  n <- length(re)
  q <- stats::qnorm(stats::ppoints(n))
  df <- data.frame(theoretical = q, sample = re)

  ggplot2::ggplot(df, ggplot2::aes(x = .data$theoretical, y = .data$sample)) +
    ggplot2::geom_point(color = "#1b9e77", alpha = 0.9) +
    ggplot2::geom_abline(
      intercept = mean(re) - stats::sd(re) * mean(q),
      slope = stats::sd(re),
      linetype = "dashed",
      color = "gray35"
    ) +
    ggplot2::labs(
      title = "Random-Effect Q-Q Plot",
      x = "Theoretical Normal Quantiles",
      y = "Estimated random intercept mode"
    ) +
    ggplot2::theme_minimal()
}

#' @keywords internal
.brsmm_autoplot_resid_group <- function(object,
                                        residual_type = c("response", "pearson"),
                                        max_groups = 25L) {
  residual_type <- match.arg(residual_type)
  res <- residuals(object, type = residual_type)
  grp <- object$group
  tab <- sort(table(grp), decreasing = TRUE)
  keep <- names(tab)[seq_len(min(max_groups, length(tab)))]

  df <- data.frame(
    residual = res,
    group = factor(as.character(grp), levels = keep)
  )
  df <- df[!is.na(df$group), , drop = FALSE]

  ggplot2::ggplot(df, ggplot2::aes(x = .data$group, y = .data$residual, fill = .data$group)) +
    ggplot2::geom_boxplot(alpha = 0.70, outlier.shape = NA) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray35") +
    ggplot2::labs(
      title = paste0(
        "Residuals by Group (top ", length(keep), " groups, ",
        residual_type, ")"
      ),
      x = "Group",
      y = "Residual"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position = "none",
      axis.text.x = ggplot2::element_text(angle = 60, hjust = 1)
    )
}
