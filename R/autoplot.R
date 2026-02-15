# ============================================================================ #
# ggplot2 autoplot methods
# ============================================================================ #

#' ggplot2 autoplot for brs models
#'
#' @description
#' Produces ggplot2 diagnostics tailored to interval-censored scale models.
#'
#' @param object A fitted \code{"brs"} object.
#' @param type Plot type:
#'   \code{"calibration"}, \code{"score_dist"}, \code{"cdf"}, or
#'   \code{"residuals_by_delta"}.
#' @param bins Number of bins used in calibration plots.
#' @param scores Optional integer vector of scores for \code{"score_dist"}.
#'   Defaults to all scores from \code{0} to \code{ncuts}.
#' @param newdata Optional data frame of covariate scenarios used by
#'   \code{type = "cdf"}.
#' @param n_grid Number of points on \eqn{(0,1)} used to draw CDF curves.
#' @param max_curves Maximum number of CDF curves shown when \code{newdata}
#'   is not provided.
#' @param residual_type Residual type passed to \code{\link{residuals.brs}}
#'   for \code{type = "residuals_by_delta"}.
#' @param ... Currently ignored.
#'
#' @return A \code{ggplot2} object.
#'
#' @details
#' \code{type = "calibration"} bins predictions and compares mean observed vs
#' mean predicted response in each bin.
#'
#' \code{type = "score_dist"} compares observed score frequencies against
#' expected frequencies implied by the fitted beta interval model.
#'
#' @references
#' Hawker, G. A., Mian, S., Kendzerska, T., and French, M. (2011).
#' Measures of adult pain: Visual Analog Scale for Pain (VAS Pain),
#' Numeric Rating Scale for Pain (NRS Pain), McGill Pain Questionnaire (MPQ),
#' Short-Form McGill Pain Questionnaire (SF-MPQ), Chronic Pain Grade Scale
#' (CPGS), Short Form-36 Bodily Pain Scale (SF-36 BPS), and Measure of
#' Intermittent and Constant Osteoarthritis Pain (ICOAP).
#' Arthritis Care and Research, 63(S11), S240-S252.
#' doi:10.1002/acr.20543.
#'
#' Hjermstad, M. J., Fayers, P. M., Haugen, D. F., et al. (2011).
#' Studies comparing Numerical Rating Scales, Verbal Rating Scales, and
#' Visual Analogue Scales for assessment of pain intensity in adults:
#' a systematic literature review.
#' Journal of Pain and Symptom Management, 41(6), 1073-1093.
#' doi:10.1016/j.jpainsymman.2010.08.016.
#'
#' @examples
#' \donttest{
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   set.seed(100)
#'   dat <- data.frame(x1 = rnorm(120), x2 = rnorm(120))
#'   sim <- brs_sim(
#'     formula = ~ x1 + x2, data = dat,
#'     beta = c(0.1, -0.3, 0.2), phi = 0.2, ncuts = 100, repar = 2
#'   )
#'   fit <- brs(y ~ x1 + x2, data = sim, repar = 2)
#'
#'   autoplot.brs(fit, type = "calibration")
#'   autoplot.brs(fit, type = "score_dist")
#' }
#' }
#'
#' @export
autoplot.brs <- function(object,
                         type = c(
                           "calibration",
                           "score_dist",
                           "cdf",
                           "residuals_by_delta"
                         ),
                         bins = 10L,
                         scores = NULL,
                         newdata = NULL,
                         n_grid = 200L,
                         max_curves = 6L,
                         residual_type = "rqr",
                         ...) {
  .check_class(object)
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for autoplot().", call. = FALSE)
  }

  type <- match.arg(type)
  bins <- as.integer(bins)
  if (!is.finite(bins) || bins < 3L) {
    stop("'bins' must be an integer >= 3.", call. = FALSE)
  }

  switch(type,
    calibration = .brs_autoplot_calibration(object, bins = bins),
    score_dist = .brs_autoplot_score_dist(object, scores = scores),
    cdf = .brs_autoplot_cdf(
      object,
      newdata = newdata,
      n_grid = as.integer(n_grid),
      max_curves = as.integer(max_curves)
    ),
    residuals_by_delta = .brs_autoplot_resid_delta(
      object,
      residual_type = residual_type
    )
  )
}

#' @keywords internal
.brs_autoplot_calibration <- function(object, bins = 10L) {
  df <- data.frame(
    observed = as.numeric(object$Y[, "yt"]),
    predicted = as.numeric(object$hatmu)
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
      title = "Calibration Plot",
      x = "Mean predicted response (bin average)",
      y = "Mean observed response (bin average)",
      size = "Bin n"
    ) +
    ggplot2::theme_minimal()
}

#' @keywords internal
.brs_autoplot_score_dist <- function(object, scores = NULL) {
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
    mu = object$hatmu,
    phi = object$hatphi,
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
      title = "Observed vs Expected Score Distribution",
      x = "Score",
      y = "Count",
      fill = ""
    ) +
    ggplot2::theme_minimal()
}

#' @keywords internal
.brs_observed_scores <- function(y, K) {
  y <- as.numeric(y)
  if (all(y >= 0 & y <= 1, na.rm = TRUE)) {
    out <- round(y * K)
  } else {
    out <- round(y)
  }
  pmin(pmax(out, 0L), K)
}

#' @keywords internal
.brs_score_prob_matrix <- function(mu, phi, repar, ncuts, lim, scores) {
  eps <- 1e-10
  shp <- brs_repar(mu = mu, phi = phi, repar = repar)

  P <- sapply(scores, function(s) {
    if (s == 0L) {
      u <- lim / ncuts
      return(stats::pbeta(u, shp$shape1, shp$shape2))
    }
    if (s == ncuts) {
      l <- (ncuts - lim) / ncuts
      return(1 - stats::pbeta(l, shp$shape1, shp$shape2))
    }
    l <- (s - lim) / ncuts
    u <- (s + lim) / ncuts
    pmax(stats::pbeta(u, shp$shape1, shp$shape2) - stats::pbeta(l, shp$shape1, shp$shape2), eps)
  })

  if (is.vector(P)) {
    P <- matrix(P, ncol = length(scores))
  }
  colnames(P) <- paste0("score_", scores)
  P
}

#' @keywords internal
.brs_autoplot_cdf <- function(object,
                              newdata = NULL,
                              n_grid = 200L,
                              max_curves = 6L) {
  n_grid <- as.integer(n_grid)
  max_curves <- as.integer(max_curves)
  if (!is.finite(n_grid) || n_grid < 20L) {
    stop("'n_grid' must be an integer >= 20.", call. = FALSE)
  }
  if (!is.finite(max_curves) || max_curves < 1L) {
    stop("'max_curves' must be an integer >= 1.", call. = FALSE)
  }

  grid <- seq(1e-4, 1 - 1e-4, length.out = n_grid)

  if (is.null(newdata)) {
    ord <- order(object$hatmu)
    idx <- unique(round(seq(1, length(ord), length.out = min(max_curves, length(ord)))))
    mu <- object$hatmu[ord[idx]]
    phi <- object$hatphi[ord[idx]]
    labels <- paste0("scenario_", seq_along(mu))
  } else {
    if (!is.data.frame(newdata)) {
      stop("'newdata' must be a data.frame.", call. = FALSE)
    }
    if (nrow(newdata) > max_curves) {
      newdata <- newdata[seq_len(max_curves), , drop = FALSE]
    }
    mu <- predict(object, newdata = newdata, type = "response")
    phi <- predict(object, newdata = newdata, type = "precision")
    labels <- paste0("scenario_", seq_along(mu))
  }

  shp <- brs_repar(mu = mu, phi = phi, repar = object$repar)
  curves <- lapply(seq_along(mu), function(i) {
    data.frame(
      y = grid,
      cdf = stats::pbeta(grid, shp$shape1[i], shp$shape2[i]),
      scenario = labels[i],
      stringsAsFactors = FALSE
    )
  })
  df <- do.call(rbind, curves)

  ggplot2::ggplot(df, ggplot2::aes(x = .data$y, y = .data$cdf, color = .data$scenario)) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::labs(
      title = "Predicted Beta CDF by Scenario",
      x = "y (unit scale)",
      y = "F(y)",
      color = ""
    ) +
    ggplot2::theme_minimal()
}

#' @keywords internal
.brs_autoplot_resid_delta <- function(object, residual_type = "rqr") {
  dlab <- c("Exact", "Left", "Right", "Interval")
  df <- data.frame(
    residual = residuals(object, type = residual_type),
    delta = factor(object$delta, levels = 0:3, labels = dlab),
    stringsAsFactors = FALSE
  )

  ggplot2::ggplot(df, ggplot2::aes(x = .data$delta, y = .data$residual, fill = .data$delta)) +
    ggplot2::geom_boxplot(alpha = 0.65, outlier.shape = NA) +
    ggplot2::geom_jitter(width = 0.18, alpha = 0.25, size = 1) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray35") +
    ggplot2::labs(
      title = paste0("Residuals by Censoring Type (", residual_type, ")"),
      x = "Censoring type",
      y = "Residual"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "none")
}
