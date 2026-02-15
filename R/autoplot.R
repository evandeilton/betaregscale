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
#'   \code{"calibration"} or \code{"score_dist"}.
#' @param bins Number of bins used in calibration plots.
#' @param scores Optional integer vector of scores for \code{"score_dist"}.
#'   Defaults to all scores from \code{0} to \code{ncuts}.
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
#'   autoplot(fit, type = "calibration")
#'   autoplot(fit, type = "score_dist")
#' }
#' }
#'
#' @export
autoplot.brs <- function(object,
                         type = c("calibration", "score_dist"),
                         bins = 10L,
                         scores = NULL,
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

  if (identical(type, "calibration")) {
    return(.brs_autoplot_calibration(object, bins = bins))
  }
  .brs_autoplot_score_dist(object, scores = scores)
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
