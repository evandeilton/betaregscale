# ============================================================================ #
# ggplot2 autoplot methods
# ============================================================================ #

#' Register S3 methods for ggplot2::autoplot
#'
#' @rawNamespace S3method(autoplot, brs)
#' @rawNamespace S3method(autoplot, brs_bootstrap)
#' @rawNamespace S3method(autoplot, brs_marginaleffects)
#' @rawNamespace S3method(autoplot, brsmm)
#' @keywords internal
"_PACKAGE"

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
#' @importFrom ggplot2 autoplot
#' @method autoplot brs
#' @export autoplot.brs
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


#' ggplot2 autoplot for bootstrap results
#'
#' @description
#' Produces visual summaries for objects returned by \code{\link{brs_bootstrap}}.
#'
#' @param object An object of class \code{"brs_bootstrap"}.
#' @param type Plot type:
#'   \code{"ci_forest"}, \code{"dist"}, \code{"qq"}, or \code{"stability"}.
#' @param parameter Optional parameter name used by \code{type = "dist"},
#'   \code{"qq"}, and \code{"stability"}. If \code{NULL}, the first parameter
#'   is used.
#' @param title Optional plot title override.
#' @param caption Optional subtitles/titles for plot types. Accepts:
#'   \itemize{
#'     \item a single string (used for the selected \code{type});
#'     \item a character vector/list with up to four entries in the order
#'       \code{ci_forest}, \code{dist}, \code{qq}, \code{stability}.
#'   }
#' @param max_parameters Maximum number of parameters shown in
#'   \code{type = "ci_forest"}.
#' @param ci_level Confidence level used in \code{type = "stability"}.
#'   Defaults to the level stored in \code{object}.
#' @param theme Optional ggplot2 theme object (e.g., \code{ggplot2::theme_bw()}).
#'   If \code{NULL}, \code{ggplot2::theme_minimal()} is used.
#' @param ... Currently ignored.
#'
#' @return A \code{ggplot2} object.
#'
#' @details
#' For \code{type = "dist"}, \code{"qq"}, and \code{"stability"},
#' bootstrap draws must be present in \code{attr(object, "boot_draws")},
#' obtained by fitting with \code{brs_bootstrap(..., keep_draws = TRUE)}.
#'
#' @method autoplot brs_bootstrap
#' @export autoplot.brs_bootstrap
autoplot.brs_bootstrap <- function(object,
                                   type = c("ci_forest", "dist", "qq", "stability"),
                                   parameter = NULL,
                                   title = NULL,
                                   caption = NULL,
                                   max_parameters = 12L,
                                   ci_level = NULL,
                                   theme = NULL,
                                   ...) {
  if (!inherits(object, "brs_bootstrap")) {
    stop("'object' must be of class 'brs_bootstrap'.", call. = FALSE)
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for autoplot().", call. = FALSE)
  }

  type <- match.arg(type)
  type_levels <- c("ci_forest", "dist", "qq", "stability")
  max_parameters <- as.integer(max_parameters)
  if (!is.finite(max_parameters) || max_parameters < 1L) {
    stop("'max_parameters' must be an integer >= 1.", call. = FALSE)
  }
  theme_obj <- .boot_resolve_theme(theme)
  cap <- .boot_pick_caption(caption = caption, type = type, type_levels = type_levels)

  draws <- attr(object, "boot_draws")

  switch(type,
    ci_forest = .boot_autoplot_ci_forest(
      object,
      max_parameters = max_parameters,
      title = title,
      caption = cap,
      theme_obj = theme_obj
    ),
    dist = .boot_autoplot_dist(
      object,
      draws = draws,
      parameter = parameter,
      title = title,
      caption = cap,
      theme_obj = theme_obj
    ),
    qq = .boot_autoplot_qq(
      object,
      draws = draws,
      parameter = parameter,
      title = title,
      caption = cap,
      theme_obj = theme_obj
    ),
    stability = .boot_autoplot_stability(
      object,
      draws = draws,
      parameter = parameter,
      ci_level = ci_level,
      title = title,
      caption = cap,
      theme_obj = theme_obj
    )
  )
}

#' @keywords internal
.boot_need_draws <- function(draws) {
  if (is.null(draws) || !is.matrix(draws) || nrow(draws) < 5L) {
    stop(
      "This plot requires bootstrap draws. Fit with ",
      "'brs_bootstrap(..., keep_draws = TRUE)'.",
      call. = FALSE
    )
  }
}

#' @keywords internal
.boot_pick_parameter <- function(object, draws = NULL, parameter = NULL) {
  params <- as.character(object$parameter)
  if (!is.null(draws) && !is.null(colnames(draws))) {
    params <- colnames(draws)
  }
  if (is.null(parameter)) {
    return(params[1L])
  }
  parameter <- as.character(parameter)[1L]
  if (!(parameter %in% params)) {
    stop(
      "'parameter' must be one of: ",
      paste(params, collapse = ", "),
      call. = FALSE
    )
  }
  parameter
}

#' @keywords internal
.boot_resolve_theme <- function(theme) {
  if (is.null(theme)) {
    return(ggplot2::theme_minimal())
  }
  if (is.function(theme)) {
    return(theme())
  }
  theme
}

#' @keywords internal
.boot_pick_caption <- function(caption, type, type_levels) {
  defaults <- c(
    ci_forest = "Parameters ordered by effect size",
    dist = "Bootstrap distribution for selected parameter",
    qq = "Normality diagnostic for bootstrap draws",
    stability = "CI bounds evolution across replicates"
  )
  if (is.null(caption)) {
    return(unname(defaults[[type]]))
  }
  cap <- as.character(unlist(caption, use.names = FALSE))
  cap <- cap[!is.na(cap)]
  if (length(cap) == 0L) {
    return(unname(defaults[[type]]))
  }
  if (length(cap) == 1L) {
    return(cap[1L])
  }
  if (length(cap) < length(type_levels)) {
    cap <- c(cap, unname(defaults[type_levels[(length(cap) + 1L):length(type_levels)]]))
  }
  cap <- cap[seq_along(type_levels)]
  cap[match(type, type_levels)]
}

#' @keywords internal
.boot_autoplot_ci_forest <- function(object, max_parameters = 12L, title = NULL, caption = NULL, theme_obj = ggplot2::theme_minimal()) {
  df <- as.data.frame(object)
  if (nrow(df) > max_parameters) {
    ord <- order(abs(df$estimate), decreasing = TRUE)
    df <- df[ord[seq_len(max_parameters)], , drop = FALSE]
  }
  df$parameter <- factor(df$parameter, levels = rev(df$parameter))

  main_title <- if (is.null(title)) {
    paste0(
      "Bootstrap CI forest (",
      round(100 * unique(df$level), 1), "%, ",
      attr(object, "ci_type"),
      ")"
    )
  } else {
    title
  }

  p <- ggplot2::ggplot(df, ggplot2::aes(y = .data$parameter)) +
    ggplot2::geom_segment(
      ggplot2::aes(x = .data$ci_lower, xend = .data$ci_upper, yend = .data$parameter),
      linewidth = 1.0, color = "#2C7FB8", alpha = 0.75
    )

  has_wald <- all(c("wald_lower", "wald_upper") %in% names(df))
  if (has_wald) {
    p <- p + ggplot2::geom_segment(
      ggplot2::aes(
        x = .data$wald_lower, xend = .data$wald_upper, yend = .data$parameter
      ),
      linewidth = 0.55, color = "gray45", alpha = 0.9
    )
  }

  if (is.null(caption) && has_wald) {
    caption <- "Blue: bootstrap CI | Gray: Wald CI"
  }

  p +
    ggplot2::geom_point(ggplot2::aes(x = .data$estimate), size = 2.2, color = "#D95F0E") +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    ggplot2::labs(
      title = main_title,
      subtitle = caption,
      x = "Parameter value",
      y = "Parameter"
    ) +
    theme_obj
}

#' @keywords internal
.boot_autoplot_dist <- function(object, draws, parameter = NULL, title = NULL, caption = NULL, theme_obj = ggplot2::theme_minimal()) {
  .boot_need_draws(draws)
  par_name <- .boot_pick_parameter(object, draws = draws, parameter = parameter)
  j <- match(par_name, colnames(draws))
  x <- draws[, j]
  x <- x[is.finite(x)]

  est <- object$estimate[match(par_name, object$parameter)]
  lo <- object$ci_lower[match(par_name, object$parameter)]
  hi <- object$ci_upper[match(par_name, object$parameter)]

  df <- data.frame(value = x)
  main_title <- if (is.null(title)) paste0("Bootstrap distribution: ", par_name) else title
  sub_cap <- if (is.null(caption)) {
    paste0(
      "estimate = ", format(round(est, 4), nsmall = 4),
      " | CI = [", format(round(lo, 4), nsmall = 4), ", ",
      format(round(hi, 4), nsmall = 4), "]"
    )
  } else {
    caption
  }

  ggplot2::ggplot(df, ggplot2::aes(x = .data$value)) +
    ggplot2::geom_histogram(
      bins = 30, fill = "#74A9CF", color = "white", alpha = 0.85
    ) +
    ggplot2::geom_vline(xintercept = est, color = "#D95F0E", linewidth = 0.9) +
    ggplot2::geom_vline(xintercept = c(lo, hi), color = "#2B8CBE", linetype = "dashed") +
    ggplot2::labs(
      title = main_title,
      subtitle = sub_cap,
      x = "Bootstrap draw",
      y = "Count"
    ) +
    theme_obj
}

#' @keywords internal
.boot_autoplot_qq <- function(object, draws, parameter = NULL, title = NULL, caption = NULL, theme_obj = ggplot2::theme_minimal()) {
  .boot_need_draws(draws)
  par_name <- .boot_pick_parameter(object, draws = draws, parameter = parameter)
  j <- match(par_name, colnames(draws))
  x <- draws[, j]
  x <- x[is.finite(x)]

  df <- data.frame(value = x)
  main_title <- if (is.null(title)) paste0("Bootstrap QQ plot: ", par_name) else title

  ggplot2::ggplot(df, ggplot2::aes(sample = .data$value)) +
    ggplot2::stat_qq(color = "#2C7FB8", alpha = 0.75, size = 1.2) +
    ggplot2::stat_qq_line(color = "#D95F0E", linewidth = 0.8) +
    ggplot2::labs(
      title = main_title,
      subtitle = caption,
      x = "Theoretical normal quantiles",
      y = "Sample quantiles"
    ) +
    theme_obj
}

#' @keywords internal
.boot_autoplot_stability <- function(object, draws, parameter = NULL, ci_level = NULL, title = NULL, caption = NULL, theme_obj = ggplot2::theme_minimal()) {
  .boot_need_draws(draws)
  par_name <- .boot_pick_parameter(object, draws = draws, parameter = parameter)
  j <- match(par_name, colnames(draws))
  x <- draws[, j]
  x <- x[is.finite(x)]

  if (is.null(ci_level)) {
    ci_level <- unique(object$level)[1L]
  }
  ci_level <- as.numeric(ci_level)[1L]
  if (!is.finite(ci_level) || ci_level <= 0 || ci_level >= 1) {
    stop("'ci_level' must be in (0, 1).", call. = FALSE)
  }
  alpha <- 1 - ci_level
  probs <- c(alpha / 2, 1 - alpha / 2)

  n <- length(x)
  idx <- seq(10L, n, by = 1L)
  lo <- vapply(idx, function(k) stats::quantile(x[seq_len(k)], probs = probs[1L], names = FALSE), numeric(1))
  hi <- vapply(idx, function(k) stats::quantile(x[seq_len(k)], probs = probs[2L], names = FALSE), numeric(1))
  df <- data.frame(
    n = rep(idx, 2L),
    bound = c(lo, hi),
    side = rep(c("Lower", "Upper"), each = length(idx))
  )

  main_title <- if (is.null(title)) {
    paste0("CI stability by number of replicates: ", par_name)
  } else {
    title
  }
  sub_cap <- if (is.null(caption)) {
    paste0(round(100 * ci_level, 1), "% interval")
  } else {
    caption
  }

  ggplot2::ggplot(df, ggplot2::aes(x = .data$n, y = .data$bound, color = .data$side)) +
    ggplot2::geom_line(linewidth = 0.8, alpha = 0.9) +
    ggplot2::scale_color_manual(values = c(Lower = "#2B8CBE", Upper = "#D95F0E")) +
    ggplot2::labs(
      title = main_title,
      subtitle = sub_cap,
      x = "Number of successful bootstrap replicates",
      y = "CI bound",
      color = ""
    ) +
    theme_obj
}


#' ggplot2 autoplot for marginal effects
#'
#' @description
#' Produces visual summaries for objects returned by
#' \code{\link{brs_marginaleffects}}.
#'
#' @param object An object of class \code{"brs_marginaleffects"}.
#' @param type Plot type: \code{"forest"}, \code{"magnitude"}, or \code{"dist"}.
#' @param variable Optional variable name for \code{type = "dist"}.
#' @param top_n Maximum number of variables shown in \code{"magnitude"}
#'   (ordered by \code{|AME|}).
#' @param title Optional plot title override.
#' @param caption Optional subtitle override.
#' @param theme Optional ggplot2 theme object. If \code{NULL},
#'   \code{ggplot2::theme_minimal()} is used.
#' @param ... Currently ignored.
#'
#' @return A \code{ggplot2} object.
#'
#' @details
#' \code{type = "dist"} requires AME simulation draws stored in
#' \code{attr(object, "ame_draws")}, which are available when marginal
#' effects are computed with \code{keep_draws = TRUE} and \code{interval = TRUE}.
#'
#' @method autoplot brs_marginaleffects
#' @export autoplot.brs_marginaleffects
autoplot.brs_marginaleffects <- function(object,
                                         type = c("forest", "magnitude", "dist"),
                                         variable = NULL,
                                         top_n = 12L,
                                         title = NULL,
                                         caption = NULL,
                                         theme = NULL,
                                         ...) {
  if (!inherits(object, "brs_marginaleffects")) {
    stop("'object' must be of class 'brs_marginaleffects'.", call. = FALSE)
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for autoplot().", call. = FALSE)
  }
  type <- match.arg(type)
  top_n <- as.integer(top_n)
  if (!is.finite(top_n) || top_n < 1L) {
    stop("'top_n' must be an integer >= 1.", call. = FALSE)
  }
  theme_obj <- .me_resolve_theme(theme)
  switch(type,
    forest = .me_autoplot_forest(object, title = title, caption = caption, theme_obj = theme_obj),
    magnitude = .me_autoplot_magnitude(object, top_n = top_n, title = title, caption = caption, theme_obj = theme_obj),
    dist = .me_autoplot_dist(object, variable = variable, title = title, caption = caption, theme_obj = theme_obj)
  )
}

#' @keywords internal
.me_resolve_theme <- function(theme) {
  if (is.null(theme)) {
    return(ggplot2::theme_minimal())
  }
  if (is.function(theme)) {
    return(theme())
  }
  theme
}

#' @keywords internal
.me_autoplot_forest <- function(object, title = NULL, caption = NULL, theme_obj = ggplot2::theme_minimal()) {
  df <- as.data.frame(object)
  df$variable <- factor(df$variable, levels = rev(df$variable))
  has_ci <- !all(is.na(df$ci.lower) | is.na(df$ci.upper))

  main_title <- if (is.null(title)) {
    paste0("Average Marginal Effects (", unique(df$model), " model, ", unique(df$type), " scale)")
  } else {
    title
  }
  sub_title <- if (is.null(caption) && has_ci) {
    paste0(round(100 * attr(object, "level"), 1), "% interval")
  } else {
    caption
  }

  p <- ggplot2::ggplot(df, ggplot2::aes(y = .data$variable))
  if (has_ci) {
    p <- p + ggplot2::geom_segment(
      ggplot2::aes(x = .data$ci.lower, xend = .data$ci.upper, yend = .data$variable),
      linewidth = 0.95, color = "#2C7FB8", alpha = 0.75
    )
  }
  p +
    ggplot2::geom_point(ggplot2::aes(x = .data$ame), size = 2.2, color = "#D95F0E") +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    ggplot2::labs(
      title = main_title,
      subtitle = sub_title,
      x = "Average marginal effect",
      y = "Variable"
    ) +
    theme_obj
}

#' @keywords internal
.me_autoplot_magnitude <- function(object, top_n = 12L, title = NULL, caption = NULL, theme_obj = ggplot2::theme_minimal()) {
  df <- as.data.frame(object)
  df$abs_ame <- abs(df$ame)
  ord <- order(df$abs_ame, decreasing = TRUE)
  df <- df[ord[seq_len(min(top_n, nrow(df)))], , drop = FALSE]
  df$variable <- factor(df$variable, levels = rev(df$variable))
  df$direction <- ifelse(df$ame >= 0, "Positive", "Negative")

  main_title <- if (is.null(title)) "AME magnitude ranking" else title
  sub_title <- if (is.null(caption)) {
    paste0("Top ", nrow(df), " variables by |AME|")
  } else {
    caption
  }

  ggplot2::ggplot(df, ggplot2::aes(x = .data$abs_ame, y = .data$variable, fill = .data$direction)) +
    ggplot2::geom_col(alpha = 0.85) +
    ggplot2::scale_fill_manual(values = c(Positive = "#1B9E77", Negative = "#D95F02")) +
    ggplot2::labs(
      title = main_title,
      subtitle = sub_title,
      x = "|AME|",
      y = "Variable",
      fill = ""
    ) +
    theme_obj
}

#' @keywords internal
.me_autoplot_dist <- function(object, variable = NULL, title = NULL, caption = NULL, theme_obj = ggplot2::theme_minimal()) {
  draws <- attr(object, "ame_draws")
  if (is.null(draws) || !is.matrix(draws) || nrow(draws) < 10L) {
    stop(
      "Distribution plot requires AME draws. Recompute with ",
      "'brs_marginaleffects(..., interval = TRUE, keep_draws = TRUE)'.",
      call. = FALSE
    )
  }
  vars <- colnames(draws)
  if (is.null(variable)) {
    variable <- vars[1L]
  } else {
    variable <- as.character(variable)[1L]
    if (!(variable %in% vars)) {
      stop(
        "'variable' must be one of: ",
        paste(vars, collapse = ", "),
        call. = FALSE
      )
    }
  }
  x <- draws[, variable]
  x <- x[is.finite(x)]
  row_i <- match(variable, object$variable)
  ame_hat <- object$ame[row_i]
  lo <- object$ci.lower[row_i]
  hi <- object$ci.upper[row_i]

  main_title <- if (is.null(title)) paste0("AME simulation distribution: ", variable) else title
  sub_title <- if (is.null(caption)) {
    paste0(
      "AME = ", format(round(ame_hat, 4), nsmall = 4),
      " | CI = [", format(round(lo, 4), nsmall = 4), ", ",
      format(round(hi, 4), nsmall = 4), "]"
    )
  } else {
    caption
  }

  ggplot2::ggplot(data.frame(ame = x), ggplot2::aes(x = .data$ame)) +
    ggplot2::geom_histogram(
      bins = 30, fill = "#74A9CF", color = "white", alpha = 0.85
    ) +
    ggplot2::geom_vline(xintercept = ame_hat, color = "#D95F0E", linewidth = 0.9) +
    ggplot2::geom_vline(xintercept = c(lo, hi), color = "#2B8CBE", linetype = "dashed") +
    ggplot2::labs(
      title = main_title,
      subtitle = sub_title,
      x = "Simulated AME draw",
      y = "Count"
    ) +
    theme_obj
}
