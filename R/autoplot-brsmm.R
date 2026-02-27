# ============================================================================ #
# ggplot2 autoplot methods for brsmm
# ============================================================================ #

#' ggplot2 autoplot for brsmm models
#'
#' @description
#' Produces ggplot2 diagnostics tailored to mixed beta interval models.
#'
#' @param object A fitted \code{"brsmm"} object.
#' @param type Plot type. One of \code{"calibration"}, \code{"score_dist"},
#'   \code{"ranef_qq"}, \code{"residuals_by_group"},
#'   \code{"ranef_caterpillar"}, \code{"ranef_density"},
#'   \code{"ranef_pairs"}, \code{"shrinkage"}, or \code{"all"}
#'   (produces all panels in a single grid).
#' @param bins Number of bins for \code{"calibration"}.
#' @param scores Optional integer vector of scores for \code{"score_dist"}.
#'   Defaults to all scores from \code{0} to \code{ncuts}.
#' @param residual_type Residual type for \code{"residuals_by_group"};
#'   passed to \code{\link{residuals.brsmm}}.
#' @param max_groups Maximum groups in \code{"residuals_by_group"}.
#' @param title Optional character: override the plot title via
#'   \code{ggplot2::labs(title = ...)}.  Ignored when \code{type = "all"}.
#' @param xlab Optional character: override the x-axis label.
#'   Ignored when \code{type = "all"}.
#' @param ylab Optional character: override the y-axis label.
#'   Ignored when \code{type = "all"}.
#' @param ncol Number of columns for the grid when \code{type = "all"}.
#'   Defaults to 2.
#' @param theme A ggplot2 theme object (e.g., \code{ggplot2::theme_bw()}) or
#'   a theme function. Applied to every panel. Defaults to
#'   \code{ggplot2::theme_minimal()}.
#' @param ... Additional arguments forwarded to \code{ggplot2::theme()} and
#'   applied on top of \code{theme}. Use named theme element arguments, e.g.
#'   \code{legend.position = "none"}.
#'
#' @return A \code{ggplot2} object, or (when \code{type = "all"}) invisibly
#'   returns the list of all panels after rendering the grid.
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
#' fit_mm <- brsmm(y ~ x1, random = ~ 1 | id, data = prep)
#' ggplot2::autoplot(fit_mm, type = "calibration", bins = 4)
#' ggplot2::autoplot(fit_mm,
#'   type = "ranef_caterpillar",
#'   title = "My title", ylab = "Mode"
#' )
#' ggplot2::autoplot(fit_mm, type = "all")
#' }
#'
#' @method autoplot brsmm
#' @importFrom rlang .data
#' @importFrom stats aggregate quantile
#' @export autoplot.brsmm
autoplot.brsmm <- function(object,
                           type = c(
                             "calibration",
                             "score_dist",
                             "ranef_qq",
                             "residuals_by_group",
                             "ranef_caterpillar",
                             "ranef_density",
                             "ranef_pairs",
                             "shrinkage",
                             "all"
                           ),
                           bins = 10L,
                           scores = NULL,
                           residual_type = c("response", "pearson"),
                           max_groups = 25L,
                           title = NULL,
                           xlab = NULL,
                           ylab = NULL,
                           ncol = 2L,
                           theme = ggplot2::theme_minimal(),
                           ...) {
  .check_class_mm(object)
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for autoplot.brsmm().", call. = FALSE)
  }

  type <- match.arg(type)
  residual_type <- match.arg(residual_type)
  bins <- as.integer(bins)
  max_groups <- as.integer(max_groups)
  ncol <- as.integer(ncol)

  if (!is.finite(bins) || bins < 3L) {
    stop("'bins' must be an integer >= 3.", call. = FALSE)
  }
  if (!is.finite(max_groups) || max_groups < 2L) {
    stop("'max_groups' must be an integer >= 2.", call. = FALSE)
  }

  # Resolve theme object (accept both function and object)
  theme_obj <- if (is.function(theme)) theme() else theme

  # Helper: apply post-hoc label overrides and ... → ggplot2::theme()
  .finalize <- function(p, ttl = title, xl = xlab, yl = ylab) {
    lbl <- Filter(Negate(is.null), list(title = ttl, x = xl, y = yl))
    if (length(lbl) > 0L) p <- p + do.call(ggplot2::labs, lbl)
    dots <- list(...)
    if (length(dots) > 0L) p <- p + do.call(ggplot2::theme, dots)
    p
  }

  # type = "all": generate every type and arrange in a grid
  if (type == "all") {
    all_types <- c(
      "calibration", "score_dist", "ranef_qq",
      "residuals_by_group", "ranef_caterpillar",
      "ranef_density", "shrinkage"
    )
    # ranef_pairs only meaningful when q_re >= 2
    if (isTRUE(object$q_re >= 2L)) all_types <- c(all_types, "ranef_pairs")

    plots <- lapply(all_types, function(tp) {
      p <- tryCatch(
        autoplot.brsmm(object,
          type = tp, theme = theme,
          bins = bins, max_groups = max_groups,
          residual_type = residual_type
        ),
        error = function(e) NULL
      )
      if (!is.null(p)) {
        dots <- list(...)
        if (length(dots) > 0L) p <- p + do.call(ggplot2::theme, dots)
      }
      p
    })
    plots <- Filter(Negate(is.null), plots)

    nrow_val <- ceiling(length(plots) / ncol)
    if (requireNamespace("gridExtra", quietly = TRUE)) {
      gridExtra::grid.arrange(grobs = plots, ncol = ncol, nrow = nrow_val)
    } else {
      message("Install 'gridExtra' for a single combined figure; printing individually.")
      for (p in plots) print(p)
    }
    return(invisible(plots))
  }

  # Single type
  p <- switch(type,
    calibration = .brsmm_autoplot_calibration(object, bins = bins, theme = theme_obj),
    score_dist = .brsmm_autoplot_score_dist(object, scores = scores, theme = theme_obj),
    ranef_qq = .brsmm_autoplot_ranef_qq(object, theme = theme_obj),
    residuals_by_group = .brsmm_autoplot_resid_group(
      object,
      residual_type = residual_type,
      max_groups = max_groups, theme = theme_obj
    ),
    ranef_caterpillar = .brsmm_autoplot_ranef_caterpillar(object, theme = theme_obj),
    ranef_density = .brsmm_autoplot_ranef_density(object, theme = theme_obj),
    ranef_pairs = .brsmm_autoplot_ranef_pairs(object, theme = theme_obj),
    shrinkage = .brsmm_autoplot_shrinkage(object, theme = theme_obj)
  )

  .finalize(p)
}


#' @keywords internal
.brsmm_autoplot_calibration <- function(object, bins = 10L,
                                        theme = ggplot2::theme_minimal()) {
  theme_obj <- if (is.function(theme)) theme() else theme
  df <- data.frame(
    observed  = as.numeric(object$Y[, "yt"]),
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

  ggplot2::ggplot(cal, ggplot2::aes(x = .data$predicted, y = .data$observed)) +
    ggplot2::geom_point(ggplot2::aes(size = .data$n), color = "#1b9e77", alpha = 0.9) +
    ggplot2::geom_line(ggplot2::aes(linewidth = .data$n), color = "#1b9e77", alpha = 0.6) +
    ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray35") +
    ggplot2::labs(
      title = "Calibration Plot (brsmm)",
      x = "Mean predicted response (bin average)",
      y = "Mean observed response (bin average)",
      size = "Bin n"
    ) +
    theme_obj
}

#' @keywords internal
.brsmm_autoplot_score_dist <- function(object, scores = NULL,
                                       theme = ggplot2::theme_minimal()) {
  theme_obj <- if (is.function(theme)) theme() else theme
  K <- as.integer(object$ncuts)
  if (is.null(scores)) scores <- 0:K
  scores <- sort(unique(as.integer(scores)))
  if (any(!is.finite(scores)) || any(scores < 0L) || any(scores > K)) {
    stop("'scores' must be integers in [0, ncuts].", call. = FALSE)
  }
  obs_scores <- .brs_observed_scores(object$Y[, "y"], K = K)
  obs_counts <- as.numeric(table(factor(obs_scores, levels = scores)))
  probs <- .brs_score_prob_matrix(
    mu = object$fitted_mu, phi = object$fitted_phi,
    repar = object$repar, ncuts = K, lim = object$lim, scores = scores
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
      x = "Score", y = "Count", fill = ""
    ) +
    theme_obj
}

#' @keywords internal
.brsmm_autoplot_ranef_qq <- function(object,
                                     theme = ggplot2::theme_minimal()) {
  theme_obj <- if (is.function(theme)) theme() else theme
  re <- .brsmm_ranef_matrix(object)
  parts <- lapply(colnames(re), function(term) {
    x <- sort(as.numeric(re[, term]))
    n <- length(x)
    q <- stats::qnorm(stats::ppoints(n))
    data.frame(
      theoretical = q, sample = x, term = term,
      intercept = mean(x) - stats::sd(x) * mean(q), slope = stats::sd(x)
    )
  })
  df <- do.call(rbind, parts)
  ggplot2::ggplot(df, ggplot2::aes(x = .data$theoretical, y = .data$sample)) +
    ggplot2::geom_point(color = "#1b9e77", alpha = 0.85, size = 1.1) +
    ggplot2::geom_abline(ggplot2::aes(intercept = .data$intercept, slope = .data$slope),
      linetype = "dashed", color = "gray35"
    ) +
    ggplot2::facet_wrap(~term, scales = "free_y") +
    ggplot2::labs(
      title = "Random-Effects Q-Q Plot",
      x = "Theoretical Normal Quantiles",
      y = "Estimated random-effect mode"
    ) +
    theme_obj
}

#' @keywords internal
.brsmm_autoplot_resid_group <- function(object,
                                        residual_type = c("response", "pearson"),
                                        max_groups = 25L,
                                        theme = ggplot2::theme_minimal()) {
  theme_obj <- if (is.function(theme)) theme() else theme
  residual_type <- match.arg(residual_type)
  res <- residuals(object, type = residual_type)
  grp <- object$group
  tab <- sort(table(grp), decreasing = TRUE)
  keep <- names(tab)[seq_len(min(max_groups, length(tab)))]
  df <- data.frame(residual = res, group = factor(as.character(grp), levels = keep))
  df <- df[!is.na(df$group), , drop = FALSE]
  ggplot2::ggplot(df, ggplot2::aes(x = .data$group, y = .data$residual, fill = .data$group)) +
    ggplot2::geom_boxplot(alpha = 0.70, outlier.shape = NA) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray35") +
    ggplot2::labs(
      title = paste0("Residuals by Group (top ", length(keep), " groups, ", residual_type, ")"),
      x = "Group", y = "Residual"
    ) +
    theme_obj +
    ggplot2::theme(
      legend.position = "none",
      axis.text.x = ggplot2::element_text(angle = 60, hjust = 1)
    )
}

#' @keywords internal
.brsmm_ranef_matrix <- function(object) {
  re <- object$random$mode_b
  if (is.matrix(re)) {
    return(re)
  }
  out <- matrix(as.numeric(re), ncol = 1L)
  rownames(out) <- names(re)
  cn <- object$random$terms
  if (is.null(cn) || length(cn) == 0L) cn <- "(Intercept)"
  colnames(out) <- cn[1L]
  out
}

#' @keywords internal
.brsmm_autoplot_ranef_caterpillar <- function(object,
                                              theme = ggplot2::theme_minimal()) {
  theme_obj <- if (is.function(theme)) theme() else theme
  re <- .brsmm_ranef_matrix(object)
  D <- object$random$D
  if (is.null(D)) {
    sd_single <- as.numeric(object$random$sd_b[1L])
    D <- matrix(sd_single^2, 1L, 1L)
  }
  model_sd <- sqrt(pmax(diag(D), 0))

  parts <- lapply(seq_along(colnames(re)), function(ci) {
    term <- colnames(re)[ci]
    v <- as.numeric(re[, ci])
    msd <- if (ci <= length(model_sd)) model_sd[ci] else model_sd[1L]
    ord <- order(v)
    data.frame(
      group = factor(rownames(re)[ord], levels = rownames(re)[ord]),
      mode  = v[ord],
      lo    = v[ord] - 1.96 * msd,
      hi    = v[ord] + 1.96 * msd,
      term  = term
    )
  })
  df <- do.call(rbind, parts)

  ggplot2::ggplot(df, ggplot2::aes(x = .data$group, y = .data$mode)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray35") +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = .data$lo, ymax = .data$hi),
      width = 0.3, color = "#7570b3", alpha = 0.5
    ) +
    ggplot2::geom_segment(
      ggplot2::aes(xend = .data$group, y = 0, yend = .data$mode),
      color = "#7570b3", alpha = 0.4
    ) +
    ggplot2::geom_point(color = "#1b9e77", size = 1.8) +
    ggplot2::facet_wrap(~term, scales = "free_y") +
    ggplot2::labs(
      title    = "Caterpillar Plot of Random Effects",
      subtitle = "Error bars: \u00b11.96 \u00d7 Model SD (marginal)",
      x        = "Group (ordered by mode)",
      y        = "Random-effect mode"
    ) +
    theme_obj +
    ggplot2::theme(axis.text.x = ggplot2::element_blank())
}

#' @keywords internal
.brsmm_autoplot_ranef_density <- function(object,
                                          theme = ggplot2::theme_minimal()) {
  theme_obj <- if (is.function(theme)) theme() else theme
  re <- .brsmm_ranef_matrix(object)
  df <- data.frame(value = as.numeric(re), term = rep(colnames(re), each = nrow(re)))
  ggplot2::ggplot(df, ggplot2::aes(x = .data$value, color = .data$term, fill = .data$term)) +
    ggplot2::geom_density(alpha = 0.20, linewidth = 0.9) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "gray35") +
    ggplot2::labs(
      title = "Density of Random Effects", x = "Random-effect mode",
      y = "Density", color = "Term", fill = "Term"
    ) +
    theme_obj
}

#' @keywords internal
.brsmm_autoplot_ranef_pairs <- function(object,
                                        theme = ggplot2::theme_minimal()) {
  theme_obj <- if (is.function(theme)) theme() else theme
  re <- .brsmm_ranef_matrix(object)
  if (ncol(re) < 2L) stop("ranef_pairs requires at least two random-effect terms.", call. = FALSE)
  df <- data.frame(x = re[, 1L], y = re[, 2L])
  xlab <- colnames(re)[1L]
  ylab <- colnames(re)[2L]
  ggplot2::ggplot(df, ggplot2::aes(x = .data$x, y = .data$y)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray35") +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "gray35") +
    ggplot2::geom_point(color = "#1b9e77", alpha = 0.8) +
    ggplot2::geom_smooth(method = "lm", se = FALSE, color = "#d95f02", linewidth = 0.8) +
    ggplot2::labs(title = "Random-Effects Pair Plot", x = xlab, y = ylab) +
    theme_obj
}

#' @keywords internal
.brsmm_autoplot_shrinkage <- function(object,
                                      theme = ggplot2::theme_minimal()) {
  theme_obj <- if (is.function(theme)) theme() else theme
  re <- .brsmm_ranef_matrix(object)
  grp <- as.character(object$group)

  # Naive per-group estimate: mean of logit-transformed midpoint response
  yt <- as.numeric(object$Y[, "yt"])
  yt_safe <- pmin(pmax(yt, 1e-6), 1 - 1e-6)
  logit_y <- log(yt_safe / (1 - yt_safe))

  grp_levels <- rownames(re)
  naive_vals <- vapply(grp_levels, function(g) {
    idx <- which(grp == g)
    if (length(idx) == 0L) {
      return(NA_real_)
    }
    mean(logit_y[idx], na.rm = TRUE)
  }, numeric(1L))

  # Global logit-mean (fixed-effects anchor)
  global_mean <- mean(logit_y, na.rm = TRUE)
  naive_re <- naive_vals - global_mean # naive RE = deviation from global

  parts <- lapply(seq_along(colnames(re)), function(ci) {
    term <- colnames(re)[ci]
    modes <- as.numeric(re[, ci])
    data.frame(
      naive  = naive_re,
      mode   = modes,
      term   = term,
      group  = grp_levels
    )
  })
  df <- do.call(rbind, parts)
  df <- df[is.finite(df$naive), , drop = FALSE]

  ggplot2::ggplot(df, ggplot2::aes(x = .data$naive, y = .data$mode)) +
    ggplot2::geom_abline(
      intercept = 0, slope = 1,
      linetype = "dashed", color = "gray40"
    ) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dotted", color = "gray60") +
    ggplot2::geom_vline(xintercept = 0, linetype = "dotted", color = "gray60") +
    ggplot2::geom_smooth(
      method = "loess", se = FALSE,
      color = "#d95f02", linewidth = 0.8, formula = y ~ x
    ) +
    ggplot2::geom_point(color = "#1b9e77", alpha = 0.8, size = 1.6) +
    ggplot2::facet_wrap(~term, scales = "free") +
    ggplot2::labs(
      title    = "Shrinkage Plot",
      subtitle = "Laplace mode vs. naive per-group logit-mean deviation",
      x        = "Naive estimate (deviation from global mean)",
      y        = "Laplace mode (random-effect mode)"
    ) +
    theme_obj
}
