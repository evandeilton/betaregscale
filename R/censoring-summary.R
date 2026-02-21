# ============================================================================ #
# Censoring summary - visual and tabular description of censoring types
# ============================================================================ #

#' Graphical and tabular censoring summary
#'
#' @description
#' Produces a visual summary of the censoring structure in a fitted
#' \code{"brs"} model or a response matrix produced by
#' \code{\link{brs_check}}.  The summary includes:
#' \enumerate{
#'   \item Bar chart of censoring type counts
#'   \item Histogram of midpoint responses colored by censoring type
#'   \item Interval plot showing \eqn{[l_i, u_i]} segments
#'   \item Proportion table of censoring types
#' }
#'
#' @param object A fitted \code{"brs"} object, a matrix
#'   returned by \code{\link{brs_check}}, or a data frame
#'   returned by \code{\link{brs_prep}} (must contain columns
#'   \code{left}, \code{right}, \code{yt}, and \code{delta}).
#' @param n_sample Integer: maximum number of observations to show
#'   in the interval plot (default 100).  If the data has more
#'   observations, a random sample is drawn.
#' @param which Integer vector selecting which panels to draw
#'   (default \code{1:4}).
#' @param caption Optional panel captions. Accepts a character vector
#'   (or list coercible to character) with up to 4 labels, in the order:
#'   burden, midpoint-by-type, width-by-type, ordered interval map.
#' @param gg     Logical: use ggplot2? (default \code{FALSE}).
#' @param title Optional global title for the plotting page.
#' @param sub.caption Optional subtitle/caption for the plotting page.
#' @param theme Optional ggplot2 theme object (e.g., \code{ggplot2::theme_bw()}).
#'   If \code{NULL}, a minimal theme is used when \code{gg = TRUE}.
#' @param palette Optional named character vector with colors for censoring
#'   types \code{Exact}, \code{Left}, \code{Right}, and \code{Interval}.
#' @param inform Logical; if \code{TRUE}, prints brief interpretation
#'   messages about boundary and interval censoring intensity.
#' @param ...    Further arguments (currently ignored).
#'
#' @return Invisibly returns a data frame with censoring counts and
#'   proportions, percentages, and interpretation flags.
#'
#' @examples
#' y <- c(0, 3, 5, 7, 10)
#' Y <- brs_check(y, ncuts = 10)
#' brs_cens(Y)
#'
#' prep <- brs_prep(data.frame(y = y), ncuts = 10)
#' brs_cens(prep)
#'
#' @references
#' Lopes, J. E. (2023). \emph{Modelos de regressao beta para dados de escala}.
#' Master's dissertation, Universidade Federal do Parana, Curitiba.
#' URI: \url{https://hdl.handle.net/1884/86624}.
#'
#' Hawker, G. A., Mian, S., Kendzerska, T., and French, M. (2011).
#' Measures of adult pain: Visual Analog Scale for Pain (VAS Pain),
#' Numeric Rating Scale for Pain (NRS Pain), McGill Pain Questionnaire (MPQ),
#' Short-Form McGill Pain Questionnaire (SF-MPQ), Chronic Pain Grade Scale
#' (CPGS), Short Form-36 Bodily Pain Scale (SF-36 BPS), and Measure of
#' Intermittent and Constant Osteoarthritis Pain (ICOAP).
#' Arthritis Care and Research, 63(S11), S240-S252.
#' \doi{10.1002/acr.20543}
#'
#' Hjermstad, M. J., Fayers, P. M., Haugen, D. F., et al. (2011).
#' Studies comparing Numerical Rating Scales, Verbal Rating Scales, and
#' Visual Analogue Scales for assessment of pain intensity in adults:
#' a systematic literature review.
#' Journal of Pain and Symptom Management, 41(6), 1073-1093.
#' \doi{10.1016/j.jpainsymman.2010.08.016}
#'
#' @importFrom graphics barplot hist par layout legend plot.new plot.window text title segments
#' @importFrom grDevices adjustcolor
#' @rdname brs_cens
#' @export
brs_cens <- function(object,
                     n_sample = 100L,
                     which = 1:4,
                     caption = NULL,
                     gg = FALSE,
                     title = "Censoring diagnostic overview",
                     sub.caption = NULL,
                     theme = NULL,
                     palette = NULL,
                     inform = FALSE,
                     ...) {
  # Extract Y and delta

  if (inherits(object, "brs")) {
    Y <- object$Y
    delta <- object$delta
  } else if (is.data.frame(object) && isTRUE(attr(object, "is_prepared")) &&
    all(c("left", "right", "yt", "delta") %in% names(object))) {
    Y <- as.matrix(object[, c("left", "right", "yt", "delta")])
    delta <- as.integer(object[["delta"]])
  } else if (is.matrix(object) && all(c("left", "right", "yt", "delta") %in% colnames(object))) {
    Y <- object
    delta <- as.integer(object[, "delta"])
  } else {
    stop(
      "object must be a 'brs' object, a matrix from brs_check(), ",
      "or a data.frame from brs_prep().",
      call. = FALSE
    )
  }

  left <- Y[, "left"]
  right <- Y[, "right"]
  yt <- Y[, "yt"]
  n <- length(delta)

  # Censoring labels and colors
  type_labels <- c("0" = "Exact", "1" = "Left", "2" = "Right", "3" = "Interval")
  type_colors <- c(
    "0" = "#009E73", # Exact
    "1" = "#D55E00", # Left
    "2" = "#0072B2", # Right
    "3" = "#CC79A7" # Interval
  )
  if (!is.null(palette)) {
    expected <- c("Exact", "Left", "Right", "Interval")
    if (!all(expected %in% names(palette))) {
      stop(
        "palette must be a named vector containing: ",
        paste(expected, collapse = ", "),
        call. = FALSE
      )
    }
    type_colors <- c(
      "0" = palette[["Exact"]],
      "1" = palette[["Left"]],
      "2" = palette[["Right"]],
      "3" = palette[["Interval"]]
    )
  }
  type_colors_lbl <- stats::setNames(
    unname(type_colors[c("0", "1", "2", "3")]),
    unname(type_labels[c("0", "1", "2", "3")])
  )

  type_factor <- factor(delta, levels = 0:3, labels = type_labels)
  width <- pmax(right - left, 0)
  which <- as.integer(which)
  which <- which[which %in% 1:4]
  if (length(which) == 0L) {
    stop("which must contain at least one value in 1:4.", call. = FALSE)
  }
  caption_default <- c(
    "Censoring burden by type",
    "Midpoint score by censoring type",
    "Interval width by censoring type",
    "Ordered interval map (by midpoint)"
  )
  if (is.null(caption)) {
    caption <- caption_default
  } else {
    caption <- as.character(unlist(caption, use.names = FALSE))
    if (length(caption) < 4L) {
      caption <- c(caption, caption_default[(length(caption) + 1L):4L])
    }
    caption <- caption[1:4]
  }
  df_obs <- data.frame(
    left = left,
    right = right,
    yt = yt,
    delta = delta,
    type = type_factor,
    width = width
  )

  # Summary table
  tab <- table(type_factor)
  props <- prop.table(tab)
  summary_df <- data.frame(
    type       = names(tab),
    count      = as.integer(tab),
    proportion = as.numeric(props),
    percentage = as.numeric(props) * 100,
    row.names  = NULL
  )
  summary_df$severity <- cut(
    summary_df$proportion,
    breaks = c(-Inf, 0.10, 0.30, Inf),
    labels = c("low", "moderate", "high"),
    right = FALSE
  )
  summary_df$interpretation <- vapply(seq_len(nrow(summary_df)), function(i) {
    tp <- summary_df$type[i]
    sev <- as.character(summary_df$severity[i])
    if (tp %in% c("Left", "Right")) {
      paste0(
        "Boundary censoring in this tail is ", sev,
        "; check possible floor/ceiling concentration."
      )
    } else if (tp == "Interval") {
      paste0(
        "Interval uncertainty is ", sev,
        "; interpret effects considering measurement granularity."
      )
    } else {
      paste0("Exact observations are ", sev, ".")
    }
  }, character(1))

  # Generic diagnostic messages for any scale-censored domain.
  prop_left <- summary_df$proportion[summary_df$type == "Left"]
  prop_right <- summary_df$proportion[summary_df$type == "Right"]
  prop_interval <- summary_df$proportion[summary_df$type == "Interval"]
  prop_boundary <- sum(prop_left, prop_right, na.rm = TRUE)
  if (isTRUE(inform) && is.finite(prop_boundary) && prop_boundary >= 0.30) {
    message(
      "Boundary censoring is substantial (",
      sprintf("%.1f%%", 100 * prop_boundary),
      "). Consider floor/ceiling-aware model interpretation."
    )
  }
  if (isTRUE(inform) && is.finite(prop_interval) && prop_interval >= 0.50) {
    message(
      "Interval-censoring is dominant (",
      sprintf("%.1f%%", 100 * prop_interval),
      "). Uncertainty width diagnostics are especially important."
    )
  }

  if (gg) {
    .brs_cens_gg(
      df_obs = df_obs,
      tab = tab,
      type_labels, type_colors, type_colors_lbl, n_sample,
      summary_df,
      which = which,
      caption = caption,
      title = title,
      sub.caption = sub.caption,
      theme = theme
    )
  } else {
    .brs_cens_base(
      df_obs = df_obs,
      tab = tab,
      type_labels, type_colors, type_colors_lbl, n_sample,
      summary_df,
      which = which,
      caption = caption,
      title = title,
      sub.caption = sub.caption
    )
  }

  invisible(summary_df)
}


# -- Base R implementation -------------------------------------------------- #

.brs_cens_base <- function(df_obs, tab, type_labels, type_colors, type_colors_lbl, n_sample,
                           summary_df, which, caption, title, sub.caption) {
  yt <- df_obs$yt
  left <- df_obs$left
  right <- df_obs$right
  type_factor <- df_obs$type
  width <- df_obs$width
  delta <- df_obs$delta
  n <- nrow(df_obs)

  nplots <- length(which)
  ncol <- min(nplots, 2L)
  nrow <- ceiling(nplots / ncol)
  op <- par(mfrow = c(nrow, ncol), mar = c(4, 4, 3, 1), oma = c(0, 0, 3, 0))
  on.exit(par(op))

  # 1. Censoring burden (count + proportion labels)
  if (1L %in% which) {
    cols <- unname(type_colors_lbl[names(tab)])
    bp <- barplot(tab,
      col = cols, main = caption[1L],
      ylab = "Count", las = 1
    )
    prop_vals <- prop.table(tab)
    text(
      x = bp, y = as.numeric(tab),
      labels = sprintf("%.1f%%", 100 * as.numeric(prop_vals)),
      pos = 3, cex = 0.8
    )
  }

  # 2. Midpoint distribution by type (boxplot + points)
  if (2L %in% which) {
    cols <- unname(type_colors_lbl[levels(type_factor)])
    graphics::boxplot(yt ~ type_factor,
      col = adjustcolor(cols, alpha.f = 0.5),
      border = cols,
      main = caption[2L],
      ylab = "Midpoint response", xlab = "", ylim = c(0, 1)
    )
    graphics::stripchart(yt ~ type_factor,
      method = "jitter",
      vertical = TRUE,
      pch = 16, cex = 0.55,
      col = adjustcolor("black", alpha.f = 0.35),
      add = TRUE
    )
  }

  # 3. Interval width by type
  if (3L %in% which) {
    cols <- unname(type_colors_lbl[levels(type_factor)])
    graphics::boxplot(width ~ type_factor,
      col = adjustcolor(cols, alpha.f = 0.5),
      border = cols,
      main = caption[3L],
      ylab = "Width (right - left)", xlab = ""
    )
    abline(h = 0, lty = 2, col = "gray45")
  }

  # 4. Ordered interval map (sampled if needed)
  if (4L %in% which) {
    if (n > n_sample) {
      samp <- sort(sample.int(n, n_sample))
    } else {
      samp <- seq_len(n)
    }
    df_int <- df_obs[samp, , drop = FALSE]
    ord <- order(df_int$yt, df_int$width)
    df_int <- df_int[ord, , drop = FALSE]
    ns <- nrow(df_int)

    cols_seg <- type_colors[as.character(df_int$delta)]
    plot(NA,
      xlim = c(0, 1), ylim = c(1, ns),
      xlab = "Unit interval", ylab = "Ordered observation",
      main = caption[4L]
    )
    segments(
      x0 = df_int$left, y0 = seq_len(ns),
      x1 = df_int$right, y1 = seq_len(ns),
      col = cols_seg, lwd = 1
    )
    graphics::points(df_int$yt, seq_len(ns), pch = 16, cex = 0.45, col = "gray20")
    legend("bottomright",
      legend = c("Exact", "Left", "Right", "Interval"),
      col = type_colors[c("0", "1", "2", "3")], lwd = 2, cex = 0.75, bty = "n"
    )
  }

  if (!is.null(sub.caption) && nzchar(sub.caption)) {
    mtext(sprintf("%s | %s", title, sub.caption), outer = TRUE, cex = 1.0, font = 2)
  } else {
    mtext(title, outer = TRUE, cex = 1.1, font = 2)
  }
}


# -- ggplot2 implementation ------------------------------------------------- #

.brs_cens_gg <- function(df_obs, tab, type_labels, type_colors, type_colors_lbl, n_sample,
                         summary_df, which, caption, title, sub.caption, theme) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for gg = TRUE.", call. = FALSE)
  }

  n <- nrow(df_obs)
  theme_obj <- .brs_cens_resolve_theme(theme)

  plots <- list()

  # 1. Censoring burden
  if (1L %in% which) {
    df_bar <- data.frame(
      type = factor(names(tab), levels = levels(df_obs$type)),
      count = as.numeric(tab),
      proportion = as.numeric(prop.table(tab))
    )
    plots[[length(plots) + 1L]] <- ggplot2::ggplot(df_bar, ggplot2::aes(x = .data$type, y = .data$count, fill = .data$type)) +
      ggplot2::geom_col(width = 0.75) +
      ggplot2::geom_text(
        ggplot2::aes(label = sprintf("%.1f%%", 100 * .data$proportion)),
        vjust = -0.35, size = 3.6
      ) +
      ggplot2::scale_fill_manual(values = type_colors_lbl, drop = FALSE) +
      ggplot2::labs(title = caption[1L], x = "", y = "Count") +
      theme_obj +
      ggplot2::theme(legend.position = "none")
  }

  # 2. Midpoint distribution by type
  if (2L %in% which) {
    plots[[length(plots) + 1L]] <- ggplot2::ggplot(df_obs, ggplot2::aes(x = .data$type, y = .data$yt, color = .data$type, fill = .data$type)) +
      ggplot2::geom_boxplot(alpha = 0.25, outlier.shape = NA, width = 0.6) +
      ggplot2::geom_jitter(width = 0.12, height = 0, alpha = 0.35, size = 1) +
      ggplot2::scale_color_manual(values = type_colors_lbl, drop = FALSE) +
      ggplot2::scale_fill_manual(values = type_colors_lbl, drop = FALSE) +
      ggplot2::coord_cartesian(ylim = c(0, 1)) +
      ggplot2::labs(
        title = caption[2L],
        x = "", y = "Midpoint response"
      ) +
      theme_obj +
      ggplot2::theme(legend.position = "none")
  }

  # 3. Interval width by type
  if (3L %in% which) {
    plots[[length(plots) + 1L]] <- ggplot2::ggplot(df_obs, ggplot2::aes(x = .data$type, y = .data$width, color = .data$type, fill = .data$type)) +
      ggplot2::geom_boxplot(alpha = 0.25, outlier.shape = NA, width = 0.6) +
      ggplot2::geom_jitter(width = 0.12, height = 0, alpha = 0.35, size = 1) +
      ggplot2::scale_color_manual(values = type_colors_lbl, drop = FALSE) +
      ggplot2::scale_fill_manual(values = type_colors_lbl, drop = FALSE) +
      ggplot2::labs(
        title = caption[3L],
        x = "", y = "Width (right - left)"
      ) +
      theme_obj +
      ggplot2::theme(legend.position = "none")
  }

  # 4. Ordered interval map
  if (4L %in% which) {
    if (n > n_sample) {
      samp <- sort(sample.int(n, n_sample))
    } else {
      samp <- seq_len(n)
    }
    df_int <- df_obs[samp, , drop = FALSE]
    ord <- order(df_int$yt, df_int$width)
    df_int <- df_int[ord, , drop = FALSE]
    df_int$idx <- seq_len(nrow(df_int))

    plots[[length(plots) + 1L]] <- ggplot2::ggplot(df_int) +
      ggplot2::geom_segment(
        ggplot2::aes(
          x = .data$left, xend = .data$right,
          y = .data$idx, yend = .data$idx,
          color = .data$type
        ),
        linewidth = 0.5, alpha = 0.85
      ) +
      ggplot2::geom_point(
        ggplot2::aes(x = .data$yt, y = .data$idx),
        color = "gray20", size = 0.85, alpha = 0.8
      ) +
      ggplot2::scale_color_manual(values = type_colors_lbl, drop = FALSE) +
      ggplot2::labs(
        title = caption[4L],
        x = "Unit interval", y = "Ordered observation", color = "Type"
      ) +
      theme_obj
  }

  np <- length(plots)
  ncol <- min(np, 2L)
  nrow <- ceiling(np / ncol)
  if (requireNamespace("gridExtra", quietly = TRUE)) {
    gridExtra::grid.arrange(
      grobs = plots, ncol = ncol, nrow = nrow,
      top = if (!is.null(title) && nzchar(title)) title else NULL,
      bottom = if (!is.null(sub.caption) && nzchar(sub.caption)) sub.caption else NULL
    )
  } else {
    for (p in plots) print(p)
  }
}

.brs_cens_resolve_theme <- function(theme) {
  if (is.null(theme)) {
    return(ggplot2::theme_minimal())
  }
  if (is.function(theme)) {
    return(theme())
  }
  theme
}
