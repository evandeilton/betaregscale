# ============================================================================ #
# Censoring summary — visual and tabular description of censoring types
# ============================================================================ #

#' Graphical and tabular censoring summary
#'
#' @description
#' Produces a visual summary of the censoring structure in a fitted
#' \code{"betaregscale"} model or a response matrix produced by
#' \code{\link{check_response}}.  The summary includes:
#' \enumerate{
#'   \item Bar chart of censoring type counts
#'   \item Histogram of midpoint responses colored by censoring type
#'   \item Interval plot showing \eqn{[l_i, u_i]} segments
#'   \item Proportion table of censoring types
#' }
#'
#' @param object A fitted \code{"betaregscale"} object, or a matrix
#'   returned by \code{\link{check_response}} (must contain columns
#'   \code{left}, \code{right}, \code{yt}, and \code{delta}).
#' @param n_sample Integer: maximum number of observations to show
#'   in the interval plot (default 100).  If the data has more
#'   observations, a random sample is drawn.
#' @param gg     Logical: use ggplot2? (default \code{FALSE}).
#' @param ...    Further arguments (currently ignored).
#'
#' @return Invisibly returns a data frame with censoring counts and
#'   proportions.
#'
#' @examples
#' \dontrun{
#' fit <- betaregscale(y ~ x1 + x2, data = sim)
#' censoring_summary(fit)
#' censoring_summary(fit, gg = TRUE)
#' }
#'
#' @importFrom graphics barplot hist par layout legend plot.new plot.window text title segments
#' @importFrom grDevices adjustcolor
#' @export
censoring_summary <- function(object, n_sample = 100L, gg = FALSE, ...) {
  # Extract Y and delta

  if (inherits(object, "betaregscale")) {
    Y <- object$Y
    delta <- object$delta
  } else if (is.matrix(object) && all(c("left", "right", "yt", "delta") %in% colnames(object))) {
    Y <- object
    delta <- as.integer(object[, "delta"])
  } else {
    stop(
      "object must be a 'betaregscale' object or a matrix from check_response().",
      call. = FALSE
    )
  }

  left <- Y[, "left"]
  right <- Y[, "right"]
  yt <- Y[, "yt"]
  n <- length(delta)

  # Censoring labels
  type_labels <- c("0" = "Exact", "1" = "Left", "2" = "Right", "3" = "Interval")
  type_colors <- c("0" = "#1b9e77", "1" = "#d95f02", "2" = "#7570b3", "3" = "#e7298a")

  type_factor <- factor(delta, levels = 0:3, labels = type_labels)

  # Summary table
  tab <- table(type_factor)
  props <- prop.table(tab)
  summary_df <- data.frame(
    type       = names(tab),
    count      = as.integer(tab),
    proportion = as.numeric(props),
    row.names  = NULL
  )

  if (gg) {
    .censoring_summary_gg(
      yt, left, right, delta, type_factor,
      type_labels, type_colors, n_sample,
      summary_df
    )
  } else {
    .censoring_summary_base(
      yt, left, right, delta, type_factor,
      type_labels, type_colors, n_sample,
      summary_df
    )
  }

  invisible(summary_df)
}


# -- Base R implementation -------------------------------------------------- #

.censoring_summary_base <- function(yt, left, right, delta, type_factor,
                                    type_labels, type_colors, n_sample,
                                    summary_df) {
  n <- length(delta)

  op <- par(mfrow = c(2, 2), mar = c(4, 4, 3, 1), oma = c(0, 0, 2, 0))
  on.exit(par(op))

  # 1. Bar chart of counts
  tab <- table(type_factor)
  cols <- type_colors[as.character(0:3)][levels(type_factor) %in% names(tab)]
  barplot(tab,
    col = cols, main = "Censoring type counts",
    ylab = "Count", las = 1
  )

  # 2. Histogram colored by type
  # Use stacked approach: overlay histograms
  br <- seq(0, 1, length.out = 30)
  all_types <- levels(type_factor)
  existing <- all_types[tab > 0]

  hist(yt,
    breaks = br, col = "gray90", border = "gray70",
    main = "Response by censoring type",
    xlab = "Midpoint response", ylab = "Frequency"
  )
  for (k in seq_along(existing)) {
    idx <- which(type_factor == existing[k])
    if (length(idx) > 0) {
      hist(yt[idx],
        breaks = br, col = adjustcolor(type_colors[as.character(which(type_labels == existing[k]) - 1)], alpha.f = 0.4),
        border = NA, add = TRUE
      )
    }
  }
  legend("topright",
    legend = existing[tab[existing] > 0],
    fill = adjustcolor(type_colors[as.character(match(existing[tab[existing] > 0], type_labels) - 1)], alpha.f = 0.4),
    cex = 0.7, bty = "n"
  )

  # 3. Interval plot
  if (n > n_sample) {
    samp <- sort(sample.int(n, n_sample))
  } else {
    samp <- seq_len(n)
  }
  ns <- length(samp)
  cols_seg <- type_colors[as.character(delta[samp])]
  plot(NA,
    xlim = c(0, 1), ylim = c(1, ns),
    xlab = "Unit interval", ylab = "Observation",
    main = "Censoring intervals"
  )
  segments(
    x0 = left[samp], y0 = seq_len(ns),
    x1 = right[samp], y1 = seq_len(ns),
    col = cols_seg, lwd = 0.8
  )

  # 4. Proportion summary as text table
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1))
  title(main = "Censoring proportions")
  y_pos <- seq(0.85, by = -0.15, length.out = nrow(summary_df))
  for (i in seq_len(nrow(summary_df))) {
    text(0.15, y_pos[i],
      sprintf(
        "%s: %d (%.1f%%)",
        summary_df$type[i],
        summary_df$count[i],
        summary_df$proportion[i] * 100
      ),
      adj = 0, cex = 1.1
    )
  }
  text(0.15, y_pos[nrow(summary_df)] - 0.15,
    sprintf("Total: %d", sum(summary_df$count)),
    adj = 0, cex = 1.1, font = 2
  )

  mtext("Censoring Summary", outer = TRUE, cex = 1.2, font = 2)
}


# -- ggplot2 implementation ------------------------------------------------- #

.censoring_summary_gg <- function(yt, left, right, delta, type_factor,
                                  type_labels, type_colors, n_sample,
                                  summary_df) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for gg = TRUE.", call. = FALSE)
  }

  n <- length(delta)

  # 1. Bar chart
  df_bar <- data.frame(type = type_factor)
  p1 <- ggplot2::ggplot(df_bar, ggplot2::aes(x = .data$type, fill = .data$type)) +
    ggplot2::geom_bar() +
    ggplot2::scale_fill_manual(values = type_colors, drop = FALSE) +
    ggplot2::labs(title = "Censoring type counts", x = "", y = "Count") +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "none")

  # 2. Histogram
  df_hist <- data.frame(yt = yt, type = type_factor)
  p2 <- ggplot2::ggplot(df_hist, ggplot2::aes(x = .data$yt, fill = .data$type)) +
    ggplot2::geom_histogram(bins = 30, alpha = 0.7, position = "identity") +
    ggplot2::scale_fill_manual(values = type_colors, drop = FALSE) +
    ggplot2::labs(
      title = "Response by censoring type",
      x = "Midpoint response", y = "Frequency", fill = "Type"
    ) +
    ggplot2::theme_minimal()

  # 3. Interval plot
  if (n > n_sample) {
    samp <- sort(sample.int(n, n_sample))
  } else {
    samp <- seq_len(n)
  }
  df_int <- data.frame(
    left = left[samp], right = right[samp],
    idx = seq_along(samp), type = type_factor[samp]
  )
  p3 <- ggplot2::ggplot(df_int, ggplot2::aes(
    y = .data$idx, xmin = .data$left, xmax = .data$right,
    color = .data$type
  )) +
    ggplot2::geom_linerange(linewidth = 0.4) +
    ggplot2::scale_color_manual(values = type_colors, drop = FALSE) +
    ggplot2::labs(
      title = "Censoring intervals",
      x = "Unit interval", y = "Observation", color = "Type"
    ) +
    ggplot2::theme_minimal()

  # 4. Table as ggplot text
  summary_df$label <- sprintf(
    "%s: %d (%.1f%%)",
    summary_df$type, summary_df$count, summary_df$proportion * 100
  )
  summary_df$y <- rev(seq_len(nrow(summary_df)))
  p4 <- ggplot2::ggplot(summary_df, ggplot2::aes(x = 0, y = .data$y, label = .data$label)) +
    ggplot2::geom_text(hjust = 0, size = 4) +
    ggplot2::annotate("text",
      x = 0, y = 0,
      label = sprintf("Total: %d", sum(summary_df$count)),
      hjust = 0, size = 4, fontface = "bold"
    ) +
    ggplot2::xlim(-0.1, 1) +
    ggplot2::ylim(-0.5, max(summary_df$y) + 0.5) +
    ggplot2::labs(title = "Censoring proportions") +
    ggplot2::theme_void()

  plots <- list(p1, p2, p3, p4)

  if (requireNamespace("gridExtra", quietly = TRUE)) {
    gridExtra::grid.arrange(
      grobs = plots, ncol = 2, nrow = 2,
      top = "Censoring Summary"
    )
  } else {
    for (p in plots) print(p)
  }
}
