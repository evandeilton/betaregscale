# ============================================================================ #
# Diagnostic plots for brsmm objects
# ============================================================================ #

#' Diagnostic plots for mixed beta interval regression
#'
#' @description
#' Produces diagnostic plots for fitted \code{"brsmm"} models:
#' residuals vs indices, Cook's distance, residuals vs linear predictor,
#' residuals vs fitted values, half-normal envelope, and predicted vs observed.
#'
#' @param x A fitted \code{"brsmm"} object.
#' @param which Integer vector selecting which panels to draw
#'   (default \code{1:4}).
#' @param type Residual type passed to \code{\link{residuals.brsmm}}
#'   (\code{"response"} or \code{"pearson"}).
#' @param nsim Number of simulations for half-normal envelope.
#' @param level Confidence level for the half-normal envelope.
#' @param caption Character vector of plot captions.
#' @param sub.caption Optional subtitle; defaults to model call.
#' @param ask Logical: prompt before each new page?
#' @param gg Logical: use ggplot2 backend?
#' @param title Optional global title for ggplot output. If \code{NULL},
#'   panel captions are used.
#' @param theme Optional ggplot2 theme object (e.g., \code{ggplot2::theme_bw()}).
#'   If \code{NULL}, a minimal theme is used.
#' @param ... Further arguments passed to base \code{plot()}.
#'
#' @return Invisibly returns \code{x}.
#'
#' @seealso \code{\link{brsmm}}, \code{\link{residuals.brsmm}},
#'   \code{\link{autoplot.brsmm}}
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
#' plot(fit, which = 1:4)
#' }
#'
#' @method plot brsmm
#' @importFrom stats fitted residuals predict quantile qqline
#' @importFrom graphics plot abline par mtext dotchart
#' @importFrom grDevices dev.interactive devAskNewPage
#' @export
plot.brsmm <- function(x,
                       which = 1:4,
                       type = c("response", "pearson"),
                       nsim = 100L,
                       level = 0.9,
                       caption = c(
                         "Residuals vs indices",
                         "Cook's distance",
                         "Residuals vs linear predictor",
                         "Residuals vs fitted values",
                         "Half-normal plot",
                         "Predicted vs observed",
                         "Random-effects Q-Q",
                         "Random-effects caterpillar"
                       ),
                       sub.caption = NULL,
                       ask = prod(par("mfcol")) < length(which) &&
                         dev.interactive(),
                       gg = FALSE,
                       title = NULL,
                       theme = NULL,
                       ...) {
  .check_class_mm(x)
  type <- match.arg(type)
  which <- as.integer(which)

  if (is.null(sub.caption)) {
    sub.caption <- deparse(x$call, width.cutoff = 80L)
    if (length(sub.caption) > 1L) {
      sub.caption <- paste(sub.caption, collapse = " ")
    }
  }

  if (gg) {
    .plot_gg_brsmm(
      x = x,
      which = which,
      type = type,
      nsim = nsim,
      level = level,
      caption = caption,
      sub.caption = sub.caption,
      title = title,
      theme = theme
    )
  } else {
    .plot_base_brsmm(
      x = x,
      which = which,
      type = type,
      nsim = nsim,
      level = level,
      caption = caption,
      sub.caption = sub.caption,
      ask = ask,
      ...
    )
  }

  invisible(x)
}

#' @keywords internal
.plot_base_brsmm <- function(x, which, type, nsim, level, caption,
                             sub.caption, ask, ...) {
  r <- residuals(x, type = type)
  mu_hat <- fitted(x, type = "mu")
  eta <- predict(x, type = "link")
  y_obs <- as.numeric(x$Y[, "yt"])

  n <- length(r)
  idx <- seq_len(n)

  p <- x$npar
  X <- x$model_matrices$X
  h <- tryCatch(
    {
      XtXinv <- solve(crossprod(X))
      rowSums((X %*% XtXinv) * X)
    },
    error = function(e) rep(1 / n, n)
  )
  cooks <- (r^2 * h) / (p * (1 - h)^2)

  show <- which
  nplots <- length(show)

  if (nplots > 1L) {
    ncol <- min(nplots, 2L)
    nrow <- ceiling(nplots / ncol)
    op <- par(mfrow = c(nrow, ncol), oma = c(0, 0, 2, 0))
    on.exit(par(op))
  }

  if (ask) {
    oask <- devAskNewPage(TRUE)
    on.exit(devAskNewPage(oask), add = TRUE)
  }

  if (1L %in% show) {
    plot(idx, r,
      xlab = "Index", ylab = "Residuals",
      main = caption[1L], pch = 20, col = "gray40", ...
    )
    abline(h = 0, lty = 2, col = "red")
  }

  if (2L %in% show) {
    plot(idx, cooks,
      type = "h", xlab = "Index", ylab = "Cook's distance",
      main = caption[2L], col = "gray40", ...
    )
    abline(h = 4 / n, lty = 2, col = "red")
  }

  if (3L %in% show) {
    plot(eta, r,
      xlab = "Linear predictor", ylab = "Residuals",
      main = caption[3L], pch = 20, col = "gray40", ...
    )
    abline(h = 0, lty = 2, col = "red")
  }

  if (4L %in% show) {
    plot(mu_hat, r,
      xlab = "Fitted values", ylab = "Residuals",
      main = caption[4L], pch = 20, col = "gray40", ...
    )
    abline(h = 0, lty = 2, col = "red")
  }

  if (5L %in% show) {
    .plot_halfnormal_base(x, r,
      nsim = nsim, level = level,
      caption = caption[5L], ...
    )
  }

  if (6L %in% show) {
    plot(mu_hat, y_obs,
      xlab = "Predicted", ylab = "Observed",
      main = caption[6L], pch = 20, col = "gray40", ...
    )
    abline(0, 1, lty = 2, col = "red")
  }

  # B1: Q-Q plot of random-effect modes
  if (7L %in% show) {
    re <- x$random$mode_b
    if (!is.matrix(re)) {
      re <- matrix(as.numeric(re),
        ncol = 1L,
        dimnames = list(names(re), x$random$terms[1L])
      )
    }
    for (ci in seq_len(ncol(re))) {
      term_label <- colnames(re)[ci]
      qqnorm(re[, ci],
        main = paste0(caption[7L], "\n(", term_label, ")"),
        pch = 20, col = "gray40"
      )
      qqline(re[, ci], lty = 2, col = "red")
    }
  }

  # B2: Caterpillar dotchart of random-effect modes
  if (8L %in% show) {
    re <- x$random$mode_b
    if (!is.matrix(re)) {
      re <- matrix(as.numeric(re),
        ncol = 1L,
        dimnames = list(names(re), x$random$terms[1L])
      )
    }
    for (ci in seq_len(ncol(re))) {
      term_label <- colnames(re)[ci]
      vals <- re[, ci]
      ord <- order(vals)
      dotchart(vals[ord],
        labels = rownames(re)[ord],
        main = paste0(caption[8L], "\n(", term_label, ")"),
        xlab = "Random-effect mode",
        pch = 20, col = "gray40"
      )
      abline(v = 0, lty = 2, col = "red")
    }
  }

  if (nplots > 1L) {
    mtext(sub.caption, outer = TRUE, cex = 0.8)
  }
}

#' @keywords internal
.plot_gg_brsmm <- function(x, which, type, nsim, level, caption, sub.caption, title, theme) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop(
      "Package 'ggplot2' is required for gg = TRUE. ",
      "Install it with install.packages('ggplot2').",
      call. = FALSE
    )
  }

  theme_obj <- .resolve_gg_theme(theme)
  has_global_title <- !is.null(title) && nzchar(title)

  r <- residuals(x, type = type)
  mu_hat <- fitted(x, type = "mu")
  eta <- predict(x, type = "link")
  y_obs <- as.numeric(x$Y[, "yt"])
  n <- length(r)
  idx <- seq_len(n)

  p <- x$npar
  X <- x$model_matrices$X
  h <- tryCatch(
    {
      XtXinv <- solve(crossprod(X))
      rowSums((X %*% XtXinv) * X)
    },
    error = function(e) rep(1 / n, n)
  )
  cooks <- (r^2 * h) / (p * (1 - h)^2)

  df <- data.frame(
    idx = idx,
    r = r,
    mu = mu_hat,
    eta = eta,
    cooks = cooks,
    y_obs = y_obs
  )

  plots <- list()
  panel_title <- function(i) {
    if (has_global_title && length(which) == 1L) title else caption[i]
  }
  panel_subtitle <- function() {
    if (length(which) == 1L && !is.null(sub.caption) && nzchar(sub.caption)) sub.caption else NULL
  }

  if (1L %in% which) {
    plots[[length(plots) + 1L]] <- ggplot2::ggplot(df, ggplot2::aes(x = .data$idx, y = .data$r)) +
      ggplot2::geom_point(color = "gray40", size = 1) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
      ggplot2::labs(title = panel_title(1L), subtitle = panel_subtitle(), x = "Index", y = "Residuals") +
      theme_obj
  }

  if (2L %in% which) {
    plots[[length(plots) + 1L]] <- ggplot2::ggplot(df, ggplot2::aes(x = .data$idx, y = .data$cooks)) +
      ggplot2::geom_segment(ggplot2::aes(xend = .data$idx, yend = 0), color = "gray40") +
      ggplot2::geom_hline(yintercept = 4 / n, linetype = "dashed", color = "red") +
      ggplot2::labs(title = panel_title(2L), subtitle = panel_subtitle(), x = "Index", y = "Cook's distance") +
      theme_obj
  }

  if (3L %in% which) {
    plots[[length(plots) + 1L]] <- ggplot2::ggplot(df, ggplot2::aes(x = .data$eta, y = .data$r)) +
      ggplot2::geom_point(color = "gray40", size = 1) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
      ggplot2::labs(title = panel_title(3L), subtitle = panel_subtitle(), x = "Linear predictor", y = "Residuals") +
      theme_obj
  }

  if (4L %in% which) {
    plots[[length(plots) + 1L]] <- ggplot2::ggplot(df, ggplot2::aes(x = .data$mu, y = .data$r)) +
      ggplot2::geom_point(color = "gray40", size = 1) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
      ggplot2::labs(title = panel_title(4L), subtitle = panel_subtitle(), x = "Fitted values", y = "Residuals") +
      theme_obj
  }

  if (5L %in% which) {
    plots[[length(plots) + 1L]] <- .plot_halfnormal_gg(
      r, nsim, level, panel_title(5L), panel_subtitle(), theme_obj
    )
  }

  if (6L %in% which) {
    plots[[length(plots) + 1L]] <- ggplot2::ggplot(df, ggplot2::aes(x = .data$mu, y = .data$y_obs)) +
      ggplot2::geom_point(color = "gray40", size = 1) +
      ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
      ggplot2::labs(title = panel_title(6L), subtitle = panel_subtitle(), x = "Predicted", y = "Observed") +
      theme_obj
  }

  np <- length(plots)
  if (np == 1L) {
    print(plots[[1L]])
  } else {
    ncol <- min(np, 2L)
    nrow <- ceiling(np / ncol)
    if (requireNamespace("gridExtra", quietly = TRUE)) {
      gridExtra::grid.arrange(
        grobs = plots, ncol = ncol, nrow = nrow,
        top = if (has_global_title) title else NULL
      )
    } else {
      for (p in plots) {
        print(p)
      }
    }
  }
}
