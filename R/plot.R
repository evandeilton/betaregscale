# ============================================================================ #
# Diagnostic plots for betaregscale objects
#
# Provides 6 standard diagnostic plots following the betareg tradition:
#   1. Residuals vs indices
#   2. Cook's distance
#   3. Residuals vs linear predictor
#   4. Residuals vs fitted values
#   5. Half-normal plot with simulated envelope
#   6. Predicted vs observed
#
# Both base R and ggplot2 backends are available (gg = TRUE).
# ============================================================================ #

#' Diagnostic plots for beta interval regression
#'
#' @description
#' Produces up to six diagnostic plots for a fitted
#' \code{"betaregscale"} model: residuals vs indices, Cook's
#' distance, residuals vs linear predictor, residuals vs fitted
#' values, a half-normal plot with simulated envelope, and
#' predicted vs observed.
#'
#' @param x      A fitted \code{"betaregscale"} object.
#' @param which  Integer vector selecting which plots to draw
#'   (default \code{1:4}).
#' @param type   Character: residual type passed to
#'   \code{\link{residuals.betaregscale}} (default \code{"rqr"}).
#' @param nsim   Integer: number of simulations for the half-normal
#'   envelope (default 100).
#' @param level  Numeric: confidence level for the envelope
#'   (default 0.9).
#' @param caption Character vector of panel captions.
#' @param sub.caption Subtitle; defaults to the model call.
#' @param ask    Logical: prompt before each page of plots?
#' @param gg     Logical: use ggplot2? (default \code{FALSE}).
#' @param ...    Further arguments passed to base \code{plot()}.
#'
#' @return Invisibly returns \code{x}.
#'
#' @method plot betaregscale
#' @importFrom stats qnorm fitted residuals hatvalues qqnorm quantile median
#' @importFrom graphics plot abline par mtext segments lines
#' @importFrom grDevices dev.interactive devAskNewPage adjustcolor
#' @export
plot.betaregscale <- function(x,
                              which = 1:4,
                              type = "rqr",
                              nsim = 100L,
                              level = 0.9,
                              caption = c(
                                "Residuals vs indices",
                                "Cook's distance",
                                "Residuals vs linear predictor",
                                "Residuals vs fitted values",
                                "Half-normal plot",
                                "Predicted vs observed"
                              ),
                              sub.caption = NULL,
                              ask = prod(par("mfcol")) < length(which) &&
                                dev.interactive(),
                              gg = FALSE,
                              ...) {
  .check_class(x)
  if (is.null(sub.caption)) {
    sub.caption <- deparse(x$call, width.cutoff = 80L)
    if (length(sub.caption) > 1L) sub.caption <- paste(sub.caption, collapse = " ")
  }

  if (gg) {
    .plot_gg(x,
      which = which, type = type, nsim = nsim, level = level,
      caption = caption, sub.caption = sub.caption
    )
  } else {
    .plot_base(x,
      which = which, type = type, nsim = nsim, level = level,
      caption = caption, sub.caption = sub.caption, ask = ask, ...
    )
  }

  invisible(x)
}


# -- Base R plotting -------------------------------------------------------- #

.plot_base <- function(x, which, type, nsim, level, caption, sub.caption,
                       ask, ...) {
  r <- residuals(x, type = type)
  mu_hat <- fitted(x, type = "mu")
  eta <- stats::make.link(x$link)$linkfun(mu_hat)
  n <- length(r)
  idx <- seq_len(n)

  # Cook's distance approximation: h_i * r_i^2 / (p * (1-h_i)^2)
  p <- x$npar
  V <- vcov(x)
  X <- x$model_matrices$X
  # Leverage via hat matrix: H = X (X'X)^{-1} X'
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

  # 1. Residuals vs indices
  if (1L %in% show) {
    plot(idx, r,
      xlab = "Index", ylab = "Residuals",
      main = caption[1L], pch = 20, col = "gray40", ...
    )
    abline(h = 0, lty = 2, col = "red")
  }

  # 2. Cook's distance
  if (2L %in% show) {
    plot(idx, cooks,
      type = "h", xlab = "Index",
      ylab = "Cook's distance", main = caption[2L],
      col = "gray40", ...
    )
    abline(h = 4 / n, lty = 2, col = "red")
  }

  # 3. Residuals vs linear predictor
  if (3L %in% show) {
    plot(eta, r,
      xlab = "Linear predictor", ylab = "Residuals",
      main = caption[3L], pch = 20, col = "gray40", ...
    )
    abline(h = 0, lty = 2, col = "red")
  }

  # 4. Residuals vs fitted
  if (4L %in% show) {
    plot(mu_hat, r,
      xlab = "Fitted values", ylab = "Residuals",
      main = caption[4L], pch = 20, col = "gray40", ...
    )
    abline(h = 0, lty = 2, col = "red")
  }

  # 5. Half-normal plot with envelope
  if (5L %in% show) {
    .plot_halfnormal_base(x, r,
      nsim = nsim, level = level,
      caption = caption[5L], ...
    )
  }

  # 6. Predicted vs observed
  if (6L %in% show) {
    y_obs <- x$Y[, "yt"]
    plot(mu_hat, y_obs,
      xlab = "Predicted", ylab = "Observed",
      main = caption[6L], pch = 20, col = "gray40", ...
    )
    abline(0, 1, lty = 2, col = "red")
  }

  if (nplots > 1L) mtext(sub.caption, outer = TRUE, cex = 0.8)
}


.plot_halfnormal_base <- function(x, r, nsim, level, caption, ...) {
  n <- length(r)
  ar <- sort(abs(r))
  qth <- stats::qnorm((seq_len(n) + n - 0.125) / (2 * n + 0.5))

  # Simulated envelope
  env <- matrix(NA_real_, nrow = n, ncol = nsim)
  for (j in seq_len(nsim)) {
    env[, j] <- sort(abs(stats::rnorm(n)))
  }
  alpha_half <- (1 - level) / 2
  lo <- apply(env, 1, quantile, probs = alpha_half)
  hi <- apply(env, 1, quantile, probs = 1 - alpha_half)
  me <- apply(env, 1, median)

  yl <- range(c(ar, lo, hi))
  plot(qth, ar,
    xlab = "Half-normal quantiles",
    ylab = "|Residuals|", main = caption,
    ylim = yl, pch = 20, col = "gray40", ...
  )
  lines(qth, me, lty = 2, col = "gray60")
  lines(qth, lo, lty = 3, col = "red")
  lines(qth, hi, lty = 3, col = "red")
}


# -- ggplot2 plotting ------------------------------------------------------- #

.plot_gg <- function(x, which, type, nsim, level, caption, sub.caption) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop(
      "Package 'ggplot2' is required for gg = TRUE. ",
      "Install it with install.packages('ggplot2').",
      call. = FALSE
    )
  }

  r <- residuals(x, type = type)
  mu_hat <- fitted(x, type = "mu")
  eta <- stats::make.link(x$link)$linkfun(mu_hat)
  n <- length(r)
  idx <- seq_len(n)
  y_obs <- x$Y[, "yt"]

  # Cook's distance
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
    idx = idx, r = r, mu = mu_hat, eta = eta,
    cooks = cooks, y_obs = y_obs
  )

  plots <- list()

  if (1L %in% which) {
    plots[[length(plots) + 1L]] <- ggplot2::ggplot(df, ggplot2::aes(x = .data$idx, y = .data$r)) +
      ggplot2::geom_point(color = "gray40", size = 1) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
      ggplot2::labs(title = caption[1L], x = "Index", y = "Residuals") +
      ggplot2::theme_minimal()
  }

  if (2L %in% which) {
    plots[[length(plots) + 1L]] <- ggplot2::ggplot(df, ggplot2::aes(x = .data$idx, y = .data$cooks)) +
      ggplot2::geom_segment(ggplot2::aes(xend = .data$idx, yend = 0),
        color = "gray40"
      ) +
      ggplot2::geom_hline(yintercept = 4 / n, linetype = "dashed", color = "red") +
      ggplot2::labs(title = caption[2L], x = "Index", y = "Cook's distance") +
      ggplot2::theme_minimal()
  }

  if (3L %in% which) {
    plots[[length(plots) + 1L]] <- ggplot2::ggplot(df, ggplot2::aes(x = .data$eta, y = .data$r)) +
      ggplot2::geom_point(color = "gray40", size = 1) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
      ggplot2::labs(title = caption[3L], x = "Linear predictor", y = "Residuals") +
      ggplot2::theme_minimal()
  }

  if (4L %in% which) {
    plots[[length(plots) + 1L]] <- ggplot2::ggplot(df, ggplot2::aes(x = .data$mu, y = .data$r)) +
      ggplot2::geom_point(color = "gray40", size = 1) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
      ggplot2::labs(title = caption[4L], x = "Fitted values", y = "Residuals") +
      ggplot2::theme_minimal()
  }

  if (5L %in% which) {
    plots[[length(plots) + 1L]] <- .plot_halfnormal_gg(r, nsim, level, caption[5L])
  }

  if (6L %in% which) {
    plots[[length(plots) + 1L]] <- ggplot2::ggplot(df, ggplot2::aes(x = .data$mu, y = .data$y_obs)) +
      ggplot2::geom_point(color = "gray40", size = 1) +
      ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
      ggplot2::labs(title = caption[6L], x = "Predicted", y = "Observed") +
      ggplot2::theme_minimal()
  }

  # Arrange plots in a grid
  np <- length(plots)
  if (np == 1L) {
    print(plots[[1L]])
  } else {
    ncol <- min(np, 2L)
    nrow <- ceiling(np / ncol)
    # Use gridExtra if available, otherwise print sequentially
    if (requireNamespace("gridExtra", quietly = TRUE)) {
      gridExtra::grid.arrange(
        grobs = plots, ncol = ncol, nrow = nrow,
        top = sub.caption
      )
    } else {
      for (p in plots) print(p)
    }
  }
}


.plot_halfnormal_gg <- function(r, nsim, level, caption) {
  n <- length(r)
  ar <- sort(abs(r))
  qth <- stats::qnorm((seq_len(n) + n - 0.125) / (2 * n + 0.5))

  env <- matrix(NA_real_, nrow = n, ncol = nsim)
  for (j in seq_len(nsim)) {
    env[, j] <- sort(abs(stats::rnorm(n)))
  }
  alpha_half <- (1 - level) / 2
  lo <- apply(env, 1, quantile, probs = alpha_half)
  hi <- apply(env, 1, quantile, probs = 1 - alpha_half)
  me <- apply(env, 1, median)

  df_hn <- data.frame(qth = qth, ar = ar, lo = lo, hi = hi, me = me)

  ggplot2::ggplot(df_hn, ggplot2::aes(x = .data$qth, y = .data$ar)) +
    ggplot2::geom_point(color = "gray40", size = 1) +
    ggplot2::geom_line(ggplot2::aes(y = .data$me),
      linetype = "dashed",
      color = "gray60"
    ) +
    ggplot2::geom_line(ggplot2::aes(y = .data$lo),
      linetype = "dotted",
      color = "red"
    ) +
    ggplot2::geom_line(ggplot2::aes(y = .data$hi),
      linetype = "dotted",
      color = "red"
    ) +
    ggplot2::labs(
      title = caption, x = "Half-normal quantiles",
      y = "|Residuals|"
    ) +
    ggplot2::theme_minimal()
}
