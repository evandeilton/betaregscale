# ============================================================================ #
# Cross-validation utilities
# ============================================================================ #

#' K-fold cross-validation for brs models
#'
#' @description
#' Performs repeated k-fold cross-validation for \code{\link{brs}} models.
#'
#' @param formula Model formula passed to \code{\link{brs}}.
#' @param data Data frame.
#' @param k Number of folds.
#' @param repeats Number of repeated k-fold runs.
#' @param ... Additional arguments forwarded to \code{\link{brs}}
#'   (e.g., \code{repar}, \code{link}, \code{method}).
#'
#' @return A data frame with one row per fold and columns:
#'   \code{repeat}, \code{fold}, \code{n_train}, \code{n_test},
#'   \code{log_score}, \code{rmse_yt}, \code{mae_yt}, \code{converged},
#'   and \code{error}. The object has class \code{"brs_cv"}.
#'
#' @details
#' The \code{log_score} is the mean log predictive contribution under the
#' complete likelihood contribution implied by each observation's
#' censoring type (\code{delta}).
#'
#' @references
#' Lopes, J. E. (2023). \emph{Modelos de regressao beta para dados de escala}.
#' Master's dissertation, Universidade Federal do Parana, Curitiba.
#' URI: https://hdl.handle.net/1884/86624.
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
#' @examples
#' \donttest{
#' dat <- data.frame(
#'   y = c(
#'     0, 5, 20, 50, 75, 90, 100, 30, 60, 45,
#'     10, 40, 55, 70, 85, 25, 35, 65, 80, 15
#'   ),
#'   x1 = rep(c(1, 2), 10),
#'   x2 = rep(c(0, 0, 1, 1), 5)
#' )
#' prep <- brs_prep(dat, ncuts = 100)
#' cv <- brs_cv(y ~ x1, data = prep, k = 3, repeats = 1)
#' cv
#' }
#'
#' @rdname brs_cv
#' @export
brs_cv <- function(formula,
                   data,
                   k = 5L,
                   repeats = 1L,
                   ...) {
  if (!is.data.frame(data)) {
    stop("'data' must be a data.frame.", call. = FALSE)
  }
  k <- as.integer(k)
  repeats <- as.integer(repeats)
  if (!is.finite(k) || k < 2L) {
    stop("'k' must be an integer >= 2.", call. = FALSE)
  }
  if (!is.finite(repeats) || repeats < 1L) {
    stop("'repeats' must be an integer >= 1.", call. = FALSE)
  }
  n <- nrow(data)
  if (k > n) {
    stop("'k' cannot exceed nrow(data).", call. = FALSE)
  }

  rows <- list()
  ii <- 1L

  for (r in seq_len(repeats)) {
    idx <- sample.int(n)
    fold_id <- rep(seq_len(k), length.out = n)
    fold_id <- fold_id[order(order(idx))]

    for (f in seq_len(k)) {
      test_idx <- which(fold_id == f)
      train_idx <- which(fold_id != f)
      train <- data[train_idx, , drop = FALSE]
      test <- data[test_idx, , drop = FALSE]

      fit <- tryCatch(
        brs(formula = formula, data = train, ...),
        error = identity
      )

      if (inherits(fit, "error")) {
        rows[[ii]] <- data.frame(
          `repeat` = r,
          fold = f,
          n_train = nrow(train),
          n_test = nrow(test),
          log_score = NA_real_,
          rmse_yt = NA_real_,
          mae_yt = NA_real_,
          converged = FALSE,
          error = conditionMessage(fit),
          stringsAsFactors = FALSE,
          check.names = FALSE
        )
        ii <- ii + 1L
        next
      }

      metrics <- tryCatch(
        .brs_cv_metrics(fit = fit, newdata = test),
        error = identity
      )

      if (inherits(metrics, "error")) {
        rows[[ii]] <- data.frame(
          `repeat` = r,
          fold = f,
          n_train = nrow(train),
          n_test = nrow(test),
          log_score = NA_real_,
          rmse_yt = NA_real_,
          mae_yt = NA_real_,
          converged = FALSE,
          error = conditionMessage(metrics),
          stringsAsFactors = FALSE,
          check.names = FALSE
        )
      } else {
        rows[[ii]] <- data.frame(
          `repeat` = r,
          fold = f,
          n_train = nrow(train),
          n_test = nrow(test),
          log_score = metrics$log_score,
          rmse_yt = metrics$rmse_yt,
          mae_yt = metrics$mae_yt,
          converged = isTRUE(fit$convergence == 0L),
          error = NA_character_,
          stringsAsFactors = FALSE,
          check.names = FALSE
        )
      }
      ii <- ii + 1L
    }
  }

  out <- do.call(rbind, rows)
  class(out) <- c("brs_cv", "data.frame")
  out
}

#' @keywords internal
.brs_cv_metrics <- function(fit, newdata) {
  mf <- stats::model.frame(fit$formula, data = newdata)
  Y <- .extract_response(
    mf = mf,
    data = newdata,
    ncuts = fit$ncuts,
    lim = fit$lim
  )

  mu <- predict(fit, newdata = newdata, type = "response")
  phi <- predict(fit, newdata = newdata, type = "precision")
  shp <- brs_repar(mu = mu, phi = phi, repar = fit$repar)

  delta <- as.integer(Y[, "delta"])
  left <- as.numeric(Y[, "left"])
  right <- as.numeric(Y[, "right"])
  yt <- as.numeric(Y[, "yt"])

  eps <- 1e-15
  p <- numeric(length(delta))

  i0 <- delta == 0L
  i1 <- delta == 1L
  i2 <- delta == 2L
  i3 <- delta == 3L

  if (any(i0)) {
    p[i0] <- stats::dbeta(yt[i0], shp$shape1[i0], shp$shape2[i0])
  }
  if (any(i1)) {
    p[i1] <- stats::pbeta(right[i1], shp$shape1[i1], shp$shape2[i1])
  }
  if (any(i2)) {
    p[i2] <- 1 - stats::pbeta(left[i2], shp$shape1[i2], shp$shape2[i2])
  }
  if (any(i3)) {
    p[i3] <- stats::pbeta(right[i3], shp$shape1[i3], shp$shape2[i3]) -
      stats::pbeta(left[i3], shp$shape1[i3], shp$shape2[i3])
  }
  p <- pmax(p, eps)

  list(
    log_score = mean(log(p)),
    rmse_yt = sqrt(mean((yt - mu)^2)),
    mae_yt = mean(abs(yt - mu))
  )
}
