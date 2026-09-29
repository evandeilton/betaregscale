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
#'   (e.g., \code{repar}, \code{link}, \code{interval}, \code{method}).
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
#' # Synthetic NRS-11 scores: 4 groups x 3 times. Simulated, not real data.
#' set.seed(3)
#' nrs <- expand.grid(id = 1:40, time = c("6h", "12h", "24h"))
#' nrs$group <- factor(paste0("g", (nrs$id - 1) %% 4 + 1))
#' eta <- -1.3 + c(0, 0.75, 0.3)[nrs$time] + c(0, -0.1, 0.05, 0.1)[nrs$group]
#' shp <- brs_repar(mu = plogis(eta), phi = 0.3, repar = 2)
#' nrs$y <- round(10 * rbeta(nrow(nrs), shp$shape1, shp$shape2))
#'
#' # 3-fold CV of two nested models, time only (m1) and time + group (m2);
#' # log_score is the mean held-out log-likelihood contribution (higher is better)
#' set.seed(5)
#' cv1 <- brs_cv(y ~ time, data = nrs, k = 3, ncuts = 10)
#' set.seed(5)
#' cv2 <- brs_cv(y ~ time + group, data = nrs, k = 3, ncuts = 10)
#' c(m1 = mean(cv1$log_score), m2 = mean(cv2$log_score))
#' cv2
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
  # Advisory lim / mixed-value warnings of the fold fits, re-emitted once after the loop
  lim_msgs <- character(0)
  keep_lim_msg <- function(w) {
    lim_msgs <<- c(lim_msgs, conditionMessage(w))
    invokeRestart("muffleWarning")
  }

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
        withCallingHandlers(brs(formula = formula, data = train, ...),
                            brs_lim_advisory = keep_lim_msg,
                            brs_mixed_advisory = keep_lim_msg),
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
  for (m in unique(lim_msgs)) warning(m, call. = FALSE)

  out <- do.call(rbind, rows)
  class(out) <- c("brs_cv", "data.frame")
  out
}

#' @keywords internal
.brs_cv_metrics <- function(fit, newdata) {
  mf <- stats::model.frame(fit$formula, data = newdata)
  # Held-out rows coarsened exactly as the training fit was
  Y <- .extract_response(
    mf = mf,
    data = newdata,
    ncuts = fit$ncuts,
    lim = fit$lim,
    interval = .brs_interval_of(fit)
  )

  # First parameter and phi for the shapes; E[Y] for the point metrics.
  pp <- .brs_predict_params(fit, newdata = newdata)
  shp <- brs_repar(mu = pp$mu, phi = pp$phi, repar = fit$repar)
  mu <- .brs_mean(pp$mu, pp$phi, fit$repar)

  delta <- as.integer(Y[, "delta"])
  left <- as.numeric(Y[, "left"])
  right <- as.numeric(Y[, "right"])
  yt <- as.numeric(Y[, "yt"])

  # Same per-observation contributions as the compiled likelihood: log-space,
  # no probability floor (a floor made held-out tail observations score a
  # constant log(1e-15) regardless of the model).
  lp <- .brs_obs_loglik(delta, left, right, yt, shp$shape1, shp$shape2)

  list(
    log_score = mean(lp),
    rmse_yt = sqrt(mean((yt - mu)^2)),
    mae_yt = mean(abs(yt - mu))
  )
}
