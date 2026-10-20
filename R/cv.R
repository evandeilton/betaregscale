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
#' With \code{cutpoints = "boundary"} in \code{...}, held-out rows use the
#' cells of the training fit's cutpoints. A training fold without one of the
#' scores 0, 1, \eqn{K - 1} or \eqn{K} cannot be fitted: its row has
#' \code{converged = FALSE} and the message in \code{error}.
#'
#' @seealso \code{\link{brsmm_cv}} for mixed models (folds by group).
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
  # Cells of the estimated boundary cutpoints (unchanged under the fixed map)
  Y <- .brs_apply_fit_cuts(Y, fit)

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

#' Group-level cross-validation for brsmm models
#'
#' @description
#' K-fold cross-validation of a \code{\link{brsmm}} model with folds formed by
#' group: all observations of a group (cluster) are in the same fold. Each
#' held-out group is scored by its marginal log-likelihood at the estimates of
#' the fit to the other folds.
#'
#' @details
#' The held-out log score of a fold is
#' \eqn{\sum_g \log \int \prod_i P(S = s_{gi} \mid b)\, \varphi(b; 0, \hat D)\, db}
#' over its groups, at the training-fold estimates
#' (\code{\link{brs_loglik_contrib}} with \code{newdata}). The random effects
#' of held-out groups are integrated out, so no random-effect estimates are
#' needed. The integral is computed by adaptive Gauss-Hermite quadrature with
#' \code{score_points} nodes per random-effect dimension, whatever
#' \code{int_method} the fits use; the default is 31 for one or two random
#' effects, 27 for three and 11 for four. For score data the log score is
#' the log probability of the observed scores (no Jacobian term), so it is on
#' the same scale as the held-out log-likelihood of any model that gives a
#' probability mass function on \eqn{0, \ldots, K}, evaluated on the same
#' folds (attribute \code{"folds"}).
#'
#' Folds: \code{k} folds of nearly equal numbers of groups, drawn at random
#' for each repeat. With \code{seed}, the folds are drawn after
#' \code{set.seed(seed)} and the RNG state is restored afterwards; without
#' it, the current RNG state is used. \code{folds} gives the assignment
#' directly (one repeat). Held-out rows use the cells of the training fit
#' (its estimated cutpoints under \code{cutpoints = "boundary"}). A fold
#' whose fit or scoring fails (for example, a factor level or a border score
#' missing from the training folds) has \code{log_score = NA} and the message
#' in \code{error}.
#'
#' @param formula,random Model and random-effects formulas passed to
#'   \code{\link{brsmm}}.
#' @param data Data frame.
#' @param k Number of folds (default 10).
#' @param repeats Number of repeated k-fold runs.
#' @param folds Optional fold assignment: a vector with one fold label per
#'   row of \code{data} (constant within each group) or a vector named by the
#'   group levels. Overrides \code{k}, \code{repeats} and \code{seed}.
#' @param seed Optional integer seed for the fold assignment.
#' @param score_points Nodes per dimension of the adaptive Gauss-Hermite rule
#'   that scores held-out groups; \code{NULL} for the default (see Details).
#' @param ... Further arguments for \code{\link{brsmm}} (e.g. \code{ncuts},
#'   \code{int_method}, \code{n_points}, \code{cutpoints}).
#'
#' @return A data frame of class \code{"brsmm_cv"} with one row per fold:
#'   \code{repeat}, \code{fold}, \code{n_groups_train}, \code{n_groups_test},
#'   \code{n_train}, \code{n_test}, \code{log_score} (sum over the held-out
#'   groups), \code{converged}, \code{error} and \code{warning} (messages of
#'   the fold fit). Attributes: \code{log_score} (total per repeat; \code{NA}
#'   if a fold failed), \code{group_scores} (held-out contribution of each
#'   group, one vector per repeat), \code{folds} (fold of each group, one
#'   vector per repeat) and \code{score_points}.
#'
#' @seealso \code{\link{brs_cv}}, \code{\link{brs_loglik_contrib}},
#'   \code{\link{brsmm}}
#'
#' @examples
#' \donttest{
#' set.seed(21)
#' nrs <- expand.grid(id = 1:30, time = c("6h", "12h", "24h"))
#' b0 <- rnorm(30, sd = 0.8)
#' eta <- -1.3 + c(0, 0.75, 0.3)[nrs$time] + b0[nrs$id]
#' shp <- brs_repar(mu = plogis(eta), phi = 0.2, repar = 2)
#' nrs$y <- round(10 * rbeta(nrow(nrs), shp$shape1, shp$shape2))
#'
#' cv <- brsmm_cv(y ~ time, random = ~ 1 | id, data = nrs, k = 5, seed = 1,
#'                ncuts = 10)
#' cv
#' attr(cv, "log_score")
#' }
#'
#' @export
brsmm_cv <- function(formula,
                     random = ~ 1 | id,
                     data,
                     k = 10L,
                     repeats = 1L,
                     folds = NULL,
                     seed = NULL,
                     score_points = NULL,
                     ...) {
  if (!is.data.frame(data)) {
    stop("'data' must be a data.frame.", call. = FALSE)
  }
  spec <- .brsmm_parse_random(random)
  gv <- spec$group_var
  if (!(gv %in% names(data))) {
    stop("Grouping variable '", gv, "' not found in data.", call. = FALSE)
  }
  if (anyNA(data[[gv]])) {
    stop("Grouping variable '", gv, "' contains missing values.", call. = FALSE)
  }
  grp <- factor(data[[gv]])
  lev <- levels(grp)
  G <- length(lev)
  q_re <- ncol(stats::model.matrix(spec$re_terms, data[1L, , drop = FALSE]))
  if (is.null(score_points)) {
    score_points <- c(31L, 31L, 27L, 11L)[min(q_re, 4L)]
  }
  score_points <- as.integer(score_points)
  if (length(score_points) != 1L || is.na(score_points) || score_points < 1L) {
    stop("'score_points' must be a positive integer.", call. = FALSE)
  }

  # Fold of each group, one vector per repeat
  fold_list <- if (!is.null(folds)) {
    list(.brsmm_cv_user_folds(folds, grp, data))
  } else {
    k <- as.integer(k)
    repeats <- as.integer(repeats)
    if (length(k) != 1L || is.na(k) || k < 2L) {
      stop("'k' must be an integer >= 2.", call. = FALSE)
    }
    if (length(repeats) != 1L || is.na(repeats) || repeats < 1L) {
      stop("'repeats' must be an integer >= 1.", call. = FALSE)
    }
    if (k > G) stop("'k' cannot exceed the number of groups (", G, ").", call. = FALSE)
    draw <- function() {
      lapply(seq_len(repeats), function(r) {
        stats::setNames(rep(seq_len(k), length.out = G)[sample.int(G)], lev)
      })
    }
    if (is.null(seed)) draw() else .brs_keep_seed({
      set.seed(seed)
      draw()
    })
  }

  rows <- list()
  group_scores <- list()
  totals <- numeric(length(fold_list))
  all_warn <- character(0)
  for (r in seq_along(fold_list)) {
    fg <- fold_list[[r]]
    fold_row <- fg[as.character(grp)]
    ids <- sort(unique(fg))
    scores_r <- stats::setNames(rep(NA_real_, G), lev)
    for (f in ids) {
      test_idx <- which(fold_row == f)
      train_idx <- which(fold_row != f)
      train <- .brs_keep_prep_attr(data[train_idx, , drop = FALSE], data)
      test <- .brs_keep_prep_attr(data[test_idx, , drop = FALSE], data)
      warn <- character(0)
      keep_warn <- function(w) {
        warn <<- c(warn, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
      res <- tryCatch(
        withCallingHandlers({
          fit <- brsmm(formula = formula, random = random, data = train, ...)
          ll <- brs_loglik_contrib(fit, newdata = test, int_method = "aghq",
                                   n_points = score_points)
          list(fit = fit, ll = ll)
        }, warning = keep_warn),
        error = identity
      )
      warn <- unique(warn)
      all_warn <- c(all_warn, warn)
      failed <- inherits(res, "error")
      if (!failed) {
        # A held-out group without complete rows contributes 0
        scores_r[as.character(unique(grp[test_idx]))] <- 0
        scores_r[names(res$ll)] <- as.numeric(res$ll)
      }
      rows[[length(rows) + 1L]] <- data.frame(
        `repeat` = r,
        fold = f,
        n_groups_train = length(unique(grp[train_idx])),
        n_groups_test = length(unique(grp[test_idx])),
        n_train = length(train_idx),
        n_test = length(test_idx),
        log_score = if (failed) NA_real_ else sum(res$ll),
        converged = !failed && isTRUE(res$fit$convergence == 0L),
        error = if (failed) conditionMessage(res) else NA_character_,
        warning = if (length(warn)) paste(warn, collapse = " | ") else NA_character_,
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    }
    group_scores[[r]] <- scores_r
    totals[r] <- sum(scores_r)
  }
  if (length(all_warn)) {
    warning("Fold fits gave warnings (see column 'warning'): ",
            paste(unique(all_warn), collapse = " | "), call. = FALSE)
  }
  out <- do.call(rbind, rows)
  attr(out, "log_score") <- totals
  attr(out, "group_scores") <- group_scores
  attr(out, "folds") <- fold_list
  attr(out, "score_points") <- score_points
  class(out) <- c("brsmm_cv", "data.frame")
  out
}

# Fold of each group from a user assignment (per row, or named by group level).
.brsmm_cv_user_folds <- function(folds, grp, data) {
  lev <- levels(grp)
  if (!is.null(names(folds)) && all(lev %in% names(folds)) && length(folds) != nrow(data)) {
    fg <- folds[lev]
  } else if (length(folds) == nrow(data)) {
    per <- tapply(as.character(folds), grp, function(v) unique(v))
    if (any(lengths(per) != 1L)) {
      stop("'folds' must be constant within each group.", call. = FALSE)
    }
    fg <- unlist(per)[lev]
  } else {
    stop("'folds' must have one value per row of 'data' or be named by the ",
         "group levels.", call. = FALSE)
  }
  if (anyNA(fg)) stop("'folds' contains missing values.", call. = FALSE)
  if (length(unique(fg)) < 2L) stop("'folds' must define at least 2 folds.", call. = FALSE)
  stats::setNames(as.vector(fg), lev)
}
