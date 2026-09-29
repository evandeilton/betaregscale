# ============================================================================ #
# Data preparation - analyst-facing pre-processing for beta interval regression
#
# brs_prep() is the bridge between raw analyst data and betaregscale().
# It validates, classifies, and rescales observations into the (0, 1) interval
# with censoring indicators compatible with the complete likelihood.
# ============================================================================ #

#' Pre-process analyst data for beta interval regression
#'
#' @description
#' Turns analyst data into the cells \eqn{[l_i, u_i]} and censoring types
#' \eqn{\delta_i} used by \code{\link{brs}} and \code{\link{brsmm}}. Scores
#' run over \eqn{0, 1, \ldots, K} with \eqn{K =} \code{ncuts} the maximum;
#' shift a scale that starts at 1 (a Likert item 1--5 becomes 0--4,
#' \code{ncuts = 4}). Four input modes are recognised per row:
#' \enumerate{
#'   \item \strong{Score only} (\code{y}): the cell of the score and
#'     \eqn{\delta} from the score, exactly as \code{\link{brs_check}}
#'     (0 \eqn{\to} 1, \eqn{K} \eqn{\to} 2, otherwise 3; a value in
#'     \eqn{(0, 1)} \eqn{\to} 0, per observation).
#'   \item \strong{Score and \code{delta}}: the analyst's censoring type,
#'     with the cell of the score.
#'   \item \strong{Bounds only} (\code{left} and/or \code{right}, \code{y}
#'     missing): an interval, or a one-sided censoring when one bound is
#'     \code{NA}.
#'   \item \strong{Score and both bounds}: the analyst's interval.
#' }
#' Covariate columns are kept unchanged.
#'
#' @details
#' A non-missing \code{delta} always wins. Otherwise the type comes from
#' the pattern of \code{left}, \code{right} and \code{y} (here
#' \eqn{K = 10}, \code{interval = "mid"}):
#' \tabular{lllll}{
#'   \code{left} \tab \code{right} \tab \code{y} \tab cell on \eqn{(0, 1)}
#'     \tab \eqn{\delta} \cr
#'   \code{NA} \tab 3 \tab \code{NA} \tab \eqn{[\epsilon, 0.3]} (below 3) \tab 1 \cr
#'   7 \tab \code{NA} \tab \code{NA} \tab \eqn{[0.7, 1 - \epsilon]} (above 7)
#'     \tab 2 \cr
#'   2 \tab 5 \tab \code{NA} \tab \eqn{[0.2, 0.5]} \tab 3 \cr
#'   4 \tab 6 \tab 5 \tab \eqn{[0.4, 0.6]} (analyst interval) \tab 3 \cr
#'   -0.5 \tab 4 \tab \code{NA} \tab \eqn{[\epsilon, 0.4]} (reaches 0) \tab 1 \cr
#'   3 \tab 10.5 \tab \code{NA} \tab \eqn{[0.3, 1 - \epsilon]} (reaches 1)
#'     \tab 2 \cr
#'   (no columns) \tab \tab 5 \tab \eqn{[0.45, 0.55]} (Mode 1) \tab 3 \cr
#'   \code{NA} \tab \code{NA} \tab 5 \tab 0.5 (exact reading) \tab 0 \cr
#'   \code{NA} \tab \code{NA} \tab 0 \tab \eqn{[\epsilon, 0.05]} \tab 1 \cr
#'   \code{NA} \tab \code{NA} \tab 10 \tab \eqn{[0.95, 1 - \epsilon]} \tab 2
#' }
#' When the data have \code{left}/\code{right} columns, a row that gives
#' only an interior score (both bounds \code{NA}) is an exact value at the
#' cell centre, with a density contribution; without those columns the same
#' score is interval-censored (Mode 1). An analyst interval that reaches 0
#' (or 1) on the unit scale is left- (or right-) censored; one that reaches
#' both borders covers the whole scale, stays \eqn{\delta = 3} and gives a
#' warning, since such a row carries no information about the parameters.
#' A row with neither a score nor a bound is an error (a \code{delta} alone
#' defines no interval).
#'
#' Score-based rows (Modes 1 and 2) use the cells of \code{interval}, as in
#' \code{\link{brs_check}}: \eqn{[s - \mathrm{lim}, s + \mathrm{lim}]/K} for
#' \code{"mid"}, \eqn{[s, s + 1]/(K + 1)} for \code{"right"}/\code{"left"};
#' a forced \eqn{\delta = 1} keeps \eqn{u_s} and a forced \eqn{\delta = 2}
#' keeps \eqn{l_s}. Analyst bounds \eqn{L} (Modes 3 and 4) are values of the
#' latent score of the chosen direction (the scale of
#' \code{predict(type = "score")}) and map to \eqn{L/K} (\code{"mid"}),
#' \eqn{L/(K + 1)} (\code{"right"}) or \eqn{(L + 1)/(K + 1)} (\code{"left"}),
#' so that the dissertation's intervals \eqn{[s - 0.5, s + 0.5]},
#' \eqn{[s, s + 1]} and \eqn{[s - 1, s]} all give the cell of score \eqn{s}.
#' They must lie in the latent range of the direction:
#' \eqn{[-0.5, K + 0.5]} (\code{"mid"}), \eqn{[0, K + 1]} (\code{"right"})
#' or \eqn{[-1, K]} (\code{"left"}); a bound outside it is an error. Rows
#' with no score get \code{y} = the latent score of their cell centre, so
#' that \code{model.frame()} keeps them.
#'
#' Endpoints are clamped to \eqn{[\epsilon, 1 - \epsilon]},
#' \eqn{\epsilon = 10^{-5}} (why this replaces the edge-effect transformation
#' of Lopes, 2023: section 'Scale change and borders' of
#' \code{\link{brs_check}}). A \eqn{\delta = 3} row with \code{left >= right}
#' is an error (zero probability; use \eqn{\delta = 0}), and so is a cell the
#' clamp squeezes to zero width (the message tells whether \code{ncuts} is
#' too large for the cells or the analyst interval lies inside the clamp).
#' As in \code{\link{brs_check}}, values in \eqn{(0, 1)} mixed with values
#' \eqn{\ge 1} give one warning. Unusual combinations, e.g. \eqn{\delta = 1}
#' with \eqn{y \neq 0}, give a warning but are kept.
#'
#' @param data Data frame with the response columns and covariates.
#' @param y Character: name of the score column (default \code{"y"}).
#' @param delta Character: name of the censoring-type column (default
#'   \code{"delta"}); values in \code{\{0, 1, 2, 3\}} or \code{NA}.
#' @param left,right Character: names of the lower and upper bound columns
#'   (defaults \code{"left"}, \code{"right"}), on the latent score scale.
#' @param ncuts Integer \eqn{K}: the maximum score (default 100); the scale is
#'   \eqn{0, 1, \ldots, K} (\eqn{K + 1} categories). Must be at least the
#'   largest \code{y}; bounds may reach the latent range given in Details.
#' @param lim Numeric in \eqn{(0, 0.5]}: half-width of the score cell under
#'   \code{interval = "mid"} (default 0.5); see \code{\link{brs_check}}.
#' @param interval Direction of the uncertainty interval, \code{"mid"}
#'   (default), \code{"right"} or \code{"left"}; see \code{\link{brs_check}}.
#'
#' @return A data frame with columns \code{left}, \code{right} (cell on
#'   \eqn{(0, 1)}), \code{yt} (cell centre), \code{y} (score, or the filled
#'   latent score for rows without one) and \code{delta}, followed by the
#'   covariates. Attributes \code{"is_prepared"} (\code{TRUE}),
#'   \code{"ncuts"}, \code{"lim"} and \code{"interval"} are reused by
#'   \code{\link{brs}} and \code{\link{brsmm}}; an explicit different value
#'   there is ignored with a warning.
#'
#' @seealso \code{\link{brs_check}} for the cells and the automatic
#'   classification; \code{\link{brs}} for fitting the model.
#'
#' @references
#' Lopes, J. E. (2023). \emph{Modelos de regressao beta para dados de escala}.
#' Master's dissertation, Universidade Federal do Parana, Curitiba.
#' URI: https://hdl.handle.net/1884/86624.
#'
#' @examples
#' # Mode 1: score only; delta from the score (0 -> 1, K -> 2, else 3)
#' d1 <- data.frame(y = c(0, 3, 10), x = c(1.2, 0.4, -0.3))
#' brs_prep(d1, ncuts = 10)
#'
#' # Mode 2: score + analyst delta (the same score read as exact and as a cell)
#' d2 <- data.frame(y = c(4, 4), delta = c(0, 3))
#' brs_prep(d2, ncuts = 10)
#'
#' # Mode 3: only left and/or right bounds (NA pattern gives delta)
#' d3 <- data.frame(left = c(NA, 7, 2), right = c(3, NA, 5))
#' brs_prep(d3, ncuts = 10)
#'
#' # Mode 4: score with analyst bounds (used as given, divided by K)
#' d4 <- data.frame(y = 5, left = 4, right = 6)
#' brs_prep(d4, ncuts = 10)
#'
#' # An analyst interval reaching a border is one-sided censoring
#' brs_prep(data.frame(left = c(-0.5, 3), right = c(4, 10.5)), ncuts = 10)
#'
#' # A Likert item 1-5: shift it to 0-4, so that ncuts = 4 is the maximum
#' lk <- data.frame(item = c(1, 2, 5, 3, 4), x = c(0.1, 0.5, 0.9, 0.3, 0.7))
#' brs_prep(data.frame(y = lk$item - 1, x = lk$x), ncuts = 4)
#'
#' # Right-direction cells [s, s + 1] / 11; the choice is stored as attributes
#' p <- brs_prep(d1, ncuts = 10, interval = "right")
#' p[, c("left", "right", "delta")]
#' attributes(p)[c("ncuts", "lim", "interval")]
#' @rdname brs_prep
#' @export
brs_prep <- function(data, y = "y", delta = "delta",
                     left = "left", right = "right",
                     ncuts = 100L, lim = 0.5,
                     interval = c("mid", "right", "left")) {
  # -- Input validation -------------------------------------------------------
  if (!is.data.frame(data)) {
    stop("'data' must be a data.frame.", call. = FALSE)
  }
  # Direction of the cells and lim rules (stored as attributes for brs())
  interval <- match.arg(interval)
  .brs_lim_check(lim, interval, warn = TRUE)

  ncuts <- as.integer(ncuts)
  eps <- 1e-5
  K <- ncuts
  n <- nrow(data)

  # Detect which columns are present
  has_y <- y %in% names(data)
  has_delta <- delta %in% names(data)
  has_left <- left %in% names(data)
  has_right <- right %in% names(data)

  # At least one usable combination must exist
  if (!has_y && !has_left && !has_right) {
    stop(
      "At least one of '", y, "', '", left, "', or '", right,
      "' must be present in 'data'.",
      call. = FALSE
    )
  }

  # Extract raw vectors (NA if column not present)
  v_y <- if (has_y) data[[y]] else rep(NA_real_, n)
  v_delta <- if (has_delta) data[[delta]] else rep(NA_integer_, n)
  v_left <- if (has_left) data[[left]] else rep(NA_real_, n)
  v_right <- if (has_right) data[[right]] else rep(NA_real_, n)
  # A column of NA only (logical in R, e.g. `right = NA`) is a numeric NA column
  all_na_num <- function(v) if (is.logical(v) && all(is.na(v))) as.numeric(v) else v
  v_y <- all_na_num(v_y)
  v_left <- all_na_num(v_left)
  v_right <- all_na_num(v_right)

  # Validate types

  if (has_y && !is.numeric(v_y)) {
    stop("Column '", y, "' must be numeric.", call. = FALSE)
  }
  if (has_y && any(!is.na(v_y) & v_y < 0)) {
    stop("Column '", y, "' must contain non-negative values.", call. = FALSE)
  }
  if (has_delta) {
    valid_delta <- v_delta[!is.na(v_delta)]
    if (length(valid_delta) > 0 && !all(valid_delta %in% 0:3)) {
      stop(
        "Column '", delta, "' must contain values in {0, 1, 2, 3}.",
        call. = FALSE
      )
    }
  }
  if (has_left && !is.numeric(v_left)) {
    stop("Column '", left, "' must be numeric.", call. = FALSE)
  }
  if (has_right && !is.numeric(v_right)) {
    stop("Column '", right, "' must be numeric.", call. = FALSE)
  }

  # Validate left <= right where both are given
  both_lr <- !is.na(v_left) & !is.na(v_right)
  if (any(both_lr & v_left > v_right, na.rm = TRUE)) {
    bad <- which(both_lr & v_left > v_right)
    stop(
      "Observation(s) ", paste(bad, collapse = ", "),
      ": left > right. Check your data.",
      call. = FALSE
    )
  }

  # Scores live on 0..K (analyst bounds are checked on the latent scale below)
  y_obs <- v_y[!is.na(v_y)]
  if (length(y_obs) > 0 && K < max(y_obs)) {
    stop(
      "'ncuts' (", K, ") must be >= the maximum observed value (",
      max(y_obs), "). Increase 'ncuts'.",
      call. = FALSE
    )
  }
  # Analyst bounds must lie on the latent scale of the direction, checked
  # before any cell is built (a bound outside it has no cell)
  rng <- .brs_latent_range(K, interval)
  for (col in c(left, right)) {
    v <- if (identical(col, left)) v_left else v_right
    out <- which(!is.na(v) & (v < rng[1] | v > rng[2]))
    if (length(out) > 0L) {
      stop(
        "Column '", col, "': observation(s) ", paste(out, collapse = ", "),
        " outside the latent scale [", rng[1], ", ", rng[2],
        "] for interval = '", interval, "'.",
        call. = FALSE
      )
    }
  }

  # -- Per-observation processing (vectorised) ---------------------------------
  # A row needs a score or a bound; a delta alone defines no interval
  all_na_rows <- is.na(v_y) & is.na(v_left) & is.na(v_right)
  if (any(all_na_rows)) {
    bad <- which(all_na_rows)
    stop("Observation(s) ", paste(bad, collapse = ", "),
         ": all relevant columns are NA (y, left and right; a delta alone ",
         "defines no interval).", call. = FALSE)
  }

  # Same per-observation rule as brs_check(): warn once on mixed (0, 1)/scores
  .brs_warn_mixed_unit(v_y)
  # Censoring type: the analyst's delta when given, else inferred per row
  auto <- is.na(v_delta)
  out_delta <- .infer_delta(v_y, v_left, v_right, K, has_lr_cols = has_left || has_right)
  out_delta[!auto] <- as.integer(v_delta[!auto])

  # An analyst interval reaching a border of (0, 1) is one-sided censoring,
  # decided before the clamp; reaching both borders (probability ~1, no
  # information) keeps delta = 3
  both <- !is.na(v_left) & !is.na(v_right)
  open_lo <- both & auto & .brs_unit_from_latent(v_left, K, interval) <= 0
  open_hi <- both & auto & .brs_unit_from_latent(v_right, K, interval) >= 1
  out_delta[open_lo & !open_hi] <- 1L
  out_delta[open_hi & !open_lo] <- 2L
  # Bounds reaching both borders: P ~ 1 whatever the parameters (no information)
  whole <- both & .brs_unit_from_latent(v_left, K, interval) <= 0 &
    .brs_unit_from_latent(v_right, K, interval) >= 1
  if (any(whole)) {
    warning("Observation(s) ", paste(which(whole), collapse = ", "),
            ": the interval covers the whole scale (both bounds at the borders), ",
            "so the row carries no information about the parameters; consider ",
            "removing it.", call. = FALSE)
  }

  ep <- .compute_endpoints(v_y, out_delta, v_left, v_right, K, lim, eps, interval)
  # Mode 3 rows have no observed score, but `y` is the formula response and
  # model.frame() would drop NA rows, silently removing the censored
  # observations from the fit. Fill with the latent score of the interval
  # midpoint; the likelihood only uses left/right/delta.
  out_y <- ifelse(!is.na(v_y), v_y, .brs_latent_score(ep$yt, K, interval))

  # Clamp to [eps, 1-eps]
  out_left <- pmin(pmax(ep$left, eps), 1 - eps)
  out_right <- pmin(pmax(ep$right, eps), 1 - eps)
  out_yt <- pmin(pmax(ep$yt, eps), 1 - eps)

  # Zero-width delta = 3 intervals: squeezed by the clamp, or empty (D2)
  .brs_stop_zero_width(out_delta == 3L & out_left >= out_right,
                       ep$right - ep$left, ep$score_based)

  # -- Build output data.frame ------------------------------------------------
  # Identify covariate columns (everything except the input y/delta/left/right)
  input_cols <- c(y, delta, left, right)
  covar_names <- setdiff(names(data), input_cols)

  result <- data.frame(
    left = out_left,
    right = out_right,
    yt = out_yt,
    y = out_y,
    delta = out_delta,
    stringsAsFactors = FALSE
  )

  if (length(covar_names) > 0) {
    result <- cbind(result, data[, covar_names, drop = FALSE])
  }

  # Fresh sequential rownames; downstream code maps model.frame rows back by
  # name, so any later subsetting/reordering of the result stays aligned.
  rownames(result) <- NULL

  # Set attributes for downstream functions
  attr(result, "is_prepared") <- TRUE
  attr(result, "ncuts") <- ncuts
  attr(result, "lim") <- lim
  attr(result, "interval") <- interval

  # Emit consistency warnings once on the final output.
  # Check against the scores the analyst actually supplied (Mode 3 fills).
  .warn_consistency(result$delta, v_y, ncuts)

  # Informative message
  tab <- table(factor(out_delta,
    levels = 0:3,
    labels = c("exact", "left", "right", "interval")
  ))
  msg_parts <- paste0(names(tab), " = ", as.integer(tab))
  message("brs_prep: n = ", n, " | ", paste(msg_parts, collapse = ", "))

  result
}


# -- Internal helpers -------------------------------------------------------- #

#' Infer censoring types from the NA pattern (vectorised)
#'
#' Used by \code{brs_prep()} for rows without an analyst \code{delta}:
#' both bounds \eqn{\to 3}; only \code{right} (no \code{y}) \eqn{\to 1};
#' only \code{left} (no \code{y}) \eqn{\to 2}; a score \eqn{\to} the rule of
#' \code{\link{brs_check}} (0 \eqn{\to} 1, \eqn{K \to} 2, \eqn{(0, 1) \to}
#' 0, else 3), except an interior score in data that have bound columns but
#' none for this row, which is exact (0).
#' @param v_y,v_left,v_right Numeric vectors (possibly \code{NA}).
#' @param K Integer: the maximum score (\code{ncuts}).
#' @param has_lr_cols Logical: whether left/right columns exist in the data.
#' @return Integer vector of censoring types.
#' @noRd
.infer_delta <- function(v_y, v_left, v_right, K, has_lr_cols = FALSE) {
  has_y <- !is.na(v_y)
  has_l <- !is.na(v_left)
  has_r <- !is.na(v_right)
  d <- rep(NA_integer_, length(v_y))
  # Scores: the brs_check() rule
  ys <- v_y[has_y]
  d[has_y] <- ifelse(ys > 0 & ys < 1, 0L, ifelse(ys == 0, 1L, ifelse(ys == K, 2L, 3L)))
  # Interior score with empty bound columns: the analyst gave an exact value
  d[has_y & has_lr_cols & !has_l & !has_r & v_y > 0 & v_y < K] <- 0L
  # Analyst bounds: both -> interval, only right -> below it, only left -> above it
  d[has_l & has_r] <- 3L
  d[!has_y & has_r & !has_l] <- 1L
  d[!has_y & has_l & !has_r] <- 2L
  d
}


#' Endpoints of every row (vectorised)
#'
#' Analyst bounds (Modes 3/4) are latent scores of \code{interval} and go
#' through \code{.brs_unit_from_latent()}; score rows (Modes 1/2, including
#' a score with a single bound, which is ignored) use the cell of the score
#' from \code{.brs_cell()}, exactly as \code{\link{brs_check}}.
#' @param v_y,v_left,v_right Numeric vectors (possibly \code{NA}).
#' @param d Integer vector of censoring types.
#' @param K Integer: the maximum score (\code{ncuts}).
#' @param lim Numeric: half-width of the uncertainty region.
#' @param eps Numeric: small constant to avoid boundary (1e-5).
#' @param interval Interval direction.
#' @return \code{list(left, right, yt, score_based)}, unclamped.
#' @noRd
.compute_endpoints <- function(v_y, d, v_left, v_right, K, lim, eps, interval = "mid") {
  has_y <- !is.na(v_y)
  has_l <- !is.na(v_left)
  has_r <- !is.na(v_right)
  both <- has_l & has_r
  only_r <- !has_y & has_r & !has_l
  only_l <- !has_y & has_l & !has_r
  score <- !both & !only_r & !only_l
  l_u <- .brs_unit_from_latent(v_left, K, interval)
  u_u <- .brs_unit_from_latent(v_right, K, interval)
  cell <- .brs_cell(v_y, K, lim, interval)
  # delta = 0 point: a proportion stays, a score uses its cell centre
  point <- ifelse(v_y > 0 & v_y < 1, v_y, cell$mid)

  # Score rows: [eps, u_s], [l_s, 1 - eps] or [l_s, u_s]; censored rows keep
  # the cell centre of the score as yt (consistent with their endpoints)
  left <- ifelse(d == 0L, point, ifelse(d == 1L, eps, cell$left))
  right <- ifelse(d == 0L, point, ifelse(d == 2L, 1 - eps, cell$right))
  yt <- ifelse(d == 0L, point, cell$mid)
  # Analyst bounds
  left[both] <- l_u[both]
  right[both] <- u_u[both]
  yt[both] <- (l_u[both] + u_u[both]) / 2
  left[only_r] <- eps
  right[only_r] <- u_u[only_r]
  yt[only_r] <- u_u[only_r] / 2
  left[only_l] <- l_u[only_l]
  right[only_l] <- 1 - eps
  yt[only_l] <- (l_u[only_l] + 1) / 2
  list(left = left, right = right, yt = yt, score_based = score)
}


#' Issue consistency warnings for unusual delta/y combinations
#'
#' These warnings are \strong{informational}, not errors.  They alert
#' the analyst when the supplied \code{delta} does not match the
#' boundary convention:
#' \itemize{
#'   \item \eqn{\delta = 1} but \eqn{y \neq 0}: left-censored on a
#'     non-zero score.  The endpoint formula adapts to the actual y
#'     (see \code{.compute_endpoints()}).
#'   \item \eqn{\delta = 2} but \eqn{y \neq K}: right-censored on a
#'     non-maximum score.  Same adaptive endpoint logic.
#'   \item \eqn{\delta = 3} but \eqn{y} at a boundary (0 or K):
#'     interval-censored at a boundary score.
#' }
#' In Monte Carlo workflows with forced delta, these warnings are
#' expected and can be suppressed with \code{suppressWarnings()}.
#' @noRd
.warn_consistency <- function(delta, y, K) {
  # delta=1 but y != 0
  idx1 <- which(delta == 1L & !is.na(y) & y != 0)
  if (length(idx1) > 0) {
    warning(
      "Observation(s) ", paste(idx1, collapse = ", "),
      ": delta = 1 (left-censored) but y != 0.",
      call. = FALSE
    )
  }
  # delta=2 but y != K
  idx2 <- which(delta == 2L & !is.na(y) & y != K)
  if (length(idx2) > 0) {
    warning(
      "Observation(s) ", paste(idx2, collapse = ", "),
      ": delta = 2 (right-censored) but y != ", K, ".",
      call. = FALSE
    )
  }
  # delta=3 but y at boundary
  idx3 <- which(delta == 3L & !is.na(y) & (y == 0 | y == K))
  if (length(idx3) > 0) {
    warning(
      "Observation(s) ", paste(idx3, collapse = ", "),
      ": delta = 3 (interval-censored) but y is at a boundary (0 or ", K, ").",
      call. = FALSE
    )
  }
}
