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
#' Validates and transforms raw data into the format required by
#' \code{\link{brs}}.
#' The analyst can supply data in several ways:
#'
#' \enumerate{
#'   \item \strong{Minimal (Mode 1)}: only the score \code{y}.
#'     Censoring is inferred automatically:
#'     \eqn{y = 0 \to \delta = 1}, \eqn{y = K \to \delta = 2},
#'     \eqn{0 < y < K \to \delta = 3},
#'     \eqn{y \in (0, 1) \to \delta = 0}.
#'   \item \strong{Classic (Mode 2)}: \code{y} + explicit
#'     \code{delta}. The analyst declares the censoring type;
#'     interval endpoints are computed using the actual \code{y}
#'     value.
#'   \item \strong{Interval (Mode 3)}: \code{left} and/or
#'     \code{right} columns (on the original scale). Censoring is
#'     inferred from the NA pattern.
#'   \item \strong{Full (Mode 4)}: \code{y}, \code{left}, and
#'     \code{right} together. The analyst's own endpoints are
#'     rescaled directly to \eqn{(0, 1)}.
#' }
#'
#' All covariate columns are preserved unchanged in the output.
#'
#' @details
#' \strong{Priority rule}: if \code{delta} is provided (non-\code{NA}),
#' it takes precedence over all automatic classification rules.
#' When \code{delta} is \code{NA}, the function infers the censoring type
#' from the pattern of \code{left}, \code{right}, and \code{y}:
#'
#' \tabular{llllll}{
#'   \code{left} \tab \code{right} \tab \code{y} \tab \code{delta}
#'   \tab Interpretation \tab Inferred \eqn{\delta} \cr
#'   \code{NA}   \tab  5  \tab \code{NA} \tab \code{NA}
#'   \tab Left-censored (below 5) \tab 1 \cr
#'   20          \tab \code{NA} \tab \code{NA} \tab \code{NA}
#'   \tab Right-censored (above 20) \tab 2 \cr
#'   30          \tab 45  \tab \code{NA} \tab \code{NA}
#'   \tab Interval-censored [30, 45] \tab 3 \cr
#'   \code{NA}   \tab \code{NA} \tab 50 \tab \code{NA}
#'   \tab Exact observation \tab 0 \cr
#'   \code{NA}   \tab \code{NA} \tab 50 \tab 3
#'   \tab Analyst says interval \tab 3 \cr
#'   \code{NA}   \tab \code{NA} \tab 0  \tab 1
#'   \tab Analyst says left-censored \tab 1 \cr
#'   \code{NA}   \tab \code{NA} \tab 99 \tab 2
#'   \tab Analyst says right-censored \tab 2 \cr
#' }
#'
#' When \code{y}, \code{left}, and \code{right} are all present for the
#' same observation, the analyst's \code{left}/\code{right} values are
#' used directly (rescaled by \eqn{K =} \code{ncuts}) and \code{delta}
#' is set to 3 (interval-censored) unless the analyst supplied
#' \code{delta} explicitly.
#'
#' \strong{Endpoint formulas for Mode 2 (y + explicit delta)}:
#'
#' When the analyst supplies \code{delta} explicitly, the endpoint
#' computation uses the actual \code{y} value to produce
#' observation-specific bounds.  This is the same logic used by
#' \code{\link{brs_check}} with a user-supplied \code{delta}
#' vector:
#'
#' \tabular{llll}{
#'   \eqn{\delta} \tab Condition \tab \eqn{l_i} (left)
#'     \tab \eqn{u_i} (right) \cr
#'   0 \tab (any) \tab \eqn{y / K} \tab \eqn{y / K} \cr
#'   1 \tab \eqn{y = 0} \tab \eqn{\epsilon}
#'     \tab \eqn{\mathrm{lim} / K} \cr
#'   1 \tab \eqn{y \neq 0} \tab \eqn{\epsilon}
#'     \tab \eqn{(y + \mathrm{lim}) / K} \cr
#'   2 \tab \eqn{y = K} \tab \eqn{(K - \mathrm{lim}) / K}
#'     \tab \eqn{1 - \epsilon} \cr
#'   2 \tab \eqn{y \neq K} \tab \eqn{(y - \mathrm{lim}) / K}
#'     \tab \eqn{1 - \epsilon} \cr
#'   3 \tab type \code{"m"} \tab \eqn{(y - \mathrm{lim}) / K}
#'     \tab \eqn{(y + \mathrm{lim}) / K} \cr
#' }
#'
#' \strong{Consistency warnings}: when the analyst supplies \code{delta}
#' values that are unusual for the given \code{y} (e.g.,
#' \eqn{\delta = 1} but \eqn{y \neq 0}), the function emits a warning
#' but proceeds.  This is by design for Monte Carlo workflows where
#' forced delta on non-boundary observations is intentional.
#'
#' All endpoints are clamped to \eqn{[\epsilon, 1 - \epsilon]} with
#' \eqn{\epsilon = 10^{-5}}.
#'
#' @param data   Data frame containing the response variable and
#'   covariates.
#' @param y      Character: name of the score column (default \code{"y"}).
#' @param delta  Character: name of the censoring indicator column
#'   (default \code{"delta"}). Values must be in \code{{0, 1, 2, 3}}.
#' @param left   Character: name of the left-endpoint column
#'   (default \code{"left"}).
#' @param right  Character: name of the right-endpoint column
#'   (default \code{"right"}).
#' @param ncuts  Integer: number of scale categories (default 100).
#' @param lim    Numeric in \eqn{(0, 0.5]}: half-width of the score cell
#'   under \code{interval = "mid"} (default 0.5). Used only when
#'   constructing intervals from \code{y} alone; see \code{\link{brs_check}}.
#' @param interval Direction of the uncertainty interval, \code{"mid"}
#'   (default), \code{"right"} or \code{"left"}; see the section 'Interval
#'   direction' of \code{\link{brs_check}}. Mode 1/2 cells follow it.
#'   Analyst endpoints \eqn{L} (Modes 3/4) are bounds on the latent score of
#'   that direction (the scale of \code{predict(type = "score")}), must lie in
#'   \eqn{[-0.5, K + 0.5]}, \eqn{[0, K + 1]} or \eqn{[-1, K]}, and map to
#'   \eqn{L / K}, \eqn{L / (K + 1)} or \eqn{(L + 1) / (K + 1)}; so
#'   \eqn{[s - 0.5, s + 0.5]}, \eqn{[s, s + 1]} and \eqn{[s - 1, s]} all give
#'   the cell of score \eqn{s}. Without an analyst \code{delta}, an interval
#'   reaching 0 (or 1) on the unit scale is left- (or right-) censored.
#'
#' @return A \code{data.frame} with the following columns appended or
#'   replaced:
#'   \describe{
#'     \item{\code{left}}{Lower endpoint on \eqn{(0, 1)}.}
#'     \item{\code{right}}{Upper endpoint on \eqn{(0, 1)}.}
#'     \item{\code{yt}}{Cell centre (point summary) on \eqn{(0, 1)}.}
#'     \item{\code{y}}{Original scale value (preserved for reference).}
#'     \item{\code{delta}}{Censoring indicator: 0 = exact, 1 = left,
#'       2 = right, 3 = interval.}
#'   }
#'   Covariate columns are preserved.
#'   The output carries attributes \code{"is_prepared"} (\code{TRUE}),
#'   \code{"ncuts"}, \code{"lim"} and \code{"interval"}, which
#'   \code{\link{brs}} and \code{\link{brsmm}} reuse (an explicit
#'   different value is ignored with a warning).
#'
#' @seealso \code{\link{brs_check}} for the automatic
#'   classification of raw scale scores;
#'   \code{\link{brs}} for fitting the model.
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
#' # --- Mode 1: y only (automatic classification, like brs_check) ---
#' d1 <- data.frame(y = c(0, 3, 5, 7, 10), x1 = rnorm(5))
#' brs_prep(d1, ncuts = 10)
#'
#' # --- Mode 2: y + explicit delta ---
#' d2 <- data.frame(
#'   y = d1$y,
#'   delta = c(0, 3, 3, 3, 0), # Force interval-censoring for 3,5,7
#'   x1 = d1$x1
#' )
#' brs_prep(d2, ncuts = 100)
#'
#' # --- Mode 3: left/right with NA patterns ---
#' d3 <- data.frame(
#'   left = c(NA, 20, 30, NA),
#'   right = c(5, NA, 45, NA),
#'   y = c(NA, NA, NA, 50),
#'   x1 = d1$x1[1:4]
#' )
#' brs_prep(d3, ncuts = 100)
#'
#' # --- Mode 4: y + left + right (analyst-supplied intervals) ---
#' d4 <- data.frame(
#'   y = c(50, 75),
#'   left = c(48, 73),
#'   right = c(52, 77),
#'   x1 = rnorm(2)
#' )
#' brs_prep(d4, ncuts = 100)
#'
#' # --- Fitting after prep ---
#' \donttest{
#' dat5 <- data.frame(
#'   y = c(
#'     0, 5, 20, 50, 75, 90, 100, 30, 60, 45,
#'     10, 40, 55, 70, 85, 25, 35, 65, 80, 15
#'   ),
#'   x1 = rep(c(1, 2), 10)
#' )
#' prep5 <- brs_prep(dat5, ncuts = 100)
#' fit5 <- brs(y ~ x1, data = prep5)
#' summary(fit5)
#' }
#'
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
#' @param K Integer: number of scale categories (ncuts).
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
#' @param K Integer: number of scale categories (ncuts).
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
