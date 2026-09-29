# ============================================================================ #
# Response checking: score -> cell on (0, 1) -> censoring type (delta).
# Cells from .brs_cell(): mid [s - lim, s + lim] / K, right/left [s, s + 1] / (K + 1).
# delta from the score: 0 -> 1, K -> 2, else 3; values in (0, 1) -> 0 (exact).
# ============================================================================ #

# Admissible interval directions (first = default).
.brs_intervals <- c("mid", "right", "left")

#' Cell of a score on (0, 1) for the chosen interval direction
#' @param s Numeric vector of scores.
#' @param K Integer: number of scale categories.
#' @param lim Half-width, used by \code{"mid"} only.
#' @param interval One of \code{"mid"}, \code{"right"}, \code{"left"}.
#' @return \code{list(left, right, mid)} on the (0, 1) scale (unclamped).
#' @keywords internal
#' @noRd
.brs_cell <- function(s, K, lim = 0.5, interval = "mid") {
  s <- as.numeric(s)
  if (identical(interval, "mid")) {
    return(list(left = (s - lim) / K, right = (s + lim) / K, mid = s / K))
  }
  # right and left share the K + 1 equal cells that partition [0, 1]
  list(left = s / (K + 1), right = (s + 1) / (K + 1), mid = (s + 0.5) / (K + 1))
}

#' Score implied by a value on (0, 1): the inverse of the cell mapping
#' @keywords internal
#' @noRd
.brs_score_from_unit <- function(y, K, interval = "mid") {
  # mid: nearest score (cells centred on s); right/left: floor over K + 1 cells
  s <- if (identical(interval, "mid")) round(y * K) else floor(y * (K + 1))
  pmin(pmax(s, 0), K)
}

#' Latent continuous score for a mean y* on (0, 1)
#' @keywords internal
#' @noRd
.brs_latent_score <- function(ey, K, interval = "mid") {
  # left reads the same cell as right one unit lower ([s - 1, s] vs [s, s + 1])
  switch(interval,
    mid = K * ey,
    right = (K + 1) * ey,
    left = (K + 1) * ey - 1
  )
}

#' Value on (0, 1) of a latent score: the inverse of .brs_latent_score()
#' @keywords internal
#' @noRd
.brs_unit_from_latent <- function(L, K, interval = "mid") {
  # Analyst endpoints (brs_prep Modes 3/4) live on this latent scale
  switch(interval,
    mid = L / K,
    right = L / (K + 1),
    left = (L + 1) / (K + 1)
  )
}

#' Admissible range of analyst bounds on the latent score scale
#' @keywords internal
#' @noRd
.brs_latent_range <- function(K, interval = "mid") {
  # The extent of the cells of scores 0..K: [-0.5, K + 0.5], [0, K + 1], [-1, K]
  switch(interval,
    mid = c(-0.5, K + 0.5),
    right = c(0, K + 1),
    left = c(-1, K)
  )
}

#' Interval direction stored in a fit ("mid" for objects without the field)
#' @keywords internal
#' @noRd
.brs_interval_of <- function(object) {
  # Fits created before `interval` existed used the mid cells
  iv <- object$interval
  if (is.null(iv)) "mid" else iv
}

#' Validate `lim` against the interval direction
#' @param warn Emit the two advisory warnings (FALSE when brs_prep() did).
#' @param generator Add the brs_sim() sentence to the lim < 0.5 warning.
#' @keywords internal
#' @noRd
.brs_lim_check <- function(lim, interval = "mid", warn = TRUE, generator = FALSE) {
  # lim is a half-width: 0.5 makes adjacent cells touch, more would overlap;
  # values within rounding error of 0.5 (e.g. 0.7 - 0.2) count as 0.5
  tol <- sqrt(.Machine$double.eps)
  if (!is.numeric(lim) || length(lim) != 1L || !is.finite(lim) ||
    lim <= 0 || lim > 0.5 + tol) {
    stop("`lim` must be a number in (0, 0.5] (half-width of a score cell).",
      call. = FALSE
    )
  }
  if (!warn) {
    return(invisible(lim))
  }
  if (identical(interval, "mid") && lim < 0.5 - tol) {
    .brs_lim_warning(
      "`lim = ", format(lim, digits = 15), "` < 0.5: the cells cover only ",
      format(2 * lim, digits = 15),
      " of each score unit, so the coarsening is partial and score ",
      "probabilities do not sum to 1.",
      if (generator) {
        " brs_sim() rounds to the nearest score, which matches the likelihood only for lim = 0.5."
      } else {
        ""
      }
    )
  }
  if (!identical(interval, "mid") && abs(lim - 0.5) > tol) {
    .brs_lim_warning(
      "`lim` is ignored for interval = \"", interval,
      "\" (cells are [s, s + 1] / (K + 1))."
    )
  }
  invisible(lim)
}

# Advisory lim warning with its own class, so internal refits can muffle it.
.brs_lim_warning <- function(...) {
  warning(structure(
    class = c("brs_lim_advisory", "warning", "condition"),
    list(message = paste0(...), call = NULL)
  ))
}

# Evaluate `expr` without the advisory lim warnings (internal refits: the
# parent call already warned once).
.brs_quiet_lim <- function(expr) {
  withCallingHandlers(expr, brs_lim_advisory = function(w) {
    invokeRestart("muffleWarning")
  })
}

# Stop on delta = 3 intervals of zero width after the 1e-5 clamp: a cell
# squeezed by the clamp (positive raw width) or an empty interval (D2).
.brs_stop_zero_width <- function(zero, raw_width, score_based) {
  if (!any(zero, na.rm = TRUE)) {
    return(invisible(NULL))
  }
  squeezed <- which(zero & raw_width > 0)
  if (length(squeezed) > 0L && any(score_based[squeezed])) {
    stop(
      "`ncuts` too large (or `lim` too small): cells narrower than the 1e-5 ",
      "border clamp near 0 or 1 (observation(s) ",
      paste(squeezed[score_based[squeezed]], collapse = ", "), ").",
      call. = FALSE
    )
  }
  if (length(squeezed) > 0L) {
    stop(
      "Observation(s) ", paste(squeezed, collapse = ", "),
      ": interval inside the 1e-5 border clamp near 0 or 1 (zero width after ",
      "clamping).",
      call. = FALSE
    )
  }
  stop(
    "Observation(s) ", paste(which(zero), collapse = ", "),
    ": delta = 3 with left >= right (zero-probability interval). ",
    "Use delta = 0 for an exact value.",
    call. = FALSE
  )
}

#' Transform and validate a scale-derived response variable
#'
#' @description
#' Maps a score on \eqn{\{0, 1, \ldots, K\}} (\eqn{K =} \code{ncuts}) to a
#' cell \eqn{[l_s, u_s]} of \eqn{(0, 1)} and to a censoring type
#' \eqn{\delta} of the complete likelihood (dissertation, eq.
#' \code{eqn_verossimilhanca_geral}): \eqn{\delta = 0} density
#' \eqn{f(y)}, \eqn{\delta = 1} \eqn{F(u)}, \eqn{\delta = 2}
#' \eqn{1 - F(l)}, \eqn{\delta = 3} \eqn{F(u) - F(l)}. A response entirely
#' in \eqn{(0, 1)} is exact (\eqn{\delta = 0}).
#'
#' @section Interval direction:
#' \code{interval} is the direction of the uncertainty interval around the
#' score (dissertation, "Mapeamento de intervalos para beta":
#' \eqn{m = [s - 0.5, s + 0.5]}, \eqn{r = [s, s + 1]}, \eqn{l = [s - 1, s]}):
#' \tabular{lll}{
#'   \code{interval} \tab cell of score \eqn{s} \tab latent score \cr
#'   \code{"mid"} \tab \eqn{[s - \mathrm{lim}, s + \mathrm{lim}] / K}
#'     \tab \eqn{K y^*} \cr
#'   \code{"right"} \tab \eqn{[s, s + 1] / (K + 1)} \tab \eqn{(K + 1) y^*} \cr
#'   \code{"left"} \tab \eqn{[s, s + 1] / (K + 1)} \tab \eqn{(K + 1) y^* - 1}
#' }
#' The \eqn{K + 1} cells of \code{"right"} and \code{"left"} are equal and
#' partition \eqn{[0, 1]}. This normalisation is a package choice that
#' differs from the dissertation, which divides \eqn{r} and \eqn{l} by
#' \eqn{K} (there \eqn{r} and \eqn{l} differ by \eqn{1/K}, and chapter 4
#' reports opposite intercept biases for them); here \code{"right"} and
#' \code{"left"} give the same likelihood and coefficients, and differ only
#' in how a fitted value is read back on the score scale (one unit), so
#' that opposite-bias signature disappears by construction. The three
#' modes are different coarsening models of the same scores: their
#' log-likelihoods are not comparable and \code{anova()} refuses to compare
#' them. \code{lim} applies to \code{"mid"} only.
#'
#' The censoring type comes from the score, before any clamping:
#' \eqn{s = 0 \to \delta = 1} with \eqn{u = u_0}, \eqn{s = K \to \delta = 2}
#' with \eqn{l = l_K}, otherwise \eqn{\delta = 3}.
#'
#' @details
#' If every value is in \eqn{(0, 1)} and \code{delta} is \code{NULL}, all
#' observations are exact. Otherwise each score gets its cell and, with
#' \code{delta = NULL}, the type above. A user-supplied \code{delta} (the
#' mechanism \code{\link{brs_sim}} uses in Monte Carlo studies) forces the
#' type per observation and keeps the cell endpoints of the score:
#' \tabular{lll}{
#'   \eqn{\delta} \tab \eqn{l_i} \tab \eqn{u_i} \cr
#'   0 \tab cell centre (or \eqn{y} when in \eqn{(0, 1)}) \tab same \cr
#'   1 \tab \eqn{\epsilon} \tab \eqn{u_s} \cr
#'   2 \tab \eqn{l_s} \tab \eqn{1 - \epsilon} \cr
#'   3 \tab \eqn{l_s} \tab \eqn{u_s}
#' }
#' Under \code{"mid"} with \code{lim = 0.5} this is \eqn{u_0 = 0.5 / K},
#' \eqn{l_K = (K - 0.5) / K} and \eqn{[l_s, u_s] = [(s - 0.5) / K,
#' (s + 0.5) / K]}. Scores outside \eqn{[0, K]} are an error (as in
#' \code{\link{brs_prep}}), and so is a \eqn{\delta = 3} observation with
#' \eqn{l_i = u_i} (zero-probability interval).
#'
#' All endpoints are clamped to \eqn{[\epsilon, 1 - \epsilon]},
#' \eqn{\epsilon = 10^{-5}}. \code{yt} is the cell centre (\eqn{s / K}
#' under \code{"mid"}, \eqn{(s + 0.5) / (K + 1)} otherwise; \eqn{y} itself
#' for exact values): it is the density argument for \eqn{\delta = 0} and
#' a point summary elsewhere; censored contributions use only
#' \code{left}/\code{right}.
#'
#' \strong{Interaction with the fitting pipeline}:
#'
#' This function is called internally by \code{.extract_response()}
#'   when the data does \emph{not} carry the \code{"is_prepared"}
#'   attribute.  If data has already been processed by
#'   \code{\link{brs_prep}} or by simulation with forced delta
#' (\code{\link{brs_sim}} with \code{delta != NULL}),
#' the pre-computed columns are used directly and
#' \code{brs_check()} is skipped.
#'
#' @param y      Numeric vector: the raw response. Can be either
#'   integer scores on the scale \eqn{\{0, 1, \ldots, K\}} or
#'   continuous values already in \eqn{(0, 1)}.
#' @param ncuts  Integer: number of scale categories \eqn{K}
#'   (default 100). Must be \eqn{\geq \max(y)}.
#' @param lim    Numeric in \eqn{(0, 0.5]}: half-width of the cell under
#'   \code{interval = "mid"} (default 0.5, adjacent cells touch). Values
#'   below 0.5 give a partial coarsening (warning); ignored, with a
#'   warning, for \code{"right"}/\code{"left"}.
#' @param delta  Integer vector or \code{NULL}. If \code{NULL}
#'   (default), censoring types are derived from the scores. If
#'   provided, must have the same length as \code{y} with elements in
#'   \code{\{0, 1, 2, 3\}}; it overrides the type per observation (see
#'   Details).
#' @param interval Direction of the uncertainty interval: \code{"mid"}
#'   (default), \code{"right"} or \code{"left"}; see the section
#'   'Interval direction'.
#'
#' @return A numeric matrix with \eqn{n} rows and 5 columns:
#' \describe{
#'   \item{\code{left}}{Lower endpoint \eqn{l_i} on \eqn{(0, 1)},
#'     clamped to \eqn{[\epsilon, 1 - \epsilon]}.}
#'   \item{\code{right}}{Upper endpoint \eqn{u_i} on \eqn{(0, 1)},
#'     clamped to \eqn{[\epsilon, 1 - \epsilon]}.}
#'   \item{\code{yt}}{Midpoint approximation \eqn{y_t} for
#'     starting-value computation. Also enters the likelihood
#'     directly as the density argument for exact observations
#'     (\eqn{\delta = 0}); for censored observations only
#'     \code{left}/\code{right} enter the likelihood.}
#'   \item{\code{y}}{Original response value (preserved unchanged).}
#'   \item{\code{delta}}{Censoring indicator: 0 = exact (density),
#'     1 = left-censored \eqn{F(u)}, 2 = right-censored
#'     \eqn{1 - F(l)}, 3 = interval-censored \eqn{F(u) - F(l)}.}
#' }
#'
#' @seealso \code{\link{brs_prep}} for the analyst-facing
#'   pre-processing function; \code{\link{brs_sim}}
#'   for simulation with forced delta.
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
#' # Scale data with boundary observations
#' y <- c(0, 3, 5, 7, 9, 10)
#' brs_check(y, ncuts = 10)
#'
#' # Right-direction intervals: cells [s, s + 1] / 11
#' brs_check(y, ncuts = 10, interval = "right")
#'
#' # Force all observations to be exact (delta = 0)
#' brs_check(y, ncuts = 10, delta = rep(0L, length(y)))
#'
#' # Force delta = 1 on non-boundary observations: u = (y + 0.5) / K
#' y2 <- c(30, 60)
#' brs_check(y2, ncuts = 100, delta = c(1L, 1L))
#' @rdname brs_check
#' @export
brs_check <- function(y, ncuts = 100L, lim = 0.5, delta = NULL,
                      interval = c("mid", "right", "left")) {
  interval <- match.arg(interval)
  n <- length(y)

  # Validate user-supplied delta
  if (!is.null(delta)) {
    delta <- as.integer(delta)
    if (length(delta) != n) {
      stop("'delta' must have the same length as 'y' (", n, ").",
        call. = FALSE
      )
    }
    if (!all(delta %in% 0:3)) {
      stop("'delta' must contain only values in {0, 1, 2, 3}.",
        call. = FALSE
      )
    }
  }
  # lim > 0.5 errors; lim < 0.5 (mid) or a non-default lim (right/left) warns
  .brs_lim_check(lim, interval, warn = TRUE)
  .brs_check_core(y, ncuts = ncuts, lim = lim, delta = delta, interval = interval)
}

# brs_check() without argument validation: the fitting path calls it after
# .validate_brs_common_args() already checked (and warned about) lim.
.brs_check_core <- function(y, ncuts, lim, delta, interval) {
  ncuts <- as.integer(ncuts)
  n <- length(y)
  eps <- 1e-5

  # Detect if data is already on (0, 1)
  is_unit <- all(y > 0 & y < 1)

  if (is_unit && is.null(delta)) {
    # Continuous data in (0, 1): treat as uncensored (delta = 0)
    message(
      "Response is already on the unit interval (0, 1); ",
      "treating as uncensored (exact) observations."
    )
    yt <- pmin(pmax(y, eps), 1 - eps)
    out <- cbind(
      left  = yt,
      right = yt,
      yt    = yt,
      y     = y,
      delta = rep(0L, n)
    )
    return(out)
  }

  # Scores outside 0..K have no cell (their clamped interval has zero width and
  # would contribute LOG_PENALTY): stop, as brs_prep() does
  if (any(y < 0, na.rm = TRUE)) {
    stop("Scores must be non-negative (minimum found: ", min(y, na.rm = TRUE), ").",
      call. = FALSE
    )
  }
  if (ncuts < max(y, na.rm = TRUE)) {
    stop(
      "Maximum response (", max(y, na.rm = TRUE), ") exceeds ncuts (", ncuts,
      "). Increase `ncuts` to at least ", ceiling(max(y, na.rm = TRUE)), ".",
      call. = FALSE
    )
  }

  # Censoring type from the score (before any clamp); a forced delta wins
  out_delta <- if (!is.null(delta)) {
    as.integer(delta)
  } else {
    ifelse(y == 0L, 1L, ifelse(y == ncuts, 2L, 3L))
  }

  # Cell of each score for the chosen direction; exact values in (0, 1) stay
  cell <- .brs_cell(y, ncuts, lim, interval)
  is_unit_obs <- (y > 0 & y < 1)
  pts_exact <- ifelse(is_unit_obs, y, cell$mid)

  # delta 1: [eps, u_s]; delta 2: [l_s, 1 - eps]; delta 3: [l_s, u_s]
  d <- out_delta
  y_left  <- ifelse(d == 0L, pts_exact, ifelse(d == 1L, eps, cell$left))
  y_right <- ifelse(d == 0L, pts_exact, ifelse(d == 2L, 1 - eps, cell$right))

  # Clamp all endpoints to (eps, 1 - eps) for numerical safety
  raw_width <- y_right - y_left
  y_left <- pmin(pmax(y_left, eps), 1 - eps)
  y_right <- pmin(pmax(y_right, eps), 1 - eps)

  # Zero-width delta = 3 cells: squeezed by the clamp, or empty (D2)
  .brs_stop_zero_width(d == 3L & y_left >= y_right, raw_width, rep(TRUE, n))

  # Cell centre (y itself for exact values), clamped like the endpoints
  yt <- pmin(pmax(pts_exact, eps), 1 - eps)

  cbind(left = y_left, right = y_right, yt = yt, y = y, delta = out_delta)
}


#' Extract the response matrix from user-supplied data
#'
#' This is the \strong{single gateway} through which every fitting,
#' log-likelihood, and starting-value function obtains the five-column
#' response matrix (\code{left}, \code{right}, \code{yt}, \code{y},
#' \code{delta}).
#'
#' \strong{Decision logic}:
#' \enumerate{
#'   \item If \code{data} carries the \code{"is_prepared"} attribute
#'     (set by \code{\link{brs_prep}} or by
#'     \code{\link{brs_sim}} with forced \code{delta}),
#'     the pre-computed columns are extracted directly. Row subsetting
#'     first tries numeric row indices from \code{rownames(mf)} and
#'     then falls back to name matching against \code{rownames(data)}
#'     when non-numeric row names are used.
#'   \item Otherwise, the raw response is extracted via
#'     \code{model.response(mf)} and passed to
#'     \code{\link{brs_check}} for automatic classification.
#' }
#'
#' This two-path design ensures that:
#' \itemize{
#'   \item Data prepared by \code{brs_prep()} (with possibly
#'     forced intervals or custom endpoints) carries the
#'     \code{"is_prepared"} attribute.
#'   \item \code{brs_sim()} also attaches the \code{"is_prepared"}
#'     attribute when \code{delta} is forced, preserving
#'     censoring indicators and observation-specific endpoints.
#'   \item Raw data without pre-processing is classified
#'     automatically by the boundary rules of
#'     \code{brs_check()}.
#' }
#'
#' @param mf Model frame (from \code{model.frame}).
#' @param data Original data frame passed by the user.  May carry
#'   the \code{"is_prepared"} attribute.
#' @param ncuts Integer: number of scale categories \eqn{K}.
#' @param lim Numeric: uncertainty half-width (only used when falling
#'   back to \code{brs_check}).
#' @param interval Interval direction (only used when falling back).
#' @return A numeric matrix with columns \code{left}, \code{right},
#'   \code{yt}, \code{y}, \code{delta} --- the same structure
#'   produced by \code{\link{brs_check}}.
#' @keywords internal
#' @noRd
.extract_response <- function(mf, data, ncuts, lim, interval = "mid") {
  if (isTRUE(attr(data, "is_prepared")) &&
    all(c("left", "right", "yt", "delta") %in% names(data))) {
    # Row names are labels, not positions: after `data[-10, ]` the row named
    # "11" is the 10th row. Always map by name (model.frame keeps the row
    # names of `data`, and data.frame row names are unique).
    rows <- match(rownames(mf), rownames(data))
    if (anyNA(rows)) {
      stop(
        "Unable to align prepared rows between `model.frame` and `data`.",
        call. = FALSE
      )
    }
    cbind(
      left  = data[["left"]][rows],
      right = data[["right"]][rows],
      yt    = data[["yt"]][rows],
      y     = stats::model.response(mf, "numeric"),
      delta = data[["delta"]][rows]
    )
  } else {
    # Raw scores: lim/interval were validated upstream, so no repeated warnings
    .brs_check_core(stats::model.response(mf, "numeric"),
      ncuts = ncuts, lim = lim, delta = NULL, interval = interval
    )
  }
}
