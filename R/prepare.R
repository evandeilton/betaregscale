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
#' @param lim    Numeric: half-width of the uncertainty region
#'   (default 0.5). Used only when constructing intervals from \code{y}
#'   alone.
#'
#' @return A \code{data.frame} with the following columns appended or
#'   replaced:
#'   \describe{
#'     \item{\code{left}}{Lower endpoint on \eqn{(0, 1)}.}
#'     \item{\code{right}}{Upper endpoint on \eqn{(0, 1)}.}
#'     \item{\code{yt}}{Midpoint approximation on \eqn{(0, 1)}.}
#'     \item{\code{y}}{Original scale value (preserved for reference).}
#'     \item{\code{delta}}{Censoring indicator: 0 = exact, 1 = left,
#'       2 = right, 3 = interval.}
#'   }
#'   Covariate columns are preserved.
#'   The output carries attributes \code{"is_prepared"} (\code{TRUE}),
#'   \code{"ncuts"} and \code{"lim"} so that
#'   \code{\link{brs}} can detect prepared data and skip the
#'   internal \code{\link{brs_check}} call.
#'
#' @seealso \code{\link{brs_check}} for the automatic
#'   classification of raw scale scores;
#'   \code{\link{brs}} for fitting the model.
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
                     ncuts = 100L, lim = 0.5) {
  # -- Input validation -------------------------------------------------------
  if (!is.data.frame(data)) {
    stop("'data' must be a data.frame.", call. = FALSE)
  }

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

  # ncuts must be >= max value
  all_vals <- c(v_y, v_left, v_right)
  all_vals <- all_vals[!is.na(all_vals)]
  if (length(all_vals) > 0) {
    max_val <- max(all_vals)
    if (K < max_val) {
      stop(
        "'ncuts' (", K, ") must be >= the maximum observed value (",
        max_val, "). Increase 'ncuts'.",
        call. = FALSE
      )
    }
  }

  # -- Per-observation processing ---------------------------------------------
  # PERF-M02: vectorize the all-NA guard first; then use vapply for the complex
  # per-row NA-pattern logic (.infer_delta / .compute_endpoints) which cannot be
  # trivially vectorized without replicating all branch conditions.

  all_na_rows <- is.na(v_y) & is.na(v_delta) & is.na(v_left) & is.na(v_right)
  if (any(all_na_rows)) {
    bad <- which(all_na_rows)
    stop("Observation(s) ", paste(bad, collapse = ", "),
         ": all relevant columns are NA.", call. = FALSE)
  }

  has_lr_cols <- has_left || has_right

  # Determine delta for each row
  out_delta <- vapply(seq_len(n), function(i) {
    di <- v_delta[i]
    if (!is.na(di) && di %in% 0:3) {
      as.integer(di)
    } else {
      .infer_delta(v_y[i], v_left[i], v_right[i], K, has_lr_cols = has_lr_cols)
    }
  }, integer(1L))

  # Compute endpoints using the determined delta
  ep_mat <- vapply(seq_len(n), function(i) {
    ep <- .compute_endpoints(
      yi = v_y[i], d = out_delta[i], li = v_left[i], ri = v_right[i],
      K = K, lim = lim, eps = eps
    )
    c(ep$left, ep$right, ep$yt)
  }, numeric(3L))

  out_left  <- ep_mat[1L, ]
  out_right <- ep_mat[2L, ]
  out_yt    <- ep_mat[3L, ]
  out_y <- ifelse(!is.na(v_y), v_y, NA_real_)

  # Clamp to [eps, 1-eps]
  out_left <- pmin(pmax(out_left, eps), 1 - eps)
  out_right <- pmin(pmax(out_right, eps), 1 - eps)
  out_yt <- pmin(pmax(out_yt, eps), 1 - eps)

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

  # Reset rownames to sequential 1:n (critical for .extract_response()

  # which indexes by as.integer(rownames(model.frame)))
  rownames(result) <- NULL

  # Set attributes for downstream functions
  attr(result, "is_prepared") <- TRUE
  attr(result, "ncuts") <- ncuts
  attr(result, "lim") <- lim

  # Emit consistency warnings once on the final output.
  .warn_consistency(result$delta, result$y, ncuts)

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

#' Infer Censoring Type from NA Pattern
#'
#' Called by \code{brs_prep()} when the analyst does not provide an
#' explicit \code{delta} value (or \code{delta} is \code{NA}) for a
#' given observation.  The inference priority is:
#'
#' \enumerate{
#'   \item Both \code{left} and \code{right} given (no \code{y})
#'     \eqn{\to \delta = 3} (interval-censored).
#'   \item Only \code{right} given
#'     \eqn{\to \delta = 1} (left-censored: value below right).
#'   \item Only \code{left} given
#'     \eqn{\to \delta = 2} (right-censored: value above left).
#'   \item \code{y} + both \code{left} + \code{right}
#'     \eqn{\to \delta = 3} (analyst-supplied interval).
#'   \item \code{y} present, left/right columns exist but both
#'     \code{NA} \eqn{\to \delta = 0} (exact observation).
#'   \item \code{y} only: boundary rules
#'     (\eqn{y = 0 \to 1}, \eqn{y = K \to 2},
#'      \eqn{y \in (0,1) \to 0}, else \eqn{\to 3}).
#' }
#'
#' @param yi Numeric scalar: the score value (or NA).
#' @param li Numeric scalar: the left endpoint (or NA).
#' @param ri Numeric scalar: the right endpoint (or NA).
#' @param K Integer: number of scale categories (ncuts).
#' @param has_lr_cols Logical: whether left/right columns exist in the
#'   original data.frame. When TRUE and both li/ri are NA but yi has a
#'   value, the observation is treated as exact (delta = 0) because the
#'   analyst explicitly chose not to provide censoring intervals.
#' @noRd
.infer_delta <- function(yi, li, ri, K, has_lr_cols = FALSE) {
  has_y <- !is.na(yi)
  has_l <- !is.na(li)
  has_r <- !is.na(ri)

  if (has_l && has_r && !has_y) {
    # Both left and right given -> interval-censored
    return(3L)
  }

  if (!has_l && has_r && !has_y) {
    # Only right given -> left-censored ("value is below right")
    return(1L)
  }
  if (has_l && !has_r && !has_y) {
    # Only left given -> right-censored ("value is above left")
    return(2L)
  }
  if (has_y && has_l && has_r) {
    # Analyst gave y AND left AND right -> interval (use analyst endpoints)
    return(3L)
  }
  if (has_y && has_lr_cols && !has_l && !has_r) {
    # Analyst provided left/right columns but left both NA for this row.
    # If y is strictly inside (0, K), treat as exact.
    # If y is 0 or K, it must be censored (boundary rule overrides "exact").
    if (yi > 0 && yi < K) {
      return(0L)
    }
    # If y is boundary, fall through to the y-only logic below (lines 417+)
    # which correctly handles 0 -> 1 and K -> 2.
  }
  if (has_y) {
    # y only (no left/right columns in data) -> classify like brs_check:
    #   y == 0 -> left-censored
    #   y == K -> right-censored
    #   0 < y < K -> interval-censored
    #   already on (0,1) -> exact
    if (yi > 0 && yi < 1) {
      return(0L)
    }
    if (yi == 0) {
      return(1L)
    }
    if (yi == K) {
      return(2L)
    }
    return(3L)
  }

  stop(
    "Cannot infer censoring type: invalid combination of y, left, right.",
    call. = FALSE
  )
}


#' Compute Endpoints from Explicit Delta
#'
#' Called by \code{brs_prep()} inside the per-observation loop.
#' Implements four endpoint-computation modes (see brs_prep docs):
#'
#' \strong{Mode 4} (analyst supplied left + right):
#'   Rescale directly: \eqn{l = l_i / K}, \eqn{u = r_i / K}.
#'
#' \strong{Mode 3} (NA-pattern inference):
#'   Left-censored (only right given): \eqn{l = \epsilon},
#'     \eqn{u = r_i / K}.
#'   Right-censored (only left given): \eqn{l = l_i / K},
#'     \eqn{u = 1 - \epsilon}.
#'
#' \strong{Modes 1 & 2} (y-based, possibly with explicit delta):
#'   Endpoint formulas depend on the delta value and whether y is
#'   at a boundary.  The key distinction:
#'   \itemize{
#'     \item \eqn{\delta = 1, y = 0}: \eqn{u = h / K} (boundary).
#'     \item \eqn{\delta = 1, y \neq 0}: \eqn{u = (y + h) / K}
#'       (forced, observation-specific).
#'     \item \eqn{\delta = 2, y = K}: \eqn{l = (K - h) / K}
#'       (boundary).
#'     \item \eqn{\delta = 2, y \neq K}: \eqn{l = (y - h) / K}
#'       (forced, observation-specific).
#'   }
#'   See \code{\link{brs_check}} for the full formula table.
#'
#' @param yi Numeric scalar: the score value (or NA).
#' @param d  Integer scalar: the censoring type (0, 1, 2, or 3).
#' @param li Numeric scalar: analyst-supplied left endpoint (or NA).
#' @param ri Numeric scalar: analyst-supplied right endpoint (or NA).
#' @param K  Integer: number of scale categories (ncuts).
#' @param lim Numeric: half-width of the uncertainty region.
#' @param eps Numeric: small constant to avoid boundary (1e-5).
#' @return A list with elements \code{left}, \code{right}, \code{yt}.
#' @noRd
.compute_endpoints <- function(yi, d, li, ri, K, lim, eps) {
  has_y <- !is.na(yi)
  has_l <- !is.na(li)
  has_r <- !is.na(ri)

  # ---- Mode 4: analyst supplied left + right (possibly with y) ----
  if (has_l && has_r) {
    left_out <- li / K
    right_out <- ri / K
    yt_out <- (li + ri) / (2 * K)
    return(list(left = left_out, right = right_out, yt = yt_out))
  }

  # ---- Mode 3: analyst supplied only left or only right (NA pattern) ----
  if (!has_y && has_r && !has_l) {
    # Left-censored: value below right
    left_out <- eps
    right_out <- ri / K
    yt_out <- ri / (2 * K)
    return(list(left = left_out, right = right_out, yt = yt_out))
  }
  if (!has_y && has_l && !has_r) {
    # Right-censored: value above left
    left_out <- li / K
    right_out <- 1 - eps
    yt_out <- (li / K + 1) / 2
    return(list(left = left_out, right = right_out, yt = yt_out))
  }

  # ---- Modes 1 & 2: y-based (possibly with explicit delta) ----
  if (!has_y) {
    stop(
      "Internal error: y is required for delta-based endpoint computation.",
      call. = FALSE
    )
  }

  switch(as.character(d),
    "0" = {
      # Exact - if y is already on (0,1), use it directly; otherwise rescale
      if (yi > 0 && yi < 1) {
        yt_out <- yi
      } else {
        yt_out <- yi / K
      }
      list(left = yt_out, right = yt_out, yt = yt_out)
    },
    "1" = {
      # Left-censored: latent value below upper bound u
      # Boundary (y=0): u = lim/K
      # Non-boundary (forced delta=1): u = (y + lim)/K
      if (yi == 0) {
        list(left = eps, right = lim / K, yt = eps)
      } else {
        list(left = eps, right = (yi + lim) / K, yt = yi / K)
      }
    },
    "2" = {
      # Right-censored: latent value above lower bound l
      # Boundary (y=K): l = (K - lim)/K
      # Non-boundary (forced delta=2): l = (y - lim)/K
      if (yi == K) {
        list(left = (K - lim) / K, right = 1 - eps, yt = 1 - eps)
      } else {
        list(left = (yi - lim) / K, right = 1 - eps, yt = yi / K)
      }
    },
    "3" = {
      # Interval-censored — midpoint geometry
      yt_out <- yi / K
      list(
        left  = (yi - lim) / K,
        right = (yi + lim) / K,
        yt    = yt_out
      )
    },
    stop("Invalid delta value: ", d, call. = FALSE)
  )
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
