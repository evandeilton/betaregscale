# ============================================================================ #
# Data preparation <U+2014> analyst-facing pre-processing for beta interval regression
#
# bs_prepare() is the bridge between raw analyst data and betaregscale().
# It validates, classifies, and rescales observations into the (0, 1) interval
# with censoring indicators compatible with the complete likelihood
# (Lopes, 2024, Eq. 2.24).
# ============================================================================ #

#' Pre-process analyst data for beta interval regression
#'
#' @description
#' Validates and transforms raw data into the format required by
#' \code{\link{betaregscale}}.
#' The analyst can supply data in several ways:
#'
#' \enumerate{
#'   \item \strong{Minimal}: only the score \code{y}. Censoring is
#'     inferred automatically (equivalent to \code{\link{check_response}}).
#'   \item \strong{Classic}: \code{y} + explicit \code{delta}. The
#'     analyst states the censoring type; interval endpoints are computed.
#'   \item \strong{Interval}: \code{left} and/or \code{right} columns
#'     (on the original scale). Censoring is inferred from the NA pattern.
#'   \item \strong{Full}: \code{y}, \code{left}, and \code{right}
#'     together. The analyst's own endpoints are rescaled directly to
#'     \eqn{(0, 1)}.
#' }
#'
#' All covariate columns are preserved unchanged in the output data frame.
#'
#' @details
#' \strong{Priority rule}: if \code{delta} is provided (non-\code{NA}),
#' it takes precedence.
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
#' All endpoints are clamped to \eqn{[\epsilon, 1 - \epsilon]} with
#' \eqn{\epsilon = 10^{-5}}.
#'
#' @param data   A \code{data.frame} containing the response and
#'   (optionally) covariates.
#' @param y      Character: name of the score column (default \code{"y"}).
#' @param delta  Character: name of the censoring indicator column
#'   (default \code{"delta"}). Values must be in \code{{0, 1, 2, 3}}.
#' @param left   Character: name of the left-endpoint column
#'   (default \code{"left"}).
#' @param right  Character: name of the right-endpoint column
#'   (default \code{"right"}).
#' @param ncuts  Integer: number of scale categories (default 100).
#' @param type   Character: interval type for interior scores when only
#'   \code{y} and \code{delta} are available.
#'   \code{"m"} = midpoint (default), \code{"l"} = left-aligned,
#'   \code{"r"} = right-aligned. See \code{\link{check_response}}.
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
#'   The output carries attributes \code{"bs_prepared"} (\code{TRUE}),
#'   \code{"ncuts"}, \code{"type"}, and \code{"lim"} so that
#'   \code{\link{betaregscale}} can detect prepared data and skip the
#'   internal \code{\link{check_response}} call.
#'
#' @seealso \code{\link{check_response}} for the automatic
#'   classification of raw scale scores;
#'   \code{\link{betaregscale}} for fitting the model.
#'
#' @examples
#' # --- Mode 1: y only (automatic classification, like check_response) ---
#' d1 <- data.frame(y = c(0, 3, 5, 7, 10), x1 = rnorm(5))
#' bs_prepare(d1, ncuts = 10)
#'
#' # --- Mode 2: y + explicit delta ---
#' d2 <- data.frame(
#'   y     = c(50, 0, 99, 50),
#'   delta = c(0, 1, 2, 3),
#'   x1    = rnorm(4)
#' )
#' bs_prepare(d2, ncuts = 100)
#'
#' # --- Mode 3: left/right with NA patterns ---
#' d3 <- data.frame(
#'   left  = c(NA, 20, 30, NA),
#'   right = c(5, NA, 45, NA),
#'   y     = c(NA, NA, NA, 50),
#'   x1    = rnorm(4)
#' )
#' bs_prepare(d3, ncuts = 100)
#'
#' # --- Mode 4: y + left + right (analyst-supplied intervals) ---
#' d4 <- data.frame(
#'   y     = c(50, 75),
#'   left  = c(48, 73),
#'   right = c(52, 77),
#'   x1    = rnorm(2)
#' )
#' bs_prepare(d4, ncuts = 100)
#'
#' # --- End-to-end workflow ---
#' \donttest{
#' set.seed(42)
#' n <- 200
#' dat <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
#' sim <- betaregscale_simulate(
#'   formula = ~ x1 + x2, data = dat,
#'   beta = c(0.2, -0.5, 0.3), phi = 1 / 5
#' )
#' prep <- bs_prepare(sim, ncuts = 100)
#' fit <- betaregscale(y ~ x1 + x2, data = prep)
#' summary(fit)
#' }
#'
#' @export
bs_prepare <- function(data, y = "y", delta = "delta",
                       left = "left", right = "right",
                       ncuts = 100L, type = "m", lim = 0.5) {
  # -- Input validation -------------------------------------------------------
  if (!is.data.frame(data)) {
    stop("'data' must be a data.frame.", call. = FALSE)
  }

  type <- match.arg(type, c("m", "l", "r"))
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
  out_left <- numeric(n)
  out_right <- numeric(n)
  out_yt <- numeric(n)
  out_delta <- integer(n)
  out_y <- numeric(n)

  for (i in seq_len(n)) {
    yi <- v_y[i]
    di <- v_delta[i]
    li <- v_left[i]
    ri <- v_right[i]

    # Check that observation is not completely NA
    all_na <- is.na(yi) && is.na(di) && is.na(li) && is.na(ri)
    if (all_na) {
      stop(
        "Observation ", i, ": all relevant columns are NA.",
        call. = FALSE
      )
    }

    # ---- Determine delta ----
    if (!is.na(di) && di %in% 0:3) {
      # Explicit delta from analyst <U+2014> use directly
      d_final <- as.integer(di)
    } else {
      # Infer delta from NA pattern
      d_final <- .infer_delta(yi, li, ri, K, has_lr_cols = has_left || has_right)
    }

    # ---- Compute endpoints ----
    endpoints <- .compute_endpoints(
      yi = yi, d = d_final, li = li, ri = ri,
      K = K, type = type, lim = lim, eps = eps
    )

    out_left[i] <- endpoints$left
    out_right[i] <- endpoints$right
    out_yt[i] <- endpoints$yt
    out_delta[i] <- d_final
    out_y[i] <- if (!is.na(yi)) yi else NA_real_
  }

  # Clamp to [eps, 1-eps]
  out_left <- pmin(pmax(out_left, eps), 1 - eps)
  out_right <- pmin(pmax(out_right, eps), 1 - eps)
  out_yt <- pmin(pmax(out_yt, eps), 1 - eps)

  # -- Consistency warnings ---------------------------------------------------
  .warn_consistency(out_delta, out_y, K)

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

  # Attach metadata
  attr(result, "bs_prepared") <- TRUE
  attr(result, "ncuts") <- K
  attr(result, "type") <- type
  attr(result, "lim") <- lim

  # Informative message
  tab <- table(factor(out_delta,
    levels = 0:3,
    labels = c("exact", "left", "right", "interval")
  ))
  msg_parts <- paste0(names(tab), " = ", as.integer(tab))
  message("bs_prepare: n = ", n, " | ", paste(msg_parts, collapse = ", "))

  result
}


# -- Internal helpers -------------------------------------------------------- #

#' Infer censoring type from NA pattern
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
    # Analyst provided left/right columns but left both NA for this row,
    # while y has a value -> exact observation (no censoring)
    return(0L)
  }
  if (has_y) {
    # y only (no left/right columns in data) -> classify like check_response:
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


#' Compute rescaled endpoints for a single observation
#' @noRd
.compute_endpoints <- function(yi, d, li, ri, K, type, lim, eps) {
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
      # Exact <U+2014> if y is already on (0,1), use it directly; otherwise rescale
      if (yi > 0 && yi < 1) {
        yt_out <- yi
      } else {
        yt_out <- yi / K
      }
      list(left = yt_out, right = yt_out, yt = yt_out)
    },
    "1" = {
      # Left-censored
      list(left = eps, right = lim / K, yt = eps)
    },
    "2" = {
      # Right-censored
      list(left = (K - lim) / K, right = 1 - eps, yt = 1 - eps)
    },
    "3" = {
      # Interval-censored <U+2014> compute via type
      yt_out <- yi / K
      switch(type,
        m = list(
          left  = (yi - lim) / K,
          right = (yi + lim) / K,
          yt    = yt_out
        ),
        l = list(
          left  = (yi - 2 * lim) / K,
          right = yi / K,
          yt    = yt_out
        ),
        r = list(
          left  = yi / K,
          right = (yi + 2 * lim) / K,
          yt    = yt_out
        )
      )
    },
    stop("Invalid delta value: ", d, call. = FALSE)
  )
}


#' Issue consistency warnings
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
