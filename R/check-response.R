# ============================================================================ #
# Response checking, interval construction, and censoring classification
#
# Implements the scale-to-unit interval transformation described in
# Lopes (2024, Sections 2.5--2.7).  Each observation is classified into
# one of four censoring types following the complete likelihood framework
# of Equation 2.24:
#
#   delta = 0 : uncensored (exact)   -- contribution f(y_i | theta)
#   delta = 1 : left-censored        -- contribution F(u_i | theta)
#   delta = 2 : right-censored       -- contribution 1 - F(l_i | theta)
#   delta = 3 : interval-censored    -- contribution F(u_i | theta) - F(l_i | theta)
#
# Three interval types are supported for interval-censored observations
# (Section 2.7, Figure 4):
#   "m" (midpoint): symmetric interval centered on the observed value
#   "l" (left):     interval shifted to the left
#   "r" (right):    interval shifted to the right
#
# ENDPOINT RULES (complete reference table)
# -----------------------------------------
# Let K = ncuts, eps = 1e-5, h = lim (half-width, default 0.5).
#
# delta | condition   | left (l_i)         | right (u_i)       | yt
# ------+-------------+--------------------+-------------------+-------
#   0   | y in (0,1)  | y                  | y                 | y
#   0   | y on scale  | y / K              | y / K             | y / K
#   1   | y == 0      | eps                | h / K             | eps
#   1   | y != 0      | eps                | (y + h) / K       | y / K
#   2   | y == K      | (K - h) / K        | 1 - eps           | 1 - eps
#   2   | y != K      | (y - h) / K        | 1 - eps           | y / K
#   3   | type = "m"  | (y - h) / K        | (y + h) / K       | y / K
#   3   | type = "l"  | (y - 2h) / K       | y / K             | y / K
#   3   | type = "r"  | y / K              | (y + 2h) / K      | y / K
#
# All final endpoints are clamped to [eps, 1 - eps].
#
# The distinction between "y == 0" and "y != 0" for delta = 1 (and
# analogously "y == K" vs "y != K" for delta = 2) is critical:
#
#   - When delta is inferred automatically (delta = NULL), y == 0 always
#     maps to delta = 1 and y == K always maps to delta = 2.  So the
#     boundary formulas are the only ones used.
#
#   - When delta is forced by the user (delta vector supplied), a
#     non-boundary observation (e.g., y = 50) can receive delta = 1 or
#     delta = 2.  In this case the endpoint formulas use the actual y
#     value to produce observation-specific bounds that preserve
#     covariate-driven variation.  Without this, forcing delta = 1 on
#     all observations would give every observation the same fixed
#     endpoint pair (eps, h/K), destroying all information about the
#     covariates and making regression impossible.
# ============================================================================ #

#' Transform and validate a scale-derived response variable
#'
#' @description
#' Takes a discrete (or continuous) response on the scale
#' \eqn{0, 1, \ldots, K} (where \eqn{K =} \code{ncuts}) and converts
#' it to a pair of interval endpoints on the open unit interval
#' \eqn{(0, 1)}.  Each observation is classified into one of four
#' censoring types following the complete likelihood of Lopes (2024,
#' Eq. 2.24):
#'
#' \describe{
#'   \item{\eqn{\delta = 0}}{Uncensored (exact): the observation is a
#'     continuous value already in \eqn{(0, 1)}.  The likelihood
#'     contribution is the density \eqn{f(y_i | \theta)}.
#'     Endpoints: \eqn{l_i = u_i = y_i} (or \eqn{y_i / K} when on
#'     the scale).}
#'   \item{\eqn{\delta = 1}}{Left-censored: the latent value is below
#'     some upper bound \eqn{u_i}.  The contribution is
#'     \eqn{F(u_i | \theta)}.
#'     When the observation is at the scale minimum (\eqn{y = 0}),
#'     the upper bound is \eqn{u_i = \mathrm{lim} / K}.
#'     When the user forces \eqn{\delta = 1} on a non-boundary
#'     observation (\eqn{y \neq 0}), the upper bound is
#'     \eqn{u_i = (y + \mathrm{lim}) / K}, preserving observation-
#'     specific variation.  In both cases \eqn{l_i = \epsilon}.}
#'   \item{\eqn{\delta = 2}}{Right-censored: the latent value is above
#'     some lower bound \eqn{l_i}.  The contribution is
#'     \eqn{1 - F(l_i | \theta)}.
#'     When the observation is at the scale maximum (\eqn{y = K}),
#'     the lower bound is \eqn{l_i = (K - \mathrm{lim}) / K}.
#'     When the user forces \eqn{\delta = 2} on a non-boundary
#'     observation (\eqn{y \neq K}), the lower bound is
#'     \eqn{l_i = (y - \mathrm{lim}) / K}, preserving observation-
#'     specific variation.  In both cases \eqn{u_i = 1 - \epsilon}.}
#'   \item{\eqn{\delta = 3}}{Interval-censored: the standard case for
#'     scale data.  The contribution is
#'     \eqn{F(u_i | \theta) - F(l_i | \theta)}.  Endpoints depend on
#'     the interval \code{type} (see Details).}
#' }
#'
#' @details
#' \strong{Automatic classification} (\code{delta = NULL}):
#'
#' If the entire input vector is already in \eqn{(0, 1)} (i.e., all
#' values satisfy \eqn{0 < y < 1}), all observations are treated as
#' uncensored (\eqn{\delta = 0}).
#'
#' Otherwise, for scale (integer) data:
#' \itemize{
#'   \item \eqn{y = 0}: left-censored (\eqn{\delta = 1}).
#'   \item \eqn{y = K}: right-censored (\eqn{\delta = 2}).
#'   \item \eqn{0 < y < K}: interval-censored (\eqn{\delta = 3}).
#' }
#'
#' \strong{User-supplied delta} (\code{delta} vector):
#'
#' When the \code{delta} argument is provided, the user-supplied
#' censoring indicators override the automatic boundary-based rules
#' on a per-observation basis.  This is the mechanism used by
#' \code{\link{betaregscale_simulate}} when the analyst forces a
#' specific censoring type in Monte Carlo studies.
#'
#' The endpoint formulas for each delta value are:
#'
#' \tabular{llll}{
#'   \eqn{\delta} \tab Condition \tab \eqn{l_i} (left) \tab \eqn{u_i} (right) \cr
#'   0 \tab \eqn{y \in (0, 1)} \tab \eqn{y}
#'     \tab \eqn{y} \cr
#'   0 \tab \eqn{y} on scale \tab \eqn{y / K}
#'     \tab \eqn{y / K} \cr
#'   1 \tab \eqn{y = 0} (boundary) \tab \eqn{\epsilon}
#'     \tab \eqn{\mathrm{lim} / K} \cr
#'   1 \tab \eqn{y \neq 0} (forced) \tab \eqn{\epsilon}
#'     \tab \eqn{(y + \mathrm{lim}) / K} \cr
#'   2 \tab \eqn{y = K} (boundary) \tab \eqn{(K - \mathrm{lim}) / K}
#'     \tab \eqn{1 - \epsilon} \cr
#'   2 \tab \eqn{y \neq K} (forced) \tab \eqn{(y - \mathrm{lim}) / K}
#'     \tab \eqn{1 - \epsilon} \cr
#'   3 \tab type \code{"m"} \tab \eqn{(y - \mathrm{lim}) / K}
#'     \tab \eqn{(y + \mathrm{lim}) / K} \cr
#'   3 \tab type \code{"l"} \tab \eqn{(y - 2\mathrm{lim}) / K}
#'     \tab \eqn{y / K} \cr
#'   3 \tab type \code{"r"} \tab \eqn{y / K}
#'     \tab \eqn{(y + 2\mathrm{lim}) / K} \cr
#' }
#'
#' All endpoints are clamped to \eqn{[\epsilon, 1 - \epsilon]} with
#' \eqn{\epsilon = 10^{-5}} to avoid boundary issues in the beta
#' likelihood.
#'
#' The midpoint approximation \code{yt} is computed as:
#' \itemize{
#'   \item \eqn{y_t = y} when \eqn{y \in (0, 1)} (continuous data).
#'   \item \eqn{y_t = y / K} when \eqn{y} is on the integer scale.
#' }
#' This value is used exclusively as an initialization aid for
#' starting-value computation and does not enter the likelihood.
#'
#' \strong{Interaction with the fitting pipeline}:
#'
#' This function is called internally by \code{.extract_response()}
#' when the data does \emph{not} carry the \code{"bs_prepared"}
#' attribute.  When data has been pre-processed by
#' \code{\link{bs_prepare}} or by simulation with forced delta
#' (\code{\link{betaregscale_simulate}} with \code{delta != NULL}),
#' the pre-computed columns are used directly and
#' \code{check_response()} is skipped.
#'
#' @param y      Numeric vector: the raw response.  Can be either
#'   integer scores on the scale \eqn{\{0, 1, \ldots, K\}} or
#'   continuous values already in \eqn{(0, 1)}.
#' @param type   \strong{Deprecated.}
#'   Character: interval type for \eqn{\delta = 3} observations.
#'   \code{"m"} = midpoint (default),
#'   \code{"l"} = left-aligned,
#'   \code{"r"} = right-aligned.
#'   Use \code{\link{bs_prepare}} to control interval geometry instead.
#' @param ncuts  Integer: number of scale categories \eqn{K}
#'   (default 100).  Must be \eqn{\geq \max(y)}.
#' @param lim    Numeric: half-width \eqn{h} of the uncertainty
#'   region (default 0.5).  Controls the width of the interval
#'   around each scale point.
#' @param delta  Integer vector or \code{NULL}.  If \code{NULL}
#'   (default), censoring types are inferred automatically from the
#'   boundary rules described above.
#'
#'   If provided, must have the same length as \code{y}, with every
#'   element in \code{\{0, 1, 2, 3\}}.  The supplied values override
#'   the automatic classification on a per-observation basis, and the
#'   endpoint formulas adapt to non-boundary observations as
#'   described in the table above.
#'
#'   This parameter is used internally by the simulation functions
#'   when the analyst forces a specific censoring type
#'   (e.g., \code{betaregscale_simulate(..., delta = 2)}).
#'
#' @return A numeric matrix with \eqn{n} rows and 5 columns:
#' \describe{
#'   \item{\code{left}}{Lower endpoint \eqn{l_i} on \eqn{(0, 1)},
#'     clamped to \eqn{[\epsilon, 1 - \epsilon]}.}
#'   \item{\code{right}}{Upper endpoint \eqn{u_i} on \eqn{(0, 1)},
#'     clamped to \eqn{[\epsilon, 1 - \epsilon]}.}
#'   \item{\code{yt}}{Midpoint approximation \eqn{y_t} for
#'     starting-value computation (does not enter the likelihood).}
#'   \item{\code{y}}{Original response value (preserved unchanged).}
#'   \item{\code{delta}}{Censoring indicator: 0 = exact (density),
#'     1 = left-censored \eqn{F(u)}, 2 = right-censored
#'     \eqn{1 - F(l)}, 3 = interval-censored \eqn{F(u) - F(l)}.}
#' }
#'
#' @seealso \code{\link{bs_prepare}} for the analyst-facing
#'   pre-processing function; \code{\link{betaregscale_simulate}}
#'   for simulation with forced delta.
#'
#' @examples
#' # Scale data with boundary observations
#' y <- c(0, 3, 5, 7, 9, 10)
#' check_response(y, ncuts = 10)
#'
#' # Force all observations to be exact (delta = 0)
#' check_response(y, ncuts = 10, delta = rep(0L, length(y)))
#'
#' # Force delta = 1 on non-boundary observations:
#' # endpoints use actual y values, preserving variation
#' y2 <- c(30, 60)
#' check_response(y2, ncuts = 100, delta = c(1L, 1L))
#' #  left = (eps, eps), right = (30.5/100, 60.5/100)
#'
#' @export
check_response <- function(y, type = "m", ncuts = 100L, lim = 0.5,
                           delta = NULL) {
  type <- match.arg(type, c("m", "l", "r"))
  ncuts <- as.integer(ncuts)
  n <- length(y)

  eps <- 1e-5

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

  # Integer/scale data: validate
  if (ncuts < max(y, na.rm = TRUE)) {
    warning(
      "Maximum response (", max(y, na.rm = TRUE), ") exceeds ncuts (",
      ncuts, "). Consider increasing ncuts.",
      call. = FALSE
    )
  }

  # Initialize output vectors
  y_left <- numeric(n)
  y_right <- numeric(n)
  out_delta <- integer(n)

  # Classify each observation
  for (i in seq_len(n)) {
    yi <- y[i]

    # Determine censoring type: user-supplied or automatic
    if (!is.null(delta)) {
      di <- delta[i]
    } else if (yi == 0) {
      di <- 1L
    } else if (yi == ncuts) {
      di <- 2L
    } else {
      di <- 3L
    }

    out_delta[i] <- di

    # Compute endpoints based on censoring type
    switch(as.character(di),
      "0" = {
        # Exact: use y directly (possibly on unit interval)
        if (yi > 0 && yi < 1) {
          y_left[i] <- yi
          y_right[i] <- yi
        } else {
          y_left[i] <- yi / ncuts
          y_right[i] <- yi / ncuts
        }
      },
      "1" = {
        # Left-censored: latent value below upper bound u
        # Boundary (y=0): u = lim/ncuts
        # Non-boundary (forced delta=1): u = (y + lim)/ncuts
        y_left[i] <- eps
        if (yi == 0) {
          y_right[i] <- lim / ncuts
        } else {
          y_right[i] <- (yi + lim) / ncuts
        }
      },
      "2" = {
        # Right-censored: latent value above lower bound l
        # Boundary (y=ncuts): l = (ncuts - lim)/ncuts
        # Non-boundary (forced delta=2): l = (y - lim)/ncuts
        y_right[i] <- 1 - eps
        if (yi == ncuts) {
          y_left[i] <- (ncuts - lim) / ncuts
        } else {
          y_left[i] <- (yi - lim) / ncuts
        }
      },
      "3" = {
        # Interval-censored
        switch(type,
          m = {
            y_left[i] <- (yi - lim) / ncuts
            y_right[i] <- (yi + lim) / ncuts
          },
          l = {
            y_left[i] <- (yi - 2 * lim) / ncuts
            y_right[i] <- yi / ncuts
          },
          r = {
            y_left[i] <- yi / ncuts
            y_right[i] <- (yi + 2 * lim) / ncuts
          }
        )
      }
    )
  }

  # Clamp all endpoints to (eps, 1 - eps) for numerical safety
  y_left <- pmin(pmax(y_left, eps), 1 - eps)
  y_right <- pmin(pmax(y_right, eps), 1 - eps)

  # Midpoint approximation for starting values
  yt <- ifelse(y > 0 & y < 1, y, y / ncuts)
  yt <- pmin(pmax(yt, eps), 1 - eps)

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
#'   \item If \code{data} carries the \code{"bs_prepared"} attribute
#'     (set by \code{\link{bs_prepare}} or by
#'     \code{\link{betaregscale_simulate}} with forced \code{delta}),
#'     the pre-computed columns are extracted directly.  Row subsetting
#'     uses \code{as.integer(rownames(mf))} to honour any subsetting
#'     performed by \code{model.frame()}.
#'   \item Otherwise, the raw response is extracted via
#'     \code{model.response(mf)} and passed to
#'     \code{\link{check_response}} for automatic classification.
#' }
#'
#' This two-path design ensures that:
#' \itemize{
#'   \item Data prepared by \code{bs_prepare()} (with possibly
#'     analyst-forced delta, custom left/right, or NA-pattern
#'     inference) is passed through unchanged.
#'   \item Simulated data with forced delta (which carries the
#'     \code{"bs_prepared"} attribute) preserves the forced
#'     censoring indicators and observation-specific endpoints.
#'   \item Raw data without pre-processing is classified
#'     automatically by the boundary rules of
#'     \code{check_response()}.
#' }
#'
#' @param mf Model frame (from \code{model.frame}).
#' @param data Original data frame passed by the user.  May carry
#'   the \code{"bs_prepared"} attribute.
#' @param ncuts Integer: number of scale categories \eqn{K}.
#' @param type Character: interval type (only used when falling back
#'   to \code{check_response}).
#' @param lim Numeric: uncertainty half-width (only used when falling
#'   back to \code{check_response}).
#' @return A numeric matrix with columns \code{left}, \code{right},
#'   \code{yt}, \code{y}, \code{delta} --- the same structure
#'   produced by \code{\link{check_response}}.
#' @keywords internal
#' @noRd
.extract_response <- function(mf, data, ncuts, type, lim) {
  if (isTRUE(attr(data, "bs_prepared")) &&
    all(c("left", "right", "yt", "delta") %in% names(data))) {
    rows <- as.integer(rownames(mf))
    cbind(
      left  = data[["left"]][rows],
      right = data[["right"]][rows],
      yt    = data[["yt"]][rows],
      y     = stats::model.response(mf, "numeric"),
      delta = data[["delta"]][rows]
    )
  } else {
    check_response(stats::model.response(mf, "numeric"),
      ncuts = ncuts, type = type, lim = lim
    )
  }
}
