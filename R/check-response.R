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
#'     contribution is the density \eqn{f(y_i | \theta)}.}
#'   \item{\eqn{\delta = 1}}{Left-censored: the observation equals the
#'     scale minimum (\eqn{y = 0}).  We only know that the latent
#'     value is below some upper bound.  The contribution is
#'     \eqn{F(u_i | \theta)}.}
#'   \item{\eqn{\delta = 2}}{Right-censored: the observation equals the
#'     scale maximum (\eqn{y = K}).  The contribution is
#'     \eqn{1 - F(l_i | \theta)}.}
#'   \item{\eqn{\delta = 3}}{Interval-censored: the standard case for
#'     scale data.  The contribution is
#'     \eqn{F(u_i | \theta) - F(l_i | \theta)}.}
#' }
#'
#' @details
#' If the input is already in \eqn{(0, 1)} (i.e., all values satisfy
#' \eqn{0 < y < 1}), observations are treated as uncensored (\eqn{\delta = 0}).
#'
#' For scale (integer) data:
#' \itemize{
#'   \item \eqn{y = 0}: left-censored. Upper bound is
#'     \eqn{u = \mathrm{lim} / K}.
#'   \item \eqn{y = K}: right-censored. Lower bound is
#'     \eqn{l = (K - \mathrm{lim}) / K}.
#'   \item \eqn{0 < y < K}: interval-censored with endpoints determined
#'     by the interval \code{type}.
#' }
#'
#' Three interval types (Section 2.7, Figure 4 of Lopes, 2024):
#' \describe{
#'   \item{\code{"m"} (midpoint)}{Symmetric: \eqn{l = (y - \mathrm{lim})/K},
#'     \eqn{u = (y + \mathrm{lim})/K}.}
#'   \item{\code{"l"} (left)}{Shifted left: \eqn{l = (y - 2\mathrm{lim})/K},
#'     \eqn{u = y/K}.}
#'   \item{\code{"r"} (right)}{Shifted right: \eqn{l = y/K},
#'     \eqn{u = (y + 2\mathrm{lim})/K}.}
#' }
#'
#' All endpoints are clamped to \eqn{[\epsilon, 1 - \epsilon]} with
#' \eqn{\epsilon = 10^{-5}} to avoid boundary issues in the beta
#' likelihood.
#'
#' @param y      Numeric vector --- the raw response.
#' @param type   Character: interval type.
#'   \code{"m"} = midpoint (default),
#'   \code{"l"} = left-aligned,
#'   \code{"r"} = right-aligned.
#' @param ncuts  Integer: number of scale categories (default 100).
#' @param lim    Numeric: half-width of the uncertainty region
#'   (default 0.5).
#'
#' @return A matrix with columns \code{left}, \code{right},
#'   \code{yt} (midpoint approximation), \code{y} (original
#'   value), and \code{delta} (censoring indicator: 0 = exact,
#'   1 = left, 2 = right, 3 = interval).
#'
#' @examples
#' # Scale data with boundary observations
#' y <- c(0, 3, 5, 7, 9, 10)
#' check_response(y, type = "m", ncuts = 10)
#'
#' # Continuous data in (0, 1) — treated as uncensored
#' y_cont <- c(0.1, 0.3, 0.5, 0.7, 0.9)
#' check_response(y_cont, type = "m", ncuts = 100)
#'
#' @export
check_response <- function(y, type = "m", ncuts = 100L, lim = 0.5) {
  type <- match.arg(type, c("m", "l", "r"))
  ncuts <- as.integer(ncuts)
  n <- length(y)

  eps <- 1e-5

  # Detect if data is already on (0, 1)
  is_unit <- all(y > 0 & y < 1)

  if (is_unit) {
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
  delta <- integer(n)

  # Classify each observation
  for (i in seq_len(n)) {
    yi <- y[i]

    if (yi == 0) {
      # Left-censored: the latent value is below the first scale point.
      # We know Y < u_i where u_i = lim / ncuts.
      # Likelihood contribution: F(u_i | theta)
      delta[i] <- 1L
      y_left[i] <- eps
      y_right[i] <- lim / ncuts
    } else if (yi == ncuts) {
      # Right-censored: the latent value exceeds the last scale point.
      # We know Y > l_i where l_i = (ncuts - lim) / ncuts.
      # Likelihood contribution: 1 - F(l_i | theta)
      delta[i] <- 2L
      y_left[i] <- (ncuts - lim) / ncuts
      y_right[i] <- 1 - eps
    } else {
      # Interval-censored: standard scale observation.
      # Likelihood contribution: F(u_i | theta) - F(l_i | theta)
      delta[i] <- 3L
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
  }

  # Clamp all endpoints to (eps, 1 - eps) for numerical safety
  y_left <- pmin(pmax(y_left, eps), 1 - eps)
  y_right <- pmin(pmax(y_right, eps), 1 - eps)

  # Midpoint approximation for starting values
  yt <- y / ncuts
  yt <- pmin(pmax(yt, eps), 1 - eps)

  cbind(left = y_left, right = y_right, yt = yt, y = y, delta = delta)
}
