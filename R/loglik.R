# ============================================================================ #
# Log-likelihood functions (R wrappers for the C++ backend)
# ============================================================================ #

#' Log-likelihood for fixed-dispersion beta interval regression
#'
#' @description
#' Computes the total log-likelihood for a beta regression model with
#' mixed-censored responses and a single (scalar) dispersion parameter.
#' The heavy computation is delegated to a compiled C++ backend.
#'
#' @details
#' The complete likelihood for observation \eqn{i} with censoring
#' indicator \eqn{\delta_i} is:
#'
#' \describe{
#'   \item{\eqn{\delta = 0} (uncensored)}{\eqn{\ell_i = \log f(y_i | a_i, b_i)}}
#'   \item{\eqn{\delta = 1} (left-censored)}{\eqn{\ell_i = \log F(u_i | a_i, b_i)}}
#'   \item{\eqn{\delta = 2} (right-censored)}{\eqn{\ell_i = \log(1 - F(l_i | a_i, b_i))}}
#'   \item{\eqn{\delta = 3} (interval-censored)}{\eqn{\ell_i = \log(F(u_i | a_i, b_i) - F(l_i | a_i, b_i))}}
#' }
#'
#' where \eqn{a_i} and \eqn{b_i} are the beta shape parameters derived
#' from the mean \eqn{\mu_i = g^{-1}(x_i'\beta)} and scalar dispersion
#' \eqn{\phi = h^{-1}(\gamma)} through the chosen reparameterization.
#'
#' @param param  Numeric vector of length \eqn{p + 1}: the first
#'   \eqn{p} elements are the regression coefficients \eqn{\beta},
#'   and the last element is the (link-scale) dispersion parameter.
#' @param formula One-sided or two-sided formula for the mean model.
#' @param data   Data frame containing the response and predictors.
#' @param link   Character: link function for the mean (default
#'   \code{"logit"}).
#' @param link_phi Character: link function for the dispersion
#'   (default \code{"logit"}).
#' @param ncuts  Integer: number of scale categories (default 100).
#' @param lim    Numeric: half-width of uncertainty region (default
#'   0.5).
#' @param repar  Integer: reparameterization scheme (0, 1, or 2;
#'   default 2).
#'
#' @return Scalar: total log-likelihood.
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
#' brs_loglik(
#'   param = c(0, 0.5, -0.2, 1 / 5),
#'   formula = y ~ x1 + x2, data = prep
#' )
#' }
#'
#' @importFrom stats model.frame model.matrix model.response terms
#' @keywords internal
#' @noRd
brs_loglik <- function(param,
                       formula,
                       data,
                       link = "logit",
                       link_phi = "logit",
                       ncuts = 100L,
                       lim = 0.5,
                       repar = 2L) {
  # Validate links
  link <- match.arg(link, .mu_links)
  link_phi <- match.arg(link_phi, .phi_links)
  repar <- as.integer(repar)

  # Build model matrices
  mf <- stats::model.frame(formula, data = data)
  Y <- .extract_response(mf, data, ncuts = ncuts, lim = lim)
  X <- stats::model.matrix(mf, data = data)

  # Dispatch to C++
  .brs_loglik_fixed_cpp(
    param         = as.numeric(param),
    X             = X,
    y_left        = as.numeric(Y[, "left"]),
    y_right       = as.numeric(Y[, "right"]),
    yt            = as.numeric(Y[, "yt"]),
    delta         = as.integer(Y[, "delta"]),
    link_mu_code  = link_to_code(link),
    link_phi_code = link_to_code(link_phi),
    repar         = repar
  )
}


#' Log-likelihood for variable-dispersion beta interval regression
#'
#' @description
#' Computes the total log-likelihood for a beta regression model with
#' mixed-censored responses and observation-specific dispersion
#' governed by a second linear predictor.  Uses the compiled C++
#' backend.
#'
#' @details
#' The formula should use the \code{\link[Formula]{Formula}} pipe
#' notation: \code{y ~ x1 + x2 | z1 + z2}, where the left-hand side
#' of \code{|} defines the mean model and the right-hand side defines
#' the dispersion model.
#'
#' @return Scalar: total log-likelihood.
#'
#' @examples
#' \donttest{
#' dat <- data.frame(
#'   y = c(
#'     0, 5, 20, 50, 75, 90, 100, 30, 60, 45,
#'     10, 40, 55, 70, 85, 25, 35, 65, 80, 15
#'   ),
#'   x1 = rep(c(1, 2), 10),
#'   x2 = rep(c(0, 0, 1, 1), 5),
#'   z1 = rep(c(0, 1), 10)
#' )
#' prep <- brs_prep(dat, ncuts = 100)
#' brs_loglik_var(
#'   param = c(0.2, -0.5, 0.3, 0.5, -0.5),
#'   formula = y ~ x1 + x2 | z1, data = prep
#' )
#' }
#'
#' @importFrom Formula as.Formula Formula
#' @importFrom stats delete.response
#' @keywords internal
#' @noRd
brs_loglik_var <- function(param,
                           formula = y ~ x1 + x2 | z1,
                           data,
                           link = "logit",
                           link_phi = "logit",
                           ncuts = 100L,
                           lim = 0.5,
                           repar = 2L) {
  # Validate
  link <- match.arg(link, .mu_links)
  link_phi <- match.arg(link_phi, .phi_links)
  repar <- as.integer(repar)

  # Parse multi-part formula
  formula <- Formula::as.Formula(formula)
  if (length(formula)[2L] < 2L) {
    formula <- Formula::as.Formula(formula(formula), ~1)
  } else if (length(formula)[2L] > 2L) {
    formula <- Formula::Formula(formula(formula, rhs = 1:2))
  }

  mf <- stats::model.frame(formula, data = data)
  mtX <- stats::terms(formula, data = data, rhs = 1L)
  mtZ <- stats::delete.response(stats::terms(formula, data = data, rhs = 2L))
  Y <- .extract_response(mf, data, ncuts = ncuts, lim = lim)
  X <- stats::model.matrix(mtX, mf)
  Z <- stats::model.matrix(mtZ, mf)

  # Dispatch to C++
  .brs_loglik_variable_cpp(
    param         = as.numeric(param),
    X             = X,
    Z             = Z,
    y_left        = as.numeric(Y[, "left"]),
    y_right       = as.numeric(Y[, "right"]),
    yt            = as.numeric(Y[, "yt"]),
    delta         = as.integer(Y[, "delta"]),
    link_mu_code  = link_to_code(link),
    link_phi_code = link_to_code(link_phi),
    repar         = repar
  )
}


#' Per-observation log-likelihood contributions (R mirror of the C++ backend)
#'
#' @description
#' Vectorised R implementation of the censored beta contributions used by
#' the compiled likelihood (\code{src/brs_common.h}, keep both in sync):
#' exact (\code{delta = 0}), left-censored (1), right-censored (2) and
#' interval-censored (3). Same rules as the C++ code: endpoints clamped to
#' \code{[1e-5, 1 - 1e-5]}, no probability floor, tail chosen by the mean
#' \eqn{a/(a+b)}, \code{pbeta()} in plain scale, and the endpoint Laplace
#' approximation below \code{1e-240}. Non-finite contributions become
#' \code{-1e6} (\code{LOG_PENALTY}).
#'
#' @param delta Integer censoring indicators.
#' @param left,right Interval endpoints on (0, 1).
#' @param yt Exact response on (0, 1) (used when \code{delta = 0}).
#' @param a,b Beta shape parameters (vectors, one per observation, or scalars).
#' @return Numeric vector of log-contributions.
#' @keywords internal
#' @noRd
.brs_obs_loglik <- function(delta, left, right, yt, a, b) {
  eps <- 1e-5
  p_tiny <- 1e-240
  n <- length(delta)
  a <- rep_len(as.numeric(a), n)
  b <- rep_len(as.numeric(b), n)
  left  <- pmin(pmax(left,  eps), 1 - eps)
  right <- pmin(pmax(right, eps), 1 - eps)
  yt    <- pmin(pmax(yt,    eps), 1 - eps)

  log1mexp <- function(z) {
    ifelse(z > log(2), log1p(-exp(-z)), log(-expm1(-z)))
  }
  # Endpoint Laplace approximation of the tail mass beyond x (vectorised):
  # log f(x) - log|g'| + log(1 + g''/g'^2) + log(1 - exp(-|g'| width)),
  # g = log f; -Inf when x is on the wrong side of the mode.
  laplace <- function(x, a, b, width, lower) {
    gp  <- (a - 1) / x - (b - 1) / (1 - x)
    ok  <- ifelse(lower, gp > 0, gp < 0)
    ok[is.na(ok)] <- FALSE
    s   <- abs(gp)
    gpp <- -(a - 1) / x^2 - (b - 1) / (1 - x)^2
    corr <- 1 + gpp / s^2
    v <- stats::dbeta(x, a, b, log = TRUE) - log(s)
    add <- is.finite(corr) & corr > 0
    v[add] <- v[add] + log(corr[add])
    fw <- is.finite(width)
    v[fw] <- v[fw] + log1mexp((s * width)[fw])
    v[!ok] <- -Inf
    v
  }
  # log(p) when p is trustworthy, else the Laplace fallback
  tail_log <- function(p, x, a, b, width, lower) {
    out <- rep(-Inf, length(p))
    use <- is.finite(p) & p >= p_tiny
    out[use] <- log(p[use])
    if (any(!use)) {
      out[!use] <- laplace(x[!use], a[!use], b[!use], width[!use], lower[!use])
    }
    out
  }

  lp <- rep(-Inf, n)
  i0 <- delta == 0L
  i1 <- delta == 1L
  i2 <- delta == 2L
  i3 <- delta == 3L
  if (any(i0)) {
    lp[i0] <- stats::dbeta(yt[i0], a[i0], b[i0], log = TRUE)
  }
  if (any(i1)) {
    lp[i1] <- tail_log(stats::pbeta(right[i1], a[i1], b[i1]),
                       right[i1], a[i1], b[i1], rep(Inf, sum(i1)),
                       rep(TRUE, sum(i1)))
  }
  if (any(i2)) {
    lp[i2] <- tail_log(stats::pbeta(left[i2], a[i2], b[i2], lower.tail = FALSE),
                       left[i2], a[i2], b[i2], rep(Inf, sum(i2)),
                       rep(FALSE, sum(i2)))
  }
  if (any(i3)) {
    l <- left[i3]; u <- right[i3]; a3 <- a[i3]; b3 <- b[i3]
    lower <- 0.5 * (l + u) <= a3 / (a3 + b3)
    lower[is.na(lower)] <- TRUE
    p1 <- ifelse(lower, stats::pbeta(u, a3, b3),
                        stats::pbeta(l, a3, b3, lower.tail = FALSE))
    p2 <- ifelse(lower, stats::pbeta(l, a3, b3),
                        stats::pbeta(u, a3, b3, lower.tail = FALSE))
    v <- rep(-Inf, length(l))
    use <- is.finite(p1) & p1 >= p_tiny
    area <- p1 - p2
    ok <- use & is.finite(area) & area > 0
    v[ok] <- log(area[ok])
    if (any(!use)) {
      v[!use] <- laplace(ifelse(lower, u, l)[!use], a3[!use], b3[!use],
                         (u - l)[!use], lower[!use])
    }
    v[!(u > l)] <- -Inf
    lp[i3] <- v
  }
  lp[!is.finite(lp)] <- -1e6
  lp
}
