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
#' indicator \eqn{\delta_i} is (Lopes, 2024, Eq. 2.24):
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
#' @param type   Character: interval type (\code{"m"}, \code{"l"},
#'   or \code{"r"}).
#' @param lim    Numeric: half-width of uncertainty region (default
#'   0.5).
#' @param repar  Integer: reparameterization scheme (0, 1, or 2;
#'   default 2).
#'
#' @return Scalar: total log-likelihood.
#'
#' @examples
#' set.seed(42)
#' n <- 100
#' dat <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
#' sim <- betaregscale_simulate(
#'   formula = ~ x1 + x2, data = dat,
#'   beta = c(0, 0.5, -0.2), phi = 1 / 5
#' )
#' betaregscale_loglik(
#'   param = c(0, 0.5, -0.2, 1 / 5),
#'   formula = y ~ x1 + x2, data = sim
#' )
#'
#' @importFrom stats model.frame model.matrix model.response terms
#' @export
betaregscale_loglik <- function(param,
                                formula,
                                data,
                                link = "logit",
                                link_phi = "logit",
                                ncuts = 100L,
                                type = "m",
                                lim = 0.5,
                                repar = 2L) {
  # Validate links
  link <- match.arg(link, .mu_links)
  link_phi <- match.arg(link_phi, .phi_links)
  repar <- as.integer(repar)

  # Build model matrices
  mf <- stats::model.frame(formula, data = data)
  Y <- .extract_response(mf, data, ncuts = ncuts, type = type, lim = lim)
  X <- stats::model.matrix(mf, data = data)

  # Dispatch to C++
  .betaregscale_loglik_fixed_cpp(
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
#' @inheritParams betaregscale_loglik
#'
#' @return Scalar: total log-likelihood.
#'
#' @examples
#' set.seed(42)
#' n <- 100
#' dat <- data.frame(
#'   x1 = rnorm(n), x2 = rnorm(n),
#'   z1 = runif(n)
#' )
#' sim <- betaregscale_simulate_z(
#'   formula_x = ~ x1 + x2, formula_z = ~z1,
#'   data = dat,
#'   beta = c(0.2, -0.5, 0.3), zeta = c(0.5, -0.5)
#' )
#' betaregscale_loglik_z(
#'   param = c(0.2, -0.5, 0.3, 0.5, -0.5),
#'   formula = y ~ x1 + x2 | z1, data = sim
#' )
#'
#' @importFrom Formula as.Formula Formula
#' @importFrom stats delete.response
#' @export
betaregscale_loglik_z <- function(param,
                                  formula = y ~ x1 + x2 | z1,
                                  data,
                                  link = "logit",
                                  link_phi = "logit",
                                  ncuts = 100L,
                                  type = "m",
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
  Y <- .extract_response(mf, data, ncuts = ncuts, type = type, lim = lim)
  X <- stats::model.matrix(mtX, mf)
  Z <- stats::model.matrix(mtZ, mf)

  # Dispatch to C++
  .betaregscale_loglik_variable_cpp(
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
