#' betaregscale: Beta Regression for Interval-Censored Scale-Derived Outcomes
#'
#' @description
#' Maximum-likelihood estimation of beta regression models for responses
#' derived from bounded rating scales. Observations are treated as
#' interval-censored on (0, 1) after a scale-to-unit transformation.
#' The complete likelihood supports mixed censoring types: uncensored
#' (exact), left-censored, right-censored, and interval-censored
#' observations (Lopes, 2024, Eq. 2.24). Both fixed- and
#' variable-dispersion submodels are supported, with flexible link
#' functions for the mean and precision components. A compiled C++
#' backend (via Rcpp and RcppArmadillo) provides numerically stable,
#' high-performance log-likelihood evaluation.
#'
#' @section Main functions:
#' \describe{
#'   \item{\code{\link{betaregscale}}}{Unified fitting interface for both
#'     fixed- and variable-dispersion models.}
#'   \item{\code{\link{betaregscale_fit}}}{Fit a fixed-dispersion model.}
#'   \item{\code{\link{betaregscale_fit_z}}}{Fit a variable-dispersion model.}
#'   \item{\code{\link{betaregscale_loglik}}}{Compute the log-likelihood
#'     (fixed dispersion).}
#'   \item{\code{\link{betaregscale_loglik_z}}}{Compute the log-likelihood
#'     (variable dispersion).}
#'   \item{\code{\link{betaregscale_simulate}}}{Simulate interval-censored
#'     data from a fixed-dispersion beta model.}
#'   \item{\code{\link{betaregscale_simulate_z}}}{Simulate data from a
#'     variable-dispersion beta model.}
#'   \item{\code{\link{censoring_summary}}}{Visual and tabular summary of
#'     censoring structure.}
#' }
#'
#' @section S3 methods:
#' Objects of class \code{"betaregscale"} support: \code{print},
#' \code{summary}, \code{coef}, \code{vcov}, \code{logLik}, \code{AIC},
#' \code{BIC}, \code{nobs}, \code{formula}, \code{model.matrix},
#' \code{fitted}, \code{residuals}, \code{predict}, \code{confint},
#' and \code{plot}.
#'
#' The \code{coef()} and \code{vcov()} methods accept a
#' \code{model = c("full", "mean", "precision")} argument following
#' the \pkg{betareg} package convention.
#'
#' @section Censoring types:
#' The complete likelihood (Lopes, 2024, Eq. 2.24) supports four
#' censoring types, classified automatically by
#' \code{\link{check_response}}:
#' \describe{
#'   \item{\eqn{\delta = 0} (exact)}{Continuous observations in (0, 1).}
#'   \item{\eqn{\delta = 1} (left-censored)}{Observations at the scale
#'     minimum (y = 0).}
#'   \item{\eqn{\delta = 2} (right-censored)}{Observations at the scale
#'     maximum (y = ncuts).}
#'   \item{\eqn{\delta = 3} (interval-censored)}{Standard scale
#'     observations between the boundaries.}
#' }
#'
#' @author Jose Eduardo Lopes \email{evandeilton@@gmail.com}
#' @references
#' Lopes, J. E. (2024). \emph{Beta Regression for Interval-Censored
#' Scale-Derived Outcomes}. MSc Dissertation, PPGMNE/UFPR.
#'
#' Ferrari, S. and Cribari-Neto, F. (2004). Beta regression for modelling
#' rates and proportions. \emph{Journal of Applied Statistics},
#' \bold{31}(7), 799--815.
#'
#' @useDynLib betaregscale, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @keywords internal
"_PACKAGE"

# -- Global variables for ggplot2 .data pronoun --------------------------------
utils::globalVariables(".data")
