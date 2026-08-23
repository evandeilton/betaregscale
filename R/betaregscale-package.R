#' betaregscale: Beta Regression for Interval-Censored Scale-Derived Outcomes
#'
#' @description
#' Maximum-likelihood estimation of beta regression models for responses
#' derived from bounded rating scales. Observations are treated as
#' interval-censored on (0, 1) after a scale-to-unit transformation.
#' The complete likelihood supports mixed censoring types: uncensored
#' (exact), left-censored, right-censored, and interval-censored
#' observations. Both fixed- and
#' variable-dispersion submodels are supported, with flexible link
#' functions for the mean and precision components. A compiled C++
#' backend (via Rcpp and RcppArmadillo) provides numerically stable,
#' high-performance log-likelihood evaluation. Standard S3 methods
#' (print(), summary(), coef(), fitted(), residuals(), predict(),
#' plot(), confint(), vcov(), logLik(), AIC(), BIC()) are available
#' for fitted objects.
#'
#' @section Main functions:
#' \describe{
#'   \item{\code{\link{brs}}}{Unified fitting interface for both
#'     fixed- and variable-dispersion models.}
#'   \item{\code{\link{brs_fit_fixed}}}{Fit a fixed-dispersion model.}
#'   \item{\code{\link{brs_fit_var}}}{Fit a variable-dispersion model.}
#'   \item{\code{\link{brsmm}}}{Fit a mixed-effects beta interval model
#'     with Gaussian random intercepts.}
#'   \item{\code{\link{brs_sim}}}{Simulate interval-censored
#'     data from fixed or variable-dispersion beta models.}
#'   \item{\code{\link{brs_bootstrap}}}{Parametric bootstrap confidence
#'     intervals for \code{brs} model parameters.}
#'   \item{\code{\link{brs_cens}}}{Visual and tabular summary of
#'     censoring structure.}
#'   \item{\code{\link{brs_prep}}}{Pre-process analyst data (validate,
#'     classify censoring, and rescale) before model fitting.}
#' }
#'
#' @section S3 methods:
#' Objects of class \code{"brs"} support: \code{print()},
#' \code{summary()}, \code{coef()}, \code{vcov()}, \code{logLik()}, \code{AIC()},
#' \code{BIC()}, \code{nobs()}, \code{formula()}, \code{model.matrix()},
#' \code{fitted()}, \code{residuals()}, \code{predict()}, \code{confint()},
#' and \code{plot()}.
#'
#' The \code{coef()} and \code{vcov()} methods accept a
#' \code{model = c("full", "mean", "precision")} argument following
#' the \pkg{betareg} package convention.
#'
#' @section Censoring types:
#' The complete likelihood supports four
#' censoring types, classified automatically by
#' \code{\link{brs_check}}:
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
#' @references
#' Lopes, J. E. (2023). \emph{Modelos de regressao beta para dados de escala}.
#' Master's dissertation, Universidade Federal do Parana, Curitiba.
#' URI: https://hdl.handle.net/1884/86624.
#'
#' Ferrari, S. L. P., and Cribari-Neto, F. (2004).
#' Beta regression for modelling rates and proportions.
#' \emph{Journal of Applied Statistics}, \bold{31}(7), 799--815.
#' \doi{10.1080/0266476042000214501}
#'
#' @useDynLib betaregscale, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @keywords internal
"_PACKAGE"

# -- Global variables for ggplot2 .data pronoun --------------------------------
utils::globalVariables(".data")
