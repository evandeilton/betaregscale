# ============================================================================ #
# Link-function utilities and beta reparameterization
# ============================================================================ #

# -- Valid link-function sets ------------------------------------------------ #

#' Valid link names for the mean submodel
#' @keywords internal
.mu_links <- c("logit", "probit", "cauchit", "cloglog")

#' Valid link names for the dispersion submodel
#' @keywords internal
.phi_links <- c(
  "logit", "probit", "cauchit", "cloglog",
  "identity", "log", "sqrt", "1/mu^2", "inverse"
)


#' Apply the inverse-link function to a linear predictor
#'
#' @description
#' Evaluates the inverse of a standard link function for a given
#' linear-predictor vector or scalar. This is a convenience wrapper
#' around \code{\link[stats]{make.link}}.
#'
#' @param eta  Numeric vector or scalar — the linear predictor
#'   \eqn{\eta = X \beta}.
#' @param link Character string naming the link function. Supported
#'   values: \code{"logit"}, \code{"probit"}, \code{"cauchit"},
#'   \code{"cloglog"}, \code{"log"}, \code{"sqrt"}, \code{"1/mu^2"},
#'   \code{"inverse"}, \code{"identity"}.
#'
#' @return Numeric vector (or scalar) of the same length as \code{eta},
#'   containing \eqn{g^{-1}(\eta)}.
#'
#' @keywords internal
apply_inv_link <- function(eta, link) {
  switch(link,
    logit = stats::make.link("logit")$linkinv(eta),
    probit = stats::make.link("probit")$linkinv(eta),
    cauchit = stats::make.link("cauchit")$linkinv(eta),
    cloglog = stats::make.link("cloglog")$linkinv(eta),
    log = stats::make.link("log")$linkinv(eta),
    sqrt = stats::make.link("sqrt")$linkinv(eta),
    "1/mu^2" = stats::make.link("1/mu^2")$linkinv(eta),
    inverse = stats::make.link("inverse")$linkinv(eta),
    identity = eta,
    stop("Unknown link function: '", link, "'.", call. = FALSE)
  )
}


#' Map link-function name to integer code for the C++ backend
#'
#' @param link Character link-function name.
#' @return Integer code consumed by the compiled likelihood.
#' @keywords internal
link_to_code <- function(link) {
  code <- match(
    link,
    c(
      "logit", "probit", "cauchit", "cloglog",
      "log", "sqrt", "inverse", "1/mu^2", "identity"
    )
  )
  if (is.na(code)) {
    stop("Unsupported link function: '", link, "'.", call. = FALSE)
  }
  code - 1L # C++ uses 0-indexed codes
}


#' Reparameterize (mu, phi) into beta shape parameters
#'
#' @description
#' Converts a mean–dispersion pair \eqn{(\mu, \phi)} to the shape
#' parameters \eqn{(a, b)} of the beta distribution under one of
#' three reparameterization schemes.
#'
#' @details
#' The three schemes are:
#' \describe{
#'   \item{\code{repar = 0}}{Direct: \eqn{a = \mu,\; b = \phi}.}
#'   \item{\code{repar = 1}}{Ferrari–Cribari-Neto:
#'     \eqn{a = \mu\phi,\; b = (1 - \mu)\phi}, where \eqn{\phi}
#'     acts as a precision parameter.}
#'   \item{\code{repar = 2}}{Mean–variance:
#'     \eqn{a = \mu(1-\phi)/\phi,\; b = (1-\mu)(1-\phi)/\phi},
#'     where \eqn{\phi \in (0,1)} is analogous to a coefficient
#'     of variation.}
#' }
#'
#' @param mu   Numeric vector of mean values in \eqn{(0, 1)}.
#' @param phi  Numeric vector (or scalar) of dispersion values.
#' @param repar Integer (0, 1, or 2) selecting the scheme.
#'
#' @return A \code{data.frame} with columns \code{shape1} and
#'   \code{shape2}.
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
#' brs_repar(mu = 0.5, phi = 0.2, repar = 2)
#'
#' @export
brs_repar <- function(mu, phi, repar = 2L) {
  repar <- as.integer(repar)
  if (!(repar %in% 0:2)) {
    stop("`repar` must be one of 0, 1, or 2.", call. = FALSE)
  }
  if (!is.numeric(mu) || !is.numeric(phi)) {
    stop("`mu` and `phi` must be numeric.", call. = FALSE)
  }

  mu <- as.numeric(mu)
  phi <- as.numeric(phi)

  if (length(phi) == 1L && length(mu) > 1L) {
    phi <- rep(phi, length(mu))
  }
  if (length(mu) == 1L && length(phi) > 1L) {
    mu <- rep(mu, length(phi))
  }
  if (length(mu) != length(phi)) {
    stop("`mu` and `phi` must have compatible lengths.", call. = FALSE)
  }
  if (any(!is.finite(mu)) || any(!is.finite(phi))) {
    stop("`mu` and `phi` must be finite.", call. = FALSE)
  }
  if (any(mu <= 0 | mu >= 1)) {
    stop("`mu` must lie in (0, 1).", call. = FALSE)
  }
  if (repar == 2L) {
    if (any(phi <= 0 | phi >= 1)) {
      stop("For `repar = 2`, `phi` must lie in (0, 1).", call. = FALSE)
    }
  } else {
    if (any(phi <= 0)) {
      stop("For `repar = 0` or `repar = 1`, `phi` must be > 0.", call. = FALSE)
    }
  }

  switch(as.character(repar),
    "0" = data.frame(
      shape1 = as.numeric(mu),
      shape2 = as.numeric(phi)
    ),
    "1" = data.frame(
      shape1 = as.numeric(mu * phi),
      shape2 = as.numeric((1 - mu) * phi)
    ),
    "2" = data.frame(
      shape1 = as.numeric(mu * ((1 - phi) / phi)),
      shape2 = as.numeric((1 - mu) * ((1 - phi) / phi))
    )
  )
}
