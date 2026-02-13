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
#' @examples
#' beta_reparam(mu = 0.5, phi = 0.2, repar = 2)
#'
#' @export
beta_reparam <- function(mu, phi, repar = 2L) {
  repar <- as.integer(repar)
  stopifnot(repar %in% 0:2)

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
