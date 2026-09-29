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


# -- Reparameterization x link compatibility -------------------------------- #

#' Link functions allowed under each reparameterization
#'
#' @description
#' Each parameter of the beta distribution is modelled through a link whose
#' inverse maps the real line onto the parameter's domain. Parameters on
#' \eqn{(0, 1)} (the mean under \code{repar = 1, 2}; the dispersion under
#' \code{repar = 2}) use \code{"logit"}, \code{"probit"}, \code{"cauchit"}
#' or \code{"cloglog"}. Parameters on \eqn{(0, \infty)} (both shapes under
#' \code{repar = 0}; the precision under \code{repar = 1}) use \code{"log"}
#' or \code{"sqrt"}. \code{"identity"}, \code{"inverse"} and \code{"1/mu^2"}
#' are not accepted for positive parameters: their inverse does not map the
#' real line onto \eqn{(0, \infty)} (negative values, a discontinuity at 0,
#' undefined for \eqn{\eta \le 0}).
#'
#' @keywords internal
.brs_link_table <- list(
  "0" = list(mu = c("log", "sqrt"), phi = c("log", "sqrt")),
  "1" = list(mu = c("logit", "probit", "cauchit", "cloglog"),
             phi = c("log", "sqrt")),
  "2" = list(mu = c("logit", "probit", "cauchit", "cloglog"),
             phi = c("logit", "probit", "cauchit", "cloglog"))
)

#' Default links implied by a reparameterization
#' @param repar Integer (0, 1, or 2).
#' @return Named character vector \code{c(link = , link_phi = )}.
#' @keywords internal
#' @noRd
.brs_default_links <- function(repar) {
  # First entry of each cell of .brs_link_table: (0,1)-links logit, positive log.
  switch(as.character(as.integer(repar)),
    "0" = c(link = "log", link_phi = "log"),
    "1" = c(link = "logit", link_phi = "log"),
    "2" = c(link = "logit", link_phi = "logit"),
    stop("`repar` must be one of 0, 1, or 2.", call. = FALSE)
  )
}

#' @keywords internal
#' @noRd
.brs_link_table_text <- function() {
  # Human-readable form of .brs_link_table for error messages.
  paste(
    "  repar = 0 (shapes p, q):       link in {log, sqrt};",
    "link_phi in {log, sqrt}",
    "\n  repar = 1 (mean, precision):   link in {logit, probit, cauchit, cloglog};",
    "link_phi in {log, sqrt}",
    "\n  repar = 2 (mean, dispersion):  link in {logit, probit, cauchit, cloglog};",
    "link_phi in {logit, probit, cauchit, cloglog}"
  )
}

#' Resolve and validate the links for a reparameterization
#'
#' @description
#' \code{NULL} links take the default implied by \code{repar}
#' (\code{0 -> log/log}, \code{1 -> logit/log}, \code{2 -> logit/logit});
#' explicit names must match \code{.brs_link_table} exactly (no partial
#' matching) or an error shows the compatibility table.
#'
#' @param link,link_phi Character or \code{NULL}.
#' @param repar Integer (0, 1, or 2).
#' @return \code{list(link = , link_phi = )}.
#' @keywords internal
#' @noRd
.resolve_links <- function(link = NULL, link_phi = NULL, repar = 2L) {
  repar <- as.integer(repar)
  if (length(repar) != 1L || is.na(repar) || !(repar %in% 0:2)) {
    stop("`repar` must be one of 0, 1, or 2.", call. = FALSE)
  }
  defaults <- .brs_default_links(repar)
  allowed <- .brs_link_table[[as.character(repar)]]
  # NULL -> default for this repar; explicit names are matched exactly.
  if (is.null(link)) link <- defaults[["link"]]
  if (is.null(link_phi)) link_phi <- defaults[["link_phi"]]
  .check_link_for_repar(link, allowed$mu, "link", repar)
  .check_link_for_repar(link_phi, allowed$phi, "link_phi", repar)
  list(link = link, link_phi = link_phi)
}

#' @keywords internal
#' @noRd
.check_link_for_repar <- function(value, allowed, arg, repar) {
  if (!is.character(value) || length(value) != 1L || is.na(value)) {
    stop("`", arg, "` must be a single link name or NULL.", call. = FALSE)
  }
  if (value %in% allowed) {
    return(invisible(value))
  }
  # Links whose inverse does not map R onto (0, Inf) get a one-line reason.
  extra <- if (value %in% c("identity", "inverse", "1/mu^2")) {
    paste0(
      "\n`identity`, `inverse` and `1/mu^2` are not accepted for positive ",
      "parameters: their inverse does not map the real line onto (0, Inf)."
    )
  } else {
    ""
  }
  stop(
    "`", arg, " = \"", value, "\"` is not compatible with `repar = ", repar,
    "`. Compatible links by reparameterization:\n",
    .brs_link_table_text(), extra,
    call. = FALSE
  )
}


# -- Mean and clamping helpers ----------------------------------------------- #

#' Mean of the beta response from the fitted parameters
#'
#' @description
#' \eqn{E[Y] = a / (a + b)}. Under \code{repar = 1, 2} the first parameter
#' is the mean itself and is returned unchanged; under \code{repar = 0} the
#' parameters are the shapes \eqn{(p, q)} and the mean is \eqn{p / (p + q)}.
#' Every user-facing "response-scale mean" goes through this helper;
#' \code{object$hatmu} and \code{brs_repar(mu = )} keep the first parameter.
#'
#' @param mu First parameter (mean, or shape \eqn{p} under \code{repar = 0}).
#' @param phi Second parameter (dispersion/precision, or shape \eqn{q}).
#' @param repar Integer (0, 1, or 2).
#' @return Numeric vector of means.
#' @keywords internal
#' @noRd
.brs_mean <- function(mu, phi, repar) {
  # E[Y] = a/(a+b): p/(p+q) under repar 0, the first parameter itself otherwise.
  if (as.integer(repar) == 0L) {
    mu <- as.numeric(mu)
    phi <- as.numeric(phi)
    return(mu / (mu + phi))
  }
  mu
}

#' Clamp the first parameter to its valid range (R mirror of the C++ code)
#'
#' @description
#' Identical to \code{clamp_mu_by_repar()} in \code{src/brs_common.h}
#' (keep both in sync), applied to the inverse-linked first parameter
#' before the shapes are formed. Under \code{repar = 1, 2} it is the mean:
#' \code{[1e-5, 1 - 1e-5]}. Under \code{repar = 0} it is the shape
#' \eqn{p > 0}: \code{[1e-5, 1e8]}, the same range
#' \code{clamp_phi_by_repar} gives the other shape. \code{+-Inf} become the
#' bounds; \code{NaN}/\code{NA} propagate (the likelihood then penalises them).
#'
#' @param mu Numeric vector.
#' @param repar Integer (0, 1, or 2).
#' @return Numeric vector.
#' @keywords internal
#' @noRd
.clamp_mu_by_repar <- function(mu, repar) {
  # Mirror of C++ clamp_mu_by_repar(): shape p in [1e-5, 1e8], mean in [1e-5, 1-1e-5];
  # +-Inf -> bounds, NaN/NA propagate (pmin/pmax keep them).
  eps_unit <- 1e-5
  hi <- if (as.integer(repar) == 0L) 1e8 else 1 - eps_unit
  pmin(pmax(as.numeric(mu), eps_unit), hi)
}

#' Clamp the second parameter to its valid range (R mirror of the C++ code)
#'
#' @description
#' Identical to \code{clamp_phi_by_repar()} in \code{src/brs_common.h}
#' (keep both in sync): \code{+-Inf} become the bounds, \code{NaN}/\code{NA}
#' propagate; \code{repar = 2} clamps to \code{[1e-5, 1 - 1e-5]}, otherwise to
#' \code{[1e-5, 1e8]}. Applied to the fitted \code{hatphi} so that the R
#' side sees the same value the compiled likelihood used.
#'
#' @param phi Numeric vector.
#' @param repar Integer (0, 1, or 2).
#' @return Numeric vector.
#' @keywords internal
#' @noRd
.clamp_phi_by_repar <- function(phi, repar) {
  # Mirror of C++ clamp_phi_by_repar(): +-Inf -> bounds, NaN/NA propagate.
  eps_unit <- 1e-5
  hi <- if (as.integer(repar) == 2L) 1 - eps_unit else 1e8
  pmin(pmax(as.numeric(phi), eps_unit), hi)
}

#' Fitted (first parameter, second parameter) for the data or for newdata
#'
#' @description
#' Returns the pair that \code{\link{brs_repar}} expects: the first
#' parameter (mean, or shape \eqn{p} under \code{repar = 0}) and the second
#' one, both on the response scale. Use it wherever shapes are needed for
#' new observations instead of feeding \code{predict(type = "response")}
#' (which is \eqn{E[Y]}) back into \code{brs_repar()}.
#'
#' @param object A \code{"brs"} or \code{"brsmm"} fit.
#' @param newdata Optional data frame.
#' @return \code{list(mu = , phi = )}, both of the same length.
#' @keywords internal
#' @noRd
.brs_predict_params <- function(object, newdata = NULL) {
  # Stored fit values, or newdata via predict(type = "link"/"precision").
  if (is.null(newdata)) {
    if (inherits(object, "brsmm")) {
      mu <- object$fitted_mu
      phi <- object$fitted_phi
    } else {
      mu <- object$hatmu
      phi <- object$hatphi
    }
  } else {
    eta <- stats::predict(object, newdata = newdata, type = "link")
    mu <- .clamp_mu_by_repar(apply_inv_link(eta, object$link), object$repar)
    phi <- stats::predict(object, newdata = newdata, type = "precision")
  }
  mu <- as.numeric(mu)
  list(mu = mu, phi = rep_len(as.numeric(phi), length(mu)))
}


#' Apply the inverse-link function to a linear predictor
#'
#' @description
#' Evaluates the inverse of a standard link function for a given
#' linear-predictor vector or scalar.
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
# BUG-H07: use direct formulas instead of make.link() which allocates a list
# of 4 closures on every call.
apply_inv_link <- function(eta, link) {
  switch(link,
    logit    = stats::plogis(eta),
    probit   = stats::pnorm(eta),
    cauchit  = 0.5 + atan(eta) / pi,
    cloglog  = -expm1(-exp(eta)),
    log      = exp(eta),
    sqrt     = pmax(eta, 0) ^ 2,
    "1/mu^2" = 1 / sqrt(pmax(eta, .Machine$double.eps)),
    inverse  = 1 / eta,
    identity = eta,
    stop("Unknown link function: '", link, "'.", call. = FALSE)
  )
}


#' Apply the forward link function to a response value
#'
#' @description
#' Evaluates the link function \eqn{g(\mu)} for a given response vector.
#' Used internally for starting-value computation.
#'
#' @param mu   Numeric vector of response values.
#' @param link Character string naming the link function (same set as
#'   \code{\link{apply_inv_link}}).
#'
#' @return Numeric vector containing \eqn{g(\mu)}.
#'
#' @keywords internal
apply_link <- function(mu, link) {
  switch(link,
    logit    = stats::qlogis(mu),
    probit   = stats::qnorm(mu),
    cauchit  = tan(pi * (mu - 0.5)),
    cloglog  = log(-log(1 - mu)),
    log      = log(mu),
    sqrt     = sqrt(mu),
    "1/mu^2" = 1 / mu^2,
    inverse  = 1 / mu,
    identity = mu,
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
#'   \item{\code{repar = 0}}{Shapes: \eqn{a = p,\; b = q} with
#'     \eqn{p, q > 0} (the first argument is \eqn{p}, the second
#'     \eqn{q}). Both parameters live on \eqn{(0, \infty)} and the
#'     mean is \eqn{E[Y] = p / (p + q)}. Regression directly on the
#'     shapes is a package extension: the dissertation presents the
#'     \eqn{(p, q)} form (its eq. \code{eqn_beta_p1}) and builds the
#'     regression models on the two reparameterizations below.}
#'   \item{\code{repar = 1}}{Ferrari–Cribari-Neto (dissertation
#'     "parametrização 1", eq. \code{eqn_beta_p2}):
#'     \eqn{a = \mu\phi,\; b = (1 - \mu)\phi}, where \eqn{\mu \in (0, 1)}
#'     is the mean and \eqn{\phi > 0} acts as a precision parameter.}
#'   \item{\code{repar = 2}}{Mean–dispersion (dissertation
#'     "parametrização 2", eq. \code{eqn_beta_p3}):
#'     \eqn{a = \mu(1-\phi)/\phi,\; b = (1-\mu)(1-\phi)/\phi},
#'     where \eqn{\mu \in (0, 1)} is the mean and \eqn{\phi \in (0,1)}
#'     is a dispersion parameter (\eqn{Var[Y] = \phi\,\mu(1-\mu)}).}
#' }
#' The admissible link functions follow from these domains; see the
#' 'Reparameterizations and links' section of \code{\link{brs}}.
#'
#' @param mu   Numeric vector: the first parameter. The mean, in
#'   \eqn{(0, 1)}, for \code{repar = 1, 2}; the shape \eqn{p > 0} for
#'   \code{repar = 0}.
#' @param phi  Numeric vector (or scalar): the second parameter
#'   (precision \eqn{\phi > 0}, dispersion \eqn{\phi \in (0, 1)}, or shape
#'   \eqn{q > 0}).
#' @param repar Integer (0, 1, or 2) selecting the scheme.
#'
#' @return A \code{data.frame} with columns \code{shape1} and
#'   \code{shape2}.
#'
#' @references
#' Lopes, J. E. (2023). \emph{Modelos de regressao beta para dados de escala}.
#' Master's dissertation, Universidade Federal do Parana, Curitiba.
#' URI: https://hdl.handle.net/1884/86624.
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
  if (repar == 0L) {
    if (any(mu <= 0)) {
      stop("For `repar = 0`, `mu` is the shape p and must be > 0.", call. = FALSE)
    }
  } else if (any(mu <= 0 | mu >= 1)) {
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
