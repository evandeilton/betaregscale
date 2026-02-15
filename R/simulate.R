# ============================================================================ #
# Starting-value computation
# ============================================================================ #

#' Compute starting values for optimization
#'
#' @description
#' Obtains rough starting values for the beta regression parameters
#' by fitting a quasi-binomial GLM on the midpoint response.  This
#' provides a reasonable initialization for the interval likelihood
#' optimizer.
#'
#' @param formula A \code{\link[Formula]{Formula}} object (possibly
#'   multi-part).
#' @param data   Data frame.
#' @param link   Mean link function name.
#' @param link_phi Dispersion link function name.
#' @param ncuts  Number of scale categories.
#' @param type   Interval type.
#' @param lim    Uncertainty half-width.
#'
#' @return Named numeric vector of starting values.
#' @keywords internal
compute_start <- function(formula, data, link = "logit",
                          link_phi = "logit", ncuts = 100L,
                          type = "m", lim = 0.5) {
  link <- match.arg(link, .mu_links)

  formula <- Formula::as.Formula(formula)
  if (length(formula)[2L] < 2L) {
    formula <- Formula::as.Formula(formula(formula), ~1)
  } else if (length(formula)[2L] > 2L) {
    formula <- Formula::Formula(formula(formula, rhs = 1:2))
  }

  mf <- stats::model.frame(formula, data = data)
  mtX <- stats::terms(formula, data = data, rhs = 1L)
  mtZ <- stats::delete.response(
    stats::terms(formula, data = data, rhs = 2L)
  )
  Y <- .extract_response(mf, data, ncuts = ncuts, type = type, lim = lim)
  x <- stats::model.matrix(mtX, mf)
  z <- stats::model.matrix(mtZ, mf)

  # Midpoint response for starting-value GLM
  y <- rowMeans(Y[, c("left", "right"), drop = FALSE], na.rm = TRUE)

  # Mean-model starting values via quasi-binomial GLM
  q <- ncol(x)
  glm_data <- if (q == 1L) {
    data.frame(y = y, x)
  } else {
    data.frame(y = y, x[, -1L, drop = FALSE])
  }
  init_beta <- stats::coef(
    stats::glm(y ~ .,
      data = glm_data,
      family = stats::quasibinomial(link = link)
    )
  )

  # Dispersion starting values
  if (is.null(z) || ncol(z) < 2L) {
    # Intercept-only dispersion: use the mean of the link-inverse
    init_phi <- mean(stats::make.link(link_phi)$linkfun(
      pmin(pmax(y, 1e-4), 1 - 1e-4)
    ), na.rm = TRUE)
    names(init_phi) <- "phi"
  } else {
    glm_data_z <- if (ncol(z) == 1L) {
      data.frame(y = y, z)
    } else {
      data.frame(y = y, z[, -1L, drop = FALSE])
    }
    init_phi <- stats::coef(
      stats::glm(y ~ .,
        data = glm_data_z,
        family = stats::quasibinomial(link = link_phi)
      )
    )
    names(init_phi) <- paste0("phi_", colnames(z))
  }

  c(init_beta, init_phi)
}


# ============================================================================ #
# Simulation functions
# ============================================================================ #

#' Simulate data from a fixed-dispersion beta interval model
#'
#' @description
#' Generates observations from a beta regression model with a single
#' (scalar) dispersion parameter.  This is useful for Monte Carlo
#' studies, power analysis, and testing.  The \code{delta} argument
#' controls the censoring type of the simulated data.
#'
#' @details
#' \strong{Data generation process}:
#' \enumerate{
#'   \item The design matrix \eqn{X} is built from \code{formula} and
#'     \code{data}.
#'   \item The linear predictor is \eqn{\eta = X \beta}, and the mean
#'     is \eqn{\mu = g^{-1}(\eta)} where \eqn{g} is the link
#'     function.
#'   \item The scalar dispersion is
#'     \eqn{\phi = h^{-1}(\texttt{phi})} where \eqn{h} is the
#'     dispersion link.
#'   \item Beta shape parameters \eqn{(a, b)} are derived from
#'     \eqn{(\mu, \phi)} via the chosen reparameterization scheme.
#'   \item Raw values \eqn{y^*_i \sim \text{Beta}(a_i, b_i)} are
#'     drawn on \eqn{(0, 1)}.
#'   \item The raw values are transformed into the response matrix
#'     by \code{.build_simulated_response()} (see below).
#' }
#'
#' \strong{Role of the \code{delta} argument}:
#'
#' When \code{delta = NULL} (default), the raw values are rounded to
#' the scale grid (\eqn{y_{\text{grid}} = \text{round}(y^* \times K)})
#' and passed to \code{\link{check_response}} for automatic
#' classification: \eqn{y = 0 \to \delta = 1}, \eqn{y = K \to
#' \delta = 2}, otherwise \eqn{\delta = 3}.  The resulting dataset
#' has a natural mix of censoring types driven by the simulated
#' values.
#'
#' When \code{delta} is an integer in \eqn{\{0, 1, 2, 3\}},
#' \strong{all} observations are forced to that censoring type, but
#' the actual simulated \eqn{y^*} values are \strong{preserved} on
#' the grid so that each observation retains its covariate-driven
#' variation.  Specifically:
#'
#' \describe{
#'   \item{\code{delta = 0} (exact)}{The continuous \eqn{y^*} values
#'     are used directly on \eqn{(0, 1)};
#'     \eqn{l_i = u_i = y_t = y^*_i}.}
#'   \item{\code{delta = 1} (left-censored)}{The grid values
#'     \eqn{y_{\text{grid}} = \text{round}(y^* K)} are kept.
#'     \code{check_response()} is called with forced
#'     \code{delta = rep(1, n)}, producing:
#'     \eqn{l_i = \epsilon},
#'     \eqn{u_i = (y_{\text{grid}} + h) / K} for non-boundary,
#'     or \eqn{u_i = h / K} when \eqn{y_{\text{grid}} = 0}.}
#'   \item{\code{delta = 2} (right-censored)}{Same logic:
#'     \eqn{u_i = 1 - \epsilon},
#'     \eqn{l_i = (y_{\text{grid}} - h) / K} for non-boundary,
#'     or \eqn{l_i = (K - h) / K} when \eqn{y_{\text{grid}} = K}.}
#'   \item{\code{delta = 3} (interval-censored)}{Grid values are
#'     clamped to \eqn{[1, K-1]} (avoiding boundaries) and
#'     \code{check_response()} is called with forced
#'     \code{delta = rep(3, n)}.}
#' }
#'
#' \strong{Attribute \code{"bs_prepared"}}:
#'
#' When \code{delta != NULL}, the returned data frame carries the
#' attribute \code{"bs_prepared" = TRUE}.  This signals to
#' \code{.extract_response()} (and thus to \code{\link{betaregscale}},
#' \code{\link{betaregscale_loglik}}, etc.) that the pre-computed
#' columns \code{left}, \code{right}, \code{yt}, and \code{delta}
#' should be used directly, bypassing the automatic classification
#' of \code{\link{check_response}}.  Without this attribute, the
#' fitting functions would re-classify the response from the \code{y}
#' column alone, which would ignore the forced delta and produce
#' incorrect censoring indicators (e.g., an observation with
#' \eqn{y = 50} and forced \eqn{\delta = 2} would be reclassified as
#' \eqn{\delta = 3} by the boundary rules).
#'
#' When \code{delta = NULL}, the attribute is \strong{not} set, so
#' the default pipeline applies.
#'
#' @param formula One-sided formula specifying the mean model
#'   predictors (e.g., \code{~ x1 + x2}).
#' @param data   Data frame containing the predictor variables.
#' @param beta   Numeric vector of regression coefficients (length
#'   must equal the number of columns of the design matrix, including
#'   the intercept).
#' @param phi    Scalar dispersion parameter (on the link scale).
#' @param link   Mean link function (default \code{"logit"}).
#'   Supported: \code{"logit"}, \code{"probit"}, \code{"cloglog"},
#'   \code{"cauchit"}, \code{"log"}.
#' @param link_phi Dispersion link function (default \code{"logit"}).
#' @param ncuts  Integer: number of scale categories \eqn{K}
#'   (default 100).
#' @param type   \strong{Deprecated.}
#'   Interval type: \code{"m"}, \code{"l"}, or \code{"r"}.
#'   Use \code{\link{bs_prepare}} to control interval geometry instead.
#' @param lim    Numeric: half-width \eqn{h} of the uncertainty
#'   region (default 0.5).
#' @param repar  Integer: reparameterization scheme (default 2).
#'   0 = direct \eqn{(a, b)}, 1 = precision/Ferrari
#'   \eqn{(\mu, \phi = a + b)}, 2 = mean-variance/Bayer
#'   \eqn{(\mu, \phi = 1/(a + b + 1))}.
#' @param delta  Integer or \code{NULL}.  If \code{NULL} (default),
#'   censoring is determined automatically by
#'   \code{\link{check_response}}.
#'
#'   If an integer in \code{\{0, 1, 2, 3\}}, \strong{all} simulated
#'   observations are forced to that censoring type.  The actual
#'   simulated values are preserved so that observation-specific
#'   endpoints reflect the underlying covariate-driven variation.
#'   See Details.
#'
#' @return A \code{data.frame} with \eqn{n} rows and columns:
#'   \code{left}, \code{right}, \code{yt}, \code{y}, \code{delta},
#'   plus the predictor columns from \code{data}.
#'   When \code{delta != NULL}, the data frame carries the attribute
#'   \code{"bs_prepared" = TRUE}.
#'
#' @seealso \code{\link{betaregscale_simulate_z}} for variable-
#'   dispersion simulation; \code{\link{check_response}} for the
#'   endpoint computation rules; \code{\link{bs_prepare}} for
#'   analyst-facing data pre-processing.
#'
#' @examples
#' set.seed(42)
#' n <- 200
#' dat <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
#' sim <- betaregscale_simulate(
#'   formula = ~ x1 + x2, data = dat,
#'   beta = c(0.2, -0.5, 0.3), phi = 1 / 5,
#'   link = "logit", link_phi = "logit"
#' )
#' head(sim)
#'
#' # Force all observations to be interval-censored
#' sim3 <- betaregscale_simulate(
#'   formula = ~ x1 + x2, data = dat,
#'   beta = c(0.2, -0.5, 0.3), phi = 1 / 5,
#'   delta = 3
#' )
#' table(sim3$delta)
#'
#' # Force right-censored: y values vary, all delta = 2
#' sim2 <- betaregscale_simulate(
#'   formula = ~ x1 + x2, data = dat,
#'   beta = c(0.2, -0.5, 0.3), phi = 1 / 5,
#'   delta = 2
#' )
#' head(sim2[, c("left", "right", "y", "delta")])
#' # Note: left varies per observation, right = 1 - eps
#'
#' @importFrom stats rbeta
#' @export
betaregscale_simulate <- function(formula, data, beta, phi = 1 / 5,
                                  link = "logit",
                                  link_phi = "logit",
                                  ncuts = 100L, type = "m",
                                  lim = 0.5, repar = 2L,
                                  delta = NULL) {
  # Deprecation warning for type
  if (!missing(type)) {
    .Deprecated(msg = paste0(
      "The 'type' argument of betaregscale_simulate() is deprecated ",
      "and will be removed in a future version. ",
      "Use bs_prepare() to control interval geometry."
    ))
  }

  # Validate delta
  if (!is.null(delta)) {
    if (length(delta) != 1L || !(delta %in% 0:3)) {
      stop("'delta' must be NULL or a single integer in {0, 1, 2, 3}.",
        call. = FALSE
      )
    }
    delta <- as.integer(delta)
  }

  link <- match.arg(link, .mu_links)
  link_phi <- match.arg(link_phi, .phi_links)
  repar <- as.integer(repar)

  mfx <- stats::model.frame(formula, data = data)
  X <- stats::model.matrix(mfx, data = data)
  n <- nrow(data)

  # Linear predictor -> mu
  eta <- X %*% beta
  mu <- apply_inv_link(eta, link)
  phi_val <- apply_inv_link(phi, link_phi)

  # Reparameterize and simulate
  pars <- beta_reparam(mu = mu, phi = phi_val, repar = repar)
  y_raw <- stats::rbeta(n, shape1 = pars$shape1, shape2 = pars$shape2)

  # Build response matrix based on delta
  out_y <- .build_simulated_response(
    y_raw = y_raw, delta = delta, ncuts = ncuts,
    type = type, lim = lim
  )

  result <- data.frame(out_y, X[, -1L, drop = FALSE])

  # When delta is forced, mark data as pre-processed so that
  # .extract_response() uses the pre-computed columns directly
  if (!is.null(delta)) {
    attr(result, "bs_prepared") <- TRUE
  }

  result
}


#' Simulate data from a variable-dispersion beta interval model
#'
#' @description
#' Generates observations from a beta regression model with
#' observation-specific dispersion governed by a second linear
#' predictor.  The \code{delta} argument controls the censoring type
#' of the simulated data.
#'
#' @details
#' The data generation is analogous to
#' \code{\link{betaregscale_simulate}}, but with observation-specific
#' dispersion:
#' \enumerate{
#'   \item Mean linear predictor: \eqn{\mu_i = g^{-1}(X_i \beta)}.
#'   \item Dispersion linear predictor:
#'     \eqn{\phi_i = h^{-1}(Z_i \zeta)}.
#'   \item Beta shape parameters \eqn{(a_i, b_i)} are derived from
#'     \eqn{(\mu_i, \phi_i)} per the reparameterization.
#'   \item \eqn{y^*_i \sim \text{Beta}(a_i, b_i)}.
#'   \item Response matrix built by
#'     \code{.build_simulated_response()}.
#' }
#'
#' The \code{delta} argument has exactly the same semantics and
#' endpoint rules as in \code{\link{betaregscale_simulate}}.  When
#' \code{delta != NULL}, the returned data frame carries
#' \code{attr(, "bs_prepared") = TRUE}.
#'
#' See \code{\link{betaregscale_simulate}} for the complete
#' documentation of the delta-forced endpoint formulas, the
#' \code{"bs_prepared"} attribute, and the interaction with the
#' fitting pipeline.
#'
#' @param formula_x One-sided formula for the mean model predictors.
#' @param formula_z One-sided formula for the dispersion model
#'   predictors.
#' @param data   Data frame containing the predictor variables.
#' @param beta   Numeric vector of mean-model coefficients.
#' @param zeta   Numeric vector of dispersion-model coefficients.
#' @param link   Mean link function (default \code{"logit"}).
#' @param link_phi Dispersion link function (default \code{"logit"}).
#' @param ncuts  Integer: number of scale categories \eqn{K}
#'   (default 100).
#' @param type   \strong{Deprecated.}
#'   Interval type (default \code{"m"}).
#'   Use \code{\link{bs_prepare}} to control interval geometry instead.
#' @param lim    Numeric: uncertainty half-width \eqn{h}
#'   (default 0.5).
#' @param repar  Integer: reparameterization scheme (default 2).
#' @param delta  Integer or \code{NULL}.  If \code{NULL} (default),
#'   censoring is determined automatically.
#'   If an integer in \code{\{0, 1, 2, 3\}}, all simulated
#'   observations are forced to that censoring type with
#'   observation-specific endpoints.
#'   See \code{\link{betaregscale_simulate}} for the full
#'   specification.
#'
#' @return A \code{data.frame} with columns: \code{left},
#'   \code{right}, \code{yt}, \code{y}, \code{delta}, plus the
#'   predictor columns from the mean and dispersion design matrices.
#'   When \code{delta != NULL}, the data frame carries the attribute
#'   \code{"bs_prepared" = TRUE}.
#'
#' @seealso \code{\link{betaregscale_simulate}} for the full
#'   documentation of delta semantics;
#'   \code{\link{check_response}} for the endpoint formulas.
#'
#' @examples
#' set.seed(42)
#' n <- 200
#' dat <- data.frame(
#'   x1 = rnorm(n), x2 = rnorm(n),
#'   z1 = runif(n), z2 = runif(n)
#' )
#' sim <- betaregscale_simulate_z(
#'   formula_x = ~ x1 + x2, formula_z = ~ z1 + z2,
#'   data = dat,
#'   beta = c(0.2, -0.5, 0.3),
#'   zeta = c(0.2, -0.4, 0.2)
#' )
#' head(sim)
#'
#' @importFrom stats rbeta
#' @export
betaregscale_simulate_z <- function(formula_x = ~ x1 + x2,
                                    formula_z = ~ z1 + z2,
                                    data,
                                    beta = c(0, 0.5, -0.2),
                                    zeta = c(1, 0.5, 0.2),
                                    link = "logit",
                                    link_phi = "logit",
                                    ncuts = 100L,
                                    type = "m",
                                    lim = 0.5,
                                    repar = 2L,
                                    delta = NULL) {
  # Deprecation warning for type
  if (!missing(type)) {
    .Deprecated(msg = paste0(
      "The 'type' argument of betaregscale_simulate_z() is deprecated ",
      "and will be removed in a future version. ",
      "Use bs_prepare() to control interval geometry."
    ))
  }

  # Validate delta
  if (!is.null(delta)) {
    if (length(delta) != 1L || !(delta %in% 0:3)) {
      stop("'delta' must be NULL or a single integer in {0, 1, 2, 3}.",
        call. = FALSE
      )
    }
    delta <- as.integer(delta)
  }

  link <- match.arg(link, .mu_links)
  link_phi <- match.arg(link_phi, .phi_links)
  repar <- as.integer(repar)

  mfx <- stats::model.frame(formula_x, data = data)
  mfz <- stats::model.frame(formula_z, data = data)
  X <- stats::model.matrix(mfx, data = data)
  Z <- stats::model.matrix(mfz, data = data)
  n <- nrow(X)

  # Linear predictors
  mu_x <- apply_inv_link(X %*% beta, link)
  mu_z <- apply_inv_link(Z %*% zeta, link_phi)

  # Reparameterize and simulate
  pars <- beta_reparam(mu = mu_x, phi = mu_z, repar = repar)
  y_raw <- stats::rbeta(n, shape1 = pars$shape1, shape2 = pars$shape2)

  # Build response matrix based on delta
  out_y <- .build_simulated_response(
    y_raw = y_raw, delta = delta, ncuts = ncuts,
    type = type, lim = lim
  )

  result <- data.frame(
    out_y,
    X[, -1L, drop = FALSE],
    Z[, -1L, drop = FALSE]
  )

  # When delta is forced, mark data as pre-processed so that
  # .extract_response() uses the pre-computed columns directly
  if (!is.null(delta)) {
    attr(result, "bs_prepared") <- TRUE
  }

  result
}


# ============================================================================ #
# Internal helper for building simulated response matrices
# ============================================================================ #

#' Build the response matrix for simulated data
#'
#' Internal helper called by \code{\link{betaregscale_simulate}} and
#' \code{\link{betaregscale_simulate_z}} to transform raw simulated
#' beta values into the five-column response matrix.
#'
#' \strong{Design principle}: for every forced censoring type, the
#' actual simulated values (rounded to the scale grid) are preserved
#' so that covariate-driven variation is maintained.  The forced
#' \code{delta} is passed as a vector to
#' \code{\link{check_response}}, which overrides the automatic
#' boundary classification and computes observation-specific
#' endpoints.
#'
#' \strong{Per-delta transformation}:
#' \describe{
#'   \item{\code{delta = NULL}}{Default path:
#'     \eqn{y_{\text{grid}} = \text{round}(y^* \times K)}, then
#'     \code{check_response(y_grid)} with automatic classification.}
#'   \item{\code{delta = 0}}{Exact observations: the continuous
#'     \eqn{y^*} values are used directly on \eqn{(0, 1)}.
#'     \code{check_response()} is \strong{not} called; the matrix
#'     is built manually with \eqn{l_i = u_i = y_t = y^*_i}.}
#'   \item{\code{delta = 1}}{Left-censored:
#'     \eqn{y_{\text{grid}} = \text{round}(y^* \times K)}, then
#'     \code{check_response(y_grid, delta = rep(1, n))}.
#'     The actual grid values are preserved; each observation gets
#'     \eqn{l_i = \epsilon},
#'     \eqn{u_i = (y_{\text{grid}} + h) / K} (non-boundary) or
#'     \eqn{u_i = h / K} (when \eqn{y_{\text{grid}} = 0}).}
#'   \item{\code{delta = 2}}{Right-censored:
#'     \eqn{y_{\text{grid}} = \text{round}(y^* \times K)}, then
#'     \code{check_response(y_grid, delta = rep(2, n))}.
#'     Each observation gets \eqn{u_i = 1 - \epsilon},
#'     \eqn{l_i = (y_{\text{grid}} - h) / K} (non-boundary) or
#'     \eqn{l_i = (K - h) / K} (when \eqn{y_{\text{grid}} = K}).}
#'   \item{\code{delta = 3}}{Interval-censored:
#'     \eqn{y_{\text{grid}} = \text{round}(y^* \times K)}, then
#'     clamped to \eqn{[1, K-1]} (boundaries avoided), then
#'     \code{check_response(y_grid, delta = rep(3, n))}.}
#' }
#'
#' @param y_raw Numeric vector of length \eqn{n}: simulated beta
#'   values on \eqn{(0, 1)}.
#' @param delta Integer scalar or \code{NULL}: forced censoring type
#'   to apply to all observations.
#' @param ncuts Integer: number of scale categories \eqn{K}.
#' @param type  Character: interval type passed to
#'   \code{check_response()} (only relevant for \code{delta = 3}
#'   or \code{delta = NULL}).
#' @param lim   Numeric: uncertainty half-width \eqn{h}.
#' @return A numeric matrix with \eqn{n} rows and columns
#'   \code{left}, \code{right}, \code{yt}, \code{y}, \code{delta}.
#' @noRd
.build_simulated_response <- function(y_raw, delta, ncuts, type, lim) {
  n <- length(y_raw)
  eps <- 1e-5

  if (is.null(delta)) {
    # Default: round to grid and classify automatically
    y_grid <- round(y_raw * ncuts, 0)
    return(check_response(y_grid, type = type, ncuts = ncuts, lim = lim))
  }

  switch(as.character(delta),
    "0" = {
      # Exact (uncensored): use continuous values directly
      yt <- pmin(pmax(y_raw, eps), 1 - eps)
      cbind(
        left  = yt,
        right = yt,
        yt    = yt,
        y     = y_raw,
        delta = rep(0L, n)
      )
    },
    "1" = {
      # Left-censored: keep actual grid values, force delta=1
      y_grid <- round(y_raw * ncuts, 0)
      check_response(y_grid,
        type = type, ncuts = ncuts, lim = lim,
        delta = rep(1L, n)
      )
    },
    "2" = {
      # Right-censored: keep actual grid values, force delta=2
      y_grid <- round(y_raw * ncuts, 0)
      check_response(y_grid,
        type = type, ncuts = ncuts, lim = lim,
        delta = rep(2L, n)
      )
    },
    "3" = {
      # Interval-censored: round to grid, avoid boundaries, force delta=3
      y_grid <- round(y_raw * ncuts, 0)
      y_grid <- pmin(pmax(y_grid, 1L), ncuts - 1L)
      check_response(y_grid,
        type = type, ncuts = ncuts, lim = lim,
        delta = rep(3L, n)
      )
    }
  )
}
