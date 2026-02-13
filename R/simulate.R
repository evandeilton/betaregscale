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
#' Generates interval-censored observations from a beta regression
#' model with a single (scalar) dispersion parameter.  This is useful
#' for Monte Carlo studies and testing.
#'
#' @param formula One-sided formula specifying the mean model
#'   predictors (e.g., \code{~ x1 + x2}).
#' @param data   Data frame containing the predictor variables.
#' @param beta   Numeric vector of regression coefficients (length
#'   must equal the number of columns of the design matrix, including
#'   the intercept).
#' @param phi    Scalar dispersion parameter (on the link scale).
#' @param link   Mean link function (default \code{"logit"}).
#' @param link_phi Dispersion link function (default \code{"logit"}).
#' @param ncuts  Number of scale categories (default 100).
#' @param type   Interval type: \code{"m"}, \code{"l"}, or \code{"r"}
#'   (default \code{"m"}).
#' @param lim    Half-width of the uncertainty region (default 0.5).
#' @param repar  Reparameterization scheme (default 2).
#'
#' @return A \code{data.frame} containing columns \code{left},
#'   \code{right}, \code{yt}, \code{y}, and the predictor variables.
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
#' @importFrom stats rbeta
#' @export
betaregscale_simulate <- function(formula, data, beta, phi = 1 / 5,
                                  link = "logit",
                                  link_phi = "logit",
                                  ncuts = 100L, type = "m",
                                  lim = 0.5, repar = 2L) {
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
  y_grid <- round(y_raw * ncuts, 0)

  out_y <- check_response(y_grid, type = type, ncuts = ncuts, lim = lim)
  data.frame(out_y, X[, -1L, drop = FALSE])
}


#' Simulate data from a variable-dispersion beta interval model
#'
#' @description
#' Generates interval-censored observations from a beta regression
#' model with observation-specific dispersion governed by a second
#' linear predictor.
#'
#' @param formula_x One-sided formula for the mean model predictors.
#' @param formula_z One-sided formula for the dispersion model
#'   predictors.
#' @param data   Data frame containing the predictor variables.
#' @param beta   Numeric vector of mean-model coefficients.
#' @param zeta   Numeric vector of dispersion-model coefficients.
#' @param link   Mean link function (default \code{"logit"}).
#' @param link_phi Dispersion link function (default \code{"logit"}).
#' @param ncuts  Number of scale categories (default 100).
#' @param type   Interval type (default \code{"m"}).
#' @param lim    Uncertainty half-width (default 0.5).
#' @param repar  Reparameterization scheme (default 2).
#'
#' @return A \code{data.frame} with interval endpoints and predictors.
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
                                    repar = 2L) {
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

  # Reparameterize and simulate (no hardcoded seed!)
  pars <- beta_reparam(mu = mu_x, phi = mu_z, repar = repar)
  y_raw <- stats::rbeta(n, shape1 = pars$shape1, shape2 = pars$shape2)
  y_grid <- round(y_raw * ncuts, 0)

  out_y <- check_response(y_grid, type = type, ncuts = ncuts, lim = lim)
  data.frame(
    out_y,
    X[, -1L, drop = FALSE],
    Z[, -1L, drop = FALSE]
  )
}
