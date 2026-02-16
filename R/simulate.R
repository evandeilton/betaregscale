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
#' @param lim    Uncertainty half-width.
#'
#' @return Named numeric vector of starting values.
#' @keywords internal
compute_start <- function(formula, data, link = "logit",
                          link_phi = "logit", ncuts = 100L,
                          lim = 0.5) {
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
  Y <- .extract_response(mf, data, ncuts = ncuts, lim = lim)
  x <- stats::model.matrix(mtX, mf)
  z <- stats::model.matrix(mtZ, mf)

  # Midpoint response for starting-value GLM
  y <- rowMeans(Y[, c("left", "right"), drop = FALSE], na.rm = TRUE)

  # Mean-model starting values via quasi-binomial GLM
  # Mean-model starting values via quasi-binomial GLM
  glm_data <- data.frame(y = y, x)
  init_beta <- stats::coef(
    stats::glm(y ~ 0 + .,
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

#' Build simulation design matrices from one- or two-part formulas
#' @keywords internal
#' @noRd
.sim_design_matrices <- function(formula, data) {
  formula <- Formula::as.Formula(formula)

  # Accept y ~ ... style by stripping the response for simulation.
  if (length(formula)[1L] > 0L) {
    formula <- Formula::Formula(
      formula(formula, lhs = 0L, rhs = seq_len(length(formula)[2L]))
    )
  }

  # Keep at most two right-hand-side parts (mean | precision).
  if (length(formula)[2L] > 2L) {
    formula <- Formula::Formula(formula(formula, rhs = 1:2))
  }

  mf <- stats::model.frame(formula, data = data)
  mtX <- stats::terms(formula, data = data, rhs = 1L)
  X <- stats::model.matrix(mtX, mf)

  if (length(formula)[2L] < 2L) {
    return(list(X = X, Z = NULL))
  }

  mtZ <- stats::delete.response(stats::terms(formula, data = data, rhs = 2L))
  Z <- stats::model.matrix(mtZ, mf)
  list(X = X, Z = Z)
}


#' Simulate data from beta interval models
#'
#' @description
#' Simulates interval-censored responses from fixed- or variable-dispersion
#' beta regression models.
#'
#' @details
#' The model structure is controlled by \code{formula} in the same style as
#' \code{\link{brs}}:
#' \itemize{
#'   \item one-part formula (\code{~ x1 + x2} or \code{y ~ x1 + x2}):
#'     fixed dispersion using scalar \code{phi}.
#'   \item two-part formula (\code{~ x1 + x2 | z1} or \code{y ~ x1 + x2 | z1}):
#'     variable dispersion using coefficient vector \code{zeta}.
#' }
#'
#' The \code{delta} argument can force a single censoring type
#' (\code{0,1,2,3}) for all observations; otherwise, censoring is classified
#' automatically from simulated scale values via \code{\link{brs_check}}.
#'
#' @param formula Model formula with one (mean) or two parts
#'   (mean \code{|} precision). A left-hand-side response is allowed but ignored.
#' @param data Data frame with predictor variables.
#' @param beta Numeric vector of mean-model coefficients.
#' @param phi Scalar dispersion parameter (link scale), used only for one-part
#'   formulas.
#' @param zeta Numeric vector of precision-model coefficients (link scale),
#'   required for two-part formulas.
#' @param link Mean link function.
#' @param link_phi Precision link function.
#' @param ncuts Number of scale categories.
#' @param lim Half-width used in interval construction.
#' @param repar Reparameterization scheme.
#' @param delta Forced censoring type (\code{0,1,2,3}) or \code{NULL}.
#'
#' @return A data frame with columns \code{left}, \code{right}, \code{yt},
#'   \code{y}, \code{delta}, plus simulated predictor columns from the model
#'   matrices. When \code{delta != NULL}, the output carries
#'   \code{attr(, "is_prepared") = TRUE}.
#'
#' @examples
#' set.seed(42)
#' n <- 200
#' dat <- data.frame(x1 = rnorm(n), x2 = rnorm(n), z1 = runif(n))
#'
#' # Fixed dispersion
#' sim_fixed <- brs_sim(
#'   formula = ~ x1 + x2, data = dat,
#'   beta = c(0.2, -0.5, 0.3), phi = 1 / 5
#' )
#'
#' # Variable dispersion
#' sim_var <- brs_sim(
#'   formula = ~ x1 + x2 | z1, data = dat,
#'   beta = c(0.2, -0.5, 0.3), zeta = c(0.5, -0.5)
#' )
#'
#' @references
#' Hawker, G. A., Mian, S., Kendzerska, T., and French, M. (2011).
#' Measures of adult pain: Visual Analog Scale for Pain (VAS Pain),
#' Numeric Rating Scale for Pain (NRS Pain), McGill Pain Questionnaire (MPQ),
#' Short-Form McGill Pain Questionnaire (SF-MPQ), Chronic Pain Grade Scale
#' (CPGS), Short Form-36 Bodily Pain Scale (SF-36 BPS), and Measure of
#' Intermittent and Constant Osteoarthritis Pain (ICOAP).
#' Arthritis Care and Research, 63(S11), S240-S252.
#' doi:10.1002/acr.20543.
#'
#' Hjermstad, M. J., Fayers, P. M., Haugen, D. F., et al. (2011).
#' Studies comparing Numerical Rating Scales, Verbal Rating Scales, and
#' Visual Analogue Scales for assessment of pain intensity in adults:
#' a systematic literature review.
#' Journal of Pain and Symptom Management, 41(6), 1073-1093.
#' doi:10.1016/j.jpainsymman.2010.08.016.
#'
#' @importFrom stats rbeta
#' @rdname brs_sim
#' @export
brs_sim <- function(formula,
                    data,
                    beta,
                    phi = 1 / 5,
                    zeta = NULL,
                    link = "logit",
                    link_phi = "logit",
                    ncuts = 100L,
                    lim = 0.5,
                    repar = 2L,
                    delta = NULL) {
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

  design <- .sim_design_matrices(formula, data)
  X <- design$X
  Z <- design$Z
  n <- nrow(X)

  if (length(beta) != ncol(X)) {
    stop(
      "'beta' length (", length(beta), ") must equal ncol(X) (",
      ncol(X), ").",
      call. = FALSE
    )
  }

  mu <- apply_inv_link(drop(X %*% beta), link)

  if (is.null(Z)) {
    if (!is.null(zeta)) {
      stop("'zeta' is only valid for formulas with '|'.", call. = FALSE)
    }
    if (length(phi) != 1L) {
      stop("'phi' must be a scalar for one-part formulas.", call. = FALSE)
    }
    phi_vec <- rep_len(apply_inv_link(phi, link_phi), n)
  } else {
    if (is.null(zeta)) {
      stop(
        "For formulas with '|', provide 'zeta' for the precision submodel.",
        call. = FALSE
      )
    }
    if (length(zeta) != ncol(Z)) {
      stop(
        "'zeta' length (", length(zeta), ") must equal ncol(Z) (",
        ncol(Z), ").",
        call. = FALSE
      )
    }
    phi_vec <- apply_inv_link(drop(Z %*% zeta), link_phi)
  }

  pars <- brs_repar(mu = mu, phi = phi_vec, repar = repar)
  y_raw <- stats::rbeta(n, shape1 = pars$shape1, shape2 = pars$shape2)

  out_y <- .build_simulated_response(
    y_raw = y_raw, delta = delta, ncuts = ncuts, lim = lim
  )

  predictors <- if (is.null(Z)) {
    X[, -1L, drop = FALSE]
  } else {
    cbind(X[, -1L, drop = FALSE], Z[, -1L, drop = FALSE])
  }

  result <- data.frame(out_y, predictors)

  if (!is.null(delta)) {
    attr(result, "is_prepared") <- TRUE
  }

  result
}


# Backward-compatibility wrapper (internal).
#' @keywords internal
#' @noRd
brs_sim_var <- function(formula_x = ~ x1 + x2,
                        formula_z = ~ z1 + z2,
                        data,
                        beta = c(0, 0.5, -0.2),
                        zeta = c(1, 0.5, 0.2),
                        link = "logit",
                        link_phi = "logit",
                        ncuts = 100L,
                        lim = 0.5,
                        repar = 2L,
                        delta = NULL) {
  brs_sim(
    formula = Formula::Formula(formula_x, formula_z),
    data = data,
    beta = beta,
    zeta = zeta,
    link = link,
    link_phi = link_phi,
    ncuts = ncuts,
    lim = lim,
    repar = repar,
    delta = delta
  )
}


# ============================================================================ #
# Internal helper for building simulated response matrices
# ============================================================================ #

#' Build the response matrix for simulated data
#'
#' Internal helper called by \code{\link{brs_sim}} to transform raw simulated
#' beta values into the five-column response matrix.
#'
#' \strong{Design principle}: for every forced censoring type, the
#' actual simulated values (rounded to the scale grid) are preserved
#' so that covariate-driven variation is maintained.  The forced
#' \code{delta} is passed as a vector to
#' \code{\link{brs_check}}, which overrides the automatic
#' boundary classification and computes observation-specific
#' endpoints.
#'
#' \strong{Per-delta transformation}:
#' \describe{
#'   \item{\code{delta = NULL}}{Default path:
#'     \eqn{y_{\text{grid}} = \text{round}(y^* \times K)}, then
#'     \code{brs_check(y_grid)} with automatic classification.}
#'   \item{\code{delta = 0}}{Exact observations: the continuous
#'     \eqn{y^*} values are used directly on \eqn{(0, 1)}.
#'     \code{brs_check()} is \strong{not} called; the matrix
#'     is built manually with \eqn{l_i = u_i = y_t = y^*_i}.}
#'   \item{\code{delta = 1}}{Left-censored:
#'     \eqn{y_{\text{grid}} = \text{round}(y^* \times K)}, then
#'     \code{brs_check(y_grid, delta = rep(1, n))}.
#'     The actual grid values are preserved; each observation gets
#'     \eqn{l_i = \epsilon},
#'     \eqn{u_i = (y_{\text{grid}} + h) / K} (non-boundary) or
#'     \eqn{u_i = h / K} (when \eqn{y_{\text{grid}} = 0}).}
#'   \item{\code{delta = 2}}{Right-censored:
#'     \eqn{y_{\text{grid}} = \text{round}(y^* \times K)}, then
#'     \code{brs_check(y_grid, delta = rep(2, n))}.
#'     Each observation gets \eqn{u_i = 1 - \epsilon},
#'     \eqn{l_i = (y_{\text{grid}} - h) / K} (non-boundary) or
#'     \eqn{l_i = (K - h) / K} (when \eqn{y_{\text{grid}} = K}).}
#'   \item{\code{delta = 3}}{Interval-censored:
#'     \eqn{y_{\text{grid}} = \text{round}(y^* \times K)}, then
#'     clamped to \eqn{[1, K-1]} (boundaries avoided), then
#'     \code{brs_check(y_grid, delta = rep(3, n))}.}
#' }
#'
#' @param y_raw Numeric vector of length \eqn{n}: simulated beta
#'   values on \eqn{(0, 1)}.
#' @param delta Integer scalar or \code{NULL}: forced censoring type
#'   to apply to all observations.
#' @param ncuts Integer: number of scale categories \eqn{K}.
#' @param lim   Numeric: uncertainty half-width \eqn{h}.
#' @return A numeric matrix with \eqn{n} rows and columns
#'   \code{left}, \code{right}, \code{yt}, \code{y}, \code{delta}.
#' @noRd
.build_simulated_response <- function(y_raw, delta, ncuts, lim) {
  n <- length(y_raw)
  eps <- 1e-5

  if (is.null(delta)) {
    # Default: round to grid and classify automatically
    y_grid <- round(y_raw * ncuts, 0)
    return(brs_check(y_grid, ncuts = ncuts, lim = lim))
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
      brs_check(y_grid,
        ncuts = ncuts, lim = lim,
        delta = rep(1L, n)
      )
    },
    "2" = {
      # Right-censored: keep actual grid values, force delta=2
      y_grid <- round(y_raw * ncuts, 0)
      brs_check(y_grid,
        ncuts = ncuts, lim = lim,
        delta = rep(2L, n)
      )
    },
    "3" = {
      # Interval-censored: round to grid, avoid boundaries, force delta=3
      y_grid <- round(y_raw * ncuts, 0)
      y_grid <- pmin(pmax(y_grid, 1L), ncuts - 1L)
      brs_check(y_grid,
        ncuts = ncuts, lim = lim,
        delta = rep(3L, n)
      )
    }
  )
}
