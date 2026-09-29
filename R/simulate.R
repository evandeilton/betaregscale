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
#' @param link   Mean link function name (\code{NULL}: default for
#'   \code{repar}).
#' @param link_phi Dispersion link function name (\code{NULL}: default for
#'   \code{repar}).
#' @param ncuts  Number of scale categories.
#' @param lim    Uncertainty half-width.
#' @param repar  Reparameterization scheme. Under \code{repar = 0} the
#'   shapes are started by the method of moments on the midpoint response
#'   (\eqn{f = m(1-m)/v - 1}, \eqn{p = m f}, \eqn{q = (1-m) f}).
#' @param interval Interval direction used when \code{data} is not prepared.
#'
#' @return Named numeric vector of starting values.
#' @keywords internal
compute_start <- function(formula, data, link = NULL,
                          link_phi = NULL, ncuts = 100L,
                          lim = 0.5, repar = 2L, interval = "mid") {
  repar <- as.integer(repar)
  # A partial interval name must not fall through to the right/left cells
  interval <- match.arg(interval, .brs_intervals)
  links <- .resolve_links(link, link_phi, repar)
  link <- links$link
  link_phi <- links$link_phi

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
  Y <- .extract_response(mf, data, ncuts = ncuts, lim = lim, interval = interval)
  x <- stats::model.matrix(mtX, mf)
  z <- stats::model.matrix(mtZ, mf)

  # Midpoint response for starting-value GLM
  y <- rowMeans(Y[, c("left", "right"), drop = FALSE], na.rm = TRUE)

  # Moment estimates on the midpoints (dispersion start; shapes under repar 0)
  y_safe <- pmin(pmax(y, 1e-4), 1 - 1e-4)
  muy    <- mean(y_safe, na.rm = TRUE)
  vy     <- stats::var(y_safe, na.rm = TRUE)
  denom  <- muy * (1 - muy)

  if (repar == 0L) {
    # Shapes by moments (f = m(1-m)/v - 1, p = m f, q = (1-m) f), intercept
    # only: the mean GLM below would start the shape p at the mean.
    f <- if (!is.na(vy) && vy > 0 && denom > 0) denom / vy - 1 else NA_real_
    if (!is.finite(f) || f <= 0) f <- 2
    init_beta <- .brs_intercept_start(x, apply_link(muy * f, link))
    init_phi <- .brs_intercept_start(z, apply_link((1 - muy) * f, link_phi))
    names(init_phi) <- if (ncol(z) < 2L) "phi" else paste0("phi_", colnames(z))
    return(c(init_beta, init_phi))
  }

  # Mean-model starting values via quasi-binomial GLM
  glm_data <- data.frame(y = y, x)
  init_beta <- stats::coef(
    stats::glm(y ~ 0 + .,
      data = glm_data,
      family = stats::quasibinomial(link = link)
    )
  )

  # Dispersion starting values
  # PERF-M04: moment-based estimate avoids a second GLM call.
  # repar=1 precision: Var(Y) = mu*(1-mu)/(phi+1) => phi = mu*(1-mu)/Var(Y) - 1
  # repar=2 mean-var:  Var(Y) = phi*mu*(1-mu)     => phi = Var(Y)/(mu*(1-mu))
  phi_start <- if (!is.na(vy) && vy > 0 && denom > 0) {
    if (repar == 1L) max(denom / vy - 1, 1.0)
    else min(max(vy / denom, 0.05), 0.9)
  } else {
    if (repar == 2L) 0.3 else 5.0
  }
  # Clamp phi_start to the valid domain of the forward link to avoid NaN.
  phi_start_clamped <- if (link_phi %in% c("logit", "probit", "cauchit", "cloglog")) {
    min(max(phi_start, 1e-6), 1 - 1e-6)
  } else {
    max(phi_start, 1e-6)
  }
  init_phi0 <- apply_link(phi_start_clamped, link_phi)
  if (!is.finite(init_phi0)) init_phi0 <- 0

  if (is.null(z) || ncol(z) < 2L) {
    init_phi <- init_phi0
    names(init_phi) <- "phi"
  } else {
    # A GLM of the mean says nothing about the dispersion/precision (all right-
    # censored data gave logit phi = 5.99): moment intercept + zero slopes.
    init_phi <- .brs_intercept_start(z, init_phi0)
    names(init_phi) <- paste0("phi_", colnames(z))
  }

  c(init_beta, init_phi)
}

# Intercept-only starting vector for a design matrix: `value` on the
# "(Intercept)" column (if any), zero elsewhere.
.brs_intercept_start <- function(M, value) {
  start <- rep(0, ncol(M))
  names(start) <- colnames(M)
  j <- match("(Intercept)", colnames(M))
  if (!is.na(j)) start[j] <- value
  start
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
#' The score is generated by the same mechanism the likelihood assumes for
#' the chosen \code{interval}: \eqn{s = \mathrm{round}(K y^*)} under
#' \code{"mid"} (cells centred on the score; exact only for
#' \code{lim = 0.5}) and \eqn{s = \lfloor (K + 1) y^* \rfloor} under
#' \code{"right"}/\code{"left"} (\eqn{K + 1} equal cells), then
#' \code{\link{brs_check}} with the same \code{interval}, so the returned
#' data are exactly what \code{\link{brs}} will see.
#'
#' @param formula Model formula with one (mean) or two parts
#'   (mean \code{|} precision). A left-hand-side response is allowed but ignored.
#' @param data Data frame with predictor variables.
#' @param beta Numeric vector of mean-model coefficients.
#' @param phi Scalar dispersion parameter (link scale), used only for one-part
#'   formulas.
#' @param zeta Numeric vector of precision-model coefficients (link scale),
#'   required for two-part formulas.
#' @param link Link for the first parameter; \code{NULL} (default) selects
#'   the link implied by \code{repar} (see \code{\link{brs}}).
#' @param link_phi Link for the second parameter; \code{NULL} (default)
#'   selects the link implied by \code{repar}.
#' @param ncuts Number of scale categories.
#' @param lim Half-width used in interval construction.
#' @param repar Reparameterization scheme. \code{beta} and \code{phi} /
#'   \code{zeta} are on the link scale of the first and second parameter of
#'   that scheme (shapes \eqn{p, q} under \code{repar = 0}).
#' @param delta Forced censoring type (\code{0,1,2,3}) or \code{NULL}.
#'   \code{delta = 1} or \code{2} censors every observation on the same side
#'   at the cell of its own value: the threshold depends on \eqn{Y}
#'   (informative censoring) and \code{\link{brs}} has no finite MLE for such
#'   data (the estimates diverge while optim may report convergence). A
#'   warning is issued, as it is whenever all simulated observations end up
#'   censored on the same side. Use it only to exercise code paths.
#' @param interval Direction of the uncertainty interval, \code{"mid"}
#'   (default), \code{"right"} or \code{"left"}; see \code{\link{brs_check}}.
#'
#' @return A data frame with columns \code{left}, \code{right}, \code{yt},
#'   \code{y}, \code{delta}, plus simulated predictor columns from the model
#'   matrices. As for \code{\link{brs_prep}}, the output carries the
#'   attributes \code{"is_prepared"} (\code{TRUE}), \code{"ncuts"},
#'   \code{"lim"} and \code{"interval"}, which \code{\link{brs}} reuses.
#'
#' @examples
#' \donttest{
#' dat <- data.frame(
#'   x1 = rep(c(1, 2), 10),
#'   x2 = rep(c(0, 0, 1, 1), 5)
#' )
#' # Fixed dispersion
#' sim_fixed <- brs_sim(
#'   formula = ~ x1 + x2, data = dat,
#'   beta = c(0.2, -0.5, 0.3), phi = 1 / 5
#' )
#' # Variable dispersion
#' sim_var <- brs_sim(
#'   formula = ~ x1 | x2, data = dat,
#'   beta = c(0.2, -0.5), zeta = c(0.5, -0.5)
#' )
#' }
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
#' @importFrom stats rbeta
#' @rdname brs_sim
#' @export
brs_sim <- function(formula,
                    data,
                    beta,
                    phi = 1 / 5,
                    zeta = NULL,
                    link = NULL,
                    link_phi = NULL,
                    ncuts = 100L,
                    lim = 0.5,
                    repar = 2L,
                    delta = NULL,
                    interval = c("mid", "right", "left")) {
  if (!is.null(delta)) {
    if (length(delta) != 1L || !(delta %in% 0:3)) {
      stop("'delta' must be NULL or a single integer in {0, 1, 2, 3}.",
        call. = FALSE
      )
    }
    delta <- as.integer(delta)
  }
  interval <- match.arg(interval)
  # Generator-specific warning: round() matches the mid likelihood only at 0.5
  .brs_lim_check(lim, interval, warn = TRUE, generator = TRUE)

  repar <- as.integer(repar)
  links <- .resolve_links(link, link_phi, repar)
  link <- links$link
  link_phi <- links$link_phi

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
    y_raw = y_raw, delta = delta, ncuts = ncuts, lim = lim, interval = interval
  )
  # All rows censored on one side: no finite MLE (forced: also informative)
  d_out <- as.integer(out_y[, "delta"])
  if (!is.null(delta) && delta %in% c(1L, 2L)) {
    warning("delta = ", delta, " censors every observation on the same side at ",
            "the cell of its own value (informative censoring): brs() has no ",
            "finite MLE for these data (estimates diverge).", call. = FALSE)
  } else if (n > 0L && (all(d_out == 1L) || all(d_out == 2L))) {
    warning("Every simulated observation is censored on the same side (delta = ",
            d_out[1L], "): brs() has no finite MLE for these data.", call. = FALSE)
  }

  # Drop the intercept by name: `0 + x` has no intercept column to drop
  no_int <- function(M) M[, colnames(M) != "(Intercept)", drop = FALSE]
  predictors <- if (is.null(Z)) no_int(X) else cbind(no_int(X), no_int(Z))

  result <- data.frame(out_y, predictors)

  # Same attributes as brs_prep(): brs() reuses the columns and ncuts/lim/interval.
  attr(result, "is_prepared") <- TRUE
  attr(result, "ncuts") <- ncuts
  attr(result, "lim") <- lim
  attr(result, "interval") <- interval

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
                        link = NULL,
                        link_phi = NULL,
                        ncuts = 100L,
                        lim = 0.5,
                        repar = 2L,
                        delta = NULL,
                        interval = "mid") {
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
    delta = delta,
    interval = interval
  )
}


# ============================================================================ #
# Internal helper for building simulated response matrices
# ============================================================================ #

#' Build the response matrix for simulated data
#'
#' Internal helper called by \code{\link{brs_sim}}: the simulated
#' \eqn{y^*} is coarsened to the score grid by the mechanism of the chosen
#' \code{interval} (\code{.brs_score_from_unit()}) and passed to the same
#' cell mapping as \code{\link{brs_check}}. A forced \code{delta} keeps the
#' grid values (covariate-driven variation) and overrides the censoring
#' type; \code{delta = 0} keeps the continuous \eqn{y^*} as exact values;
#' \code{delta = 3} avoids the border scores.
#'
#' @param y_raw Numeric vector of length \eqn{n}: simulated beta
#'   values on \eqn{(0, 1)}.
#' @param delta Integer scalar or \code{NULL}: forced censoring type
#'   to apply to all observations.
#' @param ncuts Integer: number of scale categories \eqn{K}.
#' @param lim   Numeric: uncertainty half-width \eqn{h}.
#' @param interval Interval direction.
#' @return A numeric matrix with \eqn{n} rows and columns
#'   \code{left}, \code{right}, \code{yt}, \code{y}, \code{delta}.
#' @noRd
.build_simulated_response <- function(y_raw, delta, ncuts, lim, interval = "mid") {
  n <- length(y_raw)
  eps <- 1e-5

  if (!is.null(delta) && delta == 0L) {
    # Exact (uncensored): use continuous values directly
    yt <- pmin(pmax(y_raw, eps), 1 - eps)
    return(cbind(left = yt, right = yt, yt = yt, y = y_raw, delta = rep(0L, n)))
  }

  # Score by the likelihood's own coarsening (round for mid, floor over K + 1
  # cells for right/left); brs_sim() already validated lim, so no warnings
  y_grid <- .brs_score_from_unit(y_raw, ncuts, interval)
  if (!is.null(delta) && delta == 3L) {
    # Interval-censored everywhere: keep away from the border scores
    y_grid <- pmin(pmax(y_grid, 1L), ncuts - 1L)
  }
  forced <- if (is.null(delta)) NULL else rep(delta, n)
  .brs_check_core(y_grid, ncuts = ncuts, lim = lim, delta = forced,
                  interval = interval)
}
