# ============================================================================ #
# Marginal effects
# ============================================================================ #

#' Marginal effects for brs models
#'
#' @description
#' Computes average marginal effects (AME) for numeric covariates in the
#' mean or precision submodel of a fitted \code{"brs"} object.
#'
#' @details
#' AMEs are computed by finite differences on predictions:
#' \deqn{
#' \mathrm{AME}_j = \frac{1}{n}\sum_{i=1}^{n}
#' \frac{\hat{g}_i(x_{ij} + h) - \hat{g}_i(x_{ij})}{h},
#' }
#' where \eqn{\hat{g}_i} is the selected prediction scale.
#'
#' For binary covariates coded as \code{0/1}, the effect is computed as the
#' average discrete difference \eqn{\hat{g}(x_j=1)-\hat{g}(x_j=0)}.
#'
#' If \code{interval = TRUE}, uncertainty is approximated by asymptotic
#' parameter simulation from \eqn{\mathcal{N}(\hat{\theta}, \hat{V})}.
#'
#' @param object A fitted \code{"brs"} object.
#' @param newdata Optional data frame for evaluation; defaults to the data
#'   used in fitting.
#' @param model Character; \code{"mean"} (default) or \code{"precision"}.
#' @param type Character prediction scale:
#'   \code{"response"} (default) or \code{"link"}.
#' @param variables Optional character vector of covariate names.
#'   Defaults to all numeric covariates in the selected submodel.
#' @param h Finite-difference step for non-binary numeric covariates.
#' @param interval Logical; compute interval estimates via simulation.
#' @param level Confidence level for interval estimates.
#' @param n_sim Number of parameter draws when \code{interval = TRUE}.
#' @param keep_draws Logical; if \code{TRUE} and \code{interval = TRUE},
#'   stores AME simulation draws in attribute \code{"ame_draws"}.
#'
#' @return A data frame with one row per variable and columns:
#'   \code{variable}, \code{ame}, \code{std.error}, \code{ci.lower},
#'   \code{ci.upper}, \code{model}, \code{type}, and \code{n}.
#'   The returned object has class \code{"brs_marginaleffects"} and
#'   attributes with analysis metadata.
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
#' \donttest{
#' dat <- data.frame(
#'   y = c(
#'     0, 5, 20, 50, 75, 90, 100, 30, 60, 45,
#'     10, 40, 55, 70, 85, 25, 35, 65, 80, 15
#'   ),
#'   x1 = rep(c(1, 2), 10),
#'   x2 = rep(c(0, 0, 1, 1), 5)
#' )
#' prep <- brs_prep(dat, ncuts = 100)
#' fit <- brs(y ~ x1, data = prep)
#' brs_marginaleffects(fit, model = "mean", type = "response")
#' }
#'
#' @rdname brs_marginaleffects
#' @export
brs_marginaleffects <- function(object,
                                newdata = NULL,
                                model = c("mean", "precision"),
                                type = c("response", "link"),
                                variables = NULL,
                                h = 1e-5,
                                interval = TRUE,
                                level = 0.95,
                                n_sim = 400L,
                                keep_draws = FALSE) {
  .check_class(object)
  model <- match.arg(model)
  type <- match.arg(type)
  interval <- isTRUE(interval)
  keep_draws <- isTRUE(keep_draws)

  h <- as.numeric(h)
  if (!is.finite(h) || h <= 0) {
    stop("'h' must be a positive number.", call. = FALSE)
  }
  level <- as.numeric(level)
  if (length(level) != 1L || !is.finite(level) || level <= 0 || level >= 1) {
    stop("'level' must be in (0, 1).", call. = FALSE)
  }
  n_sim <- as.integer(n_sim)
  if (interval && (!is.finite(n_sim) || n_sim < 50L)) {
    stop("'n_sim' must be >= 50 when interval = TRUE.", call. = FALSE)
  }

  eval_data <- if (is.null(newdata)) object$data else newdata
  if (!is.data.frame(eval_data)) {
    stop("'newdata' must be a data.frame.", call. = FALSE)
  }

  if (identical(model, "precision") && object$q <= 1L) {
    return(data.frame(
      variable = character(0),
      ame = numeric(0),
      std.error = numeric(0),
      ci.lower = numeric(0),
      ci.upper = numeric(0),
      model = character(0),
      type = character(0),
      n = integer(0)
    ))
  }

  vars_model <- .brs_me_model_vars(object, model = model)
  vars_numeric <- vars_model[vars_model %in% names(eval_data) &
    vapply(eval_data[vars_model], is.numeric, logical(1))]

  if (is.null(variables)) {
    variables <- vars_numeric
  } else {
    variables_in <- as.character(variables)
    variables <- intersect(variables_in, vars_numeric)
    dropped <- setdiff(variables_in, variables)
    if (length(dropped) > 0L) {
      warning(
        "Ignoring variables not available as numeric covariates in '", model, "' model: ",
        paste(dropped, collapse = ", "),
        call. = FALSE
      )
    }
  }
  if (length(variables) == 0L) {
    stop("No numeric covariates available for marginal effects.", call. = FALSE)
  }

  par_hat <- unname(object$par)
  V <- vcov(object, model = "full")
  alpha <- 1 - level

  .ame_one <- function(par, var_name) {
    base <- .brs_me_predict(object, eval_data, par = par, model = model, type = type)
    x <- eval_data[[var_name]]
    x_non_na <- stats::na.omit(x)
    if (length(x_non_na) == 0L) {
      return(NA_real_)
    }

    if (all(x_non_na %in% c(0, 1))) {
      d0 <- eval_data
      d1 <- eval_data
      d0[[var_name]] <- 0
      d1[[var_name]] <- 1
      p0 <- .brs_me_predict(object, d0, par = par, model = model, type = type)
      p1 <- .brs_me_predict(object, d1, par = par, model = model, type = type)
      return(mean(p1 - p0, na.rm = TRUE))
    }

    # Scale the perturbation to improve numerical stability across covariates.
    h_var <- h * max(stats::sd(x_non_na), 1)
    dp <- eval_data
    dm <- eval_data
    dp[[var_name]] <- dp[[var_name]] + h_var
    dm[[var_name]] <- dm[[var_name]] - h_var
    p_plus <- .brs_me_predict(object, dp, par = par, model = model, type = type)
    p_minus <- .brs_me_predict(object, dm, par = par, model = model, type = type)
    mean((p_plus - p_minus) / (2 * h_var), na.rm = TRUE)
  }

  draws <- NULL
  if (interval) {
    draws <- .brs_me_rmvnorm(n = n_sim, mu = par_hat, sigma = V)
  }

  out <- lapply(variables, function(v) {
    ame_hat <- tryCatch(
      .ame_one(par_hat, v),
      error = function(e) NA_real_
    )

    if (!interval) {
      return(data.frame(
        variable = v,
        ame = ame_hat,
        std.error = NA_real_,
        ci.lower = NA_real_,
        ci.upper = NA_real_,
        model = model,
        type = type,
        n = nrow(eval_data),
        stringsAsFactors = FALSE
      ))
    }

    ame_draws <- apply(draws, 1L, function(th) {
      tryCatch(.ame_one(th, v), error = function(e) NA_real_)
    })
    se <- stats::sd(ame_draws, na.rm = TRUE)

    data.frame(
      variable = v,
      ame = ame_hat,
      std.error = se,
      ci.lower = as.numeric(stats::quantile(ame_draws, probs = alpha / 2, na.rm = TRUE)),
      ci.upper = as.numeric(stats::quantile(ame_draws, probs = 1 - alpha / 2, na.rm = TRUE)),
      model = model,
      type = type,
      n = nrow(eval_data),
      stringsAsFactors = FALSE
    )
  })

  res <- do.call(rbind, out)
  attr(res, "level") <- level
  attr(res, "interval") <- interval
  attr(res, "n_sim") <- if (interval) n_sim else 0L
  attr(res, "model") <- model
  attr(res, "type") <- type
  attr(res, "h") <- h
  if (interval && keep_draws) {
    ame_mat <- sapply(variables, function(v) {
      apply(draws, 1L, function(th) {
        tryCatch(.ame_one(th, v), error = function(e) NA_real_)
      })
    })
    if (is.vector(ame_mat)) {
      ame_mat <- matrix(ame_mat, ncol = 1L)
      colnames(ame_mat) <- variables[1L]
    } else {
      colnames(ame_mat) <- variables
    }
    attr(res, "ame_draws") <- ame_mat
  }
  class(res) <- c("brs_marginaleffects", "data.frame")
  res
}


#' @keywords internal
.brs_me_model_vars <- function(object, model = c("mean", "precision")) {
  model <- match.arg(model)
  if (identical(model, "mean")) {
    tm <- stats::delete.response(object$terms$mean)
    return(unique(all.vars(tm)))
  }
  if (is.null(object$terms$precision)) {
    return(character(0))
  }
  unique(all.vars(object$terms$precision))
}

#' @keywords internal
.brs_me_predict <- function(object, data, par, model, type) {
  mt_mu <- stats::delete.response(object$terms$mean)
  mf_mu <- stats::model.frame(mt_mu, data = data)
  X <- stats::model.matrix(mt_mu, mf_mu)

  beta <- par[seq_len(object$p)]
  eta_mu <- as.numeric(X %*% beta)
  mu <- apply_inv_link(eta_mu, object$link)

  if (object$q > 1L && !is.null(object$terms$precision)) {
    mf_z <- stats::model.frame(object$terms$precision, data = data)
    Z <- stats::model.matrix(object$terms$precision, mf_z)
    zeta <- par[object$p + seq_len(object$q)]
    eta_phi <- as.numeric(Z %*% zeta)
  } else {
    eta_phi <- rep(par[object$p + 1L], length(mu))
  }
  phi <- apply_inv_link(eta_phi, object$link_phi)

  if (identical(model, "mean")) {
    return(if (identical(type, "response")) mu else eta_mu)
  }
  if (identical(type, "response")) phi else eta_phi
}

#' @keywords internal
.brs_me_rmvnorm <- function(n, mu, sigma) {
  p <- length(mu)
  sigma <- (sigma + t(sigma)) / 2
  ee <- eigen(sigma, symmetric = TRUE)
  vals <- pmax(ee$values, 0)
  A <- diag(sqrt(vals), nrow = p) %*% t(ee$vectors)
  Z <- matrix(stats::rnorm(n * p), nrow = n, ncol = p)
  sweep(Z %*% A, 2L, mu, FUN = "+")
}
