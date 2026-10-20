# ============================================================================ #
# Log-likelihood contributions per observation (brs) or per group (brsmm)
# ============================================================================ #

#' Log-likelihood contributions per observation or per group
#'
#' @description
#' The terms whose sum is the log-likelihood: one per observation for a
#' \code{\link{brs}} fit, and one per group (cluster) for a
#' \code{\link{brsmm}} fit, the marginal log-likelihood of the group's
#' observations with its random effects integrated out. With
#' \code{newdata}, the contributions of new observations at the fitted
#' parameters (held-out scoring).
#'
#' @details
#' For a \code{brs} fit the contribution of a score row is
#' \eqn{\log P(S = s_i)}, the log probability of its cell (the log density for
#' exact observations, \eqn{\delta = 0}), computed as in the compiled
#' likelihood. For a \code{brsmm} fit the contribution of group \eqn{g} is
#' \deqn{\log \int \prod_i P(S = s_{gi} \mid b)\, \varphi(b; 0, D)\, db,}
#' computed by the fit's \code{int_method} and \code{n_points} (or
#' \code{qmc_points}) unless \code{int_method} / \code{n_points} are given.
#' With the defaults and \code{newdata = NULL}, the contributions add up to
#' \code{logLik(object)}.
#'
#' With \code{newdata}, rows need the response, the covariates and (for
#' \code{brsmm}) the grouping variable; rows with missing values in them are
#' dropped. The groups need not be groups of the fit: their random effects
#' are integrated out, so no random-effect estimate is used. Scores are
#' coarsened with the cells of the fit (its estimated cutpoints under
#' \code{cutpoints = "boundary"}). For score data the contributions are log
#' probabilities of the observed scores (no Jacobian term), on the same scale
#' as the log-likelihood of any model that gives a probability mass function
#' on \eqn{0, \ldots, K}.
#'
#' @param object A fitted \code{"brs"} or \code{"brsmm"} object.
#' @param newdata Optional data frame of new observations (see Details).
#' @param int_method,n_points For \code{brsmm} fits: integration rule
#'   (\code{"laplace"}, \code{"aghq"} or \code{"qmc"}) and number of nodes
#'   (nodes per dimension for \code{"aghq"}, points for \code{"qmc"});
#'   \code{NULL} (default) uses those of the fit. Ignored for \code{brs} fits.
#'
#' @return Named numeric vector: by row name of the data (\code{brs}) or by
#'   group level (\code{brsmm}, with attribute \code{"n"}, the number of
#'   observations of each group).
#'
#' @seealso \code{\link{logLik.brsmm}}, \code{\link{brsmm_cv}},
#'   \code{\link{brs_marginal}}
#'
#' @examples
#' set.seed(21)
#' nrs <- expand.grid(id = 1:30, time = c("6h", "12h", "24h"))
#' b0 <- rnorm(30, sd = 0.8)
#' eta <- -1.3 + c(0, 0.75, 0.3)[nrs$time] + b0[nrs$id]
#' shp <- brs_repar(mu = plogis(eta), phi = 0.2, repar = 2)
#' nrs$y <- round(10 * rbeta(nrow(nrs), shp$shape1, shp$shape2))
#' m <- brsmm(y ~ time, random = ~ 1 | id, data = nrs, ncuts = 10)
#'
#' ll_g <- brs_loglik_contrib(m)
#' head(ll_g)
#' c(sum = sum(ll_g), logLik = as.numeric(logLik(m)))
#' # Same groups, adaptive Gauss-Hermite with 31 nodes
#' head(brs_loglik_contrib(m, int_method = "aghq", n_points = 31))
#'
#' @export
brs_loglik_contrib <- function(object, newdata = NULL, int_method = NULL,
                               n_points = NULL) {
  if (!inherits(object, c("brs", "brsmm"))) {
    stop("'object' must be a \"brs\" or \"brsmm\" fit.", call. = FALSE)
  }
  dd <- .brs_contrib_data(object, newdata)
  Y <- dd$Y
  delta <- as.integer(Y[, "delta"])
  p <- object$p
  par <- unname(object$par)

  if (!inherits(object, "brsmm")) {
    eta_mu <- as.numeric(dd$X %*% par[seq_len(p)])
    eta_phi <- if (is.null(dd$Z)) {
      rep(par[p + 1L], length(eta_mu))
    } else {
      as.numeric(dd$Z %*% par[p + seq_len(ncol(dd$Z))])
    }
    sh <- .brs_shapes_cpp(eta_mu, eta_phi, object$link, object$link_phi, object$repar)
    out <- .brs_obs_loglik(delta, as.numeric(Y[, "left"]), as.numeric(Y[, "right"]),
                           as.numeric(Y[, "yt"]), sh$a, sh$b)
    names(out) <- rownames(Y)
    return(out)
  }

  int_method <- if (is.null(int_method)) {
    object$int_method
  } else {
    match.arg(int_method, c("laplace", "aghq", "qmc"))
  }
  if (is.null(n_points)) {
    n_points <- if (int_method == "qmc") object$qmc_points else object$n_points
    if (is.null(n_points)) n_points <- if (int_method == "qmc") 1024L else 11L
  }
  n_points <- as.integer(n_points)
  if (length(n_points) != 1L || is.na(n_points) || n_points < 1L) {
    stop("'n_points' must be a positive integer.", call. = FALSE)
  }
  k0 <- p + object$q + object$k_re
  # Warm starts from other fits would only change the modes within tolerance;
  # a fresh start keeps the result reproducible
  .brsmm_reset_cache()
  on.exit(.brsmm_reset_cache(), add = TRUE)
  out <- .brsmm_loglik_groups_cpp(
    param = par[seq_len(k0)], X = dd$X, Z = dd$Z, Xr = dd$Xr,
    y_left = as.numeric(Y[, "left"]), y_right = as.numeric(Y[, "right"]),
    yt = as.numeric(Y[, "yt"]), delta = delta,
    group = as.integer(dd$group),
    link_mu = link_to_code(object$link), link_phi = link_to_code(object$link_phi),
    repar = object$repar,
    method = match(int_method, c("laplace", "aghq", "qmc")) - 1L,
    n_points = n_points
  )
  out <- as.numeric(out)
  names(out) <- levels(dd$group)
  attr(out, "n") <- as.vector(table(dd$group))
  out
}

# Response matrix, design matrices and groups of a fit (newdata = NULL) or of
# new observations coarsened with the fit's cells.
.brs_contrib_data <- function(object, newdata) {
  mixed <- inherits(object, "brsmm")
  if (is.null(newdata)) {
    mm <- object$model_matrices
    return(list(
      Y = object$Y, X = mm$X, Z = mm$Z, Xr = if (mixed) mm$Xr else NULL,
      group = if (mixed) droplevels(factor(object$group)) else NULL
    ))
  }
  if (!is.data.frame(newdata)) {
    stop("'newdata' must be a data.frame.", call. = FALSE)
  }
  vars <- all.vars(stats::formula(object$formula))
  if (mixed) {
    vars <- c(vars, all.vars(object$random$re_terms), object$random$group)
    if (!(object$random$group %in% names(newdata))) {
      stop("Grouping variable '", object$random$group, "' not found in 'newdata'.",
           call. = FALSE)
    }
  }
  vars <- intersect(unique(vars), names(newdata))
  cc <- stats::complete.cases(newdata[, vars, drop = FALSE])
  nd <- .brs_keep_prep_attr(newdata[cc, , drop = FALSE], newdata)
  if (!nrow(nd)) stop("'newdata' has no complete rows.", call. = FALSE)
  mf <- stats::model.frame(object$formula, data = nd)
  Y <- .extract_response(mf, nd, ncuts = object$ncuts, lim = object$lim,
                         interval = .brs_interval_of(object))
  # Cells of the estimated boundary cutpoints (unchanged under the fixed map)
  Y <- .brs_apply_fit_cuts(Y, object)
  des <- .brs_marg_design(object, nd)
  if (!all(des$ok) || nrow(des$X) != nrow(Y)) {
    stop("Could not align the rows of 'newdata' with the model.", call. = FALSE)
  }
  list(
    Y = Y, X = des$X, Z = des$Z, Xr = des$Xr,
    group = if (mixed) factor(nd[[object$random$group]]) else NULL
  )
}

# Copy the brs_prep() attributes of `from` to a row subset `x`.
.brs_keep_prep_attr <- function(x, from) {
  for (a in c("is_prepared", "ncuts", "lim", "interval")) {
    v <- attr(from, a, exact = TRUE)
    if (!is.null(v)) attr(x, a) <- v
  }
  x
}
