# ============================================================================ #
# ANOVA / Likelihood-ratio model comparison for brs and brsmm
# ============================================================================ #

#' Internal ANOVA table builder for brs/brsmm
#' @param models List of fitted models.
#' @param test Character test type.
#' @keywords internal
.anova_brs_family <- function(models, test = c("Chisq", "none")) {
  test <- match.arg(test)

  if (length(models) == 0L) {
    stop("At least one model must be supplied.", call. = FALSE)
  }

  ok <- vapply(models, function(m) inherits(m, "brs") || inherits(m, "brsmm"), logical(1))
  if (!all(ok)) {
    stop("All models must inherit from 'brs' or 'brsmm'.", call. = FALSE)
  }

  nobs_vec <- vapply(models, nobs, numeric(1))
  if (length(unique(nobs_vec)) != 1L) {
    stop("All models must be fitted to the same number of observations.", call. = FALSE)
  }
  # Different coarsenings (interval/ncuts/lim) are different response models:
  # their log-likelihoods are not comparable
  for (what in c("interval", "ncuts", "lim")) {
    # lim only defines the cells under "mid" (ignored for right/left)
    if (what == "lim" && !identical(.brs_interval_of(models[[1L]]), "mid")) next
    vals <- vapply(models, function(m) {
      v <- if (what == "interval") .brs_interval_of(m) else m[[what]]
      format(v)
    }, character(1))
    if (length(unique(vals)) != 1L) {
      stop(
        "All models must use the same `", what, "` (found: ",
        paste(unique(vals), collapse = ", "), "); fits with different ",
        "coarsenings of the response are not comparable.",
        call. = FALSE
      )
    }
  }

  ll_obj <- lapply(models, logLik)
  ll <- vapply(ll_obj, as.numeric, numeric(1))
  df <- vapply(ll_obj, function(x) as.numeric(attr(x, "df")), numeric(1))
  aic <- vapply(models, AIC, numeric(1))
  bic <- vapply(models, BIC, numeric(1))

  cls <- vapply(models, function(m) class(m)[1L], character(1))
  # Number of random-effect terms (0 for brs)
  nre <- vapply(models, function(m) if (inherits(m, "brsmm")) as.numeric(m$q_re) else 0, numeric(1))
  ord <- order(df, ll)
  if (!identical(ord, seq_along(models))) {
    models <- models[ord]
    ll <- ll[ord]
    df <- df[ord]
    aic <- aic[ord]
    bic <- bic[ord]
    cls <- cls[ord]
    nre <- nre[ord]
  }

  out <- data.frame(
    Model = paste0("M", seq_along(models), " (", cls, ")"),
    Df = df,
    logLik = ll,
    AIC = aic,
    BIC = bic,
    check.names = FALSE
  )

  heading <- "Likelihood-ratio comparison of brs/brsmm models"
  if (length(models) >= 2L && identical(test, "Chisq")) {
    ddf <- c(NA_real_, diff(df))
    lr <- c(NA_real_, 2 * diff(ll))
    pval <- rep(NA_real_, length(models))
    valid <- !is.na(ddf) & ddf > 0
    pval[valid] <- stats::pchisq(lr[valid], df = ddf[valid], lower.tail = FALSE)

    # One added random term: its variance is on the boundary (Self & Liang;
    # Stram & Lee), so LR ~ 1/2 chi2(d - 1) + 1/2 chi2(d)
    dre <- c(NA_real_, diff(nre))
    mix <- valid & dre == 1 & ddf >= nre
    if (any(mix)) {
      pval[mix] <- ifelse(lr[mix] <= 0, 1,
        0.5 * stats::pchisq(lr[mix], df = ddf[mix] - 1, lower.tail = FALSE) +
          0.5 * stats::pchisq(lr[mix], df = ddf[mix], lower.tail = FALSE))
      heading <- c(heading, paste0(
        "Rows ", paste0("M", which(mix), collapse = ", "), ": one added random ",
        "effect (variance on the boundary); Pr(>Chisq) from the chi-bar-square ",
        "mixture 1/2 chi2(Df - 1) + 1/2 chi2(Df)."))
    }
    other <- valid & !is.na(dre) & dre != 0 & !mix
    if (any(other)) {
      heading <- c(heading, paste0(
        "Rows ", paste0("M", which(other), collapse = ", "), ": the random-effect ",
        "structure changes by more than one term; the chi2(Df) p-value is ",
        "conservative (boundary)."))
    }

    out$Chisq <- lr
    out$`Chi Df` <- ddf
    out$`Pr(>Chisq)` <- pval
  }

  rownames(out) <- out$Model
  out$Model <- NULL
  class(out) <- c("anova", "data.frame")
  # print.anova() prints the heading (the old attr "note" was never shown)
  attr(out, "heading") <- c(heading, "")
  out
}

#' Likelihood-ratio comparison of nested beta interval models
#'
#' @description
#' Compares fitted \code{"brs"} and \code{"brsmm"} models by log-likelihood,
#' AIC, BIC and likelihood-ratio tests (Lopes, 2023, "Inferencia").
#'
#' @details
#' The models are sorted by their number of parameters. For consecutive
#' models the statistic is \eqn{LR = 2(\ell_1 - \ell_0)}, with \code{Chi Df}
#' the difference in the number of parameters, and
#' \eqn{p = P(\chi^2_{df} > LR)}. The models must be nested; this is not
#' checked. They must also describe the same response: the same observations,
#' \code{interval}, \code{ncuts} and (under \code{"mid"}) \code{lim}, or the
#' call stops, since different coarsenings are different likelihoods.
#'
#' When the larger model adds one random-effect term (a \code{"brs"} model
#' against a random-intercept \code{"brsmm"}, or one more correlated random
#' term), its variance lies on the boundary of the parameter space under
#' \eqn{H_0} and \eqn{LR} follows the mixture
#' \eqn{\frac12\chi^2_{df-1} + \frac12\chi^2_{df}} (Self and Liang, 1987;
#' Stram and Lee, 1994); for one variance component alone this is
#' \eqn{\frac12\chi^2_0 + \frac12\chi^2_1}, i.e. half the naive p-value. The
#' printed heading names the rows that use it. When more than one random term
#' is added at once, the naive \eqn{\chi^2_{df}} p-value is kept and flagged
#' as conservative.
#'
#' @param object A fitted \code{"brs"} model.
#' @param ... Further fitted \code{"brs"} and/or \code{"brsmm"} models.
#' @param test \code{"Chisq"} (default) or \code{"none"}.
#'
#' @return An object of class \code{"anova"} (a data frame) with columns
#'   \code{Df}, \code{logLik}, \code{AIC}, \code{BIC} and, for
#'   \code{test = "Chisq"}, \code{Chisq}, \code{Chi Df} and
#'   \code{Pr(>Chisq)}; the attribute \code{"heading"} explains the p-values.
#'
#' @seealso \code{\link{anova.brsmm}}, \code{\link{summary.brs}},
#'   \code{\link{logLik.brs}}
#'
#' @references
#' Lopes, J. E. (2023). \emph{Modelos de regressao beta para dados de escala}.
#' Master's dissertation, Universidade Federal do Parana, Curitiba.
#' URI: https://hdl.handle.net/1884/86624.
#'
#' Self, S. G., and Liang, K.-Y. (1987). Asymptotic properties of maximum
#' likelihood estimators and likelihood ratio tests under nonstandard
#' conditions. \emph{Journal of the American Statistical Association},
#' \bold{82}(398), 605--610. \doi{10.1080/01621459.1987.10478472}
#'
#' Stram, D. O., and Lee, J. W. (1994). Variance components testing in the
#' longitudinal mixed effects model. \emph{Biometrics}, \bold{50}(4),
#' 1171--1177. \doi{10.2307/2533455}
#'
#' @examples
#' # Synthetic NRS-11 scores: 4 groups x 3 times. Simulated, not real data.
#' set.seed(2023)
#' nrs <- expand.grid(id = 1:80, time = c("6h", "12h", "24h"))
#' nrs$group <- factor(paste0("g", (nrs$id - 1) %% 4 + 1))
#' eta <- -1.3 + c(0, 0.75, 0.3)[nrs$time] + c(0, -0.1, 0.05, 0.1)[nrs$group]
#' shp <- brs_repar(mu = plogis(eta), phi = 0.3, repar = 2)
#' nrs$y <- round(10 * rbeta(nrow(nrs), shp$shape1, shp$shape2))
#'
#' # Time only (m1) nested in time + group (m2): LR test on 3 df
#' m1 <- brs(y ~ time, data = nrs, ncuts = 10)
#' m2 <- brs(y ~ time + group, data = nrs, ncuts = 10)
#' anova(m1, m2)
#'
#' # Fits under different interval directions are different response models
#' m2_right <- brs(y ~ time + group, data = nrs, ncuts = 10, interval = "right")
#' try(anova(m2, m2_right))
#'
#' @method anova brs
#' @importFrom stats anova pchisq
#' @export
anova.brs <- function(object, ..., test = c("Chisq", "none")) {
  models <- c(list(object), list(...))
  .anova_brs_family(models = models, test = test)
}

#' Likelihood-ratio comparison involving mixed models
#'
#' @description
#' \code{anova()} for \code{"brsmm"} fits: the same table as
#' \code{\link{anova.brs}}, with the chi-bar-square mixture
#' \eqn{\frac12\chi^2_{df-1} + \frac12\chi^2_{df}} for rows that add one
#' random-effect term, whose variance is on the boundary under \eqn{H_0}.
#' This is the test to use for a variance component: the Wald statistic of
#' its log standard deviation is not meaningful (see
#' \code{\link{summary.brsmm}}).
#'
#' @param object A fitted \code{"brsmm"} model.
#' @param ... Further fitted \code{"brsmm"} and/or \code{"brs"} models.
#' @param test \code{"Chisq"} (default) or \code{"none"}.
#'
#' @return An object of class \code{"anova"}; see \code{\link{anova.brs}}.
#'
#' @seealso \code{\link{anova.brs}}, \code{\link{brsmm}},
#'   \code{\link{summary.brsmm}}
#'
#' @references
#' Self, S. G., and Liang, K.-Y. (1987). Asymptotic properties of maximum
#' likelihood estimators and likelihood ratio tests under nonstandard
#' conditions. \emph{Journal of the American Statistical Association},
#' \bold{82}(398), 605--610. \doi{10.1080/01621459.1987.10478472}
#'
#' @examples
#' set.seed(11)
#' g <- 20
#' d <- data.frame(id = factor(rep(1:g, each = 8)), x = runif(8 * g))
#' shp <- brs_repar(plogis(-0.4 + d$x + rnorm(g, sd = 0.6)[d$id]), phi = 0.25)
#' d$y <- round(10 * rbeta(nrow(d), shp$shape1, shp$shape2))
#' m0 <- brs(y ~ x, data = d, ncuts = 10)
#' m1 <- brsmm(y ~ x, random = ~ 1 | id, data = d, ncuts = 10)
#' anova(m0, m1)  # Pr(>Chisq) = half the chi2(1) tail
#'
#' @method anova brsmm
#' @importFrom stats anova pchisq
#' @export
anova.brsmm <- function(object, ..., test = c("Chisq", "none")) {
  models <- c(list(object), list(...))
  .anova_brs_family(models = models, test = test)
}
