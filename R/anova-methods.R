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

#' Model comparison by analysis of deviance (LR test) for `brs`
#'
#' @param object A fitted \code{"brs"} model.
#' @param ... Additional fitted \code{"brs"} and/or \code{"brsmm"} models to
#'   compare.
#' @param test Character; \code{"Chisq"} (default) or \code{"none"}.
#'
#' @return An object of class \code{"anova"} and \code{"data.frame"} with
#'   model-wise log-likelihood, information criteria, and (optionally) LR test
#'   columns.
#'
#' @seealso \code{\link{brs}}, \code{\link{logLik.brs}}, \code{\link{AIC.brs}},
#'   \code{\link{BIC.brs}}
#'
#' @references
#' Lopes, J. E. (2023). \emph{Modelos de regressao beta para dados de escala}.
#' Master's dissertation, Universidade Federal do Parana, Curitiba.
#' URI: https://hdl.handle.net/1884/86624.
#'
#' Ferrari, S. L. P., and Cribari-Neto, F. (2004).
#' Beta regression for modelling rates and proportions.
#' \emph{Journal of Applied Statistics}, \bold{31}(7), 799--815.
#' \doi{10.1080/0266476042000214501}
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
#' m1 <- brs(y ~ 1, data = prep)
#' m2 <- brs(y ~ x1, data = prep)
#' m3 <- brs(y ~ x1 + x2, data = prep)
#' anova(m1, m2, m3)
#' }
#'
#' @method anova brs
#' @importFrom stats anova pchisq
#' @export
anova.brs <- function(object, ..., test = c("Chisq", "none")) {
  models <- c(list(object), list(...))
  .anova_brs_family(models = models, test = test)
}

#' Model comparison by analysis of deviance (LR test) for `brsmm`
#'
#' @param object A fitted \code{"brsmm"} model.
#' @param ... Additional fitted \code{"brsmm"} and/or \code{"brs"} models to
#'   compare.
#' @param test Character; \code{"Chisq"} (default) or \code{"none"}.
#'
#' @return An object of class \code{"anova"} and \code{"data.frame"} with
#'   model-wise log-likelihood, information criteria, and (optionally) LR test
#'   columns.
#'
#' @seealso \code{\link{brsmm}}, \code{\link{logLik.brsmm}},
#'   \code{\link{AIC.brsmm}}, \code{\link{BIC.brsmm}}
#'
#' @references
#' Lopes, J. E. (2023). \emph{Modelos de regressao beta para dados de escala}.
#' Master's dissertation, Universidade Federal do Parana, Curitiba.
#' URI: https://hdl.handle.net/1884/86624.
#'
#' Ferrari, S. L. P., and Cribari-Neto, F. (2004).
#' Beta regression for modelling rates and proportions.
#' \emph{Journal of Applied Statistics}, \bold{31}(7), 799--815.
#' \doi{10.1080/0266476042000214501}
#'
#' @examples
#' \donttest{
#' dat <- data.frame(
#'   y = c(
#'     0, 5, 20, 50, 75, 90, 100, 30, 60, 45,
#'     10, 40, 55, 70, 85, 25, 35, 65, 80, 15
#'   ),
#'   x1 = rep(c(1, 2), 10),
#'   id = factor(rep(1:4, each = 5))
#' )
#' prep <- brs_prep(dat, ncuts = 100)
#' m1 <- brs(y ~ 1, data = prep)
#' m2 <- brsmm(y ~ x1, random = ~ 1 | id, data = prep)
#' anova(m1, m2)
#' }
#'
#' @method anova brsmm
#' @importFrom stats anova pchisq
#' @export
anova.brsmm <- function(object, ..., test = c("Chisq", "none")) {
  models <- c(list(object), list(...))
  .anova_brs_family(models = models, test = test)
}
