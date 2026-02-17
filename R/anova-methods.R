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

  ll_obj <- lapply(models, logLik)
  ll <- vapply(ll_obj, as.numeric, numeric(1))
  df <- vapply(ll_obj, function(x) as.numeric(attr(x, "df")), numeric(1))
  aic <- vapply(models, AIC, numeric(1))
  bic <- vapply(models, BIC, numeric(1))

  cls <- vapply(models, function(m) class(m)[1L], character(1))
  ord <- order(df, ll)
  if (!identical(ord, seq_along(models))) {
    models <- models[ord]
    ll <- ll[ord]
    df <- df[ord]
    aic <- aic[ord]
    bic <- bic[ord]
    cls <- cls[ord]
  }

  out <- data.frame(
    Model = paste0("M", seq_along(models), " (", cls, ")"),
    Df = df,
    logLik = ll,
    AIC = aic,
    BIC = bic,
    check.names = FALSE
  )

  if (length(models) >= 2L && identical(test, "Chisq")) {
    ddf <- c(NA_real_, diff(df))
    lr <- c(NA_real_, 2 * diff(ll))
    pval <- rep(NA_real_, length(models))
    valid <- !is.na(ddf) & ddf > 0
    pval[valid] <- stats::pchisq(lr[valid], df = ddf[valid], lower.tail = FALSE)

    out$Chisq <- lr
    out$`Chi Df` <- ddf
    out$`Pr(>Chisq)` <- pval
  }

  rownames(out) <- out$Model
  out$Model <- NULL
  class(out) <- c("anova", "data.frame")

  if (any(cls == "brs") && any(cls == "brsmm")) {
    attr(out, "note") <- paste(
      "Comparisons brs vs brsmm involve boundary/null random-effect conditions;",
      "use LR p-values as heuristic evidence and complement with information criteria."
    )
  }

  out
}

#' Model comparison by analysis of deviance (LR test) for `brs`
#'
#' @param object A fitted `brs` model.
#' @param ... Additional fitted `brs` and/or `brsmm` models.
#' @param test Character; `"Chisq"` (default) or `"none"`.
#'
#' @return An object of class `anova` and `data.frame` with model-wise
#'   log-likelihood, information criteria, and (optionally) LR test columns.
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
#' @param object A fitted `brsmm` model.
#' @param ... Additional fitted `brsmm` and/or `brs` models.
#' @param test Character; `"Chisq"` (default) or `"none"`.
#'
#' @return An object of class `anova` and `data.frame` with model-wise
#'   log-likelihood, information criteria, and (optionally) LR test columns.
#'
#' @method anova brsmm
#' @importFrom stats anova pchisq
#' @export
anova.brsmm <- function(object, ..., test = c("Chisq", "none")) {
  models <- c(list(object), list(...))
  .anova_brs_family(models = models, test = test)
}
