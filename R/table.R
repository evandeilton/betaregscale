# ============================================================================ #
# Model comparison tables
# ============================================================================ #

#' Compare fitted brs models in a single table
#'
#' @description
#' Builds a comparison table for one or more fitted \code{"brs"} objects,
#' summarizing fit statistics and (optionally) censoring composition.
#'
#' @param ... Fitted \code{"brs"} objects passed individually.
#' @param models Optional list of fitted \code{"brs"} objects.
#'   Use either \code{...} or \code{models}, not both.
#' @param include_censoring Logical; include censoring counts/proportions.
#'   Default is \code{TRUE}.
#' @param sort_by Character; optional sort criterion:
#'   \code{"none"} (default), \code{"AIC"}, \code{"BIC"}, or \code{"logLik"}.
#' @param decreasing Logical; sort direction when \code{sort_by != "none"}.
#' @param digits Integer number of digits used for numeric rounding.
#'
#' @return A data frame with one row per model.
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
#' m1 <- brs(y ~ 1, data = prep)
#' m2 <- brs(y ~ x1, data = prep)
#' brs_table(null = m1, x1 = m2, sort_by = "AIC")
#' }
#'
#' @rdname brs_table
#' @export
brs_table <- function(...,
                      models = NULL,
                      include_censoring = TRUE,
                      sort_by = c("none", "AIC", "BIC", "logLik"),
                      decreasing = FALSE,
                      digits = 4L) {
  dots <- list(...)
  if (length(dots) > 0L && !is.null(models)) {
    stop("Use either '...' or 'models', not both.", call. = FALSE)
  }

  if (is.null(models)) {
    models <- dots
  }
  if (inherits(models, "brs")) {
    models <- list(models)
  }
  if (!is.list(models) || length(models) == 0L) {
    stop("Provide at least one fitted 'brs' object.", call. = FALSE)
  }

  sort_by <- match.arg(sort_by)
  digits <- as.integer(digits)
  if (!is.finite(digits) || digits < 0L) {
    stop("'digits' must be a non-negative integer.", call. = FALSE)
  }

  nm <- names(models)
  if (is.null(nm)) {
    nm <- rep("", length(models))
  }
  missing_nm <- which(is.na(nm) | !nzchar(nm))
  if (length(missing_nm) > 0L) {
    nm[missing_nm] <- paste0("model_", missing_nm)
  }

  rows <- lapply(seq_along(models), function(i) {
    obj <- models[[i]]
    .check_class(obj)

    base_row <- data.frame(
      model = nm[i],
      nobs = obj$nobs,
      npar = obj$npar,
      logLik = as.numeric(logLik(obj)),
      AIC = AIC(obj),
      BIC = BIC(obj),
      pseudo_r2 = obj$pseudo.r.squared,
      stringsAsFactors = FALSE
    )

    if (!isTRUE(include_censoring)) {
      return(base_row)
    }

    delta <- obj$delta
    n <- length(delta)
    counts <- c(
      exact = sum(delta == 0L),
      left = sum(delta == 1L),
      right = sum(delta == 2L),
      interval = sum(delta == 3L)
    )

    cbind(
      base_row,
      data.frame(
        exact = counts["exact"],
        left = counts["left"],
        right = counts["right"],
        interval = counts["interval"],
        prop_exact = counts["exact"] / n,
        prop_left = counts["left"] / n,
        prop_right = counts["right"] / n,
        prop_interval = counts["interval"] / n,
        stringsAsFactors = FALSE
      )
    )
  })

  out <- do.call(rbind, rows)

  if (sort_by != "none") {
    o <- order(out[[sort_by]], decreasing = isTRUE(decreasing), na.last = TRUE)
    out <- out[o, , drop = FALSE]
  }

  num_cols <- vapply(out, is.numeric, logical(1))
  if (any(num_cols)) {
    out[num_cols] <- lapply(out[num_cols], function(x) round(x, digits = digits))
  }

  rownames(out) <- NULL
  out
}
