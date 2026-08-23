# ============================================================================ #
# Score-probability predictions
# ============================================================================ #

#' Predict score probabilities from a fitted brs model
#'
#' @description
#' Computes predicted probabilities for integer scores on the original scale
#' \eqn{\{0, 1, \ldots, K\}} implied by the fitted beta interval model.
#'
#' @details
#' For a score \eqn{s} and \eqn{K =} \code{ncuts}, probabilities are computed as:
#' \itemize{
#'   \item \eqn{P(Y=s)=F(\mathrm{lim}/K)} for \eqn{s=0},
#'   \item \eqn{P(Y=s)=1-F((K-\mathrm{lim})/K)} for \eqn{s=K},
#'   \item \eqn{P(Y=s)=F((s+\mathrm{lim})/K)-F((s-\mathrm{lim})/K)}
#'     for \eqn{s \in \{1,\ldots,K-1\}},
#' }
#' where \eqn{F} is the beta CDF under the fitted \eqn{(\mu_i,\phi_i)}.
#'
#' @param object A fitted \code{"brs"} object.
#' @param newdata Optional data frame for prediction.
#' @param scores Optional integer vector of scores to evaluate.
#'   Defaults to all scores in \code{0:ncuts}.
#' @param format Output format: \code{"matrix"} (default) or \code{"long"}.
#' @param id_col Column name for observation id when \code{format = "long"}.
#'
#' @return
#' If \code{format = "matrix"}, a numeric matrix with one row per
#' observation and one column per requested score.
#'
#' If \code{format = "long"}, a long data frame with columns
#' \code{id_col}, \code{score}, and \code{prob}.
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
#' pmat <- brs_predict_scoreprob(fit)
#' head(pmat[, 1:5])
#' plong <- brs_predict_scoreprob(fit, scores = 0:10, format = "long")
#' head(plong)
#' }
#'
#' @rdname brs_predict_scoreprob
#' @export
brs_predict_scoreprob <- function(object,
                                  newdata = NULL,
                                  scores = NULL,
                                  format = c("matrix", "long"),
                                  id_col = "id") {
  .check_class(object)
  format <- match.arg(format)

  K <- as.integer(object$ncuts)
  if (is.null(scores)) {
    scores <- 0:K
  }
  scores <- sort(unique(as.integer(scores)))
  if (any(!is.finite(scores)) || any(scores < 0L) || any(scores > K)) {
    stop("'scores' must be integers in [0, ncuts].", call. = FALSE)
  }

  if (is.null(newdata)) {
    mu <- object$hatmu
    phi <- object$hatphi
  } else {
    if (!is.data.frame(newdata)) {
      stop("'newdata' must be a data.frame.", call. = FALSE)
    }
    mu <- predict(object, newdata = newdata, type = "response")
    phi <- predict(object, newdata = newdata, type = "precision")
  }

  P <- .brs_score_prob_matrix(
    mu = mu,
    phi = phi,
    repar = object$repar,
    ncuts = K,
    lim = object$lim,
    scores = scores
  )

  if (identical(format, "matrix")) {
    return(P)
  }

  n <- nrow(P)
  data.frame(
    stats::setNames(list(rep(seq_len(n), each = length(scores))), id_col),
    score = rep(scores, times = n),
    prob = as.vector(t(P)),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}
