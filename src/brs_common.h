// brs_common.h — shared constants and inline helpers for betaregscale C++ backends
//
// Include this header from loglik.cpp and loglik_mixed_eigen.cpp ONLY.
// All definitions are inline to avoid ODR violations across translation units.
//
// BUG-H01  : eliminates code duplication between the two backends.
// BUG-C01  : inv_link case 6 (inverse) always returns positive value.
// BUG-H05  : inv_link case 5 (sqrt) guards eta < 0 to keep gradient sign correct.

#pragma once
#include <Rmath.h>
#include <cmath>
#include <vector>

// ---------------------------------------------------------------- constants --

static const double EPS_SHAPE   = 1.0e-12;
static const double MAX_SHAPE   = 1.0e8;
static const double EPS_PROB    = 1.0e-15;
static const double LOG_PENALTY = -1.0e6;
// EPS_UNIT : clamp margin for (0,1)-scale endpoints (y, phi in repar=2)
// EPS_SHAPE: minimum allowed beta shape parameter value
// EPS_PROB : tiny threshold used only by the inverse / 1/mu^2 links
// LOG_PENALTY: value returned for a non-finite contribution (never a floor)
static const double EPS_UNIT   = 1.0e-5;

// ------------------------------------------------------------------ helpers --

inline double clamp(double x, double lo, double hi) {
  return (x < lo) ? lo : ((x > hi) ? hi : x);
}

// Inverse-link functions dispatched by integer code (0-indexed, see link_to_code in R).
//   0 = logit, 1 = probit, 2 = cauchit, 3 = cloglog,
//   4 = log,   5 = sqrt,   6 = inverse, 7 = 1/mu^2,  8 = identity
inline double inv_link(double eta, int code) {
  switch (code) {
  case 0: return 1.0 / (1.0 + std::exp(-eta));
  case 1: return R::pnorm(eta, 0.0, 1.0, 1, 0);
  case 2: return 0.5 + std::atan(eta) / M_PI;
  case 3: return 1.0 - std::exp(-std::exp(eta));
  case 4: return std::exp(eta);
  case 5:
    // sqrt link: g^{-1}(eta) = eta^2 for eta >= 0, 0 below (sign-correct gradient);
    // NaN propagates (as pmax(eta, 0)^2 in R).
    return eta >= 0.0 ? eta * eta : (std::isnan(eta) ? eta : 0.0);
  case 6:
    // inverse link: g^{-1}(eta) = 1/eta.
    // BUG-C01: for |eta|~0 we must return a finite POSITIVE value (dispersion
    // is always positive). Old code returned -MAX_SHAPE for eta<0, which is wrong.
    if (std::abs(eta) < EPS_PROB) return MAX_SHAPE;
    return 1.0 / eta;
  case 7:
    // 1/mu^2 link: g^{-1}(eta) = 1/sqrt(eta). Requires eta > 0.
    if (eta <= EPS_PROB) return MAX_SHAPE;
    return 1.0 / std::sqrt(eta);
  case 8: return eta;
  default: return 1.0 / (1.0 + std::exp(-eta));
  }
}

// Clamp phi (dispersion/precision) to valid range for the chosen reparameterisation.
// repar = 1: precision phi > 0
// repar = 2: mean-variance phi in (0, 1)
// repar = 0: direct shape parameter, enforce positivity
// +-Inf map to the bounds; NaN propagates (obs_loglik then gives LOG_PENALTY).
inline double clamp_phi_by_repar(double phi, int repar) {
  if (repar == 2) return clamp(phi, EPS_UNIT, 1.0 - EPS_UNIT);
  return clamp(phi, EPS_UNIT, MAX_SHAPE);
}

// Clamp the first parameter: the mean in (0, 1) for repar 1/2, the shape p > 0
// for repar 0 (a (0, 1) clamp there capped p at 1 - EPS_UNIT). Mirror: R .clamp_mu_by_repar().
// +-Inf map to the bounds; NaN propagates (obs_loglik then gives LOG_PENALTY).
inline double clamp_mu_by_repar(double mu, int repar) {
  if (repar == 0) return clamp(mu, EPS_UNIT, MAX_SHAPE);
  return clamp(mu, EPS_UNIT, 1.0 - EPS_UNIT);
}

// Convert (mu, phi) to beta shape parameters (a, b) under the chosen reparameterisation.
// repar = 0: a = mu,            b = phi
// repar = 1: a = mu*phi,        b = (1-mu)*phi  [Ferrari & Cribari-Neto 2004]
// repar = 2: a = mu*(1-phi)/phi, b = (1-mu)*(1-phi)/phi  [mean-variance form]
inline void beta_shapes(double mu, double phi, int repar, double &a, double &b) {
  switch (repar) {
  case 0: a = mu; b = phi; break;
  case 1: a = mu * phi; b = (1.0 - mu) * phi; break;
  case 2: {
    double ratio = (1.0 - phi) / phi;
    a = mu * ratio;
    b = (1.0 - mu) * ratio;
    break;
  }
  default: a = mu; b = phi;
  }
  a = clamp(a, EPS_SHAPE, MAX_SHAPE);
  b = clamp(b, EPS_SHAPE, MAX_SHAPE);
}

// ------------------------------------------------- log-likelihood building blocks --
//
// Censored contributions are log(probabilities). Rules (mirrored 1:1 by the
// R helper .brs_obs_loglik() in R/loglik.R, keep both in sync):
//
//  * Endpoints are clamped to [EPS_UNIT, 1 - EPS_UNIT]: that is the only
//    protection of the borders of the support. There is NO probability
//    floor. The former floor (1e-15) turned every far-tail observation into
//    a constant with zero gradient, so the optimiser maximised a trimmed
//    likelihood that ignored outliers (audit 2026-09, finding C2).
//
//  * The tail is chosen by the distribution, not by the position of the
//    interval on (0, 1): an interval whose midpoint lies at or below the
//    mean a/(a+b) is measured with lower-tail CDFs, otherwise with upper
//    tails (survival). Both probabilities are then "small" quantities and
//    their difference p1 - p2 is well conditioned. The old rule (lower tail
//    iff lo + hi <= 1) subtracted two values equal to 1 to machine precision
//    whenever the interval sat far above a small mean (or below a large
//    one): catastrophic cancellation, and the floor hid it.
//
//  * pbeta() is evaluated in PLAIN scale. It is accurate to ~1e-10 relative
//    for p >= 1e-240; between underflow and ~1e-263 bratio can return
//    values wrong by up to ~10%, hence the threshold. It never warns. The
//    log-scale variant (log_p = TRUE) can emit "bpser(...) underflow to
//    -Inf" R warnings for large shapes, in a way that is not predictable
//    from log f alone, so it is not used.
//
//  * Below P_TINY the plain probability is unreliable or underflows. There
//    the endpoint Laplace approximation of the tail integral takes over
//    (log_tail_laplace): for x below the mode
//        F(x)  ~ f(x) / g'(x)  * (1 + g''(x)/g'(x)^2),   g = log f,
//    and for x above the mode the same with |g'|. Its error there is
//    <= 4e-4 log-lik units (checked against log-scale pbeta where the
//    latter is reliable), it is smooth in (a, b), and it keeps a usable
//    gradient while the optimiser is still far from the optimum. It only
//    applies on the correct side of the mode (g' > 0 for a lower tail,
//    g' < 0 for an upper tail); otherwise the contribution is -Inf, which
//    obs_loglik() maps to LOG_PENALTY, as it does for any non-finite value.

// Plain-scale probabilities below this are not trusted (see above).
static const double P_TINY = 1.0e-240;

// log f(yt | a, b)  [exact / uncensored]. The clamp keeps an exact response
// inside the open support (0, 1), where the density can diverge.
inline double log_density(double yt, double a, double b) {
  double y = clamp(yt, EPS_UNIT, 1.0 - EPS_UNIT);
  return R::dbeta(y, a, b, 1);  // log = TRUE
}

// log of the tail mass beyond x by the endpoint Laplace approximation.
//   lower_tail = true : log P(x - width < Y < x)  (x is the upper endpoint)
//   lower_tail = false: log P(x < Y < x + width)  (x is the lower endpoint)
// width = R_PosInf gives the one-sided tail log F(x) or log S(x).
// log P ~ log f(x) - log|g'| + log(1 + g''/g'^2) + log(1 - exp(-|g'| width))
inline double log_tail_laplace(double x, double a, double b, double width,
                               bool lower_tail) {
  double gp  = (a - 1.0) / x - (b - 1.0) / (1.0 - x);            // g'(x)
  if (lower_tail ? !(gp > 0.0) : !(gp < 0.0)) return R_NegInf;  // wrong side
  double s   = std::abs(gp);
  double gpp = -(a - 1.0) / (x * x) - (b - 1.0) / ((1.0 - x) * (1.0 - x));
  double corr = 1.0 + gpp / (s * s);
  double v = log_density(x, a, b) - std::log(s);
  if (corr > 0.0) v += std::log(corr);
  if (R_FINITE(width)) v += Rf_log1mexp(s * width);
  return v;
}

// log P(left < Y < right | a, b)  [interval-censored contribution]
inline double log_interval_prob(double left, double right, double a, double b) {
  double lo = clamp(left,  EPS_UNIT, 1.0 - EPS_UNIT);
  double hi = clamp(right, EPS_UNIT, 1.0 - EPS_UNIT);
  if (!(hi > lo)) return R_NegInf;
  bool lower = 0.5 * (lo + hi) <= a / (a + b);
  double p1, p2;  // p1 = mass beyond the endpoint nearest the bulk (>= p2)
  if (lower) {
    p1 = R::pbeta(hi, a, b, 1, 0);
    p2 = R::pbeta(lo, a, b, 1, 0);
  } else {
    p1 = R::pbeta(lo, a, b, 0, 0);
    p2 = R::pbeta(hi, a, b, 0, 0);
  }
  if (p1 >= P_TINY) {
    double area = p1 - p2;
    return (area > 0.0) ? std::log(area) : R_NegInf;
  }
  return log_tail_laplace(lower ? hi : lo, a, b, hi - lo, lower);
}

// log F(y | a, b)  [left-censored contribution]
inline double log_cdf(double y, double a, double b) {
  double yc = clamp(y, EPS_UNIT, 1.0 - EPS_UNIT);
  double p  = R::pbeta(yc, a, b, 1, 0);
  if (p >= P_TINY) return std::log(p);
  return log_tail_laplace(yc, a, b, R_PosInf, true);
}

// log (1 - F(y | a, b))  [right-censored contribution]
inline double log_survival(double y, double a, double b) {
  double yc = clamp(y, EPS_UNIT, 1.0 - EPS_UNIT);
  double p  = R::pbeta(yc, a, b, 0, 0);  // upper tail
  if (p >= P_TINY) return std::log(p);
  return log_tail_laplace(yc, a, b, R_PosInf, false);
}

// Per-observation log-likelihood contribution.
// delta_i : 0=exact, 1=left-censored, 2=right-censored, 3=interval-censored
inline double obs_loglik(int delta_i, double left_i, double right_i,
                         double yt_i, double a, double b) {
  double contrib;
  switch (delta_i) {
  case 0: contrib = log_density(yt_i, a, b);                   break;
  case 1: contrib = log_cdf(right_i, a, b);                    break;
  case 2: contrib = log_survival(left_i, a, b);                break;
  case 3: contrib = log_interval_prob(left_i, right_i, a, b);  break;
  default: contrib = log_interval_prob(left_i, right_i, a, b);
  }
  return std::isfinite(contrib) ? contrib : LOG_PENALTY;
}
