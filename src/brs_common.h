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
// EPS_PROB : minimum probability before applying LOG_PENALTY
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
    // sqrt link: g(phi)=sqrt(phi) => g^{-1}(eta)=eta^2.
    // BUG-H05: eta<0 gives same |value| but wrong gradient sign; clamp to 0.
    return eta >= 0.0 ? eta * eta : 0.0;
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
inline double clamp_phi_by_repar(double phi, int repar) {
  if (!std::isfinite(phi)) return (repar == 2) ? (1.0 - EPS_UNIT) : MAX_SHAPE;
  if (repar == 2) return clamp(phi, EPS_UNIT, 1.0 - EPS_UNIT);
  return clamp(phi, EPS_UNIT, MAX_SHAPE);
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

// log P(left < Y < right | a, b)
inline double log_interval_prob(double left, double right, double a, double b) {
  double lo = clamp(left,  EPS_UNIT, 1.0 - EPS_UNIT);
  double hi = clamp(right, EPS_UNIT, 1.0 - EPS_UNIT);
  double p1   = R::pbeta(lo, a, b, 1, 0);
  double p2   = R::pbeta(hi, a, b, 1, 0);
  double area = p2 - p1;
  if (area < EPS_PROB) area = EPS_PROB;
  return std::log(area);
}

// log f(yt | a, b)  [exact / uncensored]
inline double log_density(double yt, double a, double b) {
  double y = clamp(yt, EPS_UNIT, 1.0 - EPS_UNIT);
  return R::dbeta(y, a, b, 1);  // log = TRUE
}

// log F(y | a, b)  [left-censored contribution]
inline double log_cdf(double y, double a, double b) {
  double yc = clamp(y, EPS_UNIT, 1.0 - EPS_UNIT);
  double p  = R::pbeta(yc, a, b, 1, 0);
  if (p < EPS_PROB) p = EPS_PROB;
  return std::log(p);
}

// log (1 - F(y | a, b))  [right-censored contribution; uses upper tail for accuracy]
inline double log_survival(double y, double a, double b) {
  double yc = clamp(y, EPS_UNIT, 1.0 - EPS_UNIT);
  double p  = R::pbeta(yc, a, b, 0, 0);  // upper tail
  if (p < EPS_PROB) p = EPS_PROB;
  return std::log(p);
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
