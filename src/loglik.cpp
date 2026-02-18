// -------------------------------------------------------------------------- //
// betaregscale — C++ log-likelihood engine for interval-censored beta
//                regression with mixed censoring and optional variable
//                dispersion.
//
// Author : José Evandeilton Lopes
// License: MIT
//
// This file implements the core likelihood computations using
// RcppArmadillo for linear-algebra primitives and R's C-level
// pbeta / dbeta for distribution-function evaluation.
//
// Complete likelihood:
//
//   L(theta) = prod_{delta=0} f(y_i|theta)           [uncensored]
//            * prod_{delta=1} F(u_i|theta)            [left-censored]
//            * prod_{delta=2} [1 - F(l_i|theta)]      [right-censored]
//            * prod_{delta=3} [F(u_i|theta)-F(l_i|theta)] [interval]
//
// Censoring indicator delta_i:
//   0 = uncensored (exact)    -> log f(yt_i | a, b)
//   1 = left-censored         -> log F(right_i | a, b)
//   2 = right-censored        -> log(1 - F(left_i | a, b))
//   3 = interval-censored     -> log(F(right_i | a, b) - F(left_i | a, b))
//
// Numerical stability notes:
//   * Shape parameters are clamped to [EPS_SHAPE, MAX_SHAPE] before any
//     call to pbeta / dbeta.
//   * Interval probabilities  P(left < Y < right) are floored at EPS_PROB
//     to avoid log(0).
//   * Individual log-likelihood contributions that are non-finite are
//     replaced by a large negative penalty (LOG_PENALTY) so that the
//     optimizer is steered away from degenerate regions rather than
//     receiving NaN.
// -------------------------------------------------------------------------- //

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rmath.h>
#include <cmath>

// Numerical-stability constants
static const double EPS_SHAPE = 1.0e-12;
static const double MAX_SHAPE = 1.0e8;
static const double EPS_PROB = 1.0e-15;
static const double LOG_PENALTY = -1.0e6;
static const double EPS_BOUND = 1.0e-5;

// ------------------------------------------------------------------ helpers

// Clamp x into [lo, hi].
inline double clamp(double x, double lo, double hi) {
  return (x < lo) ? lo : ((x > hi) ? hi : x);
}

// Inverse-link functions.
// Dispatched by an integer code for speed inside the inner loop:
//   0 = logit, 1 = probit, 2 = cauchit, 3 = cloglog,
//   4 = log,   5 = sqrt,   6 = inverse (1/mu), 7 = 1/mu^2,
//   8 = identity
inline double inv_link(double eta, int code) {
  switch (code) {
  case 0:
    return 1.0 / (1.0 + std::exp(-eta)); // logit
  case 1:
    return R::pnorm(eta, 0.0, 1.0, 1, 0); // probit
  case 2:
    return 0.5 + std::atan(eta) / M_PI; // cauchit
  case 3:
    return 1.0 - std::exp(-std::exp(eta)); // cloglog
  case 4:
    return std::exp(eta); // log
  case 5:
    return eta * eta; // sqrt
  case 6:
    if (std::abs(eta) < EPS_PROB) {
      return (eta >= 0.0) ? MAX_SHAPE : -MAX_SHAPE;
    }
    return 1.0 / eta; // inverse
  case 7:
    if (eta <= EPS_PROB) {
      return MAX_SHAPE;
    }
    return 1.0 / std::sqrt(eta); // 1/mu^2
  case 8:
    return eta; // identity
  default:
    return 1.0 / (1.0 + std::exp(-eta));
  }
}

// Clamp phi according to the selected beta reparameterization.
// repar = 1: precision-style phi > 0
// repar = 2: mean-variance phi in (0, 1)
// repar = 0: direct shape2-like parameter, enforce positivity
inline double clamp_phi_by_repar(double phi, int repar) {
  if (!std::isfinite(phi)) {
    return (repar == 2) ? (1.0 - EPS_BOUND) : MAX_SHAPE;
  }
  if (repar == 2) {
    return clamp(phi, EPS_BOUND, 1.0 - EPS_BOUND);
  }
  return clamp(phi, EPS_BOUND, MAX_SHAPE);
}

// Beta reparameterization.
// repar = 0 : shape1 = mu,               shape2 = phi
// repar = 1 : shape1 = mu * phi,         shape2 = (1 - mu) * phi [Eq. 2.9]
// repar = 2 : shape1 = mu*(1-phi)/phi,   shape2 = (1-mu)*(1-phi)/phi [Eq. 2.11]
inline void beta_shapes(double mu, double phi, int repar, double &a,
                        double &b) {
  switch (repar) {
  case 0:
    a = mu;
    b = phi;
    break;
  case 1:
    a = mu * phi;
    b = (1.0 - mu) * phi;
    break;
  case 2: {
    double ratio = (1.0 - phi) / phi;
    a = mu * ratio;
    b = (1.0 - mu) * ratio;
    break;
  }
  default:
    a = mu;
    b = phi;
  }
  // Clamp to safe range
  a = clamp(a, EPS_SHAPE, MAX_SHAPE);
  b = clamp(b, EPS_SHAPE, MAX_SHAPE);
}

// Log of interval probability: log(P(left < Y < right | a, b))
inline double log_interval_prob(double left, double right, double a, double b) {
  double lo = clamp(left, EPS_BOUND, 1.0 - EPS_BOUND);
  double hi = clamp(right, EPS_BOUND, 1.0 - EPS_BOUND);

  double p1 = R::pbeta(lo, a, b, 1, 0);
  double p2 = R::pbeta(hi, a, b, 1, 0);

  double area = p2 - p1;
  if (area < EPS_PROB)
    area = EPS_PROB;

  return std::log(area);
}

// Log density at a point: log(f(yt | a, b))
inline double log_density(double yt, double a, double b) {
  double y = clamp(yt, EPS_BOUND, 1.0 - EPS_BOUND);
  return R::dbeta(y, a, b, 1); // log = TRUE
}

// Log CDF at a point: log(F(y | a, b))  [for left-censoring]
inline double log_cdf(double y, double a, double b) {
  double yc = clamp(y, EPS_BOUND, 1.0 - EPS_BOUND);
  double p = R::pbeta(yc, a, b, 1, 0); // lower_tail, non-log
  if (p < EPS_PROB)
    p = EPS_PROB;
  return std::log(p);
}

// Log survival at a point: log(1 - F(y | a, b))  [for right-censoring]
inline double log_survival(double y, double a, double b) {
  double yc = clamp(y, EPS_BOUND, 1.0 - EPS_BOUND);
  // Use upper tail directly for better numerical accuracy
  double p = R::pbeta(yc, a, b, 0, 0); // upper_tail, non-log
  if (p < EPS_PROB)
    p = EPS_PROB;
  return std::log(p);
}

// Compute the per-observation log-likelihood contribution given shape
// parameters (a, b), interval endpoints, midpoint, and censoring type.
//
// delta_i interpretation:
//   0 = exact / uncensored   : log f(yt_i | a, b)
//   1 = left-censored        : log F(right_i | a, b)
//   2 = right-censored       : log(1 - F(left_i | a, b))
//   3 = interval-censored    : log(F(right_i) - F(left_i))
inline double obs_loglik(int delta_i, double left_i, double right_i,
                         double yt_i, double a, double b) {
  double contrib;
  switch (delta_i) {
  case 0:
    contrib = log_density(yt_i, a, b);
    break;
  case 1:
    contrib = log_cdf(right_i, a, b);
    break;
  case 2:
    contrib = log_survival(left_i, a, b);
    break;
  case 3:
    contrib = log_interval_prob(left_i, right_i, a, b);
    break;
  default:
    contrib = log_interval_prob(left_i, right_i, a, b);
  }
  return std::isfinite(contrib) ? contrib : LOG_PENALTY;
}

// ============================================================ Fixed phi === //

//' @title C++ log-likelihood for fixed-dispersion beta interval regression
//'   with mixed censoring
//' @description Computes the total log-likelihood for a beta regression model
//'   with interval-censored responses and a single (scalar) dispersion
//'   parameter, supporting all four censoring types.
//' @param param  Numeric vector: first \code{ncol(X)} elements are beta
//'   coefficients, the last element is the scalar dispersion parameter.
//' @param X      Design matrix (n x p).
//' @param y_left  Numeric vector of left interval endpoints on (0, 1).
//' @param y_right Numeric vector of right interval endpoints on (0, 1).
//' @param yt     Numeric vector of midpoint response on (0, 1).
//' @param delta  Integer vector of censoring indicators (0,1,2,3).
//' @param link_mu_code  Integer code for the mean link function.
//' @param link_phi_code Integer code for the dispersion link function.
//' @param repar  Integer reparameterization type (0, 1, or 2).
//' @return Scalar log-likelihood value.
//' @keywords internal
// [[Rcpp::export(name = ".brs_loglik_fixed_cpp")]]
double betaregscale_loglik_fixed_cpp(const arma::vec &param, const arma::mat &X,
                                     const arma::vec &y_left,
                                     const arma::vec &y_right,
                                     const arma::vec &yt,
                                     const arma::ivec &delta, int link_mu_code,
                                     int link_phi_code, int repar) {
  int n = X.n_rows;
  int p = X.n_cols;

  arma::vec beta_vec = param.subvec(0, p - 1);
  double phi_raw = param(p);

  arma::vec eta = X * beta_vec;
  double phi = inv_link(phi_raw, link_phi_code);
  phi = clamp_phi_by_repar(phi, repar);

  double ll = 0.0;
  for (int i = 0; i < n; i++) {
    double mu_i = inv_link(eta(i), link_mu_code);
    mu_i = clamp(mu_i, EPS_BOUND, 1.0 - EPS_BOUND);

    double a, b;
    beta_shapes(mu_i, phi, repar, a, b);

    ll += obs_loglik(delta(i), y_left(i), y_right(i), yt(i), a, b);
  }

  return ll;
}

// ========================================================= Variable phi === //

//' @title C++ log-likelihood for variable-dispersion beta interval regression
//'   with mixed censoring
//' @description Computes the total log-likelihood for a beta regression model
//'   with interval-censored responses and observation-specific dispersion,
//'   supporting all four censoring types.
//' @param param Numeric vector: first \code{ncol(X)} elements are beta
//'   coefficients, next \code{ncol(Z)} elements are gamma (phi) coefficients.
//' @param X      Design matrix for the mean submodel (n x p).
//' @param Z      Design matrix for the dispersion submodel (n x q).
//' @param y_left  Numeric vector of left interval endpoints on (0, 1).
//' @param y_right Numeric vector of right interval endpoints on (0, 1).
//' @param yt     Numeric vector of midpoint response on (0, 1).
//' @param delta  Integer vector of censoring indicators (0,1,2,3).
//' @param link_mu_code  Integer code for the mean link function.
//' @param link_phi_code Integer code for the dispersion link function.
//' @param repar  Integer reparameterization type (0, 1, or 2).
//' @return Scalar log-likelihood value.
//' @keywords internal
// [[Rcpp::export(name = ".brs_loglik_variable_cpp")]]
double betaregscale_loglik_variable_cpp(
    const arma::vec &param, const arma::mat &X, const arma::mat &Z,
    const arma::vec &y_left, const arma::vec &y_right, const arma::vec &yt,
    const arma::ivec &delta, int link_mu_code, int link_phi_code, int repar) {
  int n = X.n_rows;
  int p = X.n_cols;
  int q = Z.n_cols;

  arma::vec beta_vec = param.subvec(0, p - 1);
  arma::vec gamma_vec = param.subvec(p, p + q - 1);

  arma::vec eta_mu = X * beta_vec;
  arma::vec eta_phi = Z * gamma_vec;

  double ll = 0.0;
  for (int i = 0; i < n; i++) {
    double mu_i = inv_link(eta_mu(i), link_mu_code);
    double phi_i = inv_link(eta_phi(i), link_phi_code);

    mu_i = clamp(mu_i, EPS_BOUND, 1.0 - EPS_BOUND);
    phi_i = clamp_phi_by_repar(phi_i, repar);

    double a, b;
    beta_shapes(mu_i, phi_i, repar, a, b);

    ll += obs_loglik(delta(i), y_left(i), y_right(i), yt(i), a, b);
  }

  return ll;
}

// ============================================ Numerical gradient (fixed) ===
// //

//' @title C++ gradient for fixed-dispersion log-likelihood
//' @description Returns the gradient vector of the log-likelihood with
//'   respect to all parameters (beta coefficients + scalar phi), using
//'   a central-difference numerical approximation (step = 1e-6).
//' @param param  Parameter vector (same layout as loglik function).
//' @param X      Design matrix (n x p).
//' @param y_left  Left endpoints.
//' @param y_right Right endpoints.
//' @param yt     Midpoint responses.
//' @param delta  Integer censoring indicators.
//' @param link_mu_code  Integer mean link code.
//' @param link_phi_code Integer dispersion link code.
//' @param repar  Integer reparameterization type.
//' @return Numeric gradient vector of length \code{ncol(X) + 1}.
//' @keywords internal
// [[Rcpp::export(name = ".brs_grad_fixed_cpp")]]
arma::vec
betaregscale_grad_fixed_cpp(const arma::vec &param, const arma::mat &X,
                            const arma::vec &y_left, const arma::vec &y_right,
                            const arma::vec &yt, const arma::ivec &delta,
                            int link_mu_code, int link_phi_code, int repar) {
  int npar = param.n_elem;
  arma::vec grad(npar, arma::fill::zeros);
  double h = 1.0e-6;

  for (int j = 0; j < npar; j++) {
    arma::vec p_plus = param;
    arma::vec p_minus = param;
    double step = std::max(h, h * std::abs(param(j)));
    p_plus(j) += step;
    p_minus(j) -= step;

    double f_plus =
        betaregscale_loglik_fixed_cpp(p_plus, X, y_left, y_right, yt, delta,
                                      link_mu_code, link_phi_code, repar);
    double f_minus =
        betaregscale_loglik_fixed_cpp(p_minus, X, y_left, y_right, yt, delta,
                                      link_mu_code, link_phi_code, repar);

    grad(j) = (f_plus - f_minus) / (2.0 * step);
  }

  return grad;
}

// ========================================= Numerical gradient (variable) ===
// //

//' @title C++ gradient for variable-dispersion log-likelihood
//' @description Central-difference gradient for the variable-dispersion model.
//' @param param  Parameter vector.
//' @param X      Mean design matrix.
//' @param Z      Dispersion design matrix.
//' @param y_left  Left endpoints.
//' @param y_right Right endpoints.
//' @param yt     Midpoint responses.
//' @param delta  Integer censoring indicators.
//' @param link_mu_code  Integer mean link code.
//' @param link_phi_code Integer dispersion link code.
//' @param repar  Integer reparameterization type.
//' @return Numeric gradient vector.
//' @keywords internal
// [[Rcpp::export(name = ".brs_grad_variable_cpp")]]
arma::vec betaregscale_grad_variable_cpp(
    const arma::vec &param, const arma::mat &X, const arma::mat &Z,
    const arma::vec &y_left, const arma::vec &y_right, const arma::vec &yt,
    const arma::ivec &delta, int link_mu_code, int link_phi_code, int repar) {
  int npar = param.n_elem;
  arma::vec grad(npar, arma::fill::zeros);
  double h = 1.0e-6;

  for (int j = 0; j < npar; j++) {
    arma::vec p_plus = param;
    arma::vec p_minus = param;
    double step = std::max(h, h * std::abs(param(j)));
    p_plus(j) += step;
    p_minus(j) -= step;

    double f_plus = betaregscale_loglik_variable_cpp(
        p_plus, X, Z, y_left, y_right, yt, delta, link_mu_code, link_phi_code,
        repar);
    double f_minus = betaregscale_loglik_variable_cpp(
        p_minus, X, Z, y_left, y_right, yt, delta, link_mu_code, link_phi_code,
        repar);

    grad(j) = (f_plus - f_minus) / (2.0 * step);
  }

  return grad;
}

// ============================================ Mixed model (Laplace, q=1) ===

struct GroupLaplaceResult {
  double mode_b;
  double sd_b;
  double log_li;
};

inline void build_group_index(const arma::ivec &group,
                              std::vector<std::vector<int>> &group_obs) {
  int n = group.n_elem;
  int gmax = 0;
  for (int i = 0; i < n; i++) {
    if (group(i) > gmax) {
      gmax = group(i);
    }
  }
  group_obs.assign(gmax, std::vector<int>());
  for (int i = 0; i < n; i++) {
    int gi = group(i);
    if (gi >= 1 && gi <= gmax) {
      group_obs[gi - 1].push_back(i);
    }
  }
}

inline double group_Q(
    double b,
    const std::vector<int> &idx,
    const arma::vec &eta_mu_base,
    const arma::vec &eta_phi,
    const arma::vec &y_left,
    const arma::vec &y_right,
    const arma::vec &yt,
    const arma::ivec &delta,
    double sigma_b,
    int link_mu_code,
    int link_phi_code,
    int repar) {
  double ll = 0.0;

  for (size_t k = 0; k < idx.size(); k++) {
    int ii = idx[k];
    double mu_i = inv_link(eta_mu_base(ii) + b, link_mu_code);
    double phi_i = inv_link(eta_phi(ii), link_phi_code);

    mu_i = clamp(mu_i, EPS_BOUND, 1.0 - EPS_BOUND);
    phi_i = clamp_phi_by_repar(phi_i, repar);

    double a, bb;
    beta_shapes(mu_i, phi_i, repar, a, bb);
    ll += obs_loglik(delta(ii), y_left(ii), y_right(ii), yt(ii), a, bb);
  }

  ll += R::dnorm4(b, 0.0, sigma_b, 1);
  return std::isfinite(ll) ? ll : LOG_PENALTY;
}

inline double golden_max_group(
    const std::vector<int> &idx,
    const arma::vec &eta_mu_base,
    const arma::vec &eta_phi,
    const arma::vec &y_left,
    const arma::vec &y_right,
    const arma::vec &yt,
    const arma::ivec &delta,
    double sigma_b,
    int link_mu_code,
    int link_phi_code,
    int repar,
    double lower,
    double upper,
    int maxit = 80,
    double tol = 1.0e-8) {
  const double gr = 0.6180339887498948482;
  double a = lower;
  double b = upper;
  double c = b - gr * (b - a);
  double d = a + gr * (b - a);

  double fc =
      group_Q(c, idx, eta_mu_base, eta_phi, y_left, y_right, yt, delta,
              sigma_b, link_mu_code, link_phi_code, repar);
  double fd =
      group_Q(d, idx, eta_mu_base, eta_phi, y_left, y_right, yt, delta,
              sigma_b, link_mu_code, link_phi_code, repar);

  for (int it = 0; it < maxit && std::abs(b - a) > tol; it++) {
    if (fc > fd) {
      b = d;
      d = c;
      fd = fc;
      c = b - gr * (b - a);
      fc = group_Q(c, idx, eta_mu_base, eta_phi, y_left, y_right, yt, delta,
                   sigma_b, link_mu_code, link_phi_code, repar);
    } else {
      a = c;
      c = d;
      fc = fd;
      d = a + gr * (b - a);
      fd = group_Q(d, idx, eta_mu_base, eta_phi, y_left, y_right, yt, delta,
                   sigma_b, link_mu_code, link_phi_code, repar);
    }
  }

  return (fc > fd) ? c : d;
}

inline GroupLaplaceResult laplace_group(
    const std::vector<int> &idx,
    const arma::vec &eta_mu_base,
    const arma::vec &eta_phi,
    const arma::vec &y_left,
    const arma::vec &y_right,
    const arma::vec &yt,
    const arma::ivec &delta,
    double sigma_b,
    int link_mu_code,
    int link_phi_code,
    int repar) {
  double bound = std::max(3.0, 8.0 * sigma_b + 1.0);
  double b_hat =
      golden_max_group(idx, eta_mu_base, eta_phi, y_left, y_right, yt, delta,
                       sigma_b, link_mu_code, link_phi_code, repar, -bound,
                       bound, 100, 1.0e-8);

  double h = std::max(1.0e-5, 1.0e-4 * std::max(1.0, std::abs(b_hat)));
  double q0 =
      group_Q(b_hat, idx, eta_mu_base, eta_phi, y_left, y_right, yt, delta,
              sigma_b, link_mu_code, link_phi_code, repar);
  double qp =
      group_Q(b_hat + h, idx, eta_mu_base, eta_phi, y_left, y_right, yt, delta,
              sigma_b, link_mu_code, link_phi_code, repar);
  double qm =
      group_Q(b_hat - h, idx, eta_mu_base, eta_phi, y_left, y_right, yt, delta,
              sigma_b, link_mu_code, link_phi_code, repar);

  double d2 = (qp - 2.0 * q0 + qm) / (h * h);
  double curv = -d2;
  if (!std::isfinite(curv) || curv <= EPS_PROB) {
    curv = EPS_PROB;
  }

  double log_li = q0 + 0.5 * std::log(2.0 * M_PI) - 0.5 * std::log(curv);
  if (!std::isfinite(log_li)) {
    log_li = LOG_PENALTY;
  }

  GroupLaplaceResult out;
  out.mode_b = b_hat;
  out.sd_b = 1.0 / std::sqrt(curv);
  out.log_li = log_li;
  return out;
}

//' @title C++ log-likelihood for mixed beta interval regression (Laplace)
//' @description Computes the marginal log-likelihood for a random-intercept
//'   beta interval model using a Laplace approximation of group-specific
//'   integrals over Gaussian random effects.
//' @param param Parameter vector:
//'   \code{beta (p), gamma (q), log_sigma_b (1)}.
//' @param X Mean-model design matrix.
//' @param Z Precision-model design matrix.
//' @param y_left Left interval endpoints.
//' @param y_right Right interval endpoints.
//' @param yt Midpoint responses.
//' @param delta Censoring indicator vector.
//' @param group Integer group index (1..G), one per observation.
//' @param link_mu_code Integer code for mean link.
//' @param link_phi_code Integer code for precision link.
//' @param repar Integer beta reparameterization code.
//' @return Scalar marginal log-likelihood.
//' @keywords internal
// [[Rcpp::export(name = ".brsmm_loglik_laplace_cpp")]]
double betaregscale_loglik_mixed_laplace_cpp(
    const arma::vec &param, const arma::mat &X, const arma::mat &Z,
    const arma::vec &y_left, const arma::vec &y_right, const arma::vec &yt,
    const arma::ivec &delta, const arma::ivec &group, int link_mu_code,
    int link_phi_code, int repar) {
  int p = X.n_cols;
  int q = Z.n_cols;
  if (param.n_elem != static_cast<unsigned int>(p + q + 1)) {
    return LOG_PENALTY;
  }

  arma::vec beta_vec = param.subvec(0, p - 1);
  arma::vec gamma_vec = param.subvec(p, p + q - 1);
  double sigma_b = std::exp(param(p + q));
  sigma_b = clamp(sigma_b, 1.0e-8, 1.0e8);

  arma::vec eta_mu_base = X * beta_vec;
  arma::vec eta_phi = Z * gamma_vec;

  std::vector<std::vector<int>> group_obs;
  build_group_index(group, group_obs);

  double ll = 0.0;
  for (size_t g = 0; g < group_obs.size(); g++) {
    if (group_obs[g].empty()) {
      continue;
    }
    GroupLaplaceResult res =
        laplace_group(group_obs[g], eta_mu_base, eta_phi, y_left, y_right, yt,
                      delta, sigma_b, link_mu_code, link_phi_code, repar);
    ll += res.log_li;
  }

  return ll;
}

//' @title Group-level modes for random intercept (Laplace)
//' @description Computes the posterior mode and local SD for each group random
//'   intercept under the current parameter vector.
//' @param param Parameter vector:
//'   \code{beta (p), gamma (q), log_sigma_b (1)}.
//' @param X Mean-model design matrix.
//' @param Z Precision-model design matrix.
//' @param y_left Left interval endpoints.
//' @param y_right Right interval endpoints.
//' @param yt Midpoint responses.
//' @param delta Censoring indicator vector.
//' @param group Integer group index (1..G), one per observation.
//' @param link_mu_code Integer code for mean link.
//' @param link_phi_code Integer code for precision link.
//' @param repar Integer beta reparameterization code.
//' @return Numeric matrix with columns \code{mode_b} and \code{sd_b}.
//' @keywords internal
// [[Rcpp::export(name = ".brsmm_group_modes_cpp")]]
arma::mat betaregscale_group_modes_cpp(
    const arma::vec &param, const arma::mat &X, const arma::mat &Z,
    const arma::vec &y_left, const arma::vec &y_right, const arma::vec &yt,
    const arma::ivec &delta, const arma::ivec &group, int link_mu_code,
    int link_phi_code, int repar) {
  int p = X.n_cols;
  int q = Z.n_cols;
  arma::mat out;
  if (param.n_elem != static_cast<unsigned int>(p + q + 1)) {
    return out;
  }

  arma::vec beta_vec = param.subvec(0, p - 1);
  arma::vec gamma_vec = param.subvec(p, p + q - 1);
  double sigma_b = std::exp(param(p + q));
  sigma_b = clamp(sigma_b, 1.0e-8, 1.0e8);

  arma::vec eta_mu_base = X * beta_vec;
  arma::vec eta_phi = Z * gamma_vec;

  std::vector<std::vector<int>> group_obs;
  build_group_index(group, group_obs);
  int G = static_cast<int>(group_obs.size());
  out = arma::mat(G, 2, arma::fill::zeros);

  for (int g = 0; g < G; g++) {
    if (group_obs[g].empty()) {
      out(g, 0) = 0.0;
      out(g, 1) = NA_REAL;
      continue;
    }
    GroupLaplaceResult res =
        laplace_group(group_obs[g], eta_mu_base, eta_phi, y_left, y_right, yt,
                      delta, sigma_b, link_mu_code, link_phi_code, repar);
    out(g, 0) = res.mode_b;
    out(g, 1) = res.sd_b;
  }

  return out;
}
