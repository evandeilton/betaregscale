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
//
// Constants, inverse-link, beta_shapes, and obs_loglik are in brs_common.h.
// -------------------------------------------------------------------------- //

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "brs_common.h"

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
    mu_i = clamp(mu_i, EPS_UNIT, 1.0 - EPS_UNIT);

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

    mu_i = clamp(mu_i, EPS_UNIT, 1.0 - EPS_UNIT);
    phi_i = clamp_phi_by_repar(phi_i, repar);

    double a, b;
    beta_shapes(mu_i, phi_i, repar, a, b);

    ll += obs_loglik(delta(i), y_left(i), y_right(i), yt(i), a, b);
  }

  return ll;
}

// ============================================ Numerical gradient (fixed) === //

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

// ========================================= Numerical gradient (variable) === //

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
