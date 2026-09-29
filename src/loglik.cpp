// betaregscale — compiled likelihood, gradient and Hessian for the fixed- and
// variable-dispersion beta interval regression (censoring codes 0..3, see
// obs_loglik in brs_common.h). Gradient and Hessian use the chain rule on the
// linear predictors: per-observation derivatives from brs_deriv.h, then
// X' d_mu, Z' d_phi and X' W X; the cost does not grow with p or q.

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "brs_common.h"
#include "brs_deriv.h"
#include <algorithm>

// Structural checks shared by all exports (no bounds checks under -DNDEBUG):
// sizes, param length (short or extra), NA / non-finite data, delta codes.
static void check_brs_inputs(const arma::vec &param, const arma::mat &X,
                             const arma::mat *Z, const arma::vec &y_left,
                             const arma::vec &y_right, const arma::vec &yt,
                             const arma::ivec &delta) {
  const arma::uword n = X.n_rows;
  if (n < 1) Rcpp::stop("brs: no observations.");
  if (X.n_cols < 1) Rcpp::stop("brs: X must have at least one column.");
  if (y_left.n_elem != n || y_right.n_elem != n || yt.n_elem != n ||
      delta.n_elem != n || (Z && Z->n_rows != n))
    Rcpp::stop("brs: X, Z, y_left, y_right, yt and delta must all have %d rows.",
               (int)n);
  if (!X.is_finite() || (Z && !Z->is_finite()) || !y_left.is_finite() ||
      !y_right.is_finite() || !yt.is_finite())
    Rcpp::stop("brs: X, Z, y_left, y_right and yt must not contain NA or non-finite values.");
  arma::uword k = X.n_cols + (Z ? Z->n_cols : 1);
  if (param.n_elem != k)
    Rcpp::stop("brs: param must have length %d (got %d).", (int)k,
               (int)param.n_elem);
  for (arma::uword i = 0; i < n; ++i)
    if (delta(i) < 0 || delta(i) > 3)
      Rcpp::stop("brs: delta must be in {0,1,2,3} (found %d at row %d).",
                 (int)delta(i), (int)i + 1);
}

// ============================================================ Fixed phi === //

//' @title C++ log-likelihood for fixed-dispersion beta interval regression
//' @description Total log-likelihood with a single dispersion parameter
//'   (last element of \code{param}); all four censoring types.
//' @param param  Numeric vector: \code{ncol(X)} coefficients, then phi (link scale).
//' @param X      Design matrix (n x p).
//' @param y_left,y_right Interval endpoints on (0, 1).
//' @param yt     Exact response on (0, 1) (used when \code{delta = 0}).
//' @param delta  Integer censoring indicators (0,1,2,3).
//' @param link_mu_code,link_phi_code Integer link codes (see \code{link_to_code}).
//' @param repar  Integer reparameterization type (0, 1, or 2).
//' @return Scalar log-likelihood value.
//' @keywords internal
// [[Rcpp::export(name = ".brs_loglik_fixed_cpp", rng = false)]]
double betaregscale_loglik_fixed_cpp(const arma::vec &param, const arma::mat &X,
                                     const arma::vec &y_left,
                                     const arma::vec &y_right,
                                     const arma::vec &yt,
                                     const arma::ivec &delta, int link_mu_code,
                                     int link_phi_code, int repar) {
  check_brs_inputs(param, X, nullptr, y_left, y_right, yt, delta);
  const int n = X.n_rows, p = X.n_cols;
  arma::vec eta = X * param.head(p);
  double phi = clamp_phi_by_repar(inv_link(param(p), link_phi_code), repar);
  double ll = 0.0;
  for (int i = 0; i < n; i++) {
    double mu_i = clamp_mu_by_repar(inv_link(eta(i), link_mu_code), repar);
    double a, b;
    beta_shapes(mu_i, phi, repar, a, b);
    ll += obs_loglik(delta(i), y_left(i), y_right(i), yt(i), a, b);
  }
  return ll;
}

// ========================================================= Variable phi === //

//' @title C++ log-likelihood for variable-dispersion beta interval regression
//' @description Total log-likelihood with observation-specific dispersion
//'   \code{Z gamma}; all four censoring types.
//' @param param Numeric vector: \code{ncol(X)} beta then \code{ncol(Z)} gamma.
//' @param X,Z   Design matrices of the mean (n x p) and dispersion (n x q).
//' @param y_left,y_right Interval endpoints on (0, 1).
//' @param yt     Exact response on (0, 1) (used when \code{delta = 0}).
//' @param delta  Integer censoring indicators (0,1,2,3).
//' @param link_mu_code,link_phi_code Integer link codes.
//' @param repar  Integer reparameterization type (0, 1, or 2).
//' @return Scalar log-likelihood value.
//' @keywords internal
// [[Rcpp::export(name = ".brs_loglik_variable_cpp", rng = false)]]
double betaregscale_loglik_variable_cpp(
    const arma::vec &param, const arma::mat &X, const arma::mat &Z,
    const arma::vec &y_left, const arma::vec &y_right, const arma::vec &yt,
    const arma::ivec &delta, int link_mu_code, int link_phi_code, int repar) {
  check_brs_inputs(param, X, &Z, y_left, y_right, yt, delta);
  const int n = X.n_rows, p = X.n_cols, q = Z.n_cols;
  arma::vec eta_mu = X * param.head(p);
  arma::vec eta_phi = Z * param.subvec(p, p + q - 1);
  double ll = 0.0;
  for (int i = 0; i < n; i++) {
    double mu_i  = clamp_mu_by_repar(inv_link(eta_mu(i), link_mu_code), repar);
    double phi_i = clamp_phi_by_repar(inv_link(eta_phi(i), link_phi_code), repar);
    double a, b;
    beta_shapes(mu_i, phi_i, repar, a, b);
    ll += obs_loglik(delta(i), y_left(i), y_right(i), yt(i), a, b);
  }
  return ll;
}

// ============================================================= Gradients === //

//' @title C++ gradient for the fixed-dispersion log-likelihood
//' @description Chain rule on the linear predictors:
//'   \code{crossprod(X, d_mu)} and \code{sum(d_phi)}, with per-observation
//'   derivatives by Richardson central differences (brs_deriv.h).
//' @inheritParams .brs_loglik_fixed_cpp
//' @return Numeric gradient vector of length \code{ncol(X) + 1}.
//' @keywords internal
// [[Rcpp::export(name = ".brs_grad_fixed_cpp", rng = false)]]
arma::vec
betaregscale_grad_fixed_cpp(const arma::vec &param, const arma::mat &X,
                            const arma::vec &y_left, const arma::vec &y_right,
                            const arma::vec &yt, const arma::ivec &delta,
                            int link_mu_code, int link_phi_code, int repar) {
  check_brs_inputs(param, X, nullptr, y_left, y_right, yt, delta);
  const int n = X.n_rows, p = X.n_cols;
  arma::vec eta = X * param.head(p);
  const double eta_phi = param(p);
  const LinkSpec s{link_mu_code, link_phi_code, repar};
  arma::vec dmu(n);
  double dphi = 0.0;
  for (int i = 0; i < n; ++i) {
    ObsSpec o{(int)delta(i), y_left(i), y_right(i), yt(i)};
    double d1, p1;
    obs_deriv_grad(eta(i), eta_phi, o, s, d1, p1);
    dmu(i) = d1;
    dphi += p1;
  }
  arma::vec g(p + 1);
  g.head(p) = X.t() * dmu;
  g(p) = dphi;
  return g;
}

//' @title C++ gradient for the variable-dispersion log-likelihood
//' @description Chain rule on the linear predictors:
//'   \code{crossprod(X, d_mu)} and \code{crossprod(Z, d_phi)}
//'   (see \code{.brs_grad_fixed_cpp}).
//' @inheritParams .brs_loglik_variable_cpp
//' @return Numeric gradient vector of length \code{ncol(X) + ncol(Z)}.
//' @keywords internal
// [[Rcpp::export(name = ".brs_grad_variable_cpp", rng = false)]]
arma::vec betaregscale_grad_variable_cpp(
    const arma::vec &param, const arma::mat &X, const arma::mat &Z,
    const arma::vec &y_left, const arma::vec &y_right, const arma::vec &yt,
    const arma::ivec &delta, int link_mu_code, int link_phi_code, int repar) {
  check_brs_inputs(param, X, &Z, y_left, y_right, yt, delta);
  const int n = X.n_rows, p = X.n_cols, q = Z.n_cols;
  arma::vec eta_mu = X * param.head(p);
  arma::vec eta_phi = Z * param.subvec(p, p + q - 1);
  const LinkSpec s{link_mu_code, link_phi_code, repar};
  arma::vec dmu(n), dphi(n);
  for (int i = 0; i < n; ++i) {
    ObsSpec o{(int)delta(i), y_left(i), y_right(i), yt(i)};
    obs_deriv_grad(eta_mu(i), eta_phi(i), o, s, dmu(i), dphi(i));
  }
  return arma::join_cols(X.t() * dmu, Z.t() * dphi);
}

// ============================================================== Hessians === //

//' @title C++ Hessian for the fixed-dispersion log-likelihood
//' @description Chain rule on the linear predictors:
//'   blocks \code{crossprod(X, w_mm * X)}, \code{crossprod(X, w_mp)} and
//'   \code{sum(w_pp)}, with per-observation
//'   second derivatives by Richardson central differences (17 evaluations).
//' @inheritParams .brs_loglik_fixed_cpp
//' @return Symmetric matrix of order \code{ncol(X) + 1} (log-likelihood scale).
//' @keywords internal
// [[Rcpp::export(name = ".brs_hessian_fixed_cpp", rng = false)]]
arma::mat betaregscale_hessian_fixed_cpp(
    const arma::vec &param, const arma::mat &X, const arma::vec &y_left,
    const arma::vec &y_right, const arma::vec &yt, const arma::ivec &delta,
    int link_mu_code, int link_phi_code, int repar) {
  check_brs_inputs(param, X, nullptr, y_left, y_right, yt, delta);
  const int n = X.n_rows, p = X.n_cols;
  arma::vec eta = X * param.head(p);
  const double eta_phi = param(p);
  const LinkSpec s{link_mu_code, link_phi_code, repar};
  arma::vec wmm(n), wmp(n);
  double wpp = 0.0;
  ObsDeriv r;
  for (int i = 0; i < n; ++i) {
    ObsSpec o{(int)delta(i), y_left(i), y_right(i), yt(i)};
    obs_deriv_hess(eta(i), eta_phi, o, s, r);
    wmm(i) = r.d2;
    wmp(i) = r.c11;
    wpp += r.p2;
  }
  arma::mat H(p + 1, p + 1);
  H.submat(0, 0, p - 1, p - 1) = X.t() * (X.each_col() % wmm);
  H.submat(0, p, p - 1, p) = X.t() * wmp;
  H.submat(p, 0, p, p - 1) = H.submat(0, p, p - 1, p).t();
  H(p, p) = wpp;
  return H;
}

//' @title C++ Hessian for the variable-dispersion log-likelihood
//' @description Chain rule on the linear predictors:
//'   blocks \code{crossprod(X, w_mm * X)}, \code{crossprod(X, w_mp * Z)} and
//'   \code{crossprod(Z, w_pp * Z)}.
//' @inheritParams .brs_loglik_variable_cpp
//' @return Symmetric matrix of order \code{ncol(X) + ncol(Z)} (log-likelihood scale).
//' @keywords internal
// [[Rcpp::export(name = ".brs_hessian_variable_cpp", rng = false)]]
arma::mat betaregscale_hessian_variable_cpp(
    const arma::vec &param, const arma::mat &X, const arma::mat &Z,
    const arma::vec &y_left, const arma::vec &y_right, const arma::vec &yt,
    const arma::ivec &delta, int link_mu_code, int link_phi_code, int repar) {
  check_brs_inputs(param, X, &Z, y_left, y_right, yt, delta);
  const int n = X.n_rows, p = X.n_cols, q = Z.n_cols;
  arma::vec eta_mu = X * param.head(p);
  arma::vec eta_phi = Z * param.subvec(p, p + q - 1);
  const LinkSpec s{link_mu_code, link_phi_code, repar};
  arma::vec wmm(n), wpp(n), wmp(n);
  ObsDeriv r;
  for (int i = 0; i < n; ++i) {
    ObsSpec o{(int)delta(i), y_left(i), y_right(i), yt(i)};
    obs_deriv_hess(eta_mu(i), eta_phi(i), o, s, r);
    wmm(i) = r.d2;
    wpp(i) = r.p2;
    wmp(i) = r.c11;
  }
  arma::mat H(p + q, p + q);
  H.submat(0, 0, p - 1, p - 1) = X.t() * (X.each_col() % wmm);
  H.submat(p, p, p + q - 1, p + q - 1) = Z.t() * (Z.each_col() % wpp);
  H.submat(0, p, p - 1, p + q - 1) = X.t() * (Z.each_col() % wmp);
  H.submat(p, 0, p + q - 1, p - 1) = H.submat(0, p, p - 1, p + q - 1).t();
  return H;
}
