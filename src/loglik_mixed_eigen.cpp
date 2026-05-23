// -------------------------------------------------------------------------- //
// betaregscale — Eigen backend for mixed-effects beta interval regression.
//
// Supports random-effects vectors of arbitrary dimension q_re via:
//   method = 0 : Laplace approximation
//   method = 1 : Adaptive Gauss-Hermite Quadrature (AGHQ)
//   method = 2 : Quasi-Monte Carlo with Halton sequences (QMC)
//
// Random effects are parameterised via the lower-Cholesky factor L of the
// covariance matrix D = L L'.  The packed parameter vector theta_re contains
// the lower triangle of L column-by-column, with diagonal entries stored as
// log(L_ii) to enforce positivity.
//
// Constants, inverse-link, beta_shapes, and obs_loglik are in brs_common.h.
//
// Fixes applied in this revision:
//   BUG-H01  : removed duplicate code, now uses brs_common.h.
//   BUG-C03  : find_mode_vec only accepts Newton steps that improve h.
//   BUG-C07  : find_mode_vec checks bnew.allFinite() before evaluation.
//   BUG-C05  : build_cartesian_grid guards against int overflow for large q.
//   BUG-C06  : PRIMES array extended to 50 entries + dimension guard.
//   NUM-M06  : halton() clamps r away from 0/1 before qnorm call.
//   QUAL-L06 : export names use the internal '.' prefix convention.
//   PERF-M01 : numerical_grad_hess uses a single workspace vector.
// -------------------------------------------------------------------------- //

// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include "brs_common.h"

// ------------------------------------------------------------- RE structs --

struct GroupData {
  Eigen::VectorXd delta;
  Eigen::VectorXd y_left;
  Eigen::VectorXd y_right;
  Eigen::VectorXd yt;
  Eigen::VectorXd eta_mu_fixed;
  Eigen::VectorXd eta_phi;
  Eigen::MatrixXd Zr;
};

struct RandStruct {
  Eigen::MatrixXd L;
  double logdet_D;
  int q_re;
};

// Unpack the packed lower-Cholesky theta_re into L and log|D|.
// Diagonal elements are stored as log(L_ii) in theta_re (log-parameterisation
// ensures L_ii > 0 without constraints during optimisation).
inline RandStruct unpack_re(const Eigen::VectorXd &theta_re, int q_re) {
  RandStruct rs;
  rs.q_re = q_re;
  rs.L = Eigen::MatrixXd::Zero(q_re, q_re);
  int k = 0;
  double logdet_half = 0.0;
  for (int j = 0; j < q_re; ++j) {
    for (int i = j; i < q_re; ++i) {
      double v = theta_re(k++);
      if (i == j) {
        double d = std::exp(v);
        rs.L(i, j) = d;
        logdet_half += std::log(d);
      } else {
        rs.L(i, j) = v;
      }
    }
  }
  rs.logdet_D = 2.0 * logdet_half;
  return rs;
}

// log p(b | D) = log N(b; 0, D)
inline double log_prior_b(const Eigen::VectorXd &b, const RandStruct &rs) {
  Eigen::VectorXd z = rs.L.triangularView<Eigen::Lower>().solve(b);
  double quad = z.squaredNorm();
  double cst  = rs.q_re * std::log(2.0 * M_PI);
  return -0.5 * (cst + rs.logdet_D + quad);
}

// h(b) = sum_i log_lik_i(b) + log p(b | D)
inline double h_func_vec(const Eigen::VectorXd &b, const GroupData &gd,
                         const RandStruct &rs, int link_mu_code,
                         int link_phi_code, int repar) {
  double ll = 0.0;
  int n = gd.delta.size();
  for (int i = 0; i < n; ++i) {
    double eta = gd.eta_mu_fixed(i) + gd.Zr.row(i).dot(b);
    double mu  = clamp(inv_link(eta, link_mu_code), EPS_UNIT, 1.0 - EPS_UNIT);
    double phi = clamp_phi_by_repar(inv_link(gd.eta_phi(i), link_phi_code), repar);
    double a, bb;
    beta_shapes(mu, phi, repar, a, bb);
    ll += obs_loglik((int)gd.delta(i), gd.y_left(i), gd.y_right(i), gd.yt(i),
                     a, bb);
  }
  ll += log_prior_b(b, rs);
  return std::isfinite(ll) ? ll : LOG_PENALTY;
}

// -------------------------------------------------- numerical Hessian/grad --
// PERF-M01: uses a single workspace vector bwork to avoid 4N allocations
// for off-diagonal cross terms.

inline void numerical_grad_hess(const Eigen::VectorXd &b, const GroupData &gd,
                                const RandStruct &rs, int link_mu_code,
                                int link_phi_code, int repar,
                                Eigen::VectorXd &grad, Eigen::MatrixXd &hess) {
  int q = b.size();
  grad = Eigen::VectorXd::Zero(q);
  hess = Eigen::MatrixXd::Zero(q, q);
  double f0 = h_func_vec(b, gd, rs, link_mu_code, link_phi_code, repar);

  std::vector<double> step(q, 0.0);
  for (int j = 0; j < q; ++j)
    step[j] = std::max(1e-5, 1e-4 * std::max(1.0, std::abs(b(j))));

  // Diagonal: gradient and diagonal Hessian entries.
  // Store fp and fm per dimension for potential reuse.
  std::vector<double> fp_d(q), fm_d(q);
  Eigen::VectorXd bwork = b;
  for (int j = 0; j < q; ++j) {
    bwork = b; bwork(j) += step[j];
    fp_d[j] = h_func_vec(bwork, gd, rs, link_mu_code, link_phi_code, repar);
    bwork = b; bwork(j) -= step[j];
    fm_d[j] = h_func_vec(bwork, gd, rs, link_mu_code, link_phi_code, repar);
    grad(j)    = (fp_d[j] - fm_d[j]) / (2.0 * step[j]);
    hess(j, j) = (fp_d[j] - 2.0 * f0 + fm_d[j]) / (step[j] * step[j]);
  }

  // Off-diagonal cross terms via 4-point formula.
  for (int j = 0; j < q; ++j) {
    for (int k = j + 1; k < q; ++k) {
      bwork = b; bwork(j) += step[j]; bwork(k) += step[k];
      double fpp = h_func_vec(bwork, gd, rs, link_mu_code, link_phi_code, repar);
      bwork = b; bwork(j) += step[j]; bwork(k) -= step[k];
      double fpm = h_func_vec(bwork, gd, rs, link_mu_code, link_phi_code, repar);
      bwork = b; bwork(j) -= step[j]; bwork(k) += step[k];
      double fmp = h_func_vec(bwork, gd, rs, link_mu_code, link_phi_code, repar);
      bwork = b; bwork(j) -= step[j]; bwork(k) -= step[k];
      double fmm = h_func_vec(bwork, gd, rs, link_mu_code, link_phi_code, repar);
      double hij = (fpp - fpm - fmp + fmm) / (4.0 * step[j] * step[k]);
      hess(j, k) = hij;
      hess(k, j) = hij;
    }
  }
}

// --------------------------------------------------- mode-finding (Newton) --

struct ModeResult {
  Eigen::VectorXd mode;
  Eigen::MatrixXd curvature;
  double h_at_mode;
};

// Find the mode of h(b) via Newton-Raphson with backtracking line search.
//
// BUG-C03: only accept a step if it improves h (fnew >= fcur).
//          Previous code accepted downhill steps unconditionally.
// BUG-C07: check bnew.allFinite() before calling h_func_vec.
inline ModeResult find_mode_vec(const GroupData &gd, const RandStruct &rs,
                                int link_mu_code, int link_phi_code,
                                int repar) {
  int q = rs.q_re;
  Eigen::VectorXd b = Eigen::VectorXd::Zero(q);
  double fcur = h_func_vec(b, gd, rs, link_mu_code, link_phi_code, repar);

  for (int iter = 0; iter < 50; ++iter) {
    Eigen::VectorXd g;
    Eigen::MatrixXd H;
    numerical_grad_hess(b, gd, rs, link_mu_code, link_phi_code, repar, g, H);

    // Regularise H (negative semidefinite at max) if singular.
    // Subtracting lambda*I pushes eigenvalues more negative, improving
    // conditioning of a negative-definite matrix for the LU solve.
    Eigen::FullPivLU<Eigen::MatrixXd> lu(H);
    if (!lu.isInvertible()) {
      double lambda = 1e-4 * std::max(H.cwiseAbs().maxCoeff(), 1.0);
      H -= lambda * Eigen::MatrixXd::Identity(q, q);
      lu.compute(H);
      if (!lu.isInvertible()) break;
    }

    Eigen::VectorXd step = lu.solve(g);
    if (!step.allFinite()) break;

    // BUG-C07: guard bnew before evaluation.
    Eigen::VectorXd bnew = b - step;
    if (!bnew.allFinite()) break;

    double fnew = h_func_vec(bnew, gd, rs, link_mu_code, link_phi_code, repar);

    // BUG-C03: backtracking — dampen until improvement or give up.
    double damp = 1.0;
    while (fnew < fcur && damp > 1e-4) {
      damp *= 0.5;
      bnew = b - damp * step;
      if (!bnew.allFinite()) { damp = 0.0; break; }
      fnew = h_func_vec(bnew, gd, rs, link_mu_code, link_phi_code, repar);
    }

    // BUG-C03: only accept if we actually improved.
    if (fnew >= fcur) {
      double delta_b = (bnew - b).norm();
      double delta_f = fnew - fcur;
      b    = bnew;
      fcur = fnew;
      if (delta_b < 1e-7 || delta_f < 1e-10) break;
    } else {
      break;  // no improvement possible — converged
    }
  }

  // Final Hessian at converged mode for curvature estimate.
  Eigen::VectorXd g_final;
  Eigen::MatrixXd H_final;
  numerical_grad_hess(b, gd, rs, link_mu_code, link_phi_code, repar,
                      g_final, H_final);
  Eigen::MatrixXd curv = -H_final;
  curv = 0.5 * (curv + curv.transpose());

  // Project to positive-semidefinite (eigenvalue floor).
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(curv);
  Eigen::VectorXd ev = es.eigenvalues().cwiseMax(1e-8);
  curv = es.eigenvectors() * ev.asDiagonal() * es.eigenvectors().transpose();

  ModeResult out;
  out.mode       = b;
  out.curvature  = curv;
  out.h_at_mode  = fcur;
  return out;
}

// ------------------------------------------- Gauss-Hermite quadrature rule --

void compute_gh_rule(int n, std::vector<double> &x, std::vector<double> &w) {
  Eigen::MatrixXd J = Eigen::MatrixXd::Zero(n, n);
  for (int i = 0; i < n - 1; ++i) {
    double val = std::sqrt((double)(i + 1) / 2.0);
    J(i, i + 1) = val;
    J(i + 1, i) = val;
  }
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(J);
  x.resize(n);
  w.resize(n);
  double sqrt_pi = std::sqrt(M_PI);
  for (int i = 0; i < n; ++i) {
    x[i] = es.eigenvalues()(i);
    w[i] = std::pow(es.eigenvectors()(0, i), 2) * sqrt_pi;
  }
}

// -------------------------------------------------- group data partitioning --

inline std::vector<GroupData>
build_groups(const Eigen::VectorXd &eta_mu, const Eigen::VectorXd &eta_phi,
             const Eigen::MatrixXd &Xr, const Eigen::VectorXd &y_left,
             const Eigen::VectorXd &y_right, const Eigen::VectorXd &yt,
             const Eigen::VectorXi &delta, const Eigen::VectorXi &group) {
  int G = group.maxCoeff();
  std::vector<int> counts(G, 0);
  for (int i = 0; i < group.size(); ++i)
    counts[group(i) - 1]++;

  std::vector<GroupData> groups(G);
  for (int g = 0; g < G; ++g) {
    int n = counts[g];
    groups[g].delta.resize(n);
    groups[g].y_left.resize(n);
    groups[g].y_right.resize(n);
    groups[g].yt.resize(n);
    groups[g].eta_mu_fixed.resize(n);
    groups[g].eta_phi.resize(n);
    groups[g].Zr.resize(n, Xr.cols());
  }

  std::vector<int> cur(G, 0);
  for (int i = 0; i < group.size(); ++i) {
    int g = group(i) - 1;
    int k = cur[g]++;
    groups[g].delta(k)        = delta(i);
    groups[g].y_left(k)       = y_left(i);
    groups[g].y_right(k)      = y_right(i);
    groups[g].yt(k)           = yt(i);
    groups[g].eta_mu_fixed(k) = eta_mu(i);
    groups[g].eta_phi(k)      = eta_phi(i);
    groups[g].Zr.row(k)       = Xr.row(i);
  }
  return groups;
}

// -------------------------------------------------------- Halton sequences --
// BUG-C06: extended to 50 primes; added dimension guard.
// The first 50 primes, used as bases for the Halton low-discrepancy sequence.
static const int PRIMES[50] = {
    2,   3,   5,   7,  11,  13,  17,  19,  23,  29,
   31,  37,  41,  43,  47,  53,  59,  61,  67,  71,
   73,  79,  83,  89,  97, 101, 103, 107, 109, 113,
  127, 131, 137, 139, 149, 151, 157, 163, 167, 173,
  179, 181, 191, 193, 197, 199, 211, 223, 227, 229
};
static const int N_PRIMES = 50;

// Van der Corput sequence in the given integer base.
// NUM-M06: result is clamped away from 0 and 1 before the qnorm call in the
//          caller to prevent ±Inf quantile values for extreme Halton points.
double halton(int index, int base) {
  double f = 1.0, r = 0.0;
  int i = index;
  while (i > 0) {
    f = f / base;
    r = r + f * (i % base);
    i = i / base;
  }
  return r;
}

// ------------------------------------------------- quadrature grid builders --

// BUG-C05: use int64_t to detect overflow before allocation; hard limit to
// prevent accidental allocation of multi-GB grids.
void build_cartesian_grid(const std::vector<double> &x,
                          const std::vector<double> &w, int q,
                          Eigen::MatrixXd &grid, std::vector<double> &out_w) {
  int n = (int)x.size();
  int64_t total64 = 1;
  for (int d = 0; d < q; ++d) {
    total64 *= n;
    if (total64 > 500000LL)
      Rcpp::stop(
        "AGHQ Cartesian grid would have %lld points (n_points^q_re = %d^%d). "
        "Reduce 'n_points' or use int_method='qmc'.",
        (long long)total64, n, q);
  }
  int total = (int)total64;
  grid.resize(total, q);
  out_w.resize(total);
  for (int i = 0; i < total; ++i) {
    int temp = i;
    double w_prod = 1.0;
    for (int d = 0; d < q; ++d) {
      int idx  = temp % n;
      grid(i, d) = x[idx];
      w_prod   *= w[idx];
      temp     /= n;
    }
    out_w[i] = w_prod;
  }
}

Eigen::MatrixXd build_halton_grid(int n_points, int q) {
  if (q > N_PRIMES)
    Rcpp::stop("QMC supports at most %d random-effect dimensions (q_re=%d).",
               N_PRIMES, q);
  Eigen::MatrixXd grid(n_points, q);
  for (int k = 0; k < n_points; ++k) {
    for (int d = 0; d < q; ++d) {
      int base = PRIMES[d];
      // NUM-M06: clamp before qnorm to avoid ±Inf
      double r = halton(k + 1, base);
      r = std::min(std::max(r, 1e-9), 1.0 - 1e-9);
      grid(k, d) = R::qnorm(r, 0.0, 1.0, 1, 0);
    }
  }
  return grid;
}

// ====================================================== Main exported fn === //

// Computes marginal log-likelihood using Laplace, AGHQ, or QMC.
// param: [beta, gamma, theta_re]; method: 0=Laplace, 1=AGHQ, 2=QMC
// Internal. QUAL-L06: dot-prefix naming convention.
// [[Rcpp::export(name = ".brsmm_loglik_eigen")]]
double brsmm_loglik_eigen(Eigen::VectorXd param, Eigen::MatrixXd X,
                          Eigen::MatrixXd Z, Eigen::MatrixXd Xr,
                          Eigen::VectorXd y_left, Eigen::VectorXd y_right,
                          Eigen::VectorXd yt, Eigen::VectorXi delta,
                          Eigen::VectorXi group, int link_mu, int link_phi,
                          int repar, int method, int n_points) {
  int p    = X.cols();
  int q_phi = Z.cols();
  int q_re = Xr.cols();
  int k_re = q_re * (q_re + 1) / 2;
  if (param.size() != (p + q_phi + k_re)) return LOG_PENALTY;

  Eigen::VectorXd beta     = param.head(p);
  Eigen::VectorXd gamma    = param.segment(p, q_phi);
  Eigen::VectorXd theta_re = param.tail(k_re);
  RandStruct rs = unpack_re(theta_re, q_re);

  Eigen::VectorXd eta_mu  = X * beta;
  Eigen::VectorXd eta_phi = Z * gamma;
  std::vector<GroupData> groups =
      build_groups(eta_mu, eta_phi, Xr, y_left, y_right, yt, delta, group);

  std::vector<double> gh_x, gh_w;
  Eigen::MatrixXd aghq_grid;
  std::vector<double> aghq_w;
  Eigen::MatrixXd qmc_grid;

  if (method == 1) {
    compute_gh_rule(n_points, gh_x, gh_w);
    build_cartesian_grid(gh_x, gh_w, q_re, aghq_grid, aghq_w);
  } else if (method == 2) {
    qmc_grid = build_halton_grid(n_points, q_re);
  }

  double total = 0.0;
  for (size_t g = 0; g < groups.size(); ++g) {
    if (groups[g].delta.size() == 0) continue;
    ModeResult mr = find_mode_vec(groups[g], rs, link_mu, link_phi, repar);

    if (method == 0) {
      // Laplace: log∫exp(h(b))db ≈ h(b̂) + (q/2)log(2π) - (1/2)log|curv|
      Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(mr.curvature);
      Eigen::VectorXd ev = es.eigenvalues().cwiseMax(1e-8);
      double logdet = ev.array().log().sum();
      total += mr.h_at_mode + 0.5 * q_re * std::log(2.0 * M_PI) - 0.5 * logdet;

    } else if (method == 1) {
      // AGHQ: adaptive quadrature centred at mode, scaled by curvature.
      Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(mr.curvature);
      Eigen::VectorXd ev = es.eigenvalues().cwiseMax(1e-8);
      Eigen::MatrixXd S  =
          es.eigenvectors() * ev.cwiseInverse().cwiseSqrt().asDiagonal();

      int total_pts = aghq_grid.rows();
      std::vector<double> lt(total_pts);
      for (int k = 0; k < total_pts; ++k) {
        Eigen::VectorXd z  = aghq_grid.row(k).transpose();
        Eigen::VectorXd bk = mr.mode + M_SQRT2 * S * z;
        double h = h_func_vec(bk, groups[g], rs, link_mu, link_phi, repar);
        // log-sum-exp accumulator includes the GH weight and the e^{z'z}
        // factor that removes the Gaussian weight from the quadrature rule.
        lt[k] = (aghq_w[k] > 0.0 ? std::log(aghq_w[k]) : -1e15)
                + h + z.squaredNorm();
      }
      double m = *std::max_element(lt.begin(), lt.end());
      double s = 0.0;
      for (int k = 0; k < total_pts; ++k) s += std::exp(lt[k] - m);

      double logdet_S = -0.5 * ev.array().log().sum();
      total += std::log(s) + m + q_re * std::log(M_SQRT2) + logdet_S;

    } else {
      // QMC: importance sampling with Gaussian proposal centred at mode.
      Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(mr.curvature);
      Eigen::VectorXd ev = es.eigenvalues().cwiseMax(1e-8);
      Eigen::MatrixXd S  =
          es.eigenvectors() * ev.cwiseInverse().cwiseSqrt().asDiagonal();

      int total_pts = qmc_grid.rows();
      std::vector<double> lr(total_pts);
      double m        = -1e12;
      double logdet_S = -0.5 * ev.array().log().sum();
      double cst      = -0.5 * q_re * std::log(2.0 * M_PI) - logdet_S;

      for (int k = 0; k < total_pts; ++k) {
        Eigen::VectorXd z  = qmc_grid.row(k).transpose();
        Eigen::VectorXd bk = mr.mode + S * z;
        double h   = h_func_vec(bk, groups[g], rs, link_mu, link_phi, repar);
        double lq  = cst - 0.5 * z.squaredNorm();  // log proposal density
        lr[k] = h - lq;
        if (lr[k] > m) m = lr[k];
      }
      double s = 0.0;
      for (int k = 0; k < total_pts; ++k) s += std::exp(lr[k] - m);
      total += std::log(s) - std::log((double)total_pts) + m;
    }
  }
  return total;
}

// Computes posterior mode and local SD for each group RE. Internal.
// [[Rcpp::export(name = ".brsmm_group_modes_eigen")]]
Eigen::MatrixXd brsmm_group_modes_eigen(
    Eigen::VectorXd param, Eigen::MatrixXd X, Eigen::MatrixXd Z,
    Eigen::MatrixXd Xr, Eigen::VectorXd y_left, Eigen::VectorXd y_right,
    Eigen::VectorXd yt, Eigen::VectorXi delta, Eigen::VectorXi group,
    int link_mu, int link_phi, int repar) {
  int p    = X.cols();
  int q_phi = Z.cols();
  int q_re = Xr.cols();
  int k_re = q_re * (q_re + 1) / 2;
  Eigen::MatrixXd out;
  if (param.size() != (p + q_phi + k_re)) return out;

  Eigen::VectorXd beta     = param.head(p);
  Eigen::VectorXd gamma    = param.segment(p, q_phi);
  Eigen::VectorXd theta_re = param.tail(k_re);
  RandStruct rs = unpack_re(theta_re, q_re);
  Eigen::VectorXd eta_mu  = X * beta;
  Eigen::VectorXd eta_phi = Z * gamma;
  std::vector<GroupData> groups =
      build_groups(eta_mu, eta_phi, Xr, y_left, y_right, yt, delta, group);

  int G = (int)groups.size();
  out = Eigen::MatrixXd::Zero(G, q_re);
  for (int g = 0; g < G; ++g) {
    if (groups[g].delta.size() == 0) continue;
    ModeResult mr = find_mode_vec(groups[g], rs, link_mu, link_phi, repar);
    out.row(g) = mr.mode.transpose();
  }
  return out;
}
