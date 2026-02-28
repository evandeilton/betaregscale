// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <Rmath.h>
#include <cmath>
#include <vector>

static const double EPS_SHAPE = 1.0e-12;
static const double MAX_SHAPE = 1.0e8;
static const double EPS_PROB = 1.0e-15;
static const double LOG_PENALTY = -1.0e6;
static const double EPS_BOUND = 1.0e-5;

inline double clamp(double x, double lo, double hi) {
  return (x < lo) ? lo : ((x > hi) ? hi : x);
}

inline double inv_link(double eta, int code) {
  switch (code) {
  case 0:
    return 1.0 / (1.0 + std::exp(-eta));
  case 1:
    return R::pnorm(eta, 0.0, 1.0, 1, 0);
  case 2:
    return 0.5 + std::atan(eta) / M_PI;
  case 3:
    return 1.0 - std::exp(-std::exp(eta));
  case 4:
    return std::exp(eta);
  case 5:
    return eta * eta;
  case 6:
    if (std::abs(eta) < EPS_PROB)
      return (eta >= 0.0) ? MAX_SHAPE : -MAX_SHAPE;
    return 1.0 / eta;
  case 7:
    if (eta <= EPS_PROB)
      return MAX_SHAPE;
    return 1.0 / std::sqrt(eta);
  case 8:
    return eta;
  default:
    return 1.0 / (1.0 + std::exp(-eta));
  }
}

inline double clamp_phi_by_repar(double phi, int repar) {
  if (!std::isfinite(phi)) {
    return (repar == 2) ? (1.0 - EPS_BOUND) : MAX_SHAPE;
  }
  if (repar == 2) {
    return clamp(phi, EPS_BOUND, 1.0 - EPS_BOUND);
  }
  return clamp(phi, EPS_BOUND, MAX_SHAPE);
}

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
  a = clamp(a, EPS_SHAPE, MAX_SHAPE);
  b = clamp(b, EPS_SHAPE, MAX_SHAPE);
}

inline double obs_loglik(int delta_i, double left_i, double right_i,
                         double yt_i, double a, double b) {
  double lo = clamp(left_i, EPS_BOUND, 1.0 - EPS_BOUND);
  double hi = clamp(right_i, EPS_BOUND, 1.0 - EPS_BOUND);
  double y = clamp(yt_i, EPS_BOUND, 1.0 - EPS_BOUND);

  double contrib = LOG_PENALTY;
  if (delta_i == 0) {
    contrib = R::dbeta(y, a, b, 1);
  } else if (delta_i == 1) {
    double p = R::pbeta(hi, a, b, 1, 0);
    contrib = std::log(std::max(p, EPS_PROB));
  } else if (delta_i == 2) {
    double p = R::pbeta(lo, a, b, 0, 0);
    contrib = std::log(std::max(p, EPS_PROB));
  } else {
    double p1 = R::pbeta(lo, a, b, 1, 0);
    double p2 = R::pbeta(hi, a, b, 1, 0);
    contrib = std::log(std::max(p2 - p1, EPS_PROB));
  }
  return std::isfinite(contrib) ? contrib : LOG_PENALTY;
}

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

inline double log_prior_b(const Eigen::VectorXd &b, const RandStruct &rs) {
  Eigen::VectorXd z = rs.L.triangularView<Eigen::Lower>().solve(b);
  double quad = z.squaredNorm();
  double cst = rs.q_re * std::log(2.0 * M_PI);
  return -0.5 * (cst + rs.logdet_D + quad);
}

inline double h_func_vec(const Eigen::VectorXd &b, const GroupData &gd,
                         const RandStruct &rs, int link_mu_code,
                         int link_phi_code, int repar) {
  double ll = 0.0;
  int n = gd.delta.size();
  for (int i = 0; i < n; ++i) {
    double eta = gd.eta_mu_fixed(i) + gd.Zr.row(i).dot(b);
    double mu = clamp(inv_link(eta, link_mu_code), EPS_BOUND, 1.0 - EPS_BOUND);
    double phi =
        clamp_phi_by_repar(inv_link(gd.eta_phi(i), link_phi_code), repar);
    double a, bb;
    beta_shapes(mu, phi, repar, a, bb);
    ll += obs_loglik((int)gd.delta(i), gd.y_left(i), gd.y_right(i), gd.yt(i), a,
                     bb);
  }
  ll += log_prior_b(b, rs);
  return std::isfinite(ll) ? ll : LOG_PENALTY;
}

inline void numerical_grad_hess(const Eigen::VectorXd &b, const GroupData &gd,
                                const RandStruct &rs, int link_mu_code,
                                int link_phi_code, int repar,
                                Eigen::VectorXd &grad, Eigen::MatrixXd &hess) {
  int q = b.size();
  grad = Eigen::VectorXd::Zero(q);
  hess = Eigen::MatrixXd::Zero(q, q);
  double f0 = h_func_vec(b, gd, rs, link_mu_code, link_phi_code, repar);

  std::vector<double> step(q, 0.0);
  for (int j = 0; j < q; ++j) {
    step[j] = std::max(1e-5, 1e-4 * std::max(1.0, std::abs(b(j))));
  }

  for (int j = 0; j < q; ++j) {
    Eigen::VectorXd bp = b;
    Eigen::VectorXd bm = b;
    bp(j) += step[j];
    bm(j) -= step[j];
    double fp = h_func_vec(bp, gd, rs, link_mu_code, link_phi_code, repar);
    double fm = h_func_vec(bm, gd, rs, link_mu_code, link_phi_code, repar);
    grad(j) = (fp - fm) / (2.0 * step[j]);
    hess(j, j) = (fp - 2.0 * f0 + fm) / (step[j] * step[j]);
  }

  for (int j = 0; j < q; ++j) {
    for (int k = j + 1; k < q; ++k) {
      Eigen::VectorXd bpp = b, bpm = b, bmp = b, bmm = b;
      bpp(j) += step[j];
      bpp(k) += step[k];
      bpm(j) += step[j];
      bpm(k) -= step[k];
      bmp(j) -= step[j];
      bmp(k) += step[k];
      bmm(j) -= step[j];
      bmm(k) -= step[k];

      double fpp = h_func_vec(bpp, gd, rs, link_mu_code, link_phi_code, repar);
      double fpm = h_func_vec(bpm, gd, rs, link_mu_code, link_phi_code, repar);
      double fmp = h_func_vec(bmp, gd, rs, link_mu_code, link_phi_code, repar);
      double fmm = h_func_vec(bmm, gd, rs, link_mu_code, link_phi_code, repar);
      double hij = (fpp - fpm - fmp + fmm) / (4.0 * step[j] * step[k]);
      hess(j, k) = hij;
      hess(k, j) = hij;
    }
  }
}

struct ModeResult {
  Eigen::VectorXd mode;
  Eigen::MatrixXd curvature;
  double h_at_mode;
};

inline ModeResult find_mode_vec(const GroupData &gd, const RandStruct &rs,
                                int link_mu_code, int link_phi_code,
                                int repar) {
  int q = rs.q_re;
  Eigen::VectorXd b = Eigen::VectorXd::Zero(q);
  double fcur = h_func_vec(b, gd, rs, link_mu_code, link_phi_code, repar);

  for (int iter = 0; iter < 30; ++iter) {
    Eigen::VectorXd g;
    Eigen::MatrixXd H;
    numerical_grad_hess(b, gd, rs, link_mu_code, link_phi_code, repar, g, H);

    Eigen::FullPivLU<Eigen::MatrixXd> lu(H);
    if (!lu.isInvertible()) {
      H -= 1e-6 * Eigen::MatrixXd::Identity(q, q);
      lu.compute(H);
      if (!lu.isInvertible())
        break;
    }
    Eigen::VectorXd step = lu.solve(g);
    if (!step.allFinite())
      break;

    Eigen::VectorXd bnew = b - step;
    double fnew = h_func_vec(bnew, gd, rs, link_mu_code, link_phi_code, repar);
    double damp = 1.0;
    while (fnew < fcur && damp > 1e-3) {
      damp *= 0.5;
      bnew = b - damp * step;
      fnew = h_func_vec(bnew, gd, rs, link_mu_code, link_phi_code, repar);
    }

    if (std::abs(fnew - fcur) < 1e-8 && (damp * step).norm() < 1e-6) {
      b = bnew;
      fcur = fnew;
      break;
    }
    b = bnew;
    fcur = fnew;
  }

  Eigen::VectorXd g;
  Eigen::MatrixXd H;
  numerical_grad_hess(b, gd, rs, link_mu_code, link_phi_code, repar, g, H);
  Eigen::MatrixXd curv = -H;
  curv = 0.5 * (curv + curv.transpose());
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(curv);
  Eigen::VectorXd ev = es.eigenvalues().cwiseMax(1e-8);
  curv = es.eigenvectors() * ev.asDiagonal() * es.eigenvectors().transpose();

  ModeResult out;
  out.mode = b;
  out.curvature = curv;
  out.h_at_mode = fcur;
  return out;
}

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

inline std::vector<GroupData>
build_groups(const Eigen::VectorXd &eta_mu, const Eigen::VectorXd &eta_phi,
             const Eigen::MatrixXd &Xr, const Eigen::VectorXd &y_left,
             const Eigen::VectorXd &y_right, const Eigen::VectorXd &yt,
             const Eigen::VectorXi &delta, const Eigen::VectorXi &group) {
  int G = group.maxCoeff();
  std::vector<int> counts(G, 0);
  for (int i = 0; i < group.size(); ++i) {
    counts[group(i) - 1]++;
  }

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
    groups[g].delta(k) = delta(i);
    groups[g].y_left(k) = y_left(i);
    groups[g].y_right(k) = y_right(i);
    groups[g].yt(k) = yt(i);
    groups[g].eta_mu_fixed(k) = eta_mu(i);
    groups[g].eta_phi(k) = eta_phi(i);
    groups[g].Zr.row(k) = Xr.row(i);
  }
  return groups;
}

const int PRIMES[] = {2,  3,  5,  7,  11, 13, 17, 19, 23, 29,
                      31, 37, 41, 43, 47, 53, 59, 61, 67, 71};

double halton(int index, int base) {
  double f = 1.0;
  double r = 0.0;
  int i = index;
  while (i > 0) {
    f = f / base;
    r = r + f * (i % base);
    i = i / base;
  }
  return r;
}

void build_cartesian_grid(const std::vector<double> &x,
                          const std::vector<double> &w, int q,
                          Eigen::MatrixXd &grid, std::vector<double> &out_w) {
  int n = x.size();
  int total = std::pow(n, q);
  grid.resize(total, q);
  out_w.resize(total);
  for (int i = 0; i < total; ++i) {
    int temp = i;
    double w_prod = 1.0;
    for (int d = 0; d < q; ++d) {
      int idx = temp % n;
      grid(i, d) = x[idx];
      w_prod *= w[idx];
      temp /= n;
    }
    out_w[i] = w_prod;
  }
}

Eigen::MatrixXd build_halton_grid(int n_points, int q) {
  Eigen::MatrixXd grid(n_points, q);
  for (int k = 0; k < n_points; ++k) {
    for (int d = 0; d < q; ++d) {
      int base = PRIMES[d % 20];
      double r = halton(k + 1, base);
      grid(k, d) = R::qnorm(r, 0.0, 1.0, 1, 0);
    }
  }
  return grid;
}

//' @title Mixed Model Log-Likelihood (Eigen)
//' @description Computes marginal log-likelihood using Laplace, AGHQ, or QMC.
//' @param param [beta, gamma, theta_re]
//' @param X Mean design matrix
//' @param Z Precision design matrix
//' @param Xr Random-effects design matrix
//' @param y_left,y_right,yt,delta,group Data
//' @param method 0=Laplace, 1=AGHQ, 2=QMC
//' @param n_points Number of quadrature/QMC points
//' @keywords internal
// [[Rcpp::export]]
double brsmm_loglik_eigen(Eigen::VectorXd param, Eigen::MatrixXd X,
                          Eigen::MatrixXd Z, Eigen::MatrixXd Xr,
                          Eigen::VectorXd y_left, Eigen::VectorXd y_right,
                          Eigen::VectorXd yt, Eigen::VectorXi delta,
                          Eigen::VectorXi group, int link_mu, int link_phi,
                          int repar, int method, int n_points) {
  int p = X.cols();
  int q_phi = Z.cols();
  int q_re = Xr.cols();
  int k_re = q_re * (q_re + 1) / 2;
  if (param.size() != (p + q_phi + k_re)) {
    return LOG_PENALTY;
  }

  Eigen::VectorXd beta = param.head(p);
  Eigen::VectorXd gamma = param.segment(p, q_phi);
  Eigen::VectorXd theta_re = param.tail(k_re);
  RandStruct rs = unpack_re(theta_re, q_re);

  Eigen::VectorXd eta_mu = X * beta;
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
    if (groups[g].delta.size() == 0)
      continue;
    ModeResult mr = find_mode_vec(groups[g], rs, link_mu, link_phi, repar);

    if (method == 0) {
      Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(mr.curvature);
      Eigen::VectorXd ev = es.eigenvalues().cwiseMax(1e-8);
      double logdet = ev.array().log().sum();
      total += mr.h_at_mode + 0.5 * q_re * std::log(2.0 * M_PI) - 0.5 * logdet;
    } else if (method == 1) {
      Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(mr.curvature);
      Eigen::VectorXd ev = es.eigenvalues().cwiseMax(1e-8);
      Eigen::MatrixXd S =
          es.eigenvectors() * ev.cwiseInverse().cwiseSqrt().asDiagonal();

      int total_pts = aghq_grid.rows();
      std::vector<double> lt(total_pts);
      for (int k = 0; k < total_pts; ++k) {
        Eigen::VectorXd z = aghq_grid.row(k).transpose();
        Eigen::VectorXd bk = mr.mode + M_SQRT2 * S * z;
        double h = h_func_vec(bk, groups[g], rs, link_mu, link_phi, repar);
        lt[k] = std::log(aghq_w[k]) + h + z.squaredNorm();
      }
      double m = lt[0];
      for (int k = 1; k < total_pts; ++k)
        if (lt[k] > m)
          m = lt[k];
      double s = 0.0;
      for (int k = 0; k < total_pts; ++k)
        s += std::exp(lt[k] - m);

      double logdet_S = -0.5 * ev.array().log().sum();
      total += std::log(s) + m + q_re * std::log(M_SQRT2) + logdet_S;
    } else {
      Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(mr.curvature);
      Eigen::VectorXd ev = es.eigenvalues().cwiseMax(1e-8);
      Eigen::MatrixXd S =
          es.eigenvectors() * ev.cwiseInverse().cwiseSqrt().asDiagonal();

      int total_pts = qmc_grid.rows();
      std::vector<double> lr(total_pts);
      double m = -1e12;
      double logdet_S = -0.5 * ev.array().log().sum();
      double cst = -0.5 * q_re * std::log(2.0 * M_PI) - logdet_S;

      for (int k = 0; k < total_pts; ++k) {
        Eigen::VectorXd z = qmc_grid.row(k).transpose();
        Eigen::VectorXd bk = mr.mode + S * z;
        double h = h_func_vec(bk, groups[g], rs, link_mu, link_phi, repar);
        double lq = cst - 0.5 * z.squaredNorm();
        lr[k] = h - lq;
        if (lr[k] > m)
          m = lr[k];
      }
      double s = 0.0;
      for (int k = 0; k < total_pts; ++k)
        s += std::exp(lr[k] - m);
      total += std::log(s) - std::log((double)total_pts) + m;
    }
  }
  return total;
}

//' @title Group modes for mixed model (Eigen backend)
//' @keywords internal
// [[Rcpp::export]]
Eigen::MatrixXd brsmm_group_modes_eigen(
    Eigen::VectorXd param, Eigen::MatrixXd X, Eigen::MatrixXd Z,
    Eigen::MatrixXd Xr, Eigen::VectorXd y_left, Eigen::VectorXd y_right,
    Eigen::VectorXd yt, Eigen::VectorXi delta, Eigen::VectorXi group,
    int link_mu, int link_phi, int repar) {
  int p = X.cols();
  int q_phi = Z.cols();
  int q_re = Xr.cols();
  int k_re = q_re * (q_re + 1) / 2;
  Eigen::MatrixXd out;
  if (param.size() != (p + q_phi + k_re)) {
    return out;
  }

  Eigen::VectorXd beta = param.head(p);
  Eigen::VectorXd gamma = param.segment(p, q_phi);
  Eigen::VectorXd theta_re = param.tail(k_re);
  RandStruct rs = unpack_re(theta_re, q_re);
  Eigen::VectorXd eta_mu = X * beta;
  Eigen::VectorXd eta_phi = Z * gamma;
  std::vector<GroupData> groups =
      build_groups(eta_mu, eta_phi, Xr, y_left, y_right, yt, delta, group);

  int G = (int)groups.size();
  out = Eigen::MatrixXd::Zero(G, q_re);
  for (int g = 0; g < G; ++g) {
    if (groups[g].delta.size() == 0)
      continue;
    ModeResult mr = find_mode_vec(groups[g], rs, link_mu, link_phi, repar);
    out.row(g) = mr.mode.transpose();
  }
  return out;
}
