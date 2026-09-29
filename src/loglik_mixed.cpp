// betaregscale — Armadillo backend for mixed-effects beta interval regression.
// Marginal likelihood by Laplace (method 0), AGHQ (1) or QMC (2); random effects
// b ~ N(0, LL'), theta_re packs L column-wise with log L_ii on the diagonal.
// Inner mode of h(b) = sum_i l_i + log prior: Levenberg–Marquardt Newton on
// chain-rule derivatives in the linear predictor (brs_deriv.h), warm-started
// from the previous call. Gradients differentiate each approximation itself
// (implicit-function theorem for the mode, divided differences for the
// symmetric-root quadrature scaling C^-1/2). A group whose curvature is not positive definite at the
// mode contributes LOG_PENALTY (no eigenvalue floor). The R names of the two
// original exports are kept from the Eigen backend.

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "brs_common.h"
#include "brs_deriv.h"
#include <algorithm>
#include <cfloat>
#include <cstdint>
#include <vector>

// ------------------------------------------------------------- RE structs --

// One group's observations; idx keeps the original row for the X'/Z' products.
struct GroupData {
  std::vector<ObsSpec> obs;
  std::vector<int> idx;
  arma::vec eta_mu_fixed;
  arma::vec eta_phi;
  arma::mat Zt;   // q_re x n_g: column i is the random-effects row of obs i
  int n() const { return (int)obs.size(); }
};

// L, precision P = (LL')^-1, log|D| and, per packed parameter r, dP/dtheta_r
// and dlog|D|/dtheta_r (for the gradient).
struct RandStruct {
  arma::mat L, P;
  double logdet_D;
  int q_re;
  std::vector<arma::mat> Pr;
  std::vector<double> dlogdet;
};

// Unpack theta_re; a zero or infinite L_ii makes P non-finite, so every h()
// becomes LOG_PENALTY (as in the Eigen backend).
inline RandStruct unpack_re(const arma::vec &theta_re, int q) {
  RandStruct rs;
  rs.q_re = q;
  rs.L.zeros(q, q);
  int k = 0;
  double half = 0.0;
  for (int j = 0; j < q; ++j)
    for (int i = j; i < q; ++i) {
      double v = theta_re(k++);
      if (i == j) { rs.L(i, i) = std::exp(v); half += v; }
      else rs.L(i, j) = v;
    }
  rs.logdet_D = 2.0 * half;
  // L^-1 column by column (forward substitution), then P = L^-T L^-1
  arma::mat Li(q, q, arma::fill::zeros);
  for (int c = 0; c < q; ++c)
    for (int i = 0; i < q; ++i) {
      double s = (i == c) ? 1.0 : 0.0;
      for (int j = 0; j < i; ++j) s -= rs.L(i, j) * Li(j, c);
      Li(i, c) = s / rs.L(i, i);
    }
  rs.P = Li.t() * Li;
  rs.P = 0.5 * (rs.P + rs.P.t());
  // dP/dtheta_r = -P (E L' + L E') P with E = dL/dtheta_r (a single entry)
  for (int j = 0; j < q; ++j)
    for (int i = j; i < q; ++i) {
      arma::mat E(q, q, arma::fill::zeros);
      E(i, j) = (i == j) ? rs.L(i, i) : 1.0;
      arma::mat Pr = -rs.P * (E * rs.L.t() + rs.L * E.t()) * rs.P;
      rs.Pr.push_back(0.5 * (Pr + Pr.t()));
      rs.dlogdet.push_back(i == j ? 2.0 : 0.0);
    }
  return rs;
}

// log N(b; 0, LL') = -0.5 (q log 2pi + log|D| + b'Pb).
inline double log_prior_b(const arma::vec &b, const RandStruct &rs) {
  double quad = 0.0;
  for (int i = 0; i < rs.q_re; ++i)
    for (int j = 0; j < rs.q_re; ++j) quad += b(i) * rs.P(i, j) * b(j);
  return -0.5 * (rs.q_re * std::log(2.0 * M_PI) + rs.logdet_D + quad);
}

// z_i'b for observation i of a group.
inline double zb(const GroupData &gd, int i, const arma::vec &b) {
  const double *zi = gd.Zt.colptr(i);
  double d = 0.0;
  for (arma::uword j = 0; j < b.n_elem; ++j) d += zi[j] * b(j);
  return d;
}

// h(b) = sum_i l_i(eta_i(b)) + log prior; LOG_PENALTY if not finite.
inline double h_func_vec(const arma::vec &b, const GroupData &gd,
                         const RandStruct &rs, const LinkSpec &s) {
  double ll = 0.0;
  for (int i = 0; i < gd.n(); ++i)
    ll += contrib_eta(gd.eta_mu_fixed(i) + zb(gd, i, b), gd.eta_phi(i), gd.obs[i], s);
  ll += log_prior_b(b, rs);
  return std::isfinite(ll) ? ll : LOG_PENALTY;
}

// Accumulate g += d1 z and the lower triangle of H += d2 z z' for one observation.
inline void add_obs(const double *zi, int q, double d1, double d2, arma::vec &g,
                    arma::mat &H) {
  for (int j = 0; j < q; ++j) {
    g(j) += d1 * zi[j];
    for (int l = 0; l <= j; ++l) H(j, l) += d2 * zi[j] * zi[l];
  }
}

// Mirror the lower triangle, add the prior and return h (LOG_PENALTY if not finite).
inline double finish_h(const arma::vec &b, const RandStruct &rs, double f,
                       arma::vec &g, arma::mat &H) {
  const int q = rs.q_re;
  for (int j = 0; j < q; ++j)
    for (int l = 0; l < j; ++l) H(l, j) = H(j, l);
  g -= rs.P * b;
  H -= rs.P;
  f += log_prior_b(b, rs);
  return std::isfinite(f) ? f : LOG_PENALTY;
}

// Newton iterate: h, g = sum d1_i z_i - Pb, H = sum d2_i z_i z_i' - P, exactly
// symmetric; rich = 5-point stencil (Richardson d1) near the mode, else 3-point.
inline double h_grad_hess(const arma::vec &b, const GroupData &gd,
                          const RandStruct &rs, const LinkSpec &s, arma::vec &g,
                          arma::mat &H, bool rich) {
  const int q = rs.q_re;
  g.zeros(q);
  H.zeros(q, q);
  double f = 0.0;
  for (int i = 0; i < gd.n(); ++i) {
    const double eta = gd.eta_mu_fixed(i) + zb(gd, i, b);
    double fi, d1, d2;
    if (rich) obs_deriv_mu_iter(eta, gd.eta_phi(i), gd.obs[i], s, fi, d1, d2);
    else obs_deriv_mu_rough(eta, gd.eta_phi(i), gd.obs[i], s, fi, d1, d2);
    f += fi;
    add_obs(gd.Zt.colptr(i), q, d1, d2, g, H);
  }
  return finish_h(b, rs, f, g, H);
}

// At the mode: h and H with the Richardson curvature (5 evaluations/obs).
inline double h_curv(const arma::vec &b, const GroupData &gd, const RandStruct &rs,
                     const LinkSpec &s, arma::mat &H) {
  const int q = rs.q_re;
  arma::vec g(q, arma::fill::zeros);
  H.zeros(q, q);
  double f = 0.0;
  for (int i = 0; i < gd.n(); ++i) {
    double fi, d2;
    obs_deriv_curv(gd.eta_mu_fixed(i) + zb(gd, i, b), gd.eta_phi(i), gd.obs[i], s, fi, d2);
    f += fi;
    add_obs(gd.Zt.colptr(i), q, 0.0, d2, g, H);
  }
  return finish_h(b, rs, f, g, H);
}

// Richardson gradient of h at b (diagnostics only).
inline arma::vec h_grad_rich(const arma::vec &b, const GroupData &gd,
                             const RandStruct &rs, const LinkSpec &s) {
  arma::vec g = -rs.P * b;
  ObsDeriv r;
  for (int i = 0; i < gd.n(); ++i) {
    obs_deriv_mu(gd.eta_mu_fixed(i) + zb(gd, i, b), gd.eta_phi(i), gd.obs[i], s, r);
    for (int j = 0; j < rs.q_re; ++j) g(j) += r.d1 * gd.Zt(j, i);
  }
  return g;
}

// ------------------------------------------------------ dense linear algebra --

// Solve A x = g by Cholesky (small q); false if A is not numerically PD.
inline bool chol_solve(const arma::mat &A, const arma::vec &g, arma::vec &x) {
  const int q = (int)A.n_rows;
  arma::mat L(q, q, arma::fill::zeros);
  for (int j = 0; j < q; ++j) {
    double s = A(j, j);
    for (int k = 0; k < j; ++k) s -= L(j, k) * L(j, k);
    if (!(s > 0.0) || !std::isfinite(s)) return false;
    L(j, j) = std::sqrt(s);
    for (int i = j + 1; i < q; ++i) {
      double t = A(i, j);
      for (int k = 0; k < j; ++k) t -= L(i, k) * L(j, k);
      L(i, j) = t / L(j, j);
    }
  }
  x = g;
  for (int i = 0; i < q; ++i) {          // L y = g
    for (int k = 0; k < i; ++k) x(i) -= L(i, k) * x(k);
    x(i) /= L(i, i);
  }
  for (int i = q - 1; i >= 0; --i) {     // L' x = y
    for (int k = i + 1; k < q; ++k) x(i) -= L(k, i) * x(k);
    x(i) /= L(i, i);
  }
  return x.is_finite();
}

// Symmetric eigen-decomposition with canonical eigenvector signs (largest
// |component| positive); results use only sign-free quantities (C^-1, C^-1/2).
inline bool sym_eig(const arma::mat &C, arma::vec &ev, arma::mat &V) {
  if (!C.is_finite()) return false;
  if (!arma::eig_sym(ev, V, C)) return false;
  for (arma::uword j = 0; j < V.n_cols; ++j) {
    arma::uword im = arma::index_max(arma::abs(V.col(j)));
    if (V(im, j) < 0.0) V.col(j) = -V.col(j);
  }
  return true;
}

// --------------------------------------------------- mode-finding (LM Newton) --

struct ModeResult {
  arma::vec mode;
  arma::mat C;         // curvature -H at the mode
  arma::vec ev;        // its eigenvalues (all > 0 when ok)
  arma::mat V;         // its eigenvectors (canonical signs)
  double h_at_mode;
  double grad_inf;     // |grad h|_inf at the mode (Richardson; diagnostics only)
  int iter;
  bool ok;             // false: curvature not positive definite -> LOG_PENALTY
};

// Mode of h by Levenberg–Marquardt: s solves (lambda I - H) s = g with lambda = 0
// while -H is PD, else raised until h improves. Warm start b0 is used only if
// it beats b = 0. Rough derivatives far from the mode, Richardson ones once the
// Newton step is < 1e-3; a Richardson Newton step < 1e-6 (1 + |b|) is taken from
// the quadratic model and ends the search (error ~ |s|^2). Flat h ends it at once.
inline ModeResult find_mode_vec(const GroupData &gd, const RandStruct &rs,
                                const LinkSpec &s, const arma::vec &b0,
                                bool diag = false) {
  const int q = rs.q_re;
  arma::vec b(q, arma::fill::zeros), g, step, bnew, gnew;
  arma::mat H, Hnew;
  bool rich = false;   // derivatives at b come from the 5-point stencil
  double f;
  if ((int)b0.n_elem == q && b0.is_finite() && arma::any(b0 != 0.0)) {
    f = h_grad_hess(b0, gd, rs, s, g, H, rich);
    b = b0;
    if (h_func_vec(arma::vec(q, arma::fill::zeros), gd, rs, s) > f) {
      b.zeros();
      f = h_grad_hess(b, gd, rs, s, g, H, rich);
    }
  } else {
    f = h_grad_hess(b, gd, rs, s, g, H, rich);
  }
  ModeResult out;
  out.iter = 0;
  double lambda = 0.0;
  for (int iter = 0; iter < 100; ++iter) {
    out.iter = iter + 1;
    const double scale = std::max(1.0, arma::abs(H.diag()).max());
    const double bn = arma::abs(b).max();
    double lam = lambda, fnew = f;
    bool accepted = false, converged = false, refine = false, near = false;
    for (int t = 0; t < 40; ++t) {
      arma::mat M = -H;
      M.diag() += lam;
      if (!chol_solve(M, g, step)) { lam = std::max(4.0 * lam, 1e-6 * scale); continue; }
      const double sn = arma::abs(step).max();
      if (sn < 1e-12 * (1.0 + bn) || (lam == 0.0 && sn < 1e-6 * (1.0 + bn))) {
        if (!rich) { refine = true; break; }   // redo with the Richardson gradient
        if (sn >= 1e-12 * (1.0 + bn)) b += step;
        converged = true;
        break;
      }
      near = (lam == 0.0 && sn < 1e-3 * (1.0 + bn));
      bnew = b + step;
      if (!bnew.is_finite()) { lam = std::max(4.0 * lam, 1e-6 * scale); continue; }
      fnew = h_grad_hess(bnew, gd, rs, s, gnew, Hnew, near || rich);
      if (fnew > f) { accepted = true; break; }
      lam = std::max(4.0 * lam, 1e-6 * scale);
    }
    if (refine) {
      rich = true;
      f = h_grad_hess(b, gd, rs, s, g, H, true);
      continue;
    }
    if (converged || !accepted) break;
    b = bnew; f = fnew; g = gnew; H = Hnew;
    rich = rich || near;
    lambda = (lam < 1e-4 * scale) ? 0.0 : lam / 16.0;
  }
  // value and Richardson curvature at the mode
  out.h_at_mode = h_curv(b, gd, rs, s, H);
  out.mode = b;
  out.grad_inf = diag ? arma::abs(h_grad_rich(b, gd, rs, s)).max() : NA_REAL;
  out.C = -H;
  // valid: finite prior (theta_re usable) and a positive definite curvature
  out.ok = std::isfinite(log_prior_b(b, rs)) && sym_eig(out.C, out.ev, out.V) &&
           out.ev.min() > 0.0 && out.ev.min() > 1e-14 * out.ev.max();
  return out;
}

// Warm-start cache: the previous call's modes, keyed by (G, q_re); reset per fit.
struct ModeCache {
  int G = 0, q = 0;
  arma::mat modes;
  bool valid = false;
};
static ModeCache mode_cache;

// Clear the warm-start cache. brsmm() calls it first, so a fit depends only on
// its data and start, not on earlier fits in the session.
// [[Rcpp::export(name = ".brsmm_reset_cache", rng = false)]]
void brsmm_reset_cache() {
  mode_cache = ModeCache();
}

inline arma::vec cached_start(int g, int G, int q) {
  if (mode_cache.valid && mode_cache.G == G && mode_cache.q == q)
    return mode_cache.modes.row(g).t();
  return arma::vec(q, arma::fill::zeros);
}
inline void cache_begin(int G, int q) {
  if (!(mode_cache.valid && mode_cache.G == G && mode_cache.q == q)) {
    mode_cache.G = G; mode_cache.q = q;
    mode_cache.modes.zeros(G, q);
    mode_cache.valid = true;
  }
}
inline void cache_store(int g, const arma::vec &mode) {
  if (mode.is_finite()) mode_cache.modes.row(g) = mode.t();
  else mode_cache.valid = false;   // reset on non-finite
}

// ------------------------------------------- Gauss-Hermite quadrature rule --

// Golub–Welsch: nodes = eigenvalues of the Jacobi matrix, weights sqrt(pi) V(0,i)^2.
void compute_gh_rule(int n, std::vector<double> &x, std::vector<double> &w) {
  arma::mat J(n, n, arma::fill::zeros);
  for (int i = 0; i < n - 1; ++i) {
    double val = std::sqrt((double)(i + 1) / 2.0);
    J(i, i + 1) = val;
    J(i + 1, i) = val;
  }
  arma::vec ev;
  arma::mat V;
  if (!arma::eig_sym(ev, V, J))
    Rcpp::stop("AGHQ: Gauss-Hermite rule could not be computed (n_points=%d).", n);
  x.resize(n);
  w.resize(n);
  double sqrt_pi = std::sqrt(M_PI);
  for (int i = 0; i < n; ++i) {
    x[i] = ev(i);
    w[i] = V(0, i) * V(0, i) * sqrt_pi;
  }
}

// -------------------------------------------------- group data partitioning --

// Split the observations by group (1-based codes, validated by check_mixed_inputs).
inline std::vector<GroupData>
build_groups(const arma::vec &eta_mu, const arma::vec &eta_phi,
             const arma::mat &Xr, const arma::vec &y_left,
             const arma::vec &y_right, const arma::vec &yt,
             const Rcpp::IntegerVector &delta,
             const Rcpp::IntegerVector &group) {
  const int n = (int)group.size();
  const int q = (int)Xr.n_cols;
  int G = 0;
  for (int i = 0; i < n; ++i) G = std::max(G, (int)group[i]);
  std::vector<int> counts(G, 0);
  for (int i = 0; i < n; ++i) counts[group[i] - 1]++;
  std::vector<GroupData> groups(G);
  for (int g = 0; g < G; ++g) {
    const int ng = counts[g];
    groups[g].obs.resize(ng);
    groups[g].idx.resize(ng);
    groups[g].eta_mu_fixed.set_size(ng);
    groups[g].eta_phi.set_size(ng);
    groups[g].Zt.set_size(q, ng);
  }
  std::vector<int> cur(G, 0);
  for (int i = 0; i < n; ++i) {
    const int g = group[i] - 1, k = cur[g]++;
    groups[g].obs[k] = ObsSpec{delta[i], y_left(i), y_right(i), yt(i)};
    groups[g].idx[k] = i;
    groups[g].eta_mu_fixed(k) = eta_mu(i);
    groups[g].eta_phi(k) = eta_phi(i);
    for (int j = 0; j < q; ++j) groups[g].Zt(j, k) = Xr(i, j);
  }
  return groups;
}

// -------------------------------------------------------- Halton sequences --
// The first 50 primes, used as bases for the Halton low-discrepancy sequence.
static const int PRIMES[50] = {
    2,   3,   5,   7,  11,  13,  17,  19,  23,  29,
   31,  37,  41,  43,  47,  53,  59,  61,  67,  71,
   73,  79,  83,  89,  97, 101, 103, 107, 109, 113,
  127, 131, 137, 139, 149, 151, 157, 163, 167, 173,
  179, 181, 191, 193, 197, 199, 211, 223, 227, 229
};
static const int N_PRIMES = 50;

// Van der Corput radical inverse of index in the given base (Halton coordinate).
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

// Tensor grid of the 1-D rule over q dims; int64 guard + hard limit on the size.
void build_cartesian_grid(const std::vector<double> &x,
                          const std::vector<double> &w, int q,
                          arma::mat &grid, std::vector<double> &out_w) {
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
  grid.set_size(total, q);
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

// Halton points in q dims mapped to N(0,1) by qnorm; clamped away from 0/1.
arma::mat build_halton_grid(int n_points, int q) {
  if (q > N_PRIMES)
    Rcpp::stop("QMC supports at most %d random-effect dimensions (q_re=%d).",
               N_PRIMES, q);
  arma::mat grid(n_points, q);
  for (int k = 0; k < n_points; ++k)
    for (int d = 0; d < q; ++d) {
      double r = halton(k + 1, PRIMES[d]);
      r = std::min(std::max(r, 1e-9), 1.0 - 1e-9);
      grid(k, d) = R::qnorm(r, 0.0, 1.0, 1, 0);
    }
  return grid;
}

// Structural checks once per call (no bounds checks under -DNDEBUG): sizes,
// NA / non-finite data, group and delta codes.
inline void check_mixed_inputs(const arma::mat &X, const arma::mat &Z,
                               const arma::mat &Xr, const arma::vec &y_left,
                               const arma::vec &y_right, const arma::vec &yt,
                               const Rcpp::IntegerVector &delta,
                               const Rcpp::IntegerVector &group) {
  const arma::uword n = X.n_rows;
  if (n < 1) Rcpp::stop("brsmm: no observations.");
  if (Z.n_rows != n || Xr.n_rows != n || y_left.n_elem != n ||
      y_right.n_elem != n || yt.n_elem != n ||
      (arma::uword)delta.size() != n || (arma::uword)group.size() != n) {
    Rcpp::stop("brsmm: X, Z, Xr, y_left, y_right, yt, delta and group must "
               "all have %d rows.", (int)n);
  }
  if (Xr.n_cols < 1) Rcpp::stop("brsmm: Xr must have at least one column.");
  if (!X.is_finite() || !Z.is_finite() || !Xr.is_finite() || !y_left.is_finite() ||
      !y_right.is_finite() || !yt.is_finite())
    Rcpp::stop("brsmm: X, Z, Xr, y_left, y_right and yt must not contain NA or "
               "non-finite values.");
  for (arma::uword i = 0; i < n; ++i) {
    if (group[i] == NA_INTEGER)
      Rcpp::stop("brsmm: group indices must be >= 1 (found NA at row %d).", (int)i + 1);
    if (group[i] < 1)
      Rcpp::stop("brsmm: group indices must be >= 1 (found %d at row %d).",
                 group[i], (int)i + 1);
    if ((arma::uword)group[i] > n)   // codes are 1..G with G <= n
      Rcpp::stop("brsmm: group index %d at row %d exceeds the number of rows (%d).",
                 group[i], (int)i + 1, (int)n);
    if (delta[i] == NA_INTEGER)
      Rcpp::stop("brsmm: delta must be in {0,1,2,3} (found NA at row %d).", (int)i + 1);
    if (delta[i] < 0 || delta[i] > 3)
      Rcpp::stop("brsmm: delta must be in {0,1,2,3} (found %d at row %d).",
                 delta[i], (int)i + 1);
  }
}

// ------------------------------------------------------- shared set-up --

// Everything a call needs after validation: parameter blocks, groups, grids.
struct MixedSetup {
  int p, q_phi, q_re, k_re, G;
  RandStruct rs;
  LinkSpec s;
  std::vector<GroupData> groups;
  std::vector<double> gh_x, gh_w, aghq_w;
  arma::mat aghq_grid, qmc_grid;
};

inline MixedSetup mixed_setup(const arma::vec &param, const arma::mat &X,
                              const arma::mat &Z, const arma::mat &Xr,
                              const arma::vec &y_left, const arma::vec &y_right,
                              const arma::vec &yt,
                              const Rcpp::IntegerVector &delta,
                              const Rcpp::IntegerVector &group, int link_mu,
                              int link_phi, int repar, int method, int n_points) {
  check_mixed_inputs(X, Z, Xr, y_left, y_right, yt, delta, group);
  if (method != 0 && n_points < 1) Rcpp::stop("brsmm: n_points must be >= 1.");
  MixedSetup m;
  m.p = (int)X.n_cols; m.q_phi = (int)Z.n_cols; m.q_re = (int)Xr.n_cols;
  m.k_re = m.q_re * (m.q_re + 1) / 2;
  m.s = LinkSpec{link_mu, link_phi, repar};
  if ((int)param.n_elem != m.p + m.q_phi + m.k_re)
    Rcpp::stop("brsmm: param must have length %d (got %d).", m.p + m.q_phi + m.k_re,
               (int)param.n_elem);
  arma::vec eta_mu = X * param.head(m.p);
  arma::vec eta_phi = (m.q_phi > 0) ? arma::vec(Z * param.subvec(m.p, m.p + m.q_phi - 1))
                                    : arma::vec(X.n_rows, arma::fill::zeros);
  m.rs = unpack_re(param.tail(m.k_re), m.q_re);
  m.groups = build_groups(eta_mu, eta_phi, Xr, y_left, y_right, yt, delta, group);
  m.G = (int)m.groups.size();
  if (method == 1) {
    compute_gh_rule(n_points, m.gh_x, m.gh_w);
    build_cartesian_grid(m.gh_x, m.gh_w, m.q_re, m.aghq_grid, m.aghq_w);
  } else if (method == 2) {
    m.qmc_grid = build_halton_grid(n_points, m.q_re);
  }
  cache_begin(m.G, m.q_re);
  return m;
}

// Quadrature scaling: symmetric root S = C^-1/2 = V diag(ev^-1/2) V' (unique and
// smooth in C, so independent of eigenvector signs/order); logdet_S = -0.5 sum log ev.
inline void curvature_scaling(const ModeResult &mr, int q, arma::mat &S,
                              double &logdet_S) {
  S.zeros(q, q);
  double sl = 0.0;
  for (int k = 0; k < q; ++k) {
    const double sk = std::sqrt(1.0 / mr.ev(k));
    for (int j = 0; j < q; ++j)
      for (int i = 0; i < q; ++i) S(i, j) += mr.V(i, k) * sk * mr.V(j, k);
    sl += std::log(mr.ev(k));
  }
  logdet_S = -0.5 * sl;
}

// Node b_k = mode + alpha S z_k (alpha = sqrt 2 for AGHQ, 1 for QMC) and z_k'z_k.
inline void node_point(const ModeResult &mr, const arma::mat &S,
                       const arma::mat &grid, int k, int q, double alpha,
                       arma::vec &bk, double &zz) {
  zz = 0.0;
  for (int i = 0; i < q; ++i) {
    double d = 0.0;
    for (int j = 0; j < q; ++j) d += S(i, j) * grid(k, j);
    bk(i) = mr.mode(i) + alpha * d;
    zz += grid(k, i) * grid(k, i);
  }
}

// One group's AGHQ / QMC log-integral; lt receives the log-sum-exp terms.
inline double group_quadrature(const MixedSetup &m, int g, const ModeResult &mr,
                               int method, std::vector<double> &lt) {
  const int q = m.q_re;
  arma::mat S;
  double logdet_S;
  curvature_scaling(mr, q, S, logdet_S);
  arma::vec bk(q);
  double zz;
  const bool aghq = (method == 1);
  const arma::mat &grid = aghq ? m.aghq_grid : m.qmc_grid;
  const int K = (int)grid.n_rows;
  lt.resize(K);
  const double cst = -0.5 * q * std::log(2.0 * M_PI) - logdet_S;   // QMC proposal
  for (int k = 0; k < K; ++k) {
    node_point(mr, S, grid, k, q, aghq ? M_SQRT2 : 1.0, bk, zz);
    const double h = h_func_vec(bk, m.groups[g], m.rs, m.s);
    lt[k] = aghq ? (m.aghq_w[k] > 0.0 ? std::log(m.aghq_w[k]) : -1e15) + h + zz
                 : h - (cst - 0.5 * zz);
  }
  const double mx = *std::max_element(lt.begin(), lt.end());
  double sum = 0.0;
  for (int k = 0; k < K; ++k) sum += std::exp(lt[k] - mx);
  return aghq ? std::log(sum) + mx + q * std::log(M_SQRT2) + logdet_S
              : std::log(sum) - std::log((double)K) + mx;
}

// ====================================================== Main exported fns === //

// Marginal log-likelihood by Laplace / AGHQ / QMC; param = [beta, gamma, theta_re].
// Structural errors stop; a group without a valid mode adds LOG_PENALTY.
// [[Rcpp::export(name = ".brsmm_loglik_eigen", rng = false)]]
double brsmm_loglik_eigen(const arma::vec &param, const arma::mat &X,
                          const arma::mat &Z, const arma::mat &Xr,
                          const arma::vec &y_left, const arma::vec &y_right,
                          const arma::vec &yt,
                          const Rcpp::IntegerVector &delta,
                          const Rcpp::IntegerVector &group, int link_mu,
                          int link_phi, int repar, int method, int n_points) {
  MixedSetup m = mixed_setup(param, X, Z, Xr, y_left, y_right, yt, delta, group,
                             link_mu, link_phi, repar, method, n_points);
  const int q = m.q_re;
  std::vector<double> lt;
  double total = 0.0;
  for (int g = 0; g < m.G; ++g) {
    if (m.groups[g].n() == 0) continue;
    ModeResult mr = find_mode_vec(m.groups[g], m.rs, m.s, cached_start(g, m.G, q));
    cache_store(g, mr.mode);
    if (!mr.ok) { total += LOG_PENALTY; continue; }
    if (method == 0) {
      // Laplace: h(mode) + (q/2) log 2pi - 0.5 log|C|
      double logdet = 0.0;
      for (int j = 0; j < q; ++j) logdet += std::log(mr.ev(j));
      total += mr.h_at_mode + 0.5 * q * std::log(2.0 * M_PI) - 0.5 * logdet;
    } else {
      total += group_quadrature(m, g, mr, method, lt);
    }
  }
  return total;
}

// Posterior mode of each group's random effects (G x q_re).
// [[Rcpp::export(name = ".brsmm_group_modes_eigen", rng = false)]]
arma::mat brsmm_group_modes_eigen(const arma::vec &param, const arma::mat &X,
                                  const arma::mat &Z, const arma::mat &Xr,
                                  const arma::vec &y_left,
                                  const arma::vec &y_right,
                                  const arma::vec &yt,
                                  const Rcpp::IntegerVector &delta,
                                  const Rcpp::IntegerVector &group,
                                  int link_mu, int link_phi, int repar) {
  MixedSetup m = mixed_setup(param, X, Z, Xr, y_left, y_right, yt, delta, group,
                             link_mu, link_phi, repar, 0, 1);
  arma::mat out(m.G, m.q_re, arma::fill::zeros);
  for (int g = 0; g < m.G; ++g) {
    if (m.groups[g].n() == 0) continue;
    ModeResult mr = find_mode_vec(m.groups[g], m.rs, m.s, cached_start(g, m.G, m.q_re));
    cache_store(g, mr.mode);
    out.row(g) = mr.mode.t();
  }
  return out;
}

// Inner-mode diagnostics per group (validation): mode, |grad h|_inf, smallest
// curvature eigenvalue, h at the mode, LM iterations, validity. warm = false
// clears the cache first, so every group starts from b = 0.
// [[Rcpp::export(name = ".brsmm_mode_diag_cpp", rng = false)]]
Rcpp::List brsmm_mode_diag_cpp(const arma::vec &param, const arma::mat &X,
                               const arma::mat &Z, const arma::mat &Xr,
                               const arma::vec &y_left, const arma::vec &y_right,
                               const arma::vec &yt,
                               const Rcpp::IntegerVector &delta,
                               const Rcpp::IntegerVector &group, int link_mu,
                               int link_phi, int repar, bool warm) {
  if (!warm) mode_cache.valid = false;
  MixedSetup m = mixed_setup(param, X, Z, Xr, y_left, y_right, yt, delta, group,
                             link_mu, link_phi, repar, 0, 1);
  arma::mat modes(m.G, m.q_re, arma::fill::zeros);
  arma::vec gi(m.G), emin(m.G), hm(m.G), it(m.G), ok(m.G);
  for (int g = 0; g < m.G; ++g) {
    ModeResult mr = find_mode_vec(m.groups[g], m.rs, m.s, cached_start(g, m.G, m.q_re), true);
    cache_store(g, mr.mode);
    modes.row(g) = mr.mode.t();
    gi(g) = mr.grad_inf;
    emin(g) = mr.ev.is_empty() ? NA_REAL : mr.ev.min();
    hm(g) = mr.h_at_mode;
    it(g) = mr.iter;
    ok(g) = mr.ok ? 1.0 : 0.0;
  }
  return Rcpp::List::create(Rcpp::Named("mode") = modes, Rcpp::Named("grad_inf") = gi,
                            Rcpp::Named("min_eig") = emin, Rcpp::Named("h") = hm,
                            Rcpp::Named("iter") = it, Rcpp::Named("ok") = ok);
}

// ---------------------------------------------------------------- gradient --

// dS of S = C^-1/2 under a symmetric dC (Daleckii–Krein divided differences of
// f(x) = x^-1/2): V'dS V = -A_ij s_i^2 s_j^2 / (s_i + s_j), A = V'dC V, s = ev^-1/2;
// no division by eigenvalue gaps, so it is stable for close eigenvalues.
inline arma::mat dS_of(const arma::mat &V, const arma::vec &ev, const arma::mat &dC) {
  const int q = (int)ev.n_elem;
  const arma::mat A = V.t() * dC * V;
  arma::vec sv(q);
  for (int k = 0; k < q; ++k) sv(k) = 1.0 / std::sqrt(ev(k));
  arma::mat M(q, q);
  for (int j = 0; j < q; ++j)
    for (int i = 0; i < q; ++i)
      M(i, j) = -A(i, j) * sv(i) * sv(i) * sv(j) * sv(j) / (sv(i) + sv(j));
  return V * M * V.t();
}

// Gradient in [beta, gamma, theta_re]: per perturbation t, dm = C^-1 u and
// dC = Cd - sum_j dm_j T_j (u = d grad_b h/dt, Cd = dC/dt at fixed b). Laplace:
// dh/dt - tr(C^-1 dC)/2; AGHQ/QMC: node-weighted dh/dt + G'dm + alpha tr(dS W) - tr(C^-1 dC)/2.
// [[Rcpp::export(name = ".brsmm_grad_cpp", rng = false)]]
arma::vec brsmm_grad_cpp(const arma::vec &param, const arma::mat &X,
                         const arma::mat &Z, const arma::mat &Xr,
                         const arma::vec &y_left, const arma::vec &y_right,
                         const arma::vec &yt, const Rcpp::IntegerVector &delta,
                         const Rcpp::IntegerVector &group, int link_mu,
                         int link_phi, int repar, int method, int n_points) {
  MixedSetup m = mixed_setup(param, X, Z, Xr, y_left, y_right, yt, delta, group,
                             link_mu, link_phi, repar, method, n_points);
  arma::vec grad(param.n_elem, arma::fill::zeros);
  const int q = m.q_re, n = (int)X.n_rows, kr = m.k_re;
  const double alpha = (method == 1) ? M_SQRT2 : 1.0;
  arma::vec dmu(n, arma::fill::zeros), dphi(n, arma::fill::zeros);
  arma::vec gth(kr, arma::fill::zeros);
  std::vector<double> lt;
  for (int g = 0; g < m.G; ++g) {
    const GroupData &gd = m.groups[g];
    const int ng = gd.n();
    if (ng == 0) continue;
    ModeResult mr = find_mode_vec(gd, m.rs, m.s, cached_start(g, m.G, q));
    cache_store(g, mr.mode);
    if (!mr.ok) continue;   // LOG_PENALTY group: constant, zero gradient
    // derivatives of each contribution at the mode, C^-1 and the T_j
    std::vector<ObsDeriv> od(ng);
    for (int i = 0; i < ng; ++i)
      obs_deriv_mode(gd.eta_mu_fixed(i) + zb(gd, i, mr.mode), gd.eta_phi(i),
                     gd.obs[i], m.s, od[i]);
    const arma::mat Cinv = mr.V * arma::diagmat(1.0 / mr.ev) * mr.V.t();
    std::vector<arma::mat> T(q, arma::mat(q, q, arma::fill::zeros));
    for (int i = 0; i < ng; ++i) {
      const arma::vec zi = gd.Zt.col(i);
      const arma::mat zz = zi * zi.t();
      for (int j = 0; j < q; ++j) T[j] += (od[i].d3 * zi(j)) * zz;
    }
    // direct terms dh/dt (at the mode, or node-weighted) and node statistics
    arma::vec dir_mu(ng), dir_phi(ng), dir_th(kr), G(q, arma::fill::zeros);
    arma::mat W(q, q, arma::fill::zeros);
    if (method == 0) {
      for (int i = 0; i < ng; ++i) { dir_mu(i) = od[i].d1; dir_phi(i) = od[i].p1; }
      for (int r = 0; r < kr; ++r)
        dir_th(r) = -0.5 * arma::as_scalar(mr.mode.t() * m.rs.Pr[r] * mr.mode) -
                    0.5 * m.rs.dlogdet[r];
    } else {
      group_quadrature(m, g, mr, method, lt);
      const int K = (int)lt.size();
      const double mx = *std::max_element(lt.begin(), lt.end());
      double sum = 0.0;
      for (int k = 0; k < K; ++k) sum += std::exp(lt[k] - mx);
      arma::mat S;
      double logdet_S;
      curvature_scaling(mr, q, S, logdet_S);
      const arma::mat &grid = (method == 1) ? m.aghq_grid : m.qmc_grid;
      dir_mu.zeros(); dir_phi.zeros(); dir_th.zeros();
      arma::vec bk(q), gk(q);
      double zz;
      for (int k = 0; k < K; ++k) {
        const double pk = std::exp(lt[k] - mx) / sum;
        node_point(mr, S, grid, k, q, alpha, bk, zz);
        gk = -m.rs.P * bk;
        for (int i = 0; i < ng; ++i) {
          double d1, p1;
          obs_deriv_node(gd.eta_mu_fixed(i) + zb(gd, i, bk), gd.eta_phi(i),
                         gd.obs[i], m.s, d1, p1);
          dir_mu(i) += pk * d1;
          dir_phi(i) += pk * p1;
          gk += d1 * gd.Zt.col(i);
        }
        for (int r = 0; r < kr; ++r)
          dir_th(r) += pk * (-0.5 * arma::as_scalar(bk.t() * m.rs.Pr[r] * bk) -
                             0.5 * m.rs.dlogdet[r]);
        G += pk * gk;
        W += pk * (grid.row(k).t() * gk.t());
      }
    }
    // total derivative for one perturbation (u, Cd), see the function comment
    auto total = [&](const arma::vec &u, const arma::mat &Cd) {
      const arma::vec dm = Cinv * u;
      arma::mat dC = Cd;
      for (int j = 0; j < q; ++j) dC -= dm(j) * T[j];
      double v = -0.5 * arma::accu(Cinv % dC);
      if (method != 0)
        v += arma::dot(G, dm) + alpha * arma::accu(dS_of(mr.V, mr.ev, dC) % W.t());
      return v;
    };
    for (int i = 0; i < ng; ++i) {
      const arma::vec zi = gd.Zt.col(i);
      const arma::mat zz = zi * zi.t();
      dmu(gd.idx[i]) = dir_mu(i) + total(od[i].d2 * zi, -od[i].d3 * zz);
      dphi(gd.idx[i]) = dir_phi(i) + total(od[i].c11 * zi, -od[i].c21 * zz);
    }
    for (int r = 0; r < kr; ++r)
      gth(r) += dir_th(r) + total(-m.rs.Pr[r] * mr.mode, m.rs.Pr[r]);
  }
  if (m.p > 0) grad.head(m.p) = X.t() * dmu;
  if (m.q_phi > 0) grad.subvec(m.p, m.p + m.q_phi - 1) = Z.t() * dphi;
  grad.tail(kr) = gth;
  return grad;
}

// Hessian of the marginal log-likelihood: Richardson central differences of
// the gradient (steps 1e-3 max(1,|theta_j|) and half), symmetrised. Every
// perturbed gradient starts from the modes at param, so the result does not
// depend on evaluation order.
// [[Rcpp::export(name = ".brsmm_hessian_cpp", rng = false)]]
arma::mat brsmm_hessian_cpp(const arma::vec &param, const arma::mat &X,
                            const arma::mat &Z, const arma::mat &Xr,
                            const arma::vec &y_left, const arma::vec &y_right,
                            const arma::vec &yt, const Rcpp::IntegerVector &delta,
                            const Rcpp::IntegerVector &group, int link_mu,
                            int link_phi, int repar, int method, int n_points) {
  const int npar = (int)param.n_elem;
  auto grad_at = [&](const arma::vec &p) {
    return brsmm_grad_cpp(p, X, Z, Xr, y_left, y_right, yt, delta, group,
                          link_mu, link_phi, repar, method, n_points);
  };
  grad_at(param);
  const ModeCache saved = mode_cache;
  arma::mat H(npar, npar);
  for (int j = 0; j < npar; ++j) {
    const double h = 1e-3 * std::max(1.0, std::abs(param(j)));
    const double off[4] = {h, -h, 0.5 * h, -0.5 * h};
    arma::vec gv[4];
    for (int t = 0; t < 4; ++t) {
      mode_cache = saved;
      arma::vec pp = param;
      pp(j) += off[t];
      gv[t] = grad_at(pp);
    }
    const arma::vec a = (gv[0] - gv[1]) / (2.0 * h), b = (gv[2] - gv[3]) / h;
    H.col(j) = (4.0 * b - a) / 3.0;
  }
  mode_cache = saved;
  return 0.5 * (H + H.t());
}
