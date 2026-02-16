
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <Rmath.h>
#include <cmath>
#include <vector>

// Numerical constants
static const double EPS_SHAPE = 1.0e-12;
static const double MAX_SHAPE = 1.0e8;
static const double EPS_PROB = 1.0e-15;
static const double LOG_PENALTY = -1.0e6;
static const double EPS_BOUND = 1.0e-5;

// Helper: Clamp
inline double clamp(double x, double lo, double hi) {
  return (x < lo) ? lo : ((x > hi) ? hi : x);
}

// Helper: Inverse links (scalar)
inline double inv_link(double eta, int code) {
  switch (code) {
  case 0: return 1.0 / (1.0 + std::exp(-eta)); // logit
  case 1: return R::pnorm(eta, 0.0, 1.0, 1, 0); // probit
  case 2: return 0.5 + std::atan(eta) / M_PI; // cauchit
  case 3: return 1.0 - std::exp(-std::exp(eta)); // cloglog
  case 4: return std::exp(eta); // log
  case 5: return eta * eta; // sqrt
  case 6: return 1.0 / eta; // inverse
  case 7: return 1.0 / std::sqrt(eta); // 1/mu^2
  case 8: return eta; // identity
  default: return 1.0 / (1.0 + std::exp(-eta));
  }
}

// Helper: Beta shapes
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

// Helper: Obs loglik
inline double obs_loglik(int delta_i, double left_i, double right_i,
                         double yt_i, double a, double b) {
  double contrib;
  double lo = clamp(left_i, EPS_BOUND, 1.0 - EPS_BOUND);
  double hi = clamp(right_i, EPS_BOUND, 1.0 - EPS_BOUND);
  double y  = clamp(yt_i, EPS_BOUND, 1.0 - EPS_BOUND);

  if (delta_i == 0) { // Exact
    contrib = R::dbeta(y, a, b, 1);
  } else if (delta_i == 1) { // Left-censored
    double p = R::pbeta(hi, a, b, 1, 0);
    contrib = std::log(std::max(p, EPS_PROB));
  } else if (delta_i == 2) { // Right-censored
    double p = R::pbeta(lo, a, b, 0, 0); // upper tail
    contrib = std::log(std::max(p, EPS_PROB));
  } else { // Interval-censored
    double p1 = R::pbeta(lo, a, b, 1, 0);
    double p2 = R::pbeta(hi, a, b, 1, 0);
    double area = std::max(p2 - p1, EPS_PROB);
    contrib = std::log(area);
  }
  return std::isfinite(contrib) ? contrib : LOG_PENALTY;
}

// Data structure for a group
struct GroupData {
  Eigen::VectorXd delta;
  Eigen::VectorXd y_left;
  Eigen::VectorXd y_right;
  Eigen::VectorXd yt;
  Eigen::VectorXd eta_mu_fixed; // X_i * beta
  Eigen::VectorXd eta_phi;      // Z_i * gamma
  // For q=1 random intercept, Z_random is just a vector of 1s (or implicit)
  // Generalizing slightly:
  // Eigen::MatrixXd Z_random;
};

// Function h(b) = sum log f(y|b) + log phi(b)
// For q=1 b is scalar b_val.
double h_func(double b_val, const GroupData &gd, double sigma_b,
              int link_mu_code, int link_phi_code, int repar) {
  double ll = 0.0;
  int n = gd.delta.size();
  for (int i = 0; i < n; i++) {
    double eta = gd.eta_mu_fixed(i) + b_val; // random intercept assumption for now
    double mu = inv_link(eta, link_mu_code);
    double phi = inv_link(gd.eta_phi(i), link_phi_code);
    
    mu = clamp(mu, EPS_BOUND, 1.0 - EPS_BOUND);
    phi = clamp(phi, EPS_BOUND, 1.0 - EPS_BOUND);
    
    double a, bb;
    beta_shapes(mu, phi, repar, a, bb);
    ll += obs_loglik((int)gd.delta(i), gd.y_left(i), gd.y_right(i), gd.yt(i), a, bb);
  }
  // Prior: log N(0, sigma_b^2)
  ll += R::dnorm(b_val, 0.0, sigma_b, 1);
  return ll;
}

// Numerical derivative for h'(b)
double h_grad(double b_val, const GroupData &gd, double sigma_b,
              int link_mu_code, int link_phi_code, int repar) {
  double eps = 1e-5;
  double f1 = h_func(b_val + eps, gd, sigma_b, link_mu_code, link_phi_code, repar);
  double f2 = h_func(b_val - eps, gd, sigma_b, link_mu_code, link_phi_code, repar);
  return (f1 - f2) / (2.0 * eps);
}

// Numerical 2nd derivative for h''(b)
double h_hess(double b_val, const GroupData &gd, double sigma_b,
              int link_mu_code, int link_phi_code, int repar) {
  double eps = 1e-4;
  double f0 = h_func(b_val, gd, sigma_b, link_mu_code, link_phi_code, repar);
  double f1 = h_func(b_val + eps, gd, sigma_b, link_mu_code, link_phi_code, repar);
  double f2 = h_func(b_val - eps, gd, sigma_b, link_mu_code, link_phi_code, repar);
  return (f1 - 2.0 * f0 + f2) / (eps * eps);
}

// Newton-Raphson to find mode
// Returns pair {mode, curv}
std::pair<double, double> find_mode(const GroupData &gd, double sigma_b,
                                    int link_mu_code, int link_phi_code, int repar) {
  double b = 0.0; // Start at 0
  for (int iter = 0; iter < 20; iter++) {
    double g = h_grad(b, gd, sigma_b, link_mu_code, link_phi_code, repar);
    double H = h_hess(b, gd, sigma_b, link_mu_code, link_phi_code, repar);
    if (std::abs(H) < 1e-8) break; // unsafe
    double step = g / H;
    // Simple damping or check
    if (step > 1.0) step = 1.0;
    if (step < -1.0) step = -1.0;
    b -= step;
    if (std::abs(step) < 1e-6) break;
  }
  double curv = -h_hess(b, gd, sigma_b, link_mu_code, link_phi_code, repar);
  if (curv <= 1e-8) curv = 1.0 / (sigma_b * sigma_b); // fallback to prior precision
  return {b, curv};
}

// GH Weights/Nodes generator (Golub-Welsch)
// n_points: number of points
// Returns x (nodes) and w (weights)
void compute_gh_rule(int n, std::vector<double> &x, std::vector<double> &w) {
  // For small n, hardcoded might be faster, but let's do robust eigen method
  // Construct symmetric tridiagonal matrix
  // J_i,i = 0
  // J_i,i+1 = sqrt((i+1)/2)
  Eigen::MatrixXd J = Eigen::MatrixXd::Zero(n, n);
  for (int i = 0; i < n - 1; i++) {
    double val = std::sqrt((double)(i + 1) / 2.0);
    J(i, i + 1) = val;
    J(i + 1, i) = val;
  }
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(J);
  Eigen::VectorXd nodes = es.eigenvalues();
  Eigen::MatrixXd vecs = es.eigenvectors();
  
  x.resize(n);
  w.resize(n);
  double sqrt_pi = std::sqrt(M_PI);
  for (int i = 0; i < n; i++) {
    x[i] = nodes(i);
    // w_i = v_0,i^2 * sqrt(pi) ??? 
    // Standard Golub-Welsch: w_i = (v_1,i)^2 * integral_measure
    // For Hermite, integral e^-x^2 dx is sqrt(pi).
    // The eigenvector first component squared * total mass.
    // However, typical implementations scale differently.
    // Let's rely on standard source: w_j = v_{1,j}^2 * sqrt(pi) * 2^(n-1)? No.
    // Actually w_i = v_{0,i}^2 * sqrt(pi).
    // Let's verify standard GH weights sum to sqrt(pi).
    w[i] = std::pow(vecs(0, i), 2) * sqrt_pi; 
  }
  // Weights usually for integral e^-x^2 f(x). 
  // We need for N(0,1) density. N(x) = 1/sqrt(2pi) e^-x^2/2.
  // Standard GH integrates e^-t^2. Change of variable x = t*sqrt(2).
}

// Precomputed rules cache could be added. 
// For now, compute on fly or static local for fixed N.

// -------------------------------------------------------------------------- //
// Main Exported Function
// -------------------------------------------------------------------------- //

//' @title Mixed Model Log-Likelihood (Eigen)
//' @description Computes marginal log-likelihood using Laplace, AGHQ, or QMC.
//' @param param [beta, gamma, log_sigma]
//' @param X, Z Design matrices
//' @param y_left, y_right, yt, delta Data
//' @param group Group indices
//' @param method 0=Laplace, 1=AGHQ, 2=QMC
//' @param n_points Number of quadrature/QMC points
//' @keywords internal
// [[Rcpp::export]]
double brsmm_loglik_eigen(Eigen::VectorXd param, 
                          Eigen::MatrixXd X, Eigen::MatrixXd Z,
                          Eigen::VectorXd y_left, Eigen::VectorXd y_right, 
                          Eigen::VectorXd yt, Eigen::VectorXi delta,
                          Eigen::VectorXi group,
                          int link_mu, int link_phi, int repar,
                          int method, int n_points) {
                          
  int p = X.cols();
  int q_phi = Z.cols();
  
  Eigen::VectorXd beta = param.head(p);
  Eigen::VectorXd gamma = param.segment(p, q_phi);
  double log_sigma = param(p + q_phi);
  double sigma_b = std::exp(log_sigma);
  
  Eigen::VectorXd eta_mu = X * beta;
  Eigen::VectorXd eta_phi = Z * gamma;
  
  // Identify groups
  int max_g = group.maxCoeff();
  std::vector<GroupData> groups(max_g);
  
  // Naive fill (could be optimized)
  for (int i = 0; i < group.size(); i++) {
    int g_idx = group(i) - 1; // 0-based
    // Ideally push back, but Eigen resize is costly.
    // Pre-count size logic omitted for brevity, using simplistic append or known structure.
    // Re-impl: Just iterate and track indices.
  }
  
  // Better: count group sizes first
  std::vector<int> counts(max_g, 0);
  for(int i=0; i<group.size(); ++i) counts[group(i)-1]++;
  
  // Allocate
  for(int g=0; g<max_g; ++g) {
    groups[g].delta.resize(counts[g]);
    groups[g].y_left.resize(counts[g]);
    groups[g].y_right.resize(counts[g]);
    groups[g].yt.resize(counts[g]);
    groups[g].eta_mu_fixed.resize(counts[g]);
    groups[g].eta_phi.resize(counts[g]);
  }
  
  std::vector<int> cursors(max_g, 0);
  for (int i = 0; i < group.size(); i++) {
    int g = group(i) - 1;
    int k = cursors[g]++;
    groups[g].delta(k) = delta(i);
    groups[g].y_left(k) = y_left(i);
    groups[g].y_right(k) = y_right(i);
    groups[g].yt(k) = yt(i);
    groups[g].eta_mu_fixed(k) = eta_mu(i);
    groups[g].eta_phi(k) = eta_phi(i);
  }
  
  double total_ll = 0.0;
  
  // Prepared Quadrature Rule for AGHQ
  std::vector<double> gh_x, gh_w;
  if(method == 1) { // AGHQ
    compute_gh_rule(n_points, gh_x, gh_w);
  }

  // Prepared QMC Rule
  Eigen::VectorXd qmc_pts;
  if(method == 2) { // QMC (Halton base 2)
    qmc_pts.resize(n_points);
    for(int k=0; k<n_points; ++k) {
       // Simple Halton base 2 generator
       double f = 1.0, r = 0.0;
       int i = k + 1; // offset to avoid 0
       while (i > 0) {
           f /= 2.0;
           r += f * (i % 2);
           i /= 2;
       }
       // Transform Uniform(0,1) to Normal(0,1) via Inverse CDF
       qmc_pts(k) = R::qnorm(r, 0.0, 1.0, 1, 0);
    }
  }

  // Loop Groups
  for (int g = 0; g < max_g; g++) {
    if(groups[g].delta.size() == 0) continue;
    
    // 1. Find Mode b_hat
    std::pair<double, double> opt = find_mode(groups[g], sigma_b, link_mu, link_phi, repar);
    double b_hat = opt.first;
    double curv = opt.second; // approx -H
    
    if (method == 0) { // Laplace
      // log L approx h(b_hat) - 0.5 * log(det(H)) + const
      // Actually standard formula: 
      // L approx exp(h(b_hat)) * sqrt(2pi) * H^-0.5
      // log L = h(b_hat) + 0.5*log(2pi) - 0.5*log(H)
      // Note h(b) included log prior density which has -0.5*log(2pi*sigma^2).
      // Let's check:
      // h(b) = log f(y|b) + log phi(b).
      // Integral approx exp(h(b_hat)) * sqrt(2pi/H).
      // log I approx h(b_hat) + 0.5*log(2pi) - 0.5*log(H).
      
      double val = h_func(b_hat, groups[g], sigma_b, link_mu, link_phi, repar);
      total_ll += val + 0.5 * std::log(2.0 * M_PI) - 0.5 * std::log(curv);
      
    } else if (method == 1) { // AGHQ
      // Change of variable: b = b_hat + z * s, where s = 1/sqrt(curv)
      // I = int e^h(b) db
      // approx sum w_k * e^h(b(x_k)) * Jacobian?
      // Standard GH approx: int e^-x^2 f(x) dx approx sum w_k f(x_k).
      // We want int e^(h(b) - (-b^2/2sigma^2? No, h(b) is full).
      // Let's use the transformation z ~ N(0,1).
      // b = b_hat + z / sqrt(curv).
      // We just sum over standard normal nodes.
      // E[g(z)] approx sum w_k g(x_k) / sqrt(pi)?
      // Wait, standard GH integrates e^-x^2.
      
      // Let's use the "adaptive" formulation:
      // L approx sqrt(2) / sqrt(H) * sum w_k exp( h(b_hat + sqrt(2)*x_k/sqrt(H)) + x_k^2 ) ?
      // This removes the gaussian kernel from h and puts it back?
      
      // Simpler:
      // Treat posterior as approx Normal(b_hat, 1/curv).
      // Use nodes scaled by this normal.
      // nodes_Transformed = b_hat + nodes_GH_Standard * sqrt(2) * sigma_approx ??
      
      // Let's stick to simple AGHQ logic from papers:
      // z_k are nodes of H_k(x). x_k = z_k. w_k weights.
      // These integrate exp(-x^2).
      // We assume posteror approx P(b) ~ exp( - H/2 (b-b_hat)^2 ).
      // Let's rewrite h(b) = h(b_hat) - H/2 (b-b_hat)^2 + R(b).
      // Then integral is exp(h(b_hat)) * int exp(-u^2) ...
      
      // Practical implementation:
      // nodes t_k, weights w_k from GH rule for weight e^{-t^2}.
      // b_k = b_hat + t_k * sqrt(2) / sqrt(curv).
      // log_term = h(b_k) + t_k^2. (adding back t^2 because we divide by e^-t^2 implicit in rule)
      // Integral approx: sqrt(2)/sqrt(curv) * sum w_k * exp(log_term).
      
      double sd_proxy = 1.0 / std::sqrt(curv);
      double sum_exp = 0.0;
      // Use efficient log-sum-exp if needed but direct sum is likely fine for few nodes
      
      // Need log-sum-exp for stability:
      std::vector<double> log_terms(n_points);
      
      for(int k=0; k<n_points; ++k) {
         double t_k = gh_x[k];
         // Transform
         double b_k = b_hat + t_k * M_SQRT2 * sd_proxy;
         double h_val = h_func(b_k, groups[g], sigma_b, link_mu, link_phi, repar);
         // The weight w_k is for exp(-t^2).
         // The integral is int f(t) e^-t^2 dt approx sum w_k f(t_k).
         // Our integral is int exp(h(b)) db.
         // Change var: b = b_hat + t * sqrt(2) * sd_proxy.
         // db = sqrt(2) * sd_proxy dt.
         // Int = sqrt(2)*sd * int exp(h(b(t))) dt.
         //     = sqrt(2)*sd * int exp(h(b(t)) + t^2) * e^-t^2 dt.
         // Approx = sqrt(2)*sd * sum w_k * exp(h(b(t_k)) + t_k^2).
         
         log_terms[k] = std::log(gh_w[k]) + h_val + t_k*t_k;
      }
      
      // LogSumExp
      double max_val = log_terms[0];
      for(int k=1; k<n_points; ++k) if(log_terms[k]>max_val) max_val = log_terms[k];
      double sum_s = 0;
      for(int k=0; k<n_points; ++k) sum_s += std::exp(log_terms[k] - max_val);
      
      double log_I = std::log(sum_s) + max_val + std::log(M_SQRT2 * sd_proxy);
      total_ll += log_I;
      
    } else if (method == 2) { // QMC (Importance Sampling via Shifted Normal)
       // Proposal q(b) = N(b_hat, 1/curv).
       // I = int e^h(b) db = int (e^h / q(b)) q(b) db
       //   approx 1/M sum (e^h(b_k) / q(b_k)) where b_k ~ q(b).
       // b_k generated via QMC inverse transform.
       // b_k = b_hat + sd_proxy * Z_k (Z_k from QMC Normal).
       
       double sd_proxy = 1.0 / std::sqrt(curv);
       
       std::vector<double> log_ratios(n_points);
       double max_lr = -1e9;
       
       for(int k=0; k<n_points; ++k) {
          double z = qmc_pts(k);
          double b_k = b_hat + z * sd_proxy;
          
          double h_val = h_func(b_k, groups[g], sigma_b, link_mu, link_phi, repar);
          double log_q = R::dnorm(b_k, b_hat, sd_proxy, 1);
          
          double lr = h_val - log_q;
          log_ratios[k] = lr;
          if(lr > max_lr) max_lr = lr;
       }
       
       double sum_r = 0.0;
       for(int k=0; k<n_points; ++k) sum_r += std::exp(log_ratios[k] - max_lr);
       
       double log_I = std::log(sum_r) - std::log((double)n_points) + max_lr;
       total_ll += log_I;
    }
  }
  
  return total_ll;
}
