// brs_deriv.h — derivatives of one observation's contribution l(eta_mu, eta_phi)
// in the linear predictors, by central differences (Richardson where noted).
// Shared by loglik.cpp (brs gradient / Hessian) and loglik_mixed.cpp (brsmm).
#pragma once
#include "brs_common.h"
#include <algorithm>
#include <cmath>

// One observation: censoring code and endpoints (see obs_loglik).
struct ObsSpec {
  int delta;
  double left, right, yt;
};

// Model codes shared by all observations.
struct LinkSpec {
  int link_mu, link_phi, repar;
};

// l(eta_mu, eta_phi): one observation's contribution at given predictors.
inline double contrib_eta(double eta_mu, double eta_phi, const ObsSpec &o,
                          const LinkSpec &s) {
  double mu  = clamp_mu_by_repar(inv_link(eta_mu, s.link_mu), s.repar);
  double phi = clamp_phi_by_repar(inv_link(eta_phi, s.link_phi), s.repar);
  double a, b;
  beta_shapes(mu, phi, s.repar, a, b);
  return obs_loglik(o.delta, o.left, o.right, o.yt, a, b);
}

// Derivatives of l: d = d/d eta_mu, p = d/d eta_phi, c = mixed.
struct ObsDeriv {
  double f, d1, d2, d3, p1, p2, c11, c21;
};

// Relative steps (floored at 1) from the Lote 5 numDeriv study: H1 first
// derivatives, H2 second / cross (brs Hessian), H3 third (O(h^2) stencils).
// HC: curvature at the brsmm mode; larger than H2 so its roundoff (~eps/h^2)
// does not make the Laplace value noisy (Richardson truncation stays ~h^4).
static const double FD_H1 = 1.0e-4;
static const double FD_H2 = 3.0e-4;
static const double FD_H3 = 1.0e-3;
static const double FD_HC = 1.0e-3;

inline double fd_step(double eta, double h0) {
  return h0 * std::max(1.0, std::abs(eta));
}

// Richardson first derivative from f(+-h), f(+-2h): O(h^4).
inline double rich_d1(double fm2, double fm1, double fp1, double fp2, double h) {
  double a = (fp1 - fm1) / (2.0 * h), b = (fp2 - fm2) / (4.0 * h);
  return (4.0 * a - b) / 3.0;
}

// Richardson second derivative from f(0), f(+-h), f(+-2h): O(h^4).
inline double rich_d2(double fm2, double fm1, double f0, double fp1, double fp2,
                      double h) {
  double a = (fp1 - 2.0 * f0 + fm1) / (h * h);
  double b = (fp2 - 2.0 * f0 + fm2) / (4.0 * h * h);
  return (4.0 * a - b) / 3.0;
}

// Central third derivative from f(+-h), f(+-2h): O(h^2).
inline double cent_d3(double fm2, double fm1, double fp1, double fp2, double h) {
  return (fp2 - 2.0 * fp1 + 2.0 * fm1 - fm2) / (2.0 * h * h * h);
}

// Richardson mixed second derivative at steps (h, k) and (2h, 2k): O(h^4).
inline double rich_cross(double em, double ep, double h, double k,
                         const ObsSpec &o, const LinkSpec &s) {
  double c1 = (contrib_eta(em + h, ep + k, o, s) - contrib_eta(em + h, ep - k, o, s) -
               contrib_eta(em - h, ep + k, o, s) + contrib_eta(em - h, ep - k, o, s)) /
              (4.0 * h * k);
  double c2 = (contrib_eta(em + 2.0 * h, ep + 2.0 * k, o, s) -
               contrib_eta(em + 2.0 * h, ep - 2.0 * k, o, s) -
               contrib_eta(em - 2.0 * h, ep + 2.0 * k, o, s) +
               contrib_eta(em - 2.0 * h, ep - 2.0 * k, o, s)) / (16.0 * h * k);
  return (4.0 * c1 - c2) / 3.0;
}

// d1 and p1 by Richardson (8 evaluations): brs gradient.
inline void obs_deriv_grad(double em, double ep, const ObsSpec &o,
                           const LinkSpec &s, double &d1, double &p1) {
  const double h = fd_step(em, FD_H1), k = fd_step(ep, FD_H1);
  d1 = rich_d1(contrib_eta(em - 2.0 * h, ep, o, s), contrib_eta(em - h, ep, o, s),
               contrib_eta(em + h, ep, o, s), contrib_eta(em + 2.0 * h, ep, o, s), h);
  p1 = rich_d1(contrib_eta(em, ep - 2.0 * k, o, s), contrib_eta(em, ep - k, o, s),
               contrib_eta(em, ep + k, o, s), contrib_eta(em, ep + 2.0 * k, o, s), k);
}

// d1 and p1 by plain central differences (4 evaluations): quadrature nodes.
inline void obs_deriv_node(double em, double ep, const ObsSpec &o,
                           const LinkSpec &s, double &d1, double &p1) {
  const double h = fd_step(em, FD_H1), k = fd_step(ep, FD_H1);
  d1 = (contrib_eta(em + h, ep, o, s) - contrib_eta(em - h, ep, o, s)) / (2.0 * h);
  p1 = (contrib_eta(em, ep + k, o, s) - contrib_eta(em, ep - k, o, s)) / (2.0 * k);
}

// f, d1, d2 by plain central differences at 2 H1 (3 evaluations): LM steps far
// from the mode, where only the direction matters.
inline void obs_deriv_mu_rough(double em, double ep, const ObsSpec &o,
                               const LinkSpec &s, double &f, double &d1,
                               double &d2) {
  const double h = 2.0 * fd_step(em, FD_H1);
  f = contrib_eta(em, ep, o, s);
  double m1 = contrib_eta(em - h, ep, o, s), p1 = contrib_eta(em + h, ep, o, s);
  d1 = (p1 - m1) / (2.0 * h);
  d2 = (p1 - 2.0 * f + m1) / (h * h);
}

// f, Richardson d1 and a plain d2 on the same stencil (5 evaluations): Newton
// iterations near the mode converge to the zero of the Richardson gradient.
inline void obs_deriv_mu_iter(double em, double ep, const ObsSpec &o,
                              const LinkSpec &s, double &f, double &d1,
                              double &d2) {
  const double h = fd_step(em, FD_H1);
  f = contrib_eta(em, ep, o, s);
  double m2 = contrib_eta(em - 2.0 * h, ep, o, s), m1 = contrib_eta(em - h, ep, o, s);
  double p1 = contrib_eta(em + h, ep, o, s), p2 = contrib_eta(em + 2.0 * h, ep, o, s);
  d1 = rich_d1(m2, m1, p1, p2, h);
  d2 = (p2 - 2.0 * f + m2) / (4.0 * h * h);
}

// f and Richardson d2 in eta_mu (5 evaluations): value and curvature at the mode.
inline void obs_deriv_curv(double em, double ep, const ObsSpec &o,
                           const LinkSpec &s, double &f, double &d2) {
  const double h = fd_step(em, FD_HC);
  f = contrib_eta(em, ep, o, s);
  d2 = rich_d2(contrib_eta(em - 2.0 * h, ep, o, s), contrib_eta(em - h, ep, o, s), f,
               contrib_eta(em + h, ep, o, s), contrib_eta(em + 2.0 * h, ep, o, s), h);
}

// f, d1, d2 in eta_mu by Richardson (9 evaluations): derivatives at the mode.
inline void obs_deriv_mu(double em, double ep, const ObsSpec &o,
                         const LinkSpec &s, ObsDeriv &r) {
  const double h = fd_step(em, FD_H1), h2 = fd_step(em, FD_HC);
  r.f = contrib_eta(em, ep, o, s);
  r.d1 = rich_d1(contrib_eta(em - 2.0 * h, ep, o, s), contrib_eta(em - h, ep, o, s),
                 contrib_eta(em + h, ep, o, s), contrib_eta(em + 2.0 * h, ep, o, s), h);
  r.d2 = rich_d2(contrib_eta(em - 2.0 * h2, ep, o, s), contrib_eta(em - h2, ep, o, s),
                 r.f, contrib_eta(em + h2, ep, o, s),
                 contrib_eta(em + 2.0 * h2, ep, o, s), h2);
  r.d3 = r.p1 = r.p2 = r.c11 = r.c21 = 0.0;
}

// f, d2, p2, c11 by Richardson (17 evaluations): brs Hessian.
inline void obs_deriv_hess(double em, double ep, const ObsSpec &o,
                           const LinkSpec &s, ObsDeriv &r) {
  const double h = fd_step(em, FD_H2), k = fd_step(ep, FD_H2);
  r.f = contrib_eta(em, ep, o, s);
  r.d2 = rich_d2(contrib_eta(em - 2.0 * h, ep, o, s), contrib_eta(em - h, ep, o, s),
                 r.f, contrib_eta(em + h, ep, o, s), contrib_eta(em + 2.0 * h, ep, o, s), h);
  r.p2 = rich_d2(contrib_eta(em, ep - 2.0 * k, o, s), contrib_eta(em, ep - k, o, s),
                 r.f, contrib_eta(em, ep + k, o, s), contrib_eta(em, ep + 2.0 * k, o, s), k);
  r.c11 = rich_cross(em, ep, h, k, o, s);
  r.d1 = r.d3 = r.p1 = r.c21 = 0.0;
}

// Everything the brsmm gradient needs at the mode: f, d1, d2 (as obs_deriv_mu),
// p1, c11 (Richardson), d3 and c21 at step H3 (31 evaluations).
inline void obs_deriv_mode(double em, double ep, const ObsSpec &o,
                           const LinkSpec &s, ObsDeriv &r) {
  obs_deriv_mu(em, ep, o, s, r);
  const double k1 = fd_step(ep, FD_H1);
  r.p1 = rich_d1(contrib_eta(em, ep - 2.0 * k1, o, s), contrib_eta(em, ep - k1, o, s),
                 contrib_eta(em, ep + k1, o, s), contrib_eta(em, ep + 2.0 * k1, o, s), k1);
  r.c11 = rich_cross(em, ep, fd_step(em, FD_H2), fd_step(ep, FD_H2), o, s);
  const double h3 = fd_step(em, FD_H3), k3 = fd_step(ep, FD_H3);
  r.d3 = cent_d3(contrib_eta(em - 2.0 * h3, ep, o, s), contrib_eta(em - h3, ep, o, s),
                 contrib_eta(em + h3, ep, o, s), contrib_eta(em + 2.0 * h3, ep, o, s), h3);
  double fk = contrib_eta(em, ep + k3, o, s), fmk = contrib_eta(em, ep - k3, o, s);
  double fpp = contrib_eta(em + h3, ep + k3, o, s), fpm = contrib_eta(em + h3, ep - k3, o, s);
  double fmp = contrib_eta(em - h3, ep + k3, o, s), fmm = contrib_eta(em - h3, ep - k3, o, s);
  r.c21 = ((fpp - 2.0 * fk + fmp) - (fpm - 2.0 * fmk + fmm)) / (2.0 * k3 * h3 * h3);
  r.p2 = 0.0;
}
