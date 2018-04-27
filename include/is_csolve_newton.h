#ifndef ISYS_CSOLVE_NEWTON_H
#define ISYS_CSOLVE_NEWTON_H

#include<is_cmulti.h>

// F=F(x)
void csolve_func_residual(int m, rmulti *norm, cmulti **F, cmulti **x, func_t *fF);

// satandard
int csolve_newton(int m, cmulti **x, func_t *fF, int step_max, int debug);

// adaptive
//int csolve_newton_adap(int m, cmulti **x, func_t *fF, rmulti *eps, int nmax, int debug);
int csolve_newton_adjust(int m, cmulti **x, func_t *fF, cmulti **xt, rmulti *eps, int nmax, double mu, int l, int kappa, int debug);
void csolve_newton_beta_by_recur(int n, rmulti **beta);
int csolve_newton_beta_from_mu(int n, rmulti **beta, rmulti **mu, int debug);
int csolve_newton_condition_conv(int m, cmulti **x0, cmulti **x1, cmulti **x2, rmulti *beta); // ||x2-x1||/||x1-x0||<=1/(2*beta)
int csolve_newton_condition_done(int m, cmulti **x0, cmulti **x1, rmulti *eps); // 2*||x1-x0||<=eps
int csolve_newton_gamma(rmulti *beta, rmulti *mu);
void csolve_newton_theta(int n, rmulti **theta, rmulti *theta_min, rmulti **mu, rmulti **beta);
void csolve_newton_theta_estimate(rmulti *theta_est, rmulti *theta_min);
int csolve_newton_prec_from_eta(rmulti *eta0, rmulti *eta1, int gamma_exp, int kappa, int tau);
int csolve_newton_prec_from_mu(rmulti *mu0, rmulti *eta0, int gamma_exp, int kappa, int tau);
int csolve_newton_prec_from_beta(rmulti *eta0, rmulti *beta0, rmulti *beta1, rmulti *theta0, int gamma_exp, int kappa, int tau);
void csolve_newton_e_border(rmulti *border, rmulti *eta, rmulti *beta, int gamma_exp); // gamma*eta_{n}/(8*beta_{n})
void csolve_newton_emax(rmulti *border, int b, int kappa, int tau); // 2^(-b+tau+kappa)
void csolve_newton_e_norm(rmulti *e_norm, int m, cmulti **y0, cmulti **x, func_t *fF, func_t *fJ, int prec);
void csolve_newton_map(int m, cmulti **y, cmulti **x, func_t *fF, func_t *fJ);
// krawczyk
int csolve_krawczyk(int m, cmulti **e, cmulti **x, func_t *fF, int debug);

#endif
