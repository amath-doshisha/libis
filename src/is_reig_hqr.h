#ifndef IS_REIG_HQR_H
#define IS_REIG_HQR_H

#include<is_rmulti.h>
#include<is_cmulti.h>

// QR method for computing eigenvalues of the Hessenberg-type matrix
void reig_hqr(int n, cmulti **lambda, rmulti **A, int LDA, int debug);
void reig_hqr_mt(int m, int n, rmulti **B, int LDB, cmulti **lambda, rmulti **A, int LDA, int debug);

// private functions of reig_hqr
void reig_hqr_main(int n, cmulti **lambda, rmulti **A, int LDA, int debug);
void reig_hqr_calculate_norm(rmulti *anorm, int n, rmulti **A, int LDA);
void reig_hqr_calculate_formula_11_6_23_p(rmulti *p, rmulti *r, rmulti *s, rmulti *w, rmulti *a1, rmulti *a2);
void reig_hqr_calculate_formula_11_6_23_q(rmulti *q, rmulti *a, rmulti *r, rmulti *s, rmulti *z);
void reig_hqr_calculate_formula_11_6_25(rmulti *s, rmulti *p, rmulti *q, rmulti *r);
void reig_hqr_calculate_formula_11_6_26_u(rmulti *u, rmulti *a, rmulti *q, rmulti *r);
void reig_hqr_calculate_formula_11_6_26_v(rmulti *v, rmulti *a1, rmulti *a2, rmulti *p, rmulti *z);
int reig_hqr_test_formula_11_6_26(rmulti *u, rmulti *v);
int reig_hqr_test_UDE(rmulti *s, rmulti *a);
void reig_hqr_scale(rmulti *p, rmulti *q, rmulti *r, rmulti *s);

#endif
