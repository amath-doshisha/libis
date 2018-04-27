#ifndef IS_DSVD_H
#define IS_DSVD_H

// H=[ A*v-sigma*u;  A'*u-sigma*v ]
void dsvd_residual(int m, int n, double *H, const double *A, int LDA, const double *u, const double *v, double sigma);


#endif
