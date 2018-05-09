#ifndef IS_DSVD_H
#define IS_DSVD_H

// H=[ A*v-sigma*u;  A'*u-sigma*v ]
void dsvd_residual(int m, int n, double *H, double *A, int LDA, double *u, double *v, double sigma);


#endif
