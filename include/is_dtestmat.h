#ifndef IS_DTESTMAT_H
#define IS_DTESTMAT_H

void dmat_toeplitz(int m, int n, double *A, int LDA, int k, double *a, int offset);
void dmat_frank(int m, int n, double *A, int LDA, int param);
double dfact(double x);
void dmat_ipjfact(int m, int n, double *A, int LDA, int param);

#endif
