#ifndef IS_CTESTMAT_H
#define IS_CTESTMAT_H

#include<is_cmulti.h>

void cmat_toeplitz(int m, int n, cmulti **A, int LDA, int k, double *a, int offset);
void cmat_cauchy(int m, int n, cmulti **A, int LDA);

#endif
