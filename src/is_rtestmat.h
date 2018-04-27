#ifndef IS_RTESTMAT_H
#define IS_RTESTMAT_H

#include<is_rmulti.h>

// Toeplitz matrix
void rmat_toeplitz(int m, int n, rmulti **A, int LDA, int k, double *a, int offset);

// Cauchy matrix
void rmat_cauchy(int m, int n, rmulti **A, int LDA);

#endif
