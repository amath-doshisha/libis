#ifndef ISYS_IREIG_H
#define ISYS_IREIG_H

#include<is_rmulti.h>

int ireig_krawczyk(int n, int k, rmulti **E_vec, int LDE, rmulti **E_val, rmulti **A, int LDA, rmulti **X, int LDX, rmulti **lambda, int debug);
int ireig_1pair_krawczyk(int n, rmulti **e, rmulti **A, int LDA, rmulti **x, rmulti *lambda, int debug);

// F=A*X-lambda*X
void ireig_residual(int n, rmulti **F0, rmulti **F1, rmulti **A0, int LDA0, rmulti **A1, int LDA1, rmulti **X0, rmulti **X1, rmulti *lambda0, rmulti *lambda1);

#endif




