#ifndef IS_RHPEIG_H
#define IS_RHPEIG_H

#include<is_rmulti.h>

// ニュートン型反復の最大反復回数
#define RHPEIG_STEP_MAX 100

// info type
enum { RHPEIG_NONE=0, RHPEIG_CONVERGENT, RHPEIG_DIVERGENT, RHPEIG_SINGULAR, RHPEIG_NUM };

// hyperplane constrained method
int rhpeig(int n, rmulti **X, int LDX, rmulti **Lambda, rmulti **A, int LDA, int debug);

// hyperplane constrained method for a pair
int rhpeig_1pair(int n, rmulti **x, rmulti *lambda, rmulti *e, int *Step, rmulti **A, int LDA, rmulti **z, int debug);

// JM=A-lambda*I-x*w'/C
void rhpeig_jacobi_mat(int n, rmulti **JM, int LDJM, rmulti **A, int LDA, rmulti **x, rmulti **w, rmulti *lambda, rmulti *C);

#endif
