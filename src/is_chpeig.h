#ifndef IS_CHPEIG_H
#define IS_CHPEIG_H

#include<is_cmulti.h>

// ニュートン型反復の最大反復回数
#define CHPEIG_STEP_MAX 100

// info type
enum { CHPEIG_NONE=0, CHPEIG_CONVERGENT, CHPEIG_DIVERGENT, CHPEIG_SINGULAR, CHPEIG_RECOMPUTE, CHPEIG_NUM };

// hyperplane constrained method
int chpeig(int n, cmulti **X, int LDX, cmulti **Lambda, cmulti **A, int LDA, int debug);
int chpeig_verify(int n, cmulti **X, int LDX, cmulti **Lambda, cmulti **E_vec, int LDE, cmulti **E_val, cmulti **A, int LDA, int prec_verify, int *prec, int *kprec, int debug);

// hyperplane constrained method for a pair
int chpeig_1pair(int n, cmulti **x, cmulti *lambda, rmulti *E, int *Step, cmulti **A, int LDA, cmulti **z, int debug);

// JM=A-lambda*I-x*w'/C
void chpeig_jacobi_mat(int n, cmulti **JM, int LDJM, cmulti **A, int LDA, cmulti **x, cmulti **w, cmulti *lambda, cmulti *C);

#endif
