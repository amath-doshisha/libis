#ifndef IS_ZHPEIG_H
#define IS_ZHPEIG_H

// ニュートン型反復の最大反復回数
#define ZHPEIG_STEP_MAX 35

enum { ZHPEIG_NONE=0, ZHPEIG_CONVERGENT, ZHPEIG_DIVERGENT, ZHPEIG_SINGULAR, ZHPEIG_NUM };

// hyperplane constrained method
int zhpeig(int n, dcomplex *A, int LDA, dcomplex *X, int LDX, dcomplex *Lambda, int debug);

// hyperplane constrained method for a pair
int zhpeig_1pair(int n, dcomplex *A, int LDA, dcomplex *z, dcomplex *x, dcomplex *Lambda, double *E, int *Step, int debug);

// JM=A-lambda*I-x*w'/C
void zhpeig_jacobi_mat(int n, dcomplex *JM, int LDJM, dcomplex *A, int LDA, dcomplex *X, dcomplex *W, dcomplex lambda, dcomplex C);

#endif
