#ifndef IS_DHPEIG_H
#define IS_DHPEIG_H

// ニュートン型反復の最大反復回数
#define DHPEIG_STEP_MAX 35

enum { DHPEIG_NONE=0, DHPEIG_CONVERGENT, DHPEIG_DIVERGENT, DHPEIG_SINGULAR, DHPEIG_NUM };

// hyperplane constrained method
int dhpeig(int n, const double *A, int LDA, double *X, int LDX, double *Lambda, int debug);

// hyperplane constrained method for a pair
int dhpeig_1pair(int n, const double *A, int LDA, const double *z, double *x, double *Lambda, double *E, int *Step, int debug);

// JM=A-lambda*I-x*w'/C
void dhpeig_jacobi_mat(int n, double *JM, int LDJM, const double *A, int LDA, const double *X, const double *W, double lambda, double C);

#endif
