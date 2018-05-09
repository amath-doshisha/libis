#ifndef IS_DHPSVD_H
#define IS_DHPSVD_H

// ニュートン型反復の最大反復回数
#define DHPSVD_STEP_MAX 35
// ニュートン型反復の終了条件
#define DHPSVD_EPS_1ERROR     (pow(2,-53));
#define DHPSVD_EPS_RESIDUAL   (pow(2,-53));
#define DHPSVD_EPS_REJECT     (pow(2,-50));

enum { DHPSVD_NONE=0, DHPSVD_CONVERGENT, DHPSVD_DIVERGENT, DHPSVD_SINGULAR, DHPSVD_NAN };


void dhpsvd_jacobi_mat(int m, int n, double *JM, int LDJM, double *A, int LDA, double *u, double *v, double *w, double sigma, double C);
int dhpsvd_1pair(int m, int n, double *A, int LDA, double *z, double *u, double *v, double *Sigma, double *E, int *Step, int debug);


#endif
