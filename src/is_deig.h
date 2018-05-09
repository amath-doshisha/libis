#ifndef IS_DEIG_H
#define IS_DEIG_H

void deig_sort(int n, int k, double *lambda, double *X, int LDX);
void deig_sort_index(int n, int k, double *lambda, double *X, int LDX, int *index);
void deig_sort_vector_guide(int n, int k, double *lambda, double *X, int LDX, double *X0, int LDX0);
void deig_sort_value_guide(int n, int k, double *lambda, double *X, int LDX, double *lambda0);


// F=A*X-lambda*X
void deig_residual(int n, double *F, double *A, int LDA, double *X, double lambda);
// E(k)=max(abs(A*X(:,k)-lambda(k)*X(:,k))), k=1,2,..,n
void deig_residual_norm_max(int n, double *E, double *A, int LDA, double *X, int LDX, double *lambda);

#endif
