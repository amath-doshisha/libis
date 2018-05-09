#ifndef IS_ZEIG_H
#define IS_ZEIG_H

void zeig_sort(int n, int k, dcomplex *lambda, dcomplex *X, int LDX);
void zeig_sort_index(int n, int k, dcomplex *lambda, dcomplex *X, int LDX, int *index);
void zeig_sort_vector_guide(int n, int k, dcomplex *lambda, dcomplex *X, int LDX, dcomplex *X0, int LDX0);
void zeig_sort_value_guide(int n, int k, dcomplex *lambda, dcomplex *X, int LDX, dcomplex *lambda0);

// F=A*X-lambda*X
void zeig_residual(int n, dcomplex *F, dcomplex *A, int LDA, dcomplex *X, dcomplex lambda);
// E(k)=max(abs(A*X(:,k)-lambda(k)*X(:,k))), k=1,2,..,n
void zeig_residual_norm_max(int n, double *E, dcomplex *A, int LDA, dcomplex *X, int LDX, dcomplex *lambda);


#endif
