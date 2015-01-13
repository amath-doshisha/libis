#ifndef IS_ZEIG_H
#define IS_ZEIG_H

void zeig_sort(int n, int k, dcomplex *lambda, dcomplex *X, int LDX);
void zeig_sort_index(int n, int k, dcomplex *lambda, dcomplex *X, int LDX, const int *index);
void zeig_sort_vector_guide(int n, int k, dcomplex *lambda, dcomplex *X, int LDX, dcomplex *X0, int LDX0);
void zeig_sort_value_guide(int n, int k, dcomplex *lambda, dcomplex *X, int LDX, dcomplex *lambda0);

// F=A*X-lambda*X
void zeig_residual(int n, dcomplex *F, const dcomplex *A, int LDA, const dcomplex *X, dcomplex lambda);
// E(k)=max(abs(A*X(:,k)-lambda(k)*X(:,k))), k=1,2,..,n
void zeig_residual_norm_max(int n, double *E, const dcomplex *A, int LDA, const dcomplex *X, int LDX, const dcomplex *lambda);


#endif
