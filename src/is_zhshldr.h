#ifndef IS_ZHSHLDR_H
#define IS_ZHSHLDR_H

// size(x)=n
// h[0]=0; ...; h[k-1]=0; h[k]=-s*xi; h[k+1]=x[k+1]; ...; h[n-1]=x[n-1];
void zhouseholder_vec(int n, int k, dcomplex *h, double *alpha, dcomplex *x);
void zhouseholder(int n, int k0, int nH, int k, dcomplex *h, double *alpha, dcomplex *H, int LDH, double *Alpha, dcomplex *x);

// B=A*H, H=I-alpha*h*h'
void zhouseholder_right(int m, int n, dcomplex *A, int LDA, int k, dcomplex *h, double alpha);

// B=H*A, H=I-alpha*h*h'
void zhouseholder_left(int m, int n, dcomplex *A, int LDA, int k, dcomplex *h, double alpha);

#endif
