#ifndef IS_DHSHLDR_H
#define IS_DHSHLDR_H

// size(x)=n
// h[0]=0; ...; h[k-1]=0; h[k]=-s*xi; h[k+1]=x[k+1]; ...; h[n-1]=x[n-1];
void dhouseholder_vec(int n, int k, double *h, double *alpha, const double *x);
void dhouseholder(int n, int k0, int nH, int k, double *h, double *alpha, const double *H, int LDH, const double *Alpha, const double *x);

// B=A*H, H=I-alpha*h*h'
void dhouseholder_right(int m, int n, double *A, int LDA, int k, const double *h, double alpha);

// B=H*A, H=I-alpha*h*h'
void dhouseholder_left(int m, int n, double *A, int LDA, int k, const double *h, double alpha);

#endif
