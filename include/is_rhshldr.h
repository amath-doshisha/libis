#ifndef IS_RHSHLDR_H
#define IS_RHSHLDR_H

#include<is_rmulti.h>

/**
 @file  is_rhshldr.h
*/

// size(x)=n
// h[0]=0; ...; h[k-1]=0; h[k]=-s*xi; h[k+1]=x[k+1]; ...; h[n-1]=x[n-1];
void rhouseholder_vec(int n, int k, rmulti **h, rmulti *alpha, rmulti **x);
// B=A*H, H=I-alpha*h*h'
void rhouseholder_right(int m, int n, rmulti **B, int LDB, rmulti **A, int LDA, int k, rmulti **h, rmulti *alpha); 
// B=H*A, H=I-alpha*h*h'
void rhouseholder_left(int m, int n, rmulti **B, int LDB, rmulti **A, int LDA, int k, rmulti **h, rmulti *alpha);


void rhouseholder(int n, int k0, int nH, int k, rmulti **h, rmulti *alpha, rmulti **H, int LDH, rmulti **Alpha, rmulti **x);

#endif
