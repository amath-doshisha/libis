#ifndef IS_CHSHLDR_H
#define IS_CHSHLDR_H

#include<is_rmulti.h>
#include<is_cmulti.h>

/**
 @file  is_chshldr.h
*/

// size(x)=n
// h[0]=0; ...; h[k-1]=0; h[k]=-s*xi; h[k+1]=x[k+1]; ...; h[n-1]=x[n-1];
void chouseholder_vec(int n, int k, cmulti **h, rmulti *alpha, cmulti **x);
// B=A*H, H=I-alpha*h*h'
void chouseholder_right(int m, int n, cmulti **B, int LDB, cmulti **A, int LDA, int k, cmulti **h, rmulti *alpha);
// B=H*A, H=I-alpha*h*h'
void chouseholder_left(int m, int n, cmulti **B, int LDB, cmulti **A, int LDA, int k, cmulti **h, rmulti *alpha);

void chouseholder(int n, int k0, int nH, int k, cmulti **h, rmulti *alpha, cmulti **H, int LDH, rmulti **Alpha, cmulti **x);


#endif
