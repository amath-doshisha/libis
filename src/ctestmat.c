#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"is_macros.h"
#include"is_cmulti.h"
#include"is_cvec.h"
#include"is_cmat.h"
#include"is_ctestmat.h"

void cmat_toeplitz(int m, int n, cmulti **A, int LDA, int k, double *a, int offset)
{
  int i;
  cmat_set_zeros(m,n,A,LDA);
  for(i=0; i<k; i++){
    cmat_set_diag_d(m,n,A,LDA,a[i],i-offset);
  }
}

void cmat_cauchy(int m, int n, cmulti **A, int LDA)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      cset_i(MAT(A,i,j,LDA),i+j+1);
      cinv_c(MAT(A,i,j,LDA),MAT(A,i,j,LDA));
    }
  }
}

//EOF
