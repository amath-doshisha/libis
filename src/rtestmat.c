#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"is_macros.h"
#include"is_rmulti.h"
#include"is_rvec.h"
#include"is_rmat.h"
#include"is_rtestmat.h"

void rmat_toeplitz(int m, int n, rmulti **A, int LDA, int k, double *a, int offset)
{
  int i;
  rmat_set_zeros(m,n,A,LDA);
  for(i=0; i<k; i++){
    rmat_set_diag_d(m,n,A,LDA,a[i],i-offset);
  }
}

void rmat_cauchy(int m, int n, rmulti **A, int LDA)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      rset_si(MAT(A,i,j,LDA),i+j+1);
      rinv(MAT(A,i,j,LDA),MAT(A,i,j,LDA));
    }
  }
}

//EOF
