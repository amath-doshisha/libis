#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"is_macros.h"
#include"is_ivec.h"
#include"is_zvec.h"
#include"is_zmat.h"
#include"is_ztestmat.h"

void zmat_toeplitz(int m, int n, dcomplex *A, int LDA, int k, double *a, int offset)
{
  int i;
  zmat_set_zeros(m,n,A,LDA);
  for(i=0; i<k; i++){
    zmat_diag_set_scalar_d(m,n,A,LDA,a[i],i-offset);
  }
}
