#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"is_macros.h"
#include"is_ivec.h"
#include"is_dvec.h"
#include"is_dmat.h"
#include"is_dsvd.h"

// H=[ A*v-sigma*u;  A'*u-sigma*v ]
void dsvd_residual(int m, int n, double *H, double *A, int LDA, double *u, double *v, double sigma)
{
  dvec_copy(m,&H[0],u);                // H(1:m)=u
  dvec_copy(n,&H[m],v);                // H(m+1:m+n)=u
  dvec_scale(m+n,H,(-sigma));          // H=-sigma*[u; v]
  dvec_add_lintr(m,n,&H[0],A,LDA,v);   // H(1:m)=A*v-sigma*u
  dvec_add_lintr_t(m,n,&H[m],A,LDA,u); // H(m+1:m+n)=A'*u-sigma*v
}

//EOF
