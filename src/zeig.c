#include<stdio.h>
#include<stdlib.h>
#include<stdarg.h>
#include"is_macros.h"
#include"is_ivec.h"
#include"is_zvec.h"
#include"is_dvec.h"
#include"is_zmat.h"
#include"is_zeig.h"
#include"mt19937ar.h"


// F=A*X-lambda*X
void zeig_residual(int n, dcomplex *F, dcomplex *A, int LDA, dcomplex *X, dcomplex lambda)
{
  zvec_mul_zvec_zscalar(n,F,X,lambda);     // F=lambda*X
  zvec_neg_zvec(n,F,F);                  // F=-F
  zvec_add_mul_zmat_zvec(n,n,F,A,LDA,X);      // F=A*X-lambda*X
}

// E(k)=max(abs(A*X(:,k)-lambda(k)*X(:,k))), k=1,2,..,n
void zeig_residual_norm_max(int n, double *E, dcomplex *A, int LDA, dcomplex *X, int LDX, dcomplex *lambda)
{
  int j;
  dcomplex *F=NULL;
  F=zvec_allocate(n);
  for(j=0; j<n; j++){
    zeig_residual(n,F,A,LDA,&COL(X,j,LDX),lambda[j]);
    E[j]=dnorm_max_zvec(n,F);
  }
  F=zvec_free(F);
}

/////////////////////////////////////////////////////////////////////////

void zeig_sort(int n, int k, dcomplex *lambda, dcomplex *X, int LDX)
{
  int LDY,*index=NULL;
  dcomplex *Y=NULL;
  // allocate
  index=ivec_allocate(k);
  LDY=n; Y=zmat_allocate(LDY,k);
  // sort
  zvec_sort(k,lambda,index);
  ivec_reverse(k,index);
  zvec_reverse(k,lambda);
  zmat_copy_col_index(n,k,Y,LDY,X,LDX,index);
  zmat_copy(n,k,X,LDX,Y,LDY);
  // free
  index=ivec_free(index);
  Y=zmat_free(Y);
}

/////////////////////////////////////////////////////////////////////////

void zeig_sort_index(int n, int k, dcomplex *lambda, dcomplex *X, int LDX, int *index)
{
  zvec_swap_index(k,lambda,index);
  zmat_swap_index(n,k,X,LDX,index);
}

////////////////////////////////////////////////////////

void zeig_sort_vector_guide(int n, int k, dcomplex *lambda, dcomplex *X, int LDX, dcomplex *X0, int LDX0)
{
  int i,j;
  double value,*a;
  a=dvec_allocate(k);
  for(j=0; j<k; j++){
    zmat_dist_norm_max(a,n,k,X,LDX,&COL(X0,j,LDX0));
    value=dmin_abs_dvec_index(k,a,&i);
    if(i!=j){
      zvec_swap(n,&COL(X,i,LDX),&COL(X,j,LDX));
      zvec_swap_at(lambda,i,j);
    }
  }
  a=dvec_free(a);
}

void zeig_sort_value_guide(int n, int k, dcomplex *lambda, dcomplex *X, int LDX, dcomplex *lambda0)
{
  int i,j;
  double value,*a;
  a=dvec_allocate(k);
  for(j=0; j<k; j++){
    dvec_abs_sub_zvec_zscalar(k,a,lambda,lambda0[j]);
    value=dmin_abs_dvec_index(k,a,&i);
    if(i!=j){
      zvec_swap(n,&COL(X,i,LDX),&COL(X,j,LDX));
      zvec_swap_at(lambda,i,j);
    }
  }
  a=dvec_free(a);
}


//EOF
