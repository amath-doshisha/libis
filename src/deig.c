#include<stdio.h>
#include<stdlib.h>
#include<stdarg.h>
#include"is_macros.h"
#include"is_ivec.h"
#include"is_dvec.h"
#include"is_dmat.h"
#include"is_deig.h"
#include"mt19937ar.h"


// F=A*X-lambda*X
void deig_residual(int n, double *F, double *A, int LDA, double *X, double lambda)
{
  dvec_mul_dvec_dscalar(n,F,X,-lambda); //  F=-lambda*X
  dvec_add_lintr(n,n,F,A,LDA,X);   // F=A*X-lambda*X
}

// E(k)=max(abs(A*X(:,k)-lambda(k)*X(:,k))), k=1,2,..,n
void deig_residual_norm_max(int n, double *E, double *A, int LDA, double *X, int LDX, double *lambda)
{
  int j;
  double *F=NULL;
  F=dvec_allocate(n);
  for(j=0; j<n; j++){
    deig_residual(n,F,A,LDA,&COL(X,j,LDX),lambda[j]);
    E[j]=dnorm_max_dvec(n,F);
  }
  F=dvec_free(F);
}

////////////////////////////////////////////////////////

void deig_sort(int n, int k, double *lambda, double *X, int LDX)
{
  int LDY,*index=NULL;
  double *Y=NULL;
  // allocate
  index=ivec_allocate(k);
  LDY=n; Y=dmat_allocate(LDY,k);
  // sort
  dvec_sort(k,lambda,index);
  ivec_reverse(k,index);
  dvec_reverse(k,lambda);
  dmat_copy_col_index(n,k,Y,LDY,X,LDX,index);
  dmat_copy(n,k,X,LDX,Y,LDY);
  // free
  index=ivec_free(index);
  Y=dmat_free(Y);
}

/////////////////////////////////////////////////////////////////////////

void deig_sort_index(int n, int k, double *lambda, double *X, int LDX, int *index)
{
  dvec_swap_index(k,lambda,index);
  dmat_swap_index(n,k,X,LDX,index);
}

////////////////////////////////////////////////////////

void deig_sort_vector_guide(int n, int k, double *lambda, double *X, int LDX, double *X0, int LDX0)
{
  int i,j;
  double value,*a;
  a=dvec_allocate(k);
  for(j=0; j<k; j++){
    dmat_dist_norm_max(a,n,k,X,LDX,&COL(X0,j,LDX0));
    value=dmin_abs_dvec_index(k,a,&i);
    if(i!=j){
      dvec_swap(n,&COL(X,i,LDX),&COL(X,j,LDX));
      dvec_swap_at(lambda,i,j);
    }
  }
  a=dvec_free(a);
}

void deig_sort_value_guide(int n, int k, double *lambda, double *X, int LDX, double *lambda0)
{
  int i,j;
  double value,*a;
  a=dvec_allocate(k);
  for(j=0; j<k; j++){
    dvec_abs_sub_dvec_dscalar(k,a,lambda,lambda0[j]);
    value=dmin_abs_dvec_index(k,a,&i);
    if(i!=j){
      dvec_swap(n,&COL(X,i,LDX),&COL(X,j,LDX));
      dvec_swap_at(lambda,i,j);
    }
  }
  a=dvec_free(a);
}



// EOF
