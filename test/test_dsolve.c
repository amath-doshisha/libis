#include<isys.h>

int main()
{
  int n=10,LDA;
  double *A0=NULL,*A=NULL,*b=NULL,*x=NULL,*x0=NULL,*r=NULL;

  // allocation
  LDA=n; A=dmat_allocate(LDA,n);
  LDA=n; A0=dmat_allocate(LDA,n);
  b=dvec_allocate(n);
  x=dvec_allocate(n);
  x0=dvec_allocate(n);
  r=dvec_allocate(n);

  // initialize
  init_genrand(0);
  dmat_set_rand(n,n,A,LDA,2,-1);
  dmat_copy(n,n,A0,LDA,A,LDA);
  dvec_set_ones(n,x0);
  dvec_lintr(n,n,b,A,LDA,x0);
  dvec_copy_dvec(n,x,b);


  printf("before\n");
  dmat_print(n,n,A,LDA,"A=",'f',2);
  dvec_print(n,b,"b=",'f',2);
  dvec_print(n,x0,"x0=",'f',2);

  printf("after\n");
  dsolve(n,1,x,n,A,LDA);
  dmat_print(n,n,A,LDA,"A=",'f',2);
  dvec_print(n,x,"x=",'f',2);
  dvec_print(n,x0,"x0=",'f',2);

  // residual
  dsolve_residual(n,1,r,n,A0,LDA,x,n,b,n);
  dvec_print(n,r,"r=b-A0*x=",'e',2);

  dvec_sub_dvec(n,x,x,x0);
  dvec_print(n,x,"x-x0=",'e',2);
  printf("||x-x0||=%.2e\n",dnorm_max_dvec(n,x));

  // free
  A=dmat_free(A);
  A0=dmat_free(A0);
  x=dvec_free(x);
  x0=dvec_free(x0);
  r=dvec_free(r);
  b=dvec_free(b);

  // done
  return 0;
}
