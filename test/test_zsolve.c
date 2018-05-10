#include<isys.h>

int main()
{
  int n=10,LDA;
  dcomplex *A0=NULL,*A=NULL,*b=NULL,*x=NULL,*x0=NULL,*r=NULL;

  // allocation
  LDA=n; A=zmat_allocate(LDA,n);
  LDA=n; A0=zmat_allocate(LDA,n);
  b=zvec_allocate(n);
  x=zvec_allocate(n);
  x0=zvec_allocate(n);
  r=zvec_allocate(n);

  // initialize
  init_genrand(0);
  zmat_set_rand(n,n,A,LDA,2,-1);
  zmat_copy(n,n,A0,LDA,A,LDA);
  zvec_set_ones(n,x0);
  zvec_lintr(n,n,b,A,LDA,x0);
  zvec_copy(n,x,b);

  printf("before\n");
  zmat_print(n,n,A,LDA,"A=",'f',2);
  zvec_print(n,b,"b=",'f',2);
  zvec_print(n,x0,"x0=",'f',2);

  printf("after\n");
  zsolve(n,1,x,n,A,LDA);
  zmat_print(n,n,A,LDA,"A=",'f',2);
  zvec_print(n,x,"x=",'f',2);
  zvec_print(n,x0,"x0=",'f',2);

  // residual
  zsolve_residual(n,1,r,n,A0,LDA,x,n,b,n);
  zvec_print(n,r,"r=b-A0*x=",'e',2);

  zvec_sub_zvec(n,x,x,x0);
  zvec_print(n,x,"x-x0=",'e',2);
  printf("||x-x0||=%.2e\n",zvec_norm_max(n,x));


  // free
  A0=zmat_free(A0); 
  A=zmat_free(A); 
  b=zvec_free(b); 
  x=zvec_free(x); 
  x0=zvec_free(x0); 
  r=zvec_free(r);

  // done
  return 0;
}
