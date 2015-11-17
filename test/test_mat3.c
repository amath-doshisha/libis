#include<isys.h>

int main()
{
  int M=3,N=4,L=2,LDA1=M,LDA2=N,LDB1=M,LDB2=N,LDC1=M,LDC2=N,LDX1=M,LDX2=N;
  rmulti **A=NULL,**B=NULL,**C=NULL;
  cmulti **X=NULL;

  // init
  A=rmat3_allocate(LDA1,LDA2,L);
  B=rmat3_allocate(LDB1,LDB2,L);
  C=rmat3_allocate(LDC1,LDC2,L);
  X=cmat3_allocate(LDX1,LDX2,L);

  // body
  rmat3_set_zeros(M,N,L,A,LDA1,LDA2);                   // A=0
  rmat3_set_ones(M,N,L,B,LDB1,LDB2);                    // B=1
  rmat3_neg(M,N,L,A,LDA1,LDA2,B,LDB1,LDB2);             // A=-B
  rmat3_set_rand(M,N,L,B,LDB1,LDB2,2,-1);               // B=rand(M,N,L)
  rmat3_add(M,N,L,C,LDC1,LDC2,A,LDA1,LDA2,B,LDB1,LDB2); // C=A+B
  rmat3_print(M,N,L,A,LDA1,LDA2,"A=","f",3);
  rmat3_print(M,N,L,B,LDB1,LDB2,"B=","f",3);
  rmat3_print(M,N,L,C,LDC1,LDC2,"C=","f",3);
  cmat3_print(M,N,L,X,LDX1,LDX2,"X=","f",3);

  // done
  A=rmat3_free(LDA1,LDA2,L,A);
  B=rmat3_free(LDB1,LDB2,L,B);
  C=rmat3_free(LDC1,LDC2,L,C);
  X=cmat3_free(LDX1,LDX2,L,X);
}
