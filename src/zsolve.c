#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<stdarg.h>
#include"is_macros.h"
#include"is_ivec.h"
#include"is_zmat.h"
#include"is_zsolve.h"

#define MA(I,J) (MAT(A,(I),(J),LDA))
#define MB(I,J) (MAT(B,(I),(J),LDB))
#define MC(I,J) (MAT(C,(I),(J),LDC))

// Sove A*X=B by LU decomposition
// X: input=B, output=X
// A: input=A, output=destroyed
int zsolve(int n, int NRHS, dcomplex *B, int LDB, dcomplex *A, int LDA)
{
  return zsolve_gauss_LU(n,NRHS,B,LDB,A,LDA);
}

// Sove A*X=B by LU decomposition
// X: input=B, output=X
// A: input=A, output=destroyed
int zsolve_gauss_LU(int n, int NRHS, dcomplex *B, int LDB, dcomplex *A, int LDA)
{
  int i,j,k,l,*p=NULL;
  double value,foo;
  dcomplex a,b;
  //allocate
  p=ivec_allocate(n);
  //LU分解
  for(k=0; k<n; k++){
    //pivot select
    value=Z_ABS(MA(k,k));
    for(l=k,j=k+1; j<n; j++){
      foo=Z_ABS(MA(j,k));
      if(foo>value) { l=j; value=foo; }
    }
    //エラー処理
    if(value==0) return l+1;
    p[k]=l;
    //swap A(k,k:end)<-> A(l,k:end)
    if(l!=k){
      for(j=k; j<n; j++){
	a=MA(k,j);
	MA(k,j)=MA(l,j);
	MA(l,j)=a;
      }
    }
    //L行列とU行列の作成
    for(i=k+1; i<n; i++){
      Z_SET_DIV(a,MA(i,k),MA(k,k));
      //L行列 A(i,k)=a
      MA(i,k)=a;
      //U行列 A(i,k+1:end)+=(-a)*A(k,k+1:end)
      for(j=k+1; j<n; j++){
	Z_SUB_TIMES(MA(i,j),a,MA(k,j));
      }
    }
  }
  //A*x=b -> A=L*U,L*y=b,U*x=y として解く
  //y=inv(L)*b
  for(k=0; k<n-1; k++){
    //swap b(k)<-> b(p[k])
    zmat_swap_rows(n,NRHS,B,LDB,k,p[k]);
    for(i=k+1; i<n; i++){
      for(j=0; j<NRHS; j++){
	Z_SUB_TIMES(MB(i,j),MA(i,k),MB(k,j));
      }
    }
  }
  //x=inv(U)*y
  for(k=n-1; k>=0; k--){
    for(j=0; j<NRHS; j++){
      Z_SET(a,0,0);
      for(i=k+1; i<n; i++){
	Z_ADD_TIMES(a,MA(k,i),MB(i,j));
      }
      Z_SUB(MB(k,j),a);
      Z_SET_DIV(b,MB(k,j),MA(k,k));
      MB(k,j)=b;
    }
  }
  //free
  p=ivec_free(p);
  return 0;
}

// Sove A*X=B by Gauss sweeper
// X: input=B, output=X
// A: input=A, output=destroyed
int zsolve_gauss_sweeper(int n, int NRHS, dcomplex *B, int LDB, dcomplex *A, int LDA)
{
  int i,j,k,l;
  dcomplex a,b;
  double value,foo;
  for(k=0; k<n; k++){
    //pivot select
    value=Z_ABS(MA(k,k));
    for(l=k,j=k+1; j<n; j++){
      foo=Z_ABS(MA(j,k));
      if(foo>value) { l=j; value=foo; }
    }
    if(value==0) return l+1; // error
    if(l!=k){
      zmat_swap_rows(n,NRHS,B,LDB,k,l); //swap b(k) <-> b(l)
      zmat_swap_rows(n,n,A,LDA,k,l);    //swap A(k,:) <-> A(l,:)
    }
    //軸要素を1にする
    Z_SET_INV(a,MA(k,k));
    for(j=k; j<n; j++){
      Z_SET_TIMES(b,a,MA(k,j));
      MA(k,j)=b;
    }
    for(j=0; j<NRHS; j++){
      Z_SET_TIMES(b,a,MB(k,j));
      MB(k,j)=b;
    }
    //軸要素以外が 0 になるように他の列から軸要素の列を引く
    for(i=k+1; i<n; i++){
      if(i!=k){
	a=MA(i,k);
	for(j=0; j<n; j++){
	  Z_SUB_TIMES(MA(i,j),a,MA(k,j));
	}
	for(j=0; j<NRHS; j++){
	  Z_SUB_TIMES(MB(i,j),a,MB(k,j));
	}
      }
    }
  }  
  //後退代入 x=inv(U)*b
  for(k=0; k<NRHS; k++){
    for(i=n-1; i>=0; i--){
      for(j=n-1; j>=i+1; j--){
	Z_SUB_TIMES(MB(i,k),MA(i,j),MB(j,k));
      }
    }
  }
  return 0;
}


// R=B-A*X
void zsolve_residual(int n, int NRHS, dcomplex *R, int LDR, const dcomplex *A, int LDA, const dcomplex *X, int LDX, const dcomplex *B, int LDB)
{
  zmat_copy(n,NRHS,R,LDR,B,LDB);              // R=B
  zmat_sub_prod(n,n,NRHS,R,LDR,A,LDA,X,LDX);  // R=R-A*X
}

//EOF
