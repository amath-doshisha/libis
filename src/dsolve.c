#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<stdarg.h>
#include"is_macros.h"
#include"is_ivec.h"
#include"is_dvec.h"
#include"is_dmat.h"
#include"is_dsolve.h"

#define MA(I,J) (MAT(A,(I),(J),LDA))
#define MB(I,J) (MAT(B,(I),(J),LDB))
#define MC(I,J) (MAT(C,(I),(J),LDC))

// Sove A*X=B by LU decomposition
// X: input=B, output=X
// A: input=A, output=destroyed
int dsolve(int n, int NRHS, double *B, int LDB, double *A, int LDA)
{
  return dsolve_gauss_LU(n,NRHS,B,LDB,A,LDA);
}

// Sove A*X=B by LU decomposition
// X: input=B, output=X
// A: input=A, output=destroyed
int dsolve_gauss_LU(int n, int NRHS, double *B, int LDB, double *A, int LDA)
{
  int i,j,k,l,*p=NULL;
  double h,r;
  //allocate
  p=ivec_allocate(n);
  //LU分解
  for(k=0; k<n; k++){
    //pivot select
    for(h=fabs(MA(k,k)),l=k,j=k+1; j<n; j++){
      if(fabs(MA(j,k))>h){ l=j; h=fabs(MA(j,k)); }
    }
    //エラー処理
    if(h==0) return l+1;
    p[k]=l;
    //swap A(k,k:end)<-> A(l,k:end)
    if(l!=k){
      for(j=k; j<n; j++){
	h=MA(k,j);
	MA(k,j)=MA(l,j);
	MA(l,j)=h;
      }
   }
    //L行列とU行列の作成
    for(i=k+1; i<n; i++){
      r=MA(i,k)/MA(k,k);
      //L行列 A(i,k)=r
      MA(i,k)=r;
      //U行列 A(i,k+1:end)+=(-r)*A(k,k+1:end)
      for(j=k+1; j<n; j++) MA(i,j)+=-r*MA(k,j);
    }
  }
  //A*x=b -> A=L*U,L*y=b,U*x=y として解く
  //y=inv(L)*b
  for(k=0; k<n-1; k++){
    //swap b(k)<-> b(p[k])
    dmat_swap_rows(n,NRHS,B,LDB,k,p[k]);
    for(i=k+1; i<n; i++){
      for(j=0; j<NRHS; j++) MB(i,j)-=MA(i,k)*MB(k,j);
    }
  }
  //x=inv(U)*y
  for(k=n-1; k>=0; k--){
    for(j=0; j<NRHS; j++){
      for(h=0.0,i=k+1; i<n; i++) h+=MA(k,i)*MB(i,j);
      MB(k,j)-=h;
      MB(k,j)/=MA(k,k);
    }
  }
  //free
  p=ivec_free(p);
  return 0;
}    

// Sove A*X=B by Gauss sweeper
// X: input=B, output=X
// A: input=A, output=destroyed
int dsolve_gauss_sweeper(int n, int NRHS, double *B, int LDB, double *A, int LDA)
{
  int i,j,k,l;
  double a,value;
  for(k=0; k<n; k++){
    //pivot select
    value=fabs(MA(k,k));
    for(l=k,j=k+1; j<n; j++){
      if(fabs(MA(j,k))>value) {l=j; value=fabs(MA(j,k));}
    }
    if(value==0) return l+1; //エラー処理
    if(l!=k){
      dmat_swap_rows(n,NRHS,B,LDB,k,l); //swap b(k) <-> b(l)
      dmat_swap_rows(n,n,A,LDA,k,l); //swap A(k,:) <-> A(l,:)
    }
    //ガウスの消去法
    //軸要素を1にする
    a=(1.0)/MA(k,k);
    for(j=k; j<n; j++) MA(k,j)*=a;
    for(j=0; j<NRHS; j++) MB(k,j)*=a;
    //軸要素以外が 0 になるように他の列から軸要素の列を引く
    for(i=k+1; i<n; i++){
      if(i!=k){
	a=MA(i,k);
	for(j=0; j<n; j++) MA(i,j)-=a*MA(k,j);
	for(j=0; j<NRHS; j++) MB(i,j)-=a*MB(k,j);
      }
    }
  }  
  //後退代入 x=inv(U)*b
  for(k=0; k<NRHS; k++){
    for(i=n-1; i>=0; i--){
      for(j=n-1; j>=i+1; j--) MB(i,k)-=MA(i,j)*MB(j,k);
    }
  }
  return 0;
}

// R=B-A*X
void dsolve_residual(int n, int NRHS, double *R, int LDR, double *A, int LDA, double *X, int LDX, double *B, int LDB)
{
  dmat_copy(n,NRHS,R,LDR,B,LDB);              // R=B
  dmat_sub_prod(n,n,NRHS,R,LDR,A,LDA,X,LDX);  // R=R-A*X
}

//EOF
