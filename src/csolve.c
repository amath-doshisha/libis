#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"is_macros.h"
#include"is_cmulti.h"
#include"is_ivec.h"
#include"is_cmat.h"
#include"is_csolve.h"

/**
 @file  csolve.c
 @brief 多倍長精度実数型cmultiによる線形方程式A*X=Bの解法に関する関数の定義.
 @details cmulti型の関数に関する定義は@link cmulti.c@endlinkを参照のこと.
          cmulti型のベクトルの関数に関する定義は@link cvec.c@endlinkを参照のこと.
          cmulti型の行列の関数に関する定義は@link cmat.c@endlinkを参照のこと.
 */

#define MA(I,J) (MAT(A,(I),(J),LDA))
#define MB(I,J) (MAT(B,(I),(J),LDB))
#define MC(I,J) (MAT(C,(I),(J),LDC))

/////////////////////////////////////////////////////////////////////////////////////////////////////

/** @name 線形方程式 A*X=B の解法に関する関数 */
/** @{ */

/**
 @brief 線形方程式A*X=Bの解法.
 @param[in]  B サイズ(n,NRHS)の行列B.
 @param[in]  B サイズ(n,NRHS)の解の行列X.
 @param[in]  A サイズnの正方行列で係数行列A.
 */
int csolve(int n, int NRHS, cmulti **B, int LDB, cmulti **A, int LDA, int *info)
{
  int e=0;
  cmulti **C=NULL;
  C=cmat_allocate_prec(n,n,cmat_get_prec_max(n,NRHS,B,LDB));
  e+=cmat_copy(n,n,C,n,A,LDA);
  //  e+=csolve_lu(n,NRHS,B,LDB,C,n,info);
  e+=csolve_gauss_sweeper(n,NRHS,B,LDB,C,n,info);
  C=cmat_free(n,n,C);
  return e;
}

/**
 @brief 線形方程式A*X=Bの残差 R=B-A*X
 @param[in]  A サイズ(n,n)の係数行列A.
 @param[in]  B サイズ(n,NRHS)の非同次項の行列B.
 @param[in]  X サイズ(n,NRHS)の解の行列X.
 @param[out] R サイズ(n,NRHS)の残差の行列X.
*/
int csolve_residual(int n, int NRHS, cmulti **R, int LDR, cmulti **A, int LDA, cmulti **X, int LDX, cmulti **B, int LDB)
{
  int e=0;
  e+=cmat_copy(n,NRHS,R,LDR,B,LDB);              // R=B
  e+=cmat_sub_prod(n,n,NRHS,R,LDR,A,LDA,X,LDX);  // R=R-A*X
  return e;
}

/** @} */

///////////////////////////////////////////////////////////////////////////////////////////////////////////

/** @name LU分解に関する関数 */
/** @{ */

/**
 @brief 線形方程式A*X=BのLU分解による解法.
 @param[in]  B 行列B.サイズは(n,NRHS).
 @param[out] B 解の行列X.サイズは(n,NRHS).
 @param[in]  A 係数行列A.サイズは(n,n).
 @param[out] A 破壊される.
 */
int csolve_lu(int n, int NRHS, cmulti **B, int LDB, cmulti **A, int LDA, int *info)
{
  int ret=0,*p=NULL,e=0;
  // allocate
  p=ivec_allocate(n);
  // LU decomposition
  e+=csolve_lu_decomp(n,A,LDA,p,&ret);
  // solve by A=L*U,L*y=b,U*x=y
  if(ret==0){ e+=csolve_lu_backsubs(n,NRHS,B,LDB,A,LDA,p); }
  // free
  p=ivec_free(p);
  // done
  (*info)=ret;
  return e;
}

/**
 @brief 行列AのLU分解.
 @param[in]  A サイズnの正方行列A.
 @param[out] A LU分解の結果.
 @param[in]  p サイズがnの初期化済みの配列.
 @param[out] p ピボット選択の結果.
 */
int csolve_lu_decomp(int n, cmulti **A, int LDA, int *p, int *info)
{
  int prec,ret=0,e=0,i,j,k,l;
  cmulti *a=NULL;
  rmulti *c=NULL,*value=NULL;
  // allocate
  prec=cmat_get_prec_max(n,n,A,LDA);
  a=callocate_prec(prec);
  c=rallocate_prec(prec);
  value=rallocate_prec(prec);
  // LU分解
  for(k=0; k<n; k++){
    // pivot select
    e+=cabsv(value,MA(k,k));
    for(l=k,j=k+1; j<n; j++){
      e+=cabsv(c,MA(j,k));
      if(rgt(c,value)) { l=j; e+=rcopy(value,c); } // value=c
    }
    // エラー処理
    if(ris_zero(value)) { ret=l+1; break; }
    p[k]=l;
    // swap A(k,k:end)<-> A(l,k:end)
    if(l!=k){
      for(j=k; j<n; j++){
	cswap(MA(k,j),MA(l,j));
      }
    }
    // L行列とU行列の作成
    for(i=k+1; i<n; i++){
      e+=cdiv(a,MA(i,k),MA(k,k));
      // L行列 A(i,k)=a
      e+=ccopy(MA(i,k),a);
      // U行列 A(i,k+1:end)+=(-a)*A(k,k+1:end)
      for(j=k+1; j<n; j++){
	e+=csub_mul(MA(i,j),a,MA(k,j));
      }
    }
  }
  // done
  a=cfree(a);
  c=rfree(c);
  value=rfree(value);
  (*info)=ret;
  return e;
}

/**
 @brief 線形方程式A*X=Bを後退代入 A=L*U,L*y=b,U*x=y で解く.
 @param[in]  B サイズは(n,NRHS)の行列B.
 @param[out] B 解の行列X.
 @param[in]  A サイズがnの正方行列でLU分解済みの係数行列A.
 @param[in]  p ピボット選択の要素の並び.
 */
int csolve_lu_backsubs(int n, int NRHS, cmulti **B, int LDB, cmulti **A, int LDA, int *p)
{
  int prec,e=0,i,j,k;
  cmulti *a=NULL,*b=NULL;
  // allocate
  prec=cmat_get_prec_max(n,NRHS,B,LDB);
  a=callocate_prec(prec);
  b=callocate_prec(prec);
  // y=inv(L)*b
  for(k=0; k<n-1; k++){
    // swap b(k)<-> b(p[k])
    cmat_swap_rows(n,NRHS,B,LDB,k,p[k]);
    for(i=k+1; i<n; i++){
      for(j=0; j<NRHS; j++){
	e+=csub_mul(MB(i,j),MA(i,k),MB(k,j));
      }
    }
  }
  // x=inv(U)*y
  for(k=n-1; k>=0; k--){
    for(j=0; j<NRHS; j++){
      e+=cset_d(a,0);
      for(i=k+1; i<n; i++){
	e+=cadd_mul(a,MA(k,i),MB(i,j));
      }
      e+=csub(MB(k,j),MB(k,j),a);
      e+=cdiv(b,MB(k,j),MA(k,k));
      e+=ccopy(MB(k,j),b);
    }
  }
  // done
  a=cfree(a);
  b=cfree(b);
  return e;
}

/** @} */

////////////////////////////////////////////////////////////////////////////////////////////////

/** @name ガウスの消去法に関する関数 */
/** @{ */

/**
 @brief 線形方程式A*X=Bのガウスの消去法による解法.
 @param[in]  B 行列B.サイズは(n,NRHS).
 @param[in]  B 解の行列X.サイズは(n,NRHS).
 @param[in]  A 係数行列A.サイズは(n,n).
 @param[out] A 破壊される.
 */
int csolve_gauss_sweeper(int n, int NRHS, cmulti **B, int LDB, cmulti **A, int LDA, int *info)
{
  int prec,p0,p1,ret=0,i,j,k,l,e=0;
  cmulti *a=NULL,*b=NULL;
  rmulti *value=NULL,*c=NULL;
  // allocate
  p0=cmat_get_prec_max(n,NRHS,B,LDB);
  p1=cmat_get_prec_max(n,n,A,LDA);
  prec=MAX2(p0,p1);
  a=callocate_prec(prec);
  b=callocate_prec(prec);
  c=rallocate_prec(prec);
  value=rallocate_prec(prec);
  for(k=0; k<n; k++){
    //pivot select
    e+=cabsv(value,MA(k,k));
    for(l=k,j=k+1; j<n; j++){
      e+=cabsv(c,MA(j,k));
      if(rgt(c,value)) { l=j; e+=rcopy(value,c); } // value=c
    }
    if(ris_zero(value)) { ret=l+1; break; } // error
    if(l!=k){
      cmat_swap_rows(n,NRHS,B,LDB,k,l);       //swap b(k) <-> b(l)
      if(l!=k) cmat_swap_rows(n,n,A,LDA,k,l); //swap A(k,:) <-> A(l,:)
    }
    //軸要素を1にする
    e+=cinv(a,MA(k,k));
    for(j=k; j<n; j++)    { e+=cmul(MA(k,j),a,MA(k,j)); }
    for(j=0; j<NRHS; j++) { e+=cmul(MB(k,j),a,MB(k,j)); }
    //軸要素以外が 0 になるように他の列から軸要素の列を引く
    for(i=k+1; i<n; i++){
      if(i!=k){
	e+=ccopy(a,MA(i,k));
	for(j=0; j<n; j++){
	  e+=csub_mul(MA(i,j),a,MA(k,j));
	}
	for(j=0; j<NRHS; j++){
	  e+=csub_mul(MB(i,j),a,MB(k,j));
	}
      }
    }
  }  
  //後退代入 x=inv(U)*b
  if(ret==0){
    for(k=0; k<NRHS; k++){
      for(i=n-1; i>=0; i--){
	for(j=n-1; j>=i+1; j--){
	  e+=csub_mul(MB(i,k),MA(i,j),MB(j,k));
	}
      }
    }
  }
  // free
  a=cfree(a);
  b=cfree(b);
  c=rfree(c);
  value=rfree(value);
  // done
  (*info)=ret;
  return e;
}

/** @} */

//EOF
