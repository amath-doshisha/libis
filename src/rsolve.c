#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"is_macros.h"
#include"is_rmulti.h"
#include"is_ivec.h"
#include"is_rmat.h"
#include"is_rsolve.h"
#include"is_rvec.h"

/**
 @file  rsolve.c
 @brief 多倍長精度実数型rmultiによる線形方程式A*X=Bの解法に関する関数の定義.
 @details rmulti型の関数に関する定義は@link rmulti.c@endlinkを参照のこと.
          rmulti型のベクトルの関数に関する定義は@link rvec.c@endlinkを参照のこと.
          rmulti型の行列の関数に関する定義は@link rmat.c@endlinkを参照のこと.
 */

#define MA(I,J) (MAT(A,(I),(J),LDA))
#define MB(I,J) (MAT(B,(I),(J),LDB))
#define MC(I,J) (MAT(C,(I),(J),LDC))

///////////////////////////////////////////////////////////////////////////////////////////////////////////

/** @name 線形方程式 A*X=B の解法に関する関数 */
/** @{ */

/**
 @brief 線形方程式A*X=Bの解法.
 @param[in,out]  B  [in]サイズ(n,NRHS)の行列B.[out]サイズ(n,NRHS)の解の行列X.
 @param[in]      A  サイズnの正方行列で係数行列A.
 */
void rsolve(int n, int NRHS, rmulti **B, int LDB, rmulti **A, int LDA, int *info)
{
  rmulti **C=NULL;
  C=rmat_allocate_prec(n,n,rmat_get_prec_max(n,NRHS,B,LDB));
  rmat_copy(n,n,C,n,A,LDA);
  //e+=rsolve_lu(n,NRHS,B,LDB,C,n,info);
  rsolve_gauss_sweeper(n,NRHS,B,LDB,C,n,info);
  C=rmat_free(n,n,C);
}

/**
 @brief 線形方程式A*X=Bの残差 R=B-A*X
 @param[in]  A サイズ(n,n)の係数行列A.
 @param[in]  B サイズ(n,NRHS)の非同次項の行列B.
 @param[in]  X サイズ(n,NRHS)の解の行列X.
 @param[out] R サイズ(n,NRHS)の残差の行列X.
*/
void rsolve_residual(int n, int NRHS, rmulti **R, int LDR, rmulti **A, int LDA, rmulti **X, int LDX, rmulti **B, int LDB)
{
  rmat_copy(n,NRHS,R,LDR,B,LDB);              // R=B
  rmat_sub_prod(n,n,NRHS,R,LDR,A,LDA,X,LDX);  // R=R-A*X
}

/** @} */

///////////////////////////////////////////////////////////////////////////////////////////////////////////

/** @name LU分解に関する関数 */
/** @{ */

/**
 @brief 線形方程式A*X=BのLU分解による解法.
 @param[in,out]  B   [in]行列B.サイズは(n,NRHS).[out]解の行列X.サイズは(n,NRHS).
 @param[in,out]  A   [in]係数行列A.サイズは(n,n).[out]破壊される.
 */
void rsolve_lu(int n, int NRHS, rmulti **B, int LDB, rmulti **A, int LDA, int *info)
{
  int ret=0,*p=NULL;
  // allocate
  p=ivec_allocate(n);
  // LU decomposition
  rsolve_lu_decomp(n,A,LDA,p,&ret);
  // solve by A=L*U,L*y=b,U*x=y
  if(ret==0){ rsolve_lu_backsubs(n,NRHS,B,LDB,A,LDA,p); }
  // free
  p=ivec_free(p);
  // done
  (*info)=ret;
}


/**
 @brief 行列AのLU分解.
 @param[in,out]  A   [in]サイズnの正方行列A.[out]LU分解の結果.
 @param[in,out]  p   [inサイズがnの初期化済みの配列.[out]ピボット選択の結果.
 */
void rsolve_lu_decomp(int n, rmulti **A, int LDA, int *p, int *info)
{
  int prec,ret=0,i,j,k,l;
  rmulti *a=NULL,*value=NULL;
  //allocate
  prec=rmat_get_prec_max(n,n,A,LDA);
  a=rallocate_prec(prec);
  value=rallocate_prec(prec);
  // LU分解
  for(k=0; k<n; k++){    
    // pivot select
    rabs(value,MA(k,k));
    for(l=k,j=k+1; j<n; j++){
      rabs(a,MA(j,k));
      if(rgt(a,value)) { l=j; rcopy(value,a); } // value=a
    }
    // エラー処理
    if(ris_zero(value)) { ret=l+1; break; }
    p[k]=l;
    // swap A(k,k:end)<-> A(l,k:end)
    if(l!=k){
      for(j=k; j<n; j++){
	rswap(MA(k,j),MA(l,j));
      }
    }
    // L行列とU行列の作成
    for(i=k+1; i<n; i++){
      rdiv(a,MA(i,k),MA(k,k));
      // L行列 A(i,k)=a
      rcopy(MA(i,k),a);
      // U行列 A(i,k+1:end)+=(-a)*A(k,k+1:end)
      for(j=k+1; j<n; j++){
	rsub_mul(MA(i,j),a,MA(k,j));
      }
    }
  }
  // done
  a=rfree(a);
  value=rfree(value);
  (*info)=ret;
}

/**
 @brief 線形方程式A*X=Bを後退代入 A=L*U,L*y=b,U*x=y で解く.
 @param[in,out]  B  [in]サイズは(n,NRHS)の行列B.[out]解の行列X.
 @param[in]      A  サイズがnの正方行列でLU分解済みの係数行列A.
 @param[in]      p  ピボット選択の要素の並び.
 */
void rsolve_lu_backsubs(int n, int NRHS, rmulti **B, int LDB, rmulti **A, int LDA, int *p)
{
  int prec,i,j,k;
  rmulti *a=NULL,*b=NULL;
  // allocate
  prec=rmat_get_prec_max(n,NRHS,B,LDB);
  a=rallocate_prec(prec);
  b=rallocate_prec(prec);
  // y=inv(L)*b
  for(k=0; k<n-1; k++){
    // swap b(k)<-> b(p[k])
    rmat_rows_swap(n,NRHS,B,LDB,k,p[k]);
    for(i=k+1; i<n; i++){
      for(j=0; j<NRHS; j++){
	rsub_mul(MB(i,j),MA(i,k),MB(k,j));
      }
    }
  }
  // x=inv(U)*y
  for(k=n-1; k>=0; k--){
    for(j=0; j<NRHS; j++){
      rset_d(a,0);
      for(i=k+1; i<n; i++){
	radd_mul(a,MA(k,i),MB(i,j));
      }
      rsub(MB(k,j),MB(k,j),a);
      rdiv(b,MB(k,j),MA(k,k));
      rcopy(MB(k,j),b);
    }
  }
  // done
  a=rfree(a);
  b=rfree(b);
}

/** @} */

//////////////////////////////////////////////////////////////////////////////////////////////////////////

/** @name ガウスの消去法に関する関数 */
/** @{ */

/**
 @brief 線形方程式A*X=Bのガウスの消去法による解法.
 @param[in,out]  B   [in]行列B.サイズは(n,NRHS).[out]解の行列X.サイズは(n,NRHS).
 @param[in,out]  A   [in]係数行列A.サイズは(n,n).[out]破壊される.
 */
void rsolve_gauss_sweeper(int n, int NRHS, rmulti **B, int LDB, rmulti **A, int LDA, int *info)
{
  int p0,p1,prec,ret=0,i,j,k,l;
  rmulti *a=NULL,*value=NULL;
  rmulti **rws=NULL;
  int rws_size; //ワークスペース
  // allocate
  p0=rmat_get_prec_max(n,NRHS,B,LDB);
  p1=rmat_get_prec_max(n,n,A,LDA);
  prec=MAX2(p0,p1);
  a=rallocate_prec(prec);
  value=rallocate_prec(prec);
  rws_size=1; rws=rvec_allocate_prec(rws_size,prec);
  // 上三角化
  for(k=0; k<n; k++){
    // ピボット選択
    rabs(value,MA(k,k));
    for(l=k,j=k+1; j<n; j++){
      rabs(a,MA(j,k));
      if(rgt(a,value)) { l=j; rcopy(value,a); } // value=a
    }
    if(ris_zero(value)) { ret=l+1; break; } // error
    if(l!=k){
      rmat_rows_swap(n,NRHS,B,LDB,k,l);          //swap b(k) <-> b(l)
      if(l!=k){ rmat_rows_swap(n,n,A,LDA,k,l); } //swap A(k,:) <-> A(l,:)
    }
    //軸要素を1にする
    rinv(a,MA(k,k));
    for(j=k; j<n; j++)    { rmul(MA(k,j),a,MA(k,j)); }
    for(j=0; j<NRHS; j++) { rmul(MB(k,j),a,MB(k,j)); }
    //軸要素以外が 0 になるように他の列から軸要素の列を引く
    for(i=k+1; i<n; i++){
      if(i!=k){
	rcopy(a,MA(i,k));
	for(j=k; j<n; j++){
	  rsub_mul_ws(MA(i,j),a,MA(k,j),rws_size,rws);
	}
	for(j=0; j<NRHS; j++){
	  rsub_mul_ws(MB(i,j),a,MB(k,j),rws_size,rws);	 
	}
      }
    }
  }  
  //後退代入 x=inv(U)*b
  if(ret==0){
    for(k=0; k<NRHS; k++){
      for(i=n-1; i>=0; i--){
	for(j=n-1; j>=i+1; j--){
	  rsub_mul_ws(MB(i,k),MA(i,j),MB(j,k),rws_size,rws);
	}
      }
    }
  }
  // free
  a=rfree(a);
  value=rfree(value);
  rws=rvec_free(rws_size,rws);
  // done
  (*info)=ret;
}

/** @} */

//EOF
