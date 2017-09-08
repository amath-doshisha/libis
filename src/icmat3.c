#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<stdarg.h>
#include"is_macros.h"
#include"is_dcomplex.h"
#include"is_rmulti.h"
#include"is_cmulti.h"
#include"is_irmulti.h"
#include"is_icmulti.h"
#include"is_irmat3.h"
#include"is_icmat3.h"

/**
 @file  icmat.c
 @brief 多倍長精度複素数型cmultiの行列に関する関数の定義
 @details cmulti型の関数に関する定義は@link cmulti.c@endlinkを参照のこと.
          cmulti型のベクトルの関数に関する定義は@link cvec.c@endlinkを参照のこと.
 */

/////////////////////////////////////////////////////

/** @name cmulti型の行列の初期化に関する関数 */
/** @{ */


/** @} */

/////////////////////////////////////////////////////

/** @name cmulti型の行列のメンバ変数に関する関数 */
/** @{ */

/** @} */

////////////////////////////////////////////////////////////////////////

/** @name cmulti型の行列の入出力に関する関数 */
/** @{ */


/** @} */

/////////////////////////////////////////////////////////

/** @name cmulti型の行列の値の設定に関する関数 */
/** @{ */

/**
 @brief cmulti型の3次元配列の値のコピー [B0,B1]=[A0,A1]
*/
int icmat3_copy(int m, int n, int l, cmulti **B0, cmulti **B1, int LDB1, int LDB2, cmulti **A0, cmulti **A1, int LDA1, int LDA2)
{
  int i,j,k,e=0;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	e+=iccopy(MAT3(B0,i,j,k,LDB1,LDB2),MAT3(B1,i,j,k,LDB1,LDB2),MAT3(A0,i,j,k,LDA1,LDA2),MAT3(A1,i,j,k,LDA1,LDA2));
      }
    }
  }
  return e;
}

/**
 @brief cmulti型の3次元配列の値を倍精度実数から設定 [B0,B1]=cmulti(A).
*/
int icmat3_set_z(int m, int n, int l, cmulti **B0, cmulti **B1, int LDB1, int LDB2, dcomplex *A, int LDA1, int LDA2)
{
  int i,j,k,e=0;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	e+=icset_z(MAT3(B0,i,j,k,LDB1,LDB2),MAT3(B1,i,j,k,LDB1,LDB2),MAT3(A,i,j,k,LDA1,LDA2));
      }
    }
  }
  return e;
}

/**
 @brief cmulti型の3次元配列の値を倍精度実数から設定 [B0,B1]=cmulti(A).
*/
int icmat3_set_d(int m, int n, int l, cmulti **B0, cmulti **B1, int LDB1, int LDB2, double *A, int LDA1, int LDA2)
{
  int i,j,k,e=0;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	e+=icset_d(MAT3(B0,i,j,k,LDB1,LDB2),MAT3(B1,i,j,k,LDB1,LDB2),MAT3(A,i,j,k,LDA1,LDA2));
      }
    }
  }
  return e;
}

/** @} */




////////////////////////////////////////////////////////////////////////

/** @name その他 */
/** @{ */


/**
 @brief 区間の中心 [m-r,m+r]=[A0,A1]
*/
int icmat3_mid(int m, int n, int l, cmulti **mid, int LD1mid, int LD2mid, cmulti **A0, cmulti **A1, int LDA1, int LDA2)
{
  int i,j,k,e=0;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	e+=icmid(MAT3(mid,i,j,k,LD1mid,LD2mid),MAT3(A0,i,j,k,LDA1,LDA2),MAT3(A1,i,j,k,LDA1,LDA2));
      }
    }
  }
  return e;
}

/**
 @brief 区間の半径 [m-r,m+r]=[A0,A1]
*/
int icmat3_rad(int m, int n, int l, cmulti **rad, int LD1rad, int LD2rad, cmulti **A0, cmulti **A1, int LDA1, int LDA2)
{
  int i,j,k,e=0;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	e+=icrad(MAT3(rad,i,j,k,LD1rad,LD2rad),MAT3(A0,i,j,k,LDA1,LDA2),MAT3(A1,i,j,k,LDA1,LDA2));
      }
    }
  }
  return e;
}

/**
 @brief 区間の中心と半径 [m-r,m+r]=[A0,A1]
*/
int icmat3_mr(int m, int n, int l, cmulti **mid, int LD1mid, int LD2mid, cmulti **rad, int LD1rad, int LD2rad, cmulti **A0, cmulti **A1, int LDA1, int LDA2)
{
  int i,j,k,e=0;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	e+=icmr(MAT3(mid,i,j,k,LD1mid,LD2mid),MAT3(rad,i,j,k,LD1rad,LD2rad),MAT3(A0,i,j,k,LDA1,LDA2),MAT3(A1,i,j,k,LDA1,LDA2));
      }
    }
  }
  return e;
}

/** @} */

//EOF
