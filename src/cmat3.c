#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<stdarg.h>
#include"is_macros.h"
#include"is_dcomplex.h"
#include"is_rmulti.h"
#include"is_cmulti.h"
#include"is_cmat.h"
#include"is_rmat.h"
#include"is_rmat3.h"
#include"is_zmat.h"
#include"is_dmat.h"
#include"is_cvec.h"
#include"is_cmat.h"
#include"is_cmat3.h"
#include"is_func.h"

/**
 @file  cmat.c
 @brief 多倍長精度複素数型cmultiの3次元配列に関する関数の定義
 @details cmulti型の関数に関する定義は@link cmulti.c@endlinkを参照のこと.
          cmulti型のベクトルの関数に関する定義は@link cvec.c@endlinkを参照のこと.
          cmulti型の行列の関数に関する定義は@link cmat.c@endlinkを参照のこと.
 */

/////////////////////////////////////////////////////

#define FILE_NAME_LENGTH_MAX 100000

/////////////////////////////////////////////////////

/** @name cmulti型の3次元配列の初期化に関する関数 */
/** @{ */


/**
 @brief cmulti型の3次元配列(LDA1,LDA2,l)の新規生成.
*/
cmulti **cmat3_allocate(int LDA1, int LDA2, int l)
{
  int i,j,k;
  cmulti **A=NULL;
  A=(cmulti**)malloc(sizeof(cmulti*)*LDA1*LDA2*l);
  for(k=0; k<l; k++){
    for(j=0; j<LDA2; j++){
      for(i=0; i<LDA1; i++){
	MAT3(A,i,j,k,LDA1,LDA2)=callocate();
      }
    }
  }
  return A;
}

/**
 @brief cmulti型の3次元配列A(LDA1,LDA2,l)の終了処理.
*/
cmulti **cmat3_free(int LDA1, int LDA2, int l, cmulti **A)
{
  int i,j,k;
  if(A==NULL) return NULL;
  for(k=0; k<l; k++){
    for(j=0; j<LDA2; j++){
      for(i=0; i<LDA1; i++){
	MAT3(A,i,j,k,LDA1,LDA2)=cfree(MAT3(A,i,j,k,LDA1,LDA2));
      }
    }
  }
  free(A);
  return A=NULL;
}

/** @} */

/////////////////////////////////////////////////////

/** @name cmulti型の3次元配列のメンバ変数に関する関数 */
/** @{ */

/** @} */

////////////////////////////////////////////////////////////////////////

/** @name cmulti型の3次元配列の入出力に関する関数 */
/** @{ */

/**
 @brief cmulti型の行列の表示.
*/
void cmat3_print(int m, int n, int l, cmulti **A, int LDA1, int LDA2, const char *name, const char *f, int digits)
{
  int i,j,k,p;
  char format[128];
  if(STR_EQ(f,"f") || STR_EQ(f,"F")){
    p=3;
    sprintf(format,"%%%d.%dR%s %%%d.%dR%s  ",digits+p,digits,f,digits+p,digits,f);
  }
  else if((strcmp(f,"e")==0) || (strcmp(f,"E")==0)){
    p=7;
    sprintf(format,"%%%d.%dR%s %%%d.%dR%s  ",digits+p,digits,f,digits+p,digits,f);
  }
  else{
    p=3;
    sprintf(format,"(%%%d.%dR%s, %%%d.%dR%s) ",digits+p,digits,f,digits+p,digits,f);
  }
  if(A==NULL){
    if(name!=NULL){ printf("%s=NULL\n",name); }
    else          { printf("NULL\n"); }
    return;
  }
  if(name!=NULL){ printf("%s\n",name); }
  if(A==NULL || m<=0 || n<=0) return;
  for(k=0; k<l; k++){
    for(i=0; i<m; i++){
      for(j=0; j<n; j++){
	mpfr_printf(format,C_R(MAT3(A,i,j,k,LDA1,LDA2)),C_I(MAT3(A,i,j,k,LDA1,LDA2)));
      }
      printf("\n");
    }
  }
}


/** @} */

/////////////////////////////////////////////////////////

/** @name cmulti型の3次元配列の値の設定に関する関数 */

/** @} */

////////////////////////////////////////////////////////////////////////

/** @name cmulti型の3次元配列の型変換に関する関数 */
/** @{ */

/** @} */

////////////////////////////////////////////////////////////////////////

/** @name cmulti型の3次元配列の要素の配置に関する関数 */
/** @{ */

/** @} */

////////////////////////////////////////////////////////////////////////

/** @name cmulti型の3次元配列の自動精度調整モードが機能する関数 */
/** @{ */

/** @} */

///////////////////////////////////

/** @name cmulti型の3次元配列の数学関数に関する関数 */
/** @{ */

/** @} */

////////////////////////////////////////////////////////////////////////

/** @name cmulti型の3次元配列の要素の比較に関する関数 */
/** @{ */

/** @} */

////////////////////////////////////////////////////////////////////

/** @name cmulti型の写像に関する関数 */
/** @{ */

/** @} */

//EOF
