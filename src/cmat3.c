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

/**
 @brief cmulti型の3次元配列の値を倍精度浮動小数点数から設定.
*/
int cmat3_set_all_z(int m, int n, int l, cmulti **A, int LDA1, int LDA2, dcomplex a)
{
  int i,j,k,e=0;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	e+=cset_z(MAT3(A,i,j,k,LDA1,LDA2),a);
      }
    }
  }
  return e;
}

/**
 @brief cmulti型の3次元配列の値を倍精度浮動小数点数から設定.
*/
int cmat3_set_all_dd(int m, int n, int l, cmulti **A, int LDA1, int LDA2, double a_r, double a_i)
{
  int i,j,k,e=0;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	e+=cset_dd(MAT3(A,i,j,k,LDA1,LDA2),a_r,a_i);
      }
    }
  }
  return e;
}

/**
 @brief cmulti型の3次元配列の値を倍精度浮動小数点数から設定.
*/
int cmat3_set_all_d(int m, int n, int l, cmulti **A, int LDA1, int LDA2, double a)
{
  int i,j,k,e=0;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	e+=cset_d(MAT3(A,i,j,k,LDA1,LDA2),a);
      }
    }
  }
  return e;
}

/**
 @brief cmulti型の3次元配列の値を零に設定.
*/
int cmat3_set_zeros(int m, int n, int l, cmulti **A, int LDA1, int LDA2)
{
  return cmat3_set_all_d(m,n,l,A,LDA1,LDA2,0);
}

/**
 @brief cmulti型の3次元配列の値を1に設定.
*/
int cmat3_set_ones(int m, int n, int l, cmulti **A, int LDA1, int LDA2)
{
  return cmat3_set_all_d(m,n,l,A,LDA1,LDA2,1);
}

/**
 @brief cmulti型の3次元配列の値を区間(b,a+b)の疑似乱数値を設定.
*/
void cmat3_set_rand(int m, int n, int l, cmulti **A, int LDA1, int LDA2, double a, double b)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	cset_rand(MAT3(A,i,j,k,LDA1,LDA2));
	cmul_d(MAT3(A,i,j,k,LDA1,LDA2),MAT3(A,i,j,k,LDA1,LDA2),a);
	cadd_d(MAT3(A,i,j,k,LDA1,LDA2),MAT3(A,i,j,k,LDA1,LDA2),b);
      }
    }
  }
}


/** @} */

////////////////////////////////////////////////////////////////////////

/** @name cmulti型の3次元配列の型変換に関する関数 */
/** @{ */

/**
 @brief rmulti型の3次元配列を倍精度型に変換 B=dcomplex(A).
 */
void cmat3_get_z(int m, int n, int l, dcomplex *B, int LDB1, int LDB2, cmulti **A, int LDA1, int LDA2)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	MAT3(B,i,j,k,LDB1,LDB2)=cget_z(MAT3(A,i,j,k,LDA1,LDA2));
      }
    }
  }
  return;
}


/** @} */

////////////////////////////////////////////////////////////////////////

/** @name cmulti型の3次元配列の要素の配置に関する関数 */
/** @{ */

/** @} */

////////////////////////////////////////////////////////////////////////

/** @name cmulti型の3次元配列の自動精度調整モードが機能する関数 */
/** @{ */

/**
 @brief cmulti型の3次元配列の値のコピー B=A
*/
int cmat3_copy(int m, int n, int l, cmulti **B, int LDB1, int LDB2, cmulti **A, int LDA1, int LDA2)
{
  int i,j,k,e=0;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	e+=ccopy(MAT3(B,i,j,k,LDB1,LDB2),MAT3(A,i,j,k,LDA1,LDA2));
      }
    }
  }
  return e;
}

/**
 @brief cmulti型の3次元配列の符号反転 B=-A.
*/
int cmat3_neg(int m, int n, int l, cmulti **B, int LDB1, int LDB2, cmulti **A, int LDA1, int LDA2)
{
  int i,j,k,e=0;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	e+=cneg(MAT3(B,i,j,k,LDB1,LDB2),MAT3(A,i,j,k,LDA1,LDA2));
      }
    }
  }
  return e;
}

/**
 @brief cmulti型の3次元配列の足し算 C=A+B.
*/
int cmat3_add(int m, int n, int l, cmulti **C, int LDC1, int LDC2, cmulti **A, int LDA1, int LDA2, cmulti **B, int LDB1, int LDB2)
{
  int i,j,k,e=0;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	e+=cadd(MAT3(C,i,j,k,LDC1,LDC2),MAT3(A,i,j,k,LDA1,LDA2),MAT3(B,i,j,k,LDB1,LDB2));
      }
    }
  }
  return e;
}

/**
 @brief cmulti型の3次元配列の足し算 C=A+B.
*/
int cmat3_add_rmat3(int m, int n, int l, cmulti **C, int LDC1, int LDC2, cmulti **A, int LDA1, int LDA2, rmulti **B, int LDB1, int LDB2)
{
  int i,j,k,e=0;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	e+=cadd_r(MAT3(C,i,j,k,LDC1,LDC2),MAT3(A,i,j,k,LDA1,LDA2),MAT3(B,i,j,k,LDB1,LDB2));
      }
    }
  }
  return e;
}

/**
 @brief cmulti型の3次元配列の足し算 C=A+b.
*/
int cmat3_add_c(int m, int n, int l, cmulti **C, int LDC1, int LDC2, cmulti **A, int LDA1, int LDA2, cmulti *b)
{
  int i,j,k,e=0;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	e+=cadd(MAT3(C,i,j,k,LDC1,LDC2),MAT3(A,i,j,k,LDA1,LDA2),b);
      }
    }
  }
  return e;
}

/**
 @brief cmulti型の3次元配列の足し算 C=A+b.
*/
int cmat3_add_r(int m, int n, int l, cmulti **C, int LDC1, int LDC2, cmulti **A, int LDA1, int LDA2, rmulti *b)
{
  int i,j,k,e=0;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	e+=cadd_r(MAT3(C,i,j,k,LDC1,LDC2),MAT3(A,i,j,k,LDA1,LDA2),b);
      }
    }
  }
  return e;
}

/**
 @brief cmulti型の3次元配列の足し算 C=A+b.
*/
int rmat3_add_c(int m, int n, int l, cmulti **C, int LDC1, int LDC2, rmulti **A, int LDA1, int LDA2, cmulti *b)
{
  int i,j,k,e=0;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	e+=cadd_r(MAT3(C,i,j,k,LDC1,LDC2),b,MAT3(A,i,j,k,LDA1,LDA2));
      }
    }
  }
  return e;
}


/** @} */

///////////////////////////////////

/** @name cmulti型の3次元配列の数学関数に関する関数 */
/** @{ */

/** @} */

////////////////////////////////////////////////////////////////////////

/** @name cmulti型の3次元配列の要素の比較に関する関数 */
/** @{ */


/**
 @brief cmulti型の3次元配列の比較　C=(A==B).
*/
void cmat3_eq(int m, int n, int l, int *C, int LDC1, int LDC2, cmulti **A, int LDA1, int LDA2, cmulti **B, int LDB1, int LDB2)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	MAT3(C,i,j,k,LDC1,LDC2)=ceq(MAT3(A,i,j,k,LDA1,LDA2),MAT3(B,i,j,k,LDB1,LDB2));
      }
    }
  }
  return;
}

/**
 @brief cmulti型の3次元配列の比較　C=(A==B).
*/
void cmat3_eq_rmat3(int m, int n, int l, int *C, int LDC1, int LDC2, cmulti **A, int LDA1, int LDA2, rmulti **B, int LDB1, int LDB2)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	MAT3(C,i,j,k,LDC1,LDC2)=ceq_r(MAT3(A,i,j,k,LDA1,LDA2),MAT3(B,i,j,k,LDB1,LDB2));
      }
    }
  }
  return;
}

/**
 @brief cmulti型の3次元配列の比較　C=(A==b).
*/
void cmat3_eq_c(int m, int n, int l, int *C, int LDC1, int LDC2, cmulti **A, int LDA1, int LDA2, cmulti *b)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	MAT3(C,i,j,k,LDC1,LDC2)=ceq(MAT3(A,i,j,k,LDA1,LDA2),b);
      }
    }
  }
  return;
}

/**
 @brief cmulti型の3次元配列の比較　C=(A==b).
*/
void cmat3_eq_r(int m, int n, int l, int *C, int LDC1, int LDC2, cmulti **A, int LDA1, int LDA2, rmulti *b)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	MAT3(C,i,j,k,LDC1,LDC2)=ceq_r(MAT3(A,i,j,k,LDA1,LDA2),b);
      }
    }
  }
  return;
}

/**
 @brief cmulti型の3次元配列の比較　C=(A==b).
*/
void rmat3_eq_c(int m, int n, int l, int *C, int LDC1, int LDC2, rmulti **A, int LDA1, int LDA2, cmulti *b)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	MAT3(C,i,j,k,LDC1,LDC2)=ceq_r(b,MAT3(A,i,j,k,LDA1,LDA2));
      }
    }
  }
  return;
}


/** @} */

////////////////////////////////////////////////////////////////////

/** @name cmulti型の写像に関する関数 */
/** @{ */

/** @} */

//EOF
