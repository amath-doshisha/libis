#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<stdarg.h>
#include"is_macros.h"
#include"is_rmulti.h"
#include"is_rmat.h"
#include"is_dmat.h"
#include"is_rvec.h"
#include"is_rmat3.h"

/**
 @file  rmat.c
 @brief 多倍長精度実数型rmultiの3次元配列に関する関数の定義
 @details rmulti型の関数に関する定義は@link rmulti.c@endlinkを参照のこと.
          rmulti型のベクトルの関数に関する定義は@link rvec.c@endlinkを参照のこと.
          rmulti型の行列の関数に関する定義は@link rmat.c@endlinkを参照のこと.
 */

/////////////////////////////////////////////////////

#define FILE_NAME_LENGTH_MAX 100000

/////////////////////////////////////////////////////

/** @name rmulti型の3次元配列の初期化に関する関数 */
/** @{ */


/**
 @brief rmulti型の3次元配列(LDA1,LDA2,l)の新規生成.
*/
rmulti **rmat3_allocate(int LDA1, int LDA2, int l)
{
  int i,j,k;
  rmulti **A=NULL;
  A=(rmulti**)malloc(sizeof(rmulti*)*LDA1*LDA2*l);
  for(k=0; k<l; k++){
    for(j=0; j<LDA2; j++){
      for(i=0; i<LDA1; i++){
	MAT3(A,i,j,k,LDA1,LDA2)=rallocate();
      }
    }
  }
  return A;
}

/**
 @brief rmulti型の3次元配列(LDA1,LDA2,l)の精度を指定しての新規生成.
*/
rmulti **rmat3_allocate_prec(int LDA1, int LDA2, int l, int prec)
{
  int i,j,k;
  rmulti **A=NULL;
  A=(rmulti**)malloc(sizeof(rmulti*)*LDA1*LDA2*l);
  for(k=0; k<l; k++){
    for(j=0; j<LDA2; j++){
      for(i=0; i<LDA1; i++){
	MAT3(A,i,j,k,LDA1,LDA2)=rallocate_prec(prec);
      }
    }
  }
  return A;
}


/**
 @brief rmulti型の3次元配列の複製の生成 A=B．
*/
rmulti **rmat3_allocate_clone(int LDB1, int LDB2, int l, rmulti **B)
{
  int i,j,k,LDA1,LDA2;
  rmulti **A=NULL;
  LDA1=LDB1;
  LDA2=LDB2;
  A=(rmulti**)malloc(sizeof(rmulti*)*LDA1*LDA2*l);
  for(k=0; k<l; k++){
    for(j=0; j<LDA2; j++){
      for(i=0; i<LDA1; i++){
	MAT3(A,i,j,k,LDA1,LDA2)=rallocate_clone(MAT3(B,i,j,k,LDB1,LDB2));
      }
    }
  }
  return A;
}

/**
 @brief rmulti型の3次元配列の複製を添字を指定して生成 A=B(I,J,K)．
*/
rmulti **rmat3_allocate_clone_index(int LDA1, int LDA2, int l, rmulti **B, int LDB1, int LDB2, int *I, int *J, int *K)
{
  int i,j,k;
  rmulti **A=NULL;
  A=(rmulti**)malloc(sizeof(rmulti*)*LDA1*LDA2*l);
  for(k=0; k<l; k++){
    for(j=0; j<LDA2; j++){
      for(i=0; i<LDA1; i++){
	MAT3(A,i,j,k,LDA1,LDA2)=rallocate_clone(MAT3(B,I[i],J[j],K[k],LDB1,LDB2));
      }
    }
  }
  return A;
}

/**
 @brief rmulti型の3次元配列A(LDA1,LDA2,l)の終了処理.
*/
rmulti **rmat3_free(int LDA1, int LDA2, int l, rmulti **A)
{
  int i,j,k;
  if(A==NULL){ return NULL; }
  for(k=0; k<l; k++){
    for(j=0; j<LDA2; j++){
      for(i=0; i<LDA1; i++){
	MAT3(A,i,j,k,LDA1,LDA2)=rfree(MAT3(A,i,j,k,LDA1,LDA2));
      }
    }
  }
  free(A);
  return A=NULL;
}


/**
 @brief 初期化済みのrmulti型の3次元配列の浮動小数点数の精度(ビット数)を変更し再初期化.
*/
int rmat3_round(int m, int n, int l, rmulti **A, int LDA1, int LDA2, int prec)
{
  int i,j,k,e=0;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	e+=rround(MAT3(A,i,j,k,LDA1,LDA2),prec);
      }
    }
  }
  return e;
}


/**
 @brief rmulti型の3次元配列を複成 B=A.
*/
int rmat3_clone(int m, int n, int l, rmulti **B, int LDB1, int LDB2, rmulti **A, int LDA1, int LDA2)
{
  int i,j,k,e=0;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	e+=rclone(MAT3(B,i,j,k,LDB1,LDB2),MAT3(A,i,j,k,LDA1,LDA2));
      }
    }
  }
  return e;
}

/**
 @brief rmulti型の3次元配列を添字を指定して複成 B=A(I,J,K).
*/
int rmat3_clone_index(int m, int n, int l, rmulti **B, int LDB1, int LDB2, rmulti **A, int LDA1, int LDA2, int *I, int *J, int *K)
{
  int i,j,k,e=0;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	e+=rclone(MAT3(B,i,j,k,LDB1,LDB2),MAT3(A,I[i],J[j],K[k],LDA1,LDA2));
      }
    }
  }
  return e;
}

/**
 @brief rmulti型の3次元配列の値の交換 A<=>B.
*/
void rmat3_swap(int m, int n, int l, rmulti **A, int LDA1, int LDA2, rmulti **B, int LDB1, int LDB2)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	rswap(MAT3(A,i,j,k,LDA1,LDA2),MAT3(B,i,j,k,LDB1,LDB2));
      }
    }
  }
}

/** @} */


/////////////////////////////////////////////////////

/** @name rmulti型の3次元配列のメンバ変数に関する関数 */
/** @{ */


/**
 @brief rmulti型の3次元配列の浮動小数点数の精度(ビット数)を取得.
*/
void rmat3_get_prec(int m, int n, int l, int *P, int LDP1, int LDP2, rmulti **A, int LDA1, int LDA2)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	MAT3(P,i,j,k,LDP1,LDP2)=rget_prec(MAT3(A,i,j,k,LDA1,LDA2));
      }
    }
  }
}

/**
 @brief rmulti型行列の精度(ビット数)の最大値の取得.
*/
int rmat3_get_prec_max(int m, int n, int l, rmulti **A, int LDA1, int LDA2)
{
  int i,j,k,p,value=0;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	p=rget_prec(MAT3(A,i,j,k,LDA1,LDA2));
	if(p>value){ value=p; }
      }
    }
  }
  return value;
}

/**
 @brief rmulti型の3次元配列の浮動小数点数の指数部を取得.
*/
void rmat3_get_exp(int m, int n, int l, int *P, int LDP1, int LDP2, rmulti **A, int LDA1, int LDA2)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	MAT3(P,i,j,k,LDP1,LDP2)=rget_exp(MAT3(A,i,j,k,LDA1,LDA2));
      }
    }
  }
}

/**
 @brief rmulti型の3次元配列の浮動小数点数の符号部を取得.
*/
void rmat3_get_sgn(int m, int n, int l, int *P, int LDP1, int LDP2, rmulti **A, int LDA1, int LDA2)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	MAT3(P,i,j,k,LDP1,LDP2)=rget_sgn(MAT3(A,i,j,k,LDA1,LDA2));
      }
    }
  }
}


/** @} */


////////////////////////////////////////////////////////////////////////

/** @name rmulti型の3次元配列の入出力に関する関数 */
/** @{ */

/**
 @brief rmulti型の3次元配列の表示.
*/
void rmat3_print(int m, int n, int l, rmulti **A, int LDA1, int LDA2, const char *name, const char *f, int digits)
{
  int i,j,k,p;
  char format[128];
  if(STR_EQ(f,"f") || STR_EQ(f,"F")){ p=3; }
  else if((strcmp(f,"e")==0) || (strcmp(f,"E")==0)){ p=7; }
  else{ p=3; }
  sprintf(format,"%%%d.%dR%s ",digits+p,digits,f);
  if(A==NULL){
    if(name!=NULL){ printf("%s=NULL\n",name); }
    else          { printf("NULL\n"); }
    return;
  }
  if(name!=NULL){ printf("%s\n",name); }
  if(A==NULL) return;
  for(k=0; k<l; k++){
    for(i=0; i<m; i++){
      for(j=0; j<n; j++){
	mpfr_printf(format,MAT3(A,i,j,k,LDA1,LDA2));
      }
      printf("\n");
    }
    printf("\n");
  }
}

/** @} */

/////////////////////////////////////////////////////////

/** @name rmulti型の3次元配列の値の設定に関する関数 */

/**
 @brief rmulti型の3次元配列の値を倍精度実数から設定 B=rmulti(A).
*/
int rmat3_set_d(int m, int n, int l, rmulti **B, int LDB1, int LDB2, double *A, int LDA1, int LDA2)
{
  int i,j,k,e=0;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	e+=rset_d(MAT3(B,i,j,k,LDB1,LDB2),MAT3(A,i,j,k,LDA1,LDA2));
      }
    }
  }
  return e;
}

/**
 @brief rmulti型の3次元配列の値を文字列型から設定 B=rmulti(A).
*/
void rmat3_set_s(int m, int n, int l, rmulti **B, int LDB1, int LDB2, char **A, int LDA1, int LDA2)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	rset_s(MAT3(B,i,j,k,LDB1,LDB2),MAT3(A,i,j,k,LDA1,LDA2));
      }
    }
  }
}

/**
 @brief rmulti型の3次元配列の値を倍精度浮動小数点数から設定.
*/
int rmat3_set_all_d(int m, int n, int l, rmulti **A, int LDA1, int LDA2, double a)
{
  int i,j,k,e=0;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	e+=rset_d(MAT3(A,i,j,k,LDA1,LDA2),a);
      }
    }
  }
  return e;
}

/**
 @brief rmulti型の3次元配列の値を零に設定.
*/
int rmat3_set_zeros(int m, int n, int l, rmulti **A, int LDA1, int LDA2)
{
  return rmat3_set_all_d(m,n,l,A,LDA1,LDA2,0);
}

/**
 @brief rmulti型の3次元配列の値を1に設定.
*/
int rmat3_set_ones(int m, int n, int l, rmulti **A, int LDA1, int LDA2)
{
  return rmat3_set_all_d(m,n,l,A,LDA1,LDA2,1);
}

/**
 @brief rmulti型の3次元配列の値をnanに設定.
*/
int rmat3_set_nan(int m, int n, int l, rmulti **A, int LDA1, int LDA2)
{
  int i,j,k,e=0;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	rset_nan(MAT3(A,i,j,k,LDA1,LDA2));
      }
    }
  }
  return e;
}

/**
 @brief rmulti型の3次元配列の値をinfに設定.
*/
int rmat3_set_inf(int m, int n, int l, rmulti **A, int LDA1, int LDA2)
{
  int i,j,k,e=0;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	rset_inf(MAT3(A,i,j,k,LDA1,LDA2),1);
      }
    }
  }
  return e;
}

/**
 @brief rmulti型の3次元配列の値を区間(b,a+b)の疑似乱数値を設定.
*/
void rmat3_set_rand(int m, int n, int l, rmulti **A, int LDA1, int LDA2, double a, double b)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	rset_rand(MAT3(A,i,j,k,LDA1,LDA2));
	rmul_d(MAT3(A,i,j,k,LDA1,LDA2),MAT3(A,i,j,k,LDA1,LDA2),a);
	radd_d(MAT3(A,i,j,k,LDA1,LDA2),MAT3(A,i,j,k,LDA1,LDA2),b);
      }
    }
  }
}


/** @} */


////////////////////////////////////////////////////////////////////////


/** @name rmulti型の3次元配列の型変換に関する関数 */
/** @{ */

/**
 @brief rmulti型の3次元配列を倍精度型に変換 B=double(A).
 */
void rmat3_get_d(int m, int n, int l, double *B, int LDB1, int LDB2, rmulti **A, int LDA1, int LDA2)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	MAT3(B,i,j,k,LDB1,LDB2)=rget_d(MAT3(A,i,j,k,LDA1,LDA2));
      }
    }
  }
  return;
}



/** @} */

////////////////////////////////////////////////////////////////////////

/** @name rmulti型の3次元配列の要素の配置に関する関数 */
/** @{ */

/** @} */

////////////////////////////////////////////////////////////////////////

/** @name rmulti型の3次元配列の自動精度調整モードが機能する関数 */
/** @{ */

/**
 @brief rmulti型の3次元配列の値のコピー B=A
*/
int rmat3_copy(int m, int n, int l, rmulti **B, int LDB1, int LDB2, rmulti **A, int LDA1, int LDA2)
{
  int i,j,k,e=0;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	e+=rcopy(MAT3(B,i,j,k,LDB1,LDB2),MAT3(A,i,j,k,LDA1,LDA2));
      }
    }
  }
  return e;
}

/**
 @brief rmulti型の3次元配列を添字を指定してコピー B=A(I,J,K).
*/
int rmat3_copy_index(int m, int n, int l, rmulti **B, int LDB1, int LDB2, rmulti **A, int LDA1, int LDA2, int *I, int *J, int *K)
{
  int i,j,k,e=0;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	e+=rcopy(MAT3(B,i,j,k,LDB1,LDB2),MAT3(A,I[i],J[j],K[k],LDA1,LDA2));
      }
    }
  }
  return e;
}

/**
 @brief rmulti型の3次元配列を添字を指定してコピー B(I,J,K)=A.
*/
int rmat3_index_copy(int m, int n, int l, rmulti **B, int LDB1, int LDB2, rmulti **A, int LDA1, int LDA2, int *I, int *J, int *K)
{
  int i,j,k,e=0;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	e+=rcopy(MAT3(B,I[i],J[j],K[k],LDB1,LDB2),MAT3(A,i,j,k,LDA1,LDA2));
      }
    }
  }
  return e;
}

/**
 @brief rmulti型の3次元配列の符号反転 B=-A.
*/
int rmat3_neg(int m, int n, int l, rmulti **B, int LDB1, int LDB2, rmulti **A, int LDA1, int LDA2)
{
  int i,j,k,e=0;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	e+=rneg(MAT3(B,i,j,k,LDB1,LDB2),MAT3(A,i,j,k,LDA1,LDA2));
      }
    }
  }
  return e;
}

/**
 @brief rmulti型の3次元配列の足し算 C=A+B.
*/
int rmat3_add(int m, int n, int l, rmulti **C, int LDC1, int LDC2, rmulti **A, int LDA1, int LDA2, rmulti **B, int LDB1, int LDB2)
{
  int i,j,k,e=0;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	e+=radd(MAT3(C,i,j,k,LDC1,LDC2),MAT3(A,i,j,k,LDA1,LDA2),MAT3(B,i,j,k,LDB1,LDB2));
      }
    }
  }
  return e;
}

/**
 @brief rmulti型の3次元配列の足し算 C=A+b.
*/
int rmat3_add_r(int m, int n, int l, rmulti **C, int LDC1, int LDC2, rmulti **A, int LDA1, int LDA2, rmulti *b)
{
  int i,j,k,e=0;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	e+=radd(MAT3(C,i,j,k,LDC1,LDC2),MAT3(A,i,j,k,LDA1,LDA2),b);
      }
    }
  }
  return e;
}

/**
 @brief rmulti型の3次元配列の引き算 C=A-B.
*/
int rmat3_sub(int m, int n, int l, rmulti **C, int LDC1, int LDC2, rmulti **A, int LDA1, int LDA2, rmulti **B, int LDB1, int LDB2)
{
  int i,j,k,e=0;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	e+=rsub(MAT3(C,i,j,k,LDC1,LDC2),MAT3(A,i,j,k,LDA1,LDA2),MAT3(B,i,j,k,LDB1,LDB2));
      }
    }
  }
  return e;
}

/**
 @brief rmulti型の3次元配列の引き算 C=a-B.
*/
int rmat3_sub_r1(int m, int n, int l, rmulti **C, int LDC1, int LDC2, rmulti *a, rmulti **B, int LDB1, int LDB2)
{
  int i,j,k,e=0;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	e+=rsub(MAT3(C,i,j,k,LDC1,LDC2),a,MAT3(B,i,j,k,LDB1,LDB2));
      }
    }
  }
  return e;
}

/**
 @brief rmulti型の3次元配列の引き算 C=A-b.
*/
int rmat3_sub_r2(int m, int n, int l, rmulti **C, int LDC1, int LDC2, rmulti **A, int LDA1, int LDA2, rmulti *b)
{
  int i,j,k,e=0;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	e+=rsub(MAT3(C,i,j,k,LDC1,LDC2),MAT3(A,i,j,k,LDA1,LDA2),b);
      }
    }
  }
  return e;
}

//追加

/**
 @brief rmulti型の3次元配列の引き算 C=A-b.
*/
int rmat3_sub_d2(int m, int n, int l, rmulti **C, int LDC1, int LDC2, rmulti **A, int LDA1, int LDA2, double b)
{
  int i,j,k,e=0;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	e+=rsub_d2(MAT3(C,i,j,k,LDC1,LDC2),MAT3(A,i,j,k,LDA1,LDA2),b);
      }
    }
  }
  return e;
}

//ここまで

/**
 @brief rmulti型の3次元配列の掛け算 C=A.*B.
*/
int rmat3_mul(int m, int n, int l, rmulti **C, int LDC1, int LDC2, rmulti **A, int LDA1, int LDA2, rmulti **B, int LDB1, int LDB2)
{
  int i,j,k,e=0;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	e+=rmul(MAT3(C,i,j,k,LDC1,LDC2),MAT3(A,i,j,k,LDA1,LDA2),MAT3(B,i,j,k,LDB1,LDB2));
      }
    }
  }
  return e;
}

/**
 @brief rmulti型の3次元配列の掛け算 C=A.*b.
*/
int rmat3_mul_r(int m, int n, int l, rmulti **C, int LDC1, int LDC2, rmulti **A, int LDA1, int LDA2, rmulti *b)
{
  int i,j,k,e=0;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	e+=rmul(MAT3(C,i,j,k,LDC1,LDC2),MAT3(A,i,j,k,LDA1,LDA2),b);
      }
    }
  }
  return e;
}

/** @} */

////////////////////////////////////////////////////////////////////////

/** @name rmulti型の3次元配列の数学関数に関する関数 */
/** @{ */

//追加

/**
 @brief rmulti型の3次元配列の整数への丸め C=floor(A).
*/
int rmat3_floor(int m, int n, int l, rmulti **C, int LDC1, int LDC2, rmulti **A, int LDA1, int LDA2)
{
  int i,j,k,e=0;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	e+=rfloor(MAT3(C,i,j,k,LDC1,LDC2),MAT3(A,i,j,k,LDA1,LDA2));
      }
    }
  }
  return e;
}

//ここまで

/**
 @brief rmulti型の3次元配列の絶対値 C=abs(A).
*/
int rmat3_abs(int m, int n, int l, rmulti **C, int LDC1, int LDC2, rmulti **A, int LDA1, int LDA2)
{
  int i,j,k,e=0;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	e+=rabs(MAT3(C,i,j,k,LDC1,LDC2),MAT3(A,i,j,k,LDA1,LDA2));
      }
    }
  }
  return e;
}

//編集済み

/**
 @brief rmulti型の3次元配列のlog10を取った値 C=log10(A).
*/
int rmat3_log10(int m, int n, int l, rmulti **C, int LDC1, int LDC2, rmulti **A, int LDA1, int LDA2)
{
  int i,j,k,e=0;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	e+=rlog10(MAT3(C,i,j,k,LDC1,LDC2),MAT3(A,i,j,k,LDA1,LDA2));
      }
    }
  }
  return e;
}

/**
 @brief 3次元配列の平方根 C=sqrt(A).
*/
int rmat3_sqrt(int m, int n, int l, rmulti **C, int LDC1, int LDC2, rmulti **A, int LDA1, int LDA2)
{
  int i,j,k,e=0;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	e+=rsqrt(MAT3(C,i,j,k,LDC1,LDC2),MAT3(A,i,j,k,LDA1,LDA2));
      }
    }
  }
  return e;
}

/**
 @brief rmulti型の3次元配列の和 B=sum(A).
*/
int rmat3_sum(int m, int n, int l, rmulti **B, int LDB1, int LDB2, rmulti **A, int LDA1, int LDA2)
{
  int k,e=0;
  for(k=0; k<l; k++){
    e+=rvec_sum_rmat(m,n,&MAT3(B,0,0,k,LDB1,LDB2),&MAT3(A,0,0,k,LDA1,LDA2),LDA1);
  }
  return e;
}

//ここまで

/**
 @brief rmulti型の3次元配列の最大値 B=max(A).
*/
int rmat3_max(int m, int n, int l, rmulti **B, int LDB1, int LDB2, rmulti **A, int LDA1, int LDA2)
{
  int k,e=0;
  for(k=0; k<l; k++){
    e+=rvec_max_rmat(m,n,&MAT3(B,0,0,k,LDB1,LDB2),&MAT3(A,0,0,k,LDA1,LDA2),LDA1);
  }
  return e;
}

//追加

/**
 @brief rmulti型の3次元配列A,Bの最大要素 C=max(A,B).
*/
int rmat3_max2(int m, int n, int l, rmulti **C, int LDC1, int LDC2, rmulti **A, int LDA1, int LDA2, rmulti **B, int LDB1, int LDB2)
{
  int i,j,k,e=0;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	if(rgt(MAT3(A,i,j,k,LDA1,LDA2),MAT3(B,i,j,k,LDB1,LDB2))){
	  //A>B
	  e+=rcopy(MAT3(C,i,j,k,LDC1,LDC2),MAT3(A,i,j,k,LDA1,LDA2));
	}else{
	  //A<=B
	  e+=rcopy(MAT3(C,i,j,k,LDC1,LDC2),MAT3(B,i,j,k,LDB1,LDB2));
	}
      }
    }
  }
  return e;
}
//ここまで

/**
 @brief rmulti型の3次元配列の最小値 B=min(A).
*/
int rmat3_min(int m, int n, int l, rmulti **B, int LDB1, int LDB2, rmulti **A, int LDA1, int LDA2)
{
  int k,e=0;
  for(k=0; k<l; k++){
    e+=rvec_min_rmat(m,n,&MAT3(B,0,0,k,LDB1,LDB2),&MAT3(A,0,0,k,LDA1,LDA2),LDA1);
  }
  return e;
}

/**
 @brief rmulti型の3次元配列の割り算 C=A./B.
*/
int rmat3_div(int m, int n, int l, rmulti **C, int LDC1, int LDC2, rmulti **A, int LDA1, int LDA2, rmulti **B, int LDB1, int LDB2)
{
  int i,j,k,e=0;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	e+=rdiv(MAT3(C,i,j,k,LDC1,LDC2),MAT3(A,i,j,k,LDA1,LDA2),MAT3(B,i,j,k,LDB1,LDB2));
      }
    }
  }
  return e;
}

/**
 @brief rmulti型の3次元配列の割り算 C=a./B.
*/
int rmat3_div_r1(int m, int n, int l, rmulti **C, int LDC1, int LDC2, rmulti *a, rmulti **B, int LDB1, int LDB2)
{
  int i,j,k,e=0;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	e+=rdiv(MAT3(C,i,j,k,LDC1,LDC2),a,MAT3(B,i,j,k,LDB1,LDB2));
      }
    }
  }
  return e;
}

/**
 @brief rmulti型の3次元配列の割り算 C=A./b.
*/
int rmat3_div_r2(int m, int n, int l, rmulti **C, int LDC1, int LDC2, rmulti **A, int LDA1, int LDA2, rmulti *b)
{
  int i,j,k,e=0;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	e+=rdiv(MAT3(C,i,j,k,LDC1,LDC2),MAT3(A,i,j,k,LDA1,LDA2),b);
      }
    }
  }
  return e;
}

/**
 @brief rmulti型の3次元配列のべき乗 C=A.^B.
*/
int rmat3_pow(int m, int n, int l, rmulti **C, int LDC1, int LDC2, rmulti **A, int LDA1, int LDA2, rmulti **B, int LDB1, int LDB2)
{
  int i,j,k,e=0;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	e+=rpow(MAT3(C,i,j,k,LDC1,LDC2),MAT3(A,i,j,k,LDA1,LDA2),MAT3(B,i,j,k,LDB1,LDB2));
      }
    }
  }
  return e;
}

/**
 @brief rmulti型の3次元配列のべき乗 C=a.^B.
*/
int rmat3_pow_r1(int m, int n, int l, rmulti **C, int LDC1, int LDC2, rmulti *a, rmulti **B, int LDB1, int LDB2)
{
  int i,j,k,e=0;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	e+=rpow(MAT3(C,i,j,k,LDC1,LDC2),a,MAT3(B,i,j,k,LDB1,LDB2));
      }
    }
  }
  return e;
}

/**
 @brief rmulti型の3次元配列のべき乗 C=A.^b.
*/
int rmat3_pow_r2(int m, int n, int l, rmulti **C, int LDC1, int LDC2, rmulti **A, int LDA1, int LDA2, rmulti *b)
{
  int i,j,k,e=0;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	e+=rpow(MAT3(C,i,j,k,LDC1,LDC2),MAT3(A,i,j,k,LDA1,LDA2),b);
      }
    }
  }
  return e;
}

//追加

/**
 @brief rmulti型の3次元配列の指数部で評価 C=get_exp10(x,offset)
*/
int rmat3_get_exp10(int m, int n, int l, rmulti **C, int LDC1, int LDC2, rmulti **A, int LDA1, int LDA2, double b)
{
  int i,j,k,e=0;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	e+=rget_exp10(MAT3(C,i,j,k,LDC1,LDC2),MAT3(A,i,j,k,LDA1,LDA2),b);
      }
    }
  }
  return e;
}

/**
 @brief rmulti型の3次元配列の指数部で評価 C=get_exp2(x,offset)
*/
int rmat3_get_exp2(int m, int n, int l, rmulti **C, int LDC1, int LDC2, rmulti **A, int LDA1, int LDA2, double b)
{
  int i,j,k,e=0;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	e+=rget_exp2(MAT3(C,i,j,k,LDC1,LDC2),MAT3(A,i,j,k,LDA1,LDA2),b);
      }
    }
  }
  return e;
}

//ここまで
/** @} */


////////////////////////////////////////////////////////////////////////

/** @name rmulti型の3次元配列の要素の比較に関する関数 */
/** @{ */

/**
 @brief rmulti型の3次元配列の比較 B=isnan(A).
 */
void rmat3_isnan(int m, int n, int l, int *B, int LDB1, int LDB2, rmulti **A, int LDA1, int LDA2)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	MAT3(B,i,j,k,LDB1,LDB2)=ris_nan(MAT3(A,i,j,k,LDA1,LDA2));
      }
    }
  }
  return;
}

/**
 @brief rmulti型の3次元配列の比較 B=isinf(A).
*/
void rmat3_isinf(int m, int n, int l, int *B, int LDB1, int LDB2, rmulti **A, int LDA1, int LDA2)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	MAT3(B,i,j,k,LDB1,LDB2)=ris_inf(MAT3(A,i,j,k,LDA1,LDA2));
      }
    }
  }
  return;
}

/**
 @brief rmulti型の3次元配列の比較 C=(A==B).
*/
void rmat3_eq(int m, int n, int l, int *C, int LDC1, int LDC2, rmulti **A, int LDA1, int LDA2, rmulti **B, int LDB1, int LDB2)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	MAT3(C,i,j,k,LDC1,LDC2)=req(MAT3(A,i,j,k,LDA1,LDA2),MAT3(B,i,j,k,LDB1,LDB2));
      }
    }
  }
  return;
}

/**
 @brief rmulti型の3次元配列の比較 C=(A==b).
*/
void rmat3_eq_r(int m, int n, int l, int *C, int LDC1, int LDC2, rmulti **A, int LDA1, int LDA2, rmulti *b)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	MAT3(C,i,j,k,LDC1,LDC2)=req(MAT3(A,i,j,k,LDA1,LDA2),b);
      }
    }
  }
  return;
}

/**
 @brief rmulti型の3次元配列の比較 C=(A!=B).
*/
void rmat3_ne(int m, int n, int l, int *C, int LDC1, int LDC2, rmulti **A, int LDA1, int LDA2, rmulti **B, int LDB1, int LDB2)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	MAT3(C,i,j,k,LDC1,LDC2)=rne(MAT3(A,i,j,k,LDA1,LDA2),MAT3(B,i,j,k,LDB1,LDB2));
      }
    }
  }
  return;
}

/**
 @brief rmulti型の3次元配列の比較 C=(A!=b).
*/
void rmat3_ne_r(int m, int n, int l, int *C, int LDC1, int LDC2, rmulti **A, int LDA1, int LDA2, rmulti *b)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	MAT3(C,i,j,k,LDC1,LDC2)=rne(MAT3(A,i,j,k,LDA1,LDA2),b);
      }
    }
  }
  return;
}

/**
 @brief rmulti型の3次元配列の比較 C=(A>=B).
*/
void rmat3_ge(int m, int n, int l, int *C, int LDC1, int LDC2, rmulti **A, int LDA1, int LDA2, rmulti **B, int LDB1, int LDB2)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	MAT3(C,i,j,k,LDC1,LDC2)=rge(MAT3(A,i,j,k,LDA1,LDA2),MAT3(B,i,j,k,LDB1,LDB2));
      }
    }
  }
  return;
}

/**
 @brief rmulti型の3次元配列の比較 C=(A>=b).
*/
void rmat3_ge_r(int m, int n, int l, int *C, int LDC1, int LDC2, rmulti **A, int LDA1, int LDA2, rmulti *b)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	MAT3(C,i,j,k,LDC1,LDC2)=rge(MAT3(A,i,j,k,LDA1,LDA2),b);
      }
    }
  }
  return;
}

/**
 @brief rmulti型の3次元配列の比較 C=(A>B).
*/
void rmat3_gt(int m, int n, int l, int *C, int LDC1, int LDC2, rmulti **A, int LDA1, int LDA2, rmulti **B, int LDB1, int LDB2)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	MAT3(C,i,j,k,LDC1,LDC2)=rgt(MAT3(A,i,j,k,LDA1,LDA2),MAT3(B,i,j,k,LDB1,LDB2));
      }
    }
  }
  return;
}

/**
 @brief rmulti型の3次元配列の比較 C=(A>b).
*/
void rmat3_gt_r(int m, int n, int l, int *C, int LDC1, int LDC2, rmulti **A, int LDA1, int LDA2, rmulti *b)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	MAT3(C,i,j,k,LDC1,LDC2)=rgt(MAT3(A,i,j,k,LDA1,LDA2),b);
      }
    }
  }
  return;
}

/**
 @brief rmulti型の3次元配列の比較 C=(A<=B).
*/
void rmat3_le(int m, int n, int l, int *C, int LDC1, int LDC2, rmulti **A, int LDA1, int LDA2, rmulti **B, int LDB1, int LDB2)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	MAT3(C,i,j,k,LDC1,LDC2)=rle(MAT3(A,i,j,k,LDA1,LDA2),MAT3(B,i,j,k,LDB1,LDB2));
      }
    }
  }
  return;
}

/**
 @brief rmulti型の3次元配列の比較 C=(A<=b).
*/
void rmat3_le_r(int m, int n, int l, int *C, int LDC1, int LDC2, rmulti **A, int LDA1, int LDA2, rmulti *b)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	MAT3(C,i,j,k,LDC1,LDC2)=rle(MAT3(A,i,j,k,LDA1,LDA2),b);
      }
    }
  }
  return;
}

/**
 @brief rmulti型の3次元配列の比較 C=(A<B).
*/
void rmat3_lt(int m, int n, int l, int *C, int LDC1, int LDC2, rmulti **A, int LDA1, int LDA2, rmulti **B, int LDB1, int LDB2)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	MAT3(C,i,j,k,LDC1,LDC2)=rlt(MAT3(A,i,j,k,LDA1,LDA2),MAT3(B,i,j,k,LDB1,LDB2));
      }
    }
  }
  return;
}

/**
 @brief rmulti型の3次元配列の比較 C=(A<b).
*/
void rmat3_lt_r(int m, int n, int l, int *C, int LDC1, int LDC2, rmulti **A, int LDA1, int LDA2, rmulti *b)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	MAT3(C,i,j,k,LDC1,LDC2)=rlt(MAT3(A,i,j,k,LDA1,LDA2),b);
      }
    }
  }
  return;
}

/** @} */

////////////////////////////////////////////////////////////////////

/** @name rmulti型の写像に関する関数 */
/** @{ */

/** @} */

//EOF
