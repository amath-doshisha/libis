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
 @file  cmat3.c
 @brief 多倍長精度複素数型cmultiの3次元配列に関する関数の定義
 @details cmulti型の関数に関する定義は@link cmulti.c@endlinkを参照のこと.
          cmulti型のベクトルの関数に関する定義は@link cvec.c@endlinkを参照のこと.
          cmulti型の行列の関数に関する定義は@link cmat.c@endlinkを参照のこと.
 */

/////////////////////////////////////////////////////

#define FILE_NAME_LENGTH_MAX 100000

/////////////////////////////////////////////////////


/** @name icmulti型の3次元配列値の設定に関する関数 */
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
 @brief cmulti型の3次元配列の複製を添字を指定して生成 A=B(I,J,K)．
*/
cmulti **cmat3_allocate_clone_index(int LDA1, int LDA2, int l, cmulti **B, int LDB1, int LDB2, int *I, int *J, int *K)
{
  int i,j,k;
  cmulti **A=NULL;
  A=(cmulti**)malloc(sizeof(cmulti*)*LDA1*LDA2*l);
  for(k=0; k<l; k++){
    for(j=0; j<LDA2; j++){
      for(i=0; i<LDA1; i++){
	MAT3(A,i,j,k,LDA1,LDA2)=callocate_clone(MAT3(B,I[i],J[j],K[k],LDB1,LDB2));
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
	MAT3(A,i,j,k,LDA1,LDA2)=cmfree(MAT3(A,i,j,k,LDA1,LDA2));
      }
    }
  }
  free(A);
  return A=NULL;
}

/**
 @brief cmulti型の3次元配列の値の複製 B=A
*/
void cmat3_clone(int m, int n, int l, cmulti **B, int LDB1, int LDB2, cmulti **A, int LDA1, int LDA2)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	cclone_c(MAT3(B,i,j,k,LDB1,LDB2),MAT3(A,i,j,k,LDA1,LDA2));
      }
    }
  }
}

/**
 @brief cmulti型の3次元配列の値を添字を指定して複製 B=A(I,J,K)
*/
void cmat3_clone_index(int m, int n, int l, cmulti **B, int LDB1, int LDB2, cmulti **A, int LDA1, int LDA2, int *I, int *J, int *K)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	cclone_c(MAT3(B,i,j,k,LDB1,LDB2),MAT3(A,I[i],J[j],K[k],LDA1,LDA2));
      }
    }
  }
}

/** @} */

/////////////////////////////////////////////////////

/** @name cmulti型の3次元配列のメンバ変数に関する関数 */
/** @{ */

/**
 @brief cmulti型の3次元配列の浮動小数点数の精度(ビット数)を取得.
*/
void cmat3_get_prec(int m, int n, int l, int *P, int LDP1, int LDP2, cmulti **A, int LDA1, int LDA2)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	MAT3(P,i,j,k,LDP1,LDP2)=cget_prec(MAT3(A,i,j,k,LDA1,LDA2));
      }
    }
  }
}

int cmat3_get_prec_max(int m, int n, int l, cmulti **A, int LDA1, int LDA2)
{
  int value,p,i=0,j=0,k=0;
  value=cget_prec(MAT3(A,i,j,k,LDA1,LDA2));
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	p=cget_prec(MAT3(A,i,j,k,LDA1,LDA2));
	if(p>value){ value=p; }
      }
    }
  }
  return value;
}

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
 @brief cmulti型の3次元配列の値を倍精度実数から設定 B=cmulti(A).
*/
void cmat3_set_z(int m, int n, int l, cmulti **B, int LDB1, int LDB2, dcomplex *A, int LDA1, int LDA2)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	cset_z(MAT3(B,i,j,k,LDB1,LDB2),MAT3(A,i,j,k,LDA1,LDA2));
      }
    }
  }
}

/**
 @brief cmulti型の3次元配列の値を倍精度実数から設定 B=cmulti(A).
*/
void cmat3_set_d(int m, int n, int l, cmulti **B, int LDB1, int LDB2, double *A, int LDA1, int LDA2)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	cset_d(MAT3(B,i,j,k,LDB1,LDB2),MAT3(A,i,j,k,LDA1,LDA2));
      }
    }
  }
}

/**
 @brief cmulti型の3次元配列の値を文字列型から設定 B=cmulti(A).
*/
void cmat3_set_s(int m, int n, int l, cmulti **B, int LDB1, int LDB2, char **A, int LDA1, int LDA2)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	cset_s(MAT3(B,i,j,k,LDB1,LDB2),MAT3(A,i,j,k,LDA1,LDA2));
      }
    }
  }
}


/**
 @brief cmulti型の3次元配列の値を倍精度浮動小数点数から設定.
*/
void cmat3_set_all_z(int m, int n, int l, cmulti **A, int LDA1, int LDA2, dcomplex a)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	cset_z(MAT3(A,i,j,k,LDA1,LDA2),a);
      }
    }
  }
}

/**
 @brief cmulti型の3次元配列の値を倍精度浮動小数点数から設定.
*/
void cmat3_set_all_dd(int m, int n, int l, cmulti **A, int LDA1, int LDA2, double a_r, double a_i)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	cset_dd(MAT3(A,i,j,k,LDA1,LDA2),a_r,a_i);
      }
    }
  }
}

/**
 @brief cmulti型の3次元配列の値を倍精度浮動小数点数から設定.
*/
void cmat3_set_all_d(int m, int n, int l, cmulti **A, int LDA1, int LDA2, double a)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	cset_d(MAT3(A,i,j,k,LDA1,LDA2),a);
      }
    }
  }
}

/**
 @brief cmulti型の3次元配列の値を零に設定.
*/
void cmat3_set_zeros(int m, int n, int l, cmulti **A, int LDA1, int LDA2)
{
  return cmat3_set_all_d(m,n,l,A,LDA1,LDA2,0);
}

/**
 @brief cmulti型の3次元配列の値を1に設定.
*/
void cmat3_set_ones(int m, int n, int l, cmulti **A, int LDA1, int LDA2)
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
	cmul_cd(MAT3(A,i,j,k,LDA1,LDA2),MAT3(A,i,j,k,LDA1,LDA2),a);
	cadd_cd(MAT3(A,i,j,k,LDA1,LDA2),MAT3(A,i,j,k,LDA1,LDA2),b);
      }
    }
  }
}

/**
 @brief cmulti型の3次元配列の値をnanに設定.
*/
void cmat3_set_nan(int m, int n, int l, cmulti **A, int LDA1, int LDA2)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	cset_nan(MAT3(A,i,j,k,LDA1,LDA2));
      }
    }
  }
}

/**
 @brief cmulti型の3次元配列の値をinfに設定.
*/
void cmat3_set_inf(int m, int n, int l, cmulti **A, int LDA1, int LDA2)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	cset_inf(MAT3(A,i,j,k,LDA1,LDA2),1,1);
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
void cmat3_copy(int m, int n, int l, cmulti **B, int LDB1, int LDB2, cmulti **A, int LDA1, int LDA2)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	cset_c(MAT3(B,i,j,k,LDB1,LDB2),MAT3(A,i,j,k,LDA1,LDA2));
      }
    }
  }
}

/**
 @brief cmulti型の3次元配列の値のコピー B=A
*/
void cmat3_copy_rmat3(int m, int n, int l, cmulti **B, int LDB1, int LDB2, rmulti **A, int LDA1, int LDA2)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	cset_r(MAT3(B,i,j,k,LDB1,LDB2),MAT3(A,i,j,k,LDA1,LDA2));
      }
    }
  }
}

//追加

/**
 @brief cmulti型の3次元配列の値のコピー C=A+Bi
*/
void cmat3_copy_rmat3_rmat3(int m, int n, int l, cmulti **C, int LDC1, int LDC2, rmulti **A, int LDA1, int LDA2, rmulti **B, int LDB1, int LDB2)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	cset_rr(MAT3(C,i,j,k,LDC1,LDC2),MAT3(A,i,j,k,LDA1,LDA2),MAT3(B,i,j,k,LDB1,LDB2));
      }
    }
  }
}

/**
 @brief cmulti型の3次元配列の値のコピー C=a+Bi
*/
void cmat3_copy_r_rmat3(int m, int n, int l, cmulti **C, int LDC1, int LDC2, rmulti *a, rmulti **B, int LDB1, int LDB2)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	cset_rr(MAT3(C,i,j,k,LDC1,LDC2),a,MAT3(B,i,j,k,LDB1,LDB2));
      }
    }
  }
}

/**
 @brief cmulti型の3次元配列の値のコピー C=A+bi
*/
void cmat3_copy_rmat3_r(int m, int n, int l, cmulti **C, int LDC1, int LDC2, rmulti **A, int LDA1, int LDA2, rmulti *b)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	cset_rr(MAT3(C,i,j,k,LDC1,LDC2),MAT3(A,i,j,k,LDA1,LDA2),b);
      }
    }
  }
}
//ここまで

/**
 @brief cmulti型の3次元配列の値を添字を指定してコピー B=A(I,J,K)
*/
void cmat3_copy_index(int m, int n, int l, cmulti **B, int LDB1, int LDB2, cmulti **A, int LDA1, int LDA2, int *I, int *J, int *K)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	cset_c(MAT3(B,i,j,k,LDB1,LDB2),MAT3(A,I[i],J[j],K[k],LDA1,LDA2));
      }
    }
  }
}

/**
 @brief cmulti型の3次元配列の値を添字を指定してコピー B(I,J,K)=A
*/
void cmat3_index_copy(int m, int n, int l, cmulti **B, int LDB1, int LDB2, cmulti **A, int LDA1, int LDA2, int *I, int *J, int *K)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	cset_c(MAT3(B,I[i],J[j],K[k],LDB1,LDB2),MAT3(A,i,j,k,LDA1,LDA2));
      }
    }
  }
}

/**
 @brief cmulti型の3次元配列の値を添字を指定してコピー B(I,J,K)=A
*/
void cmat3_index_copy_rmat3(int m, int n, int l, cmulti **B, int LDB1, int LDB2, rmulti **A, int LDA1, int LDA2, int *I, int *J, int *K)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	cset_r(MAT3(B,I[i],J[j],K[k],LDB1,LDB2),MAT3(A,i,j,k,LDA1,LDA2));
      }
    }
  }
}


/**
 @brief cmulti型の3次元配列の実部 B=real(A)
*/
void cmat3_real(int m, int n, int l, rmulti **B, int LDB1, int LDB2, cmulti **A, int LDA1, int LDA2)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	rset_r(MAT3(B,i,j,k,LDB1,LDB2),C_R(MAT3(A,i,j,k,LDA1,LDA2)));
      }
    }
  }
}

/**
 @brief cmulti型の3次元配列の虚部 B=imag(A)
*/
void cmat3_imag(int m, int n, int l, rmulti **B, int LDB1, int LDB2, cmulti **A, int LDA1, int LDA2)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	rset_r(MAT3(B,i,j,k,LDB1,LDB2),C_I(MAT3(A,i,j,k,LDA1,LDA2)));
      }
    }
  }
}

/**
 @brief cmulti型の3次元配列の共役 B=conj(A)
*/
void cmat3_conj(int m, int n, int l, cmulti **B, int LDB1, int LDB2, cmulti **A, int LDA1, int LDA2)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	cconj_c(MAT3(B,i,j,k,LDB1,LDB2),MAT3(A,i,j,k,LDA1,LDA2));
      }
    }
  }
}

/**
 @brief cmulti型の3次元配列の符号反転 B=-A.
*/
void cmat3_neg(int m, int n, int l, cmulti **B, int LDB1, int LDB2, cmulti **A, int LDA1, int LDA2)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	cneg_c(MAT3(B,i,j,k,LDB1,LDB2),MAT3(A,i,j,k,LDA1,LDA2));
      }
    }
  }
}

/**
 @brief cmulti型の3次元配列の足し算 C=A+B.
*/
void cmat3_add(int m, int n, int l, cmulti **C, int LDC1, int LDC2, cmulti **A, int LDA1, int LDA2, cmulti **B, int LDB1, int LDB2)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	cadd_cc(MAT3(C,i,j,k,LDC1,LDC2),MAT3(A,i,j,k,LDA1,LDA2),MAT3(B,i,j,k,LDB1,LDB2));
      }
    }
  }
}

/**
 @brief cmulti型の3次元配列の足し算 C=A+B.
*/
void cmat3_add_rmat3(int m, int n, int l, cmulti **C, int LDC1, int LDC2, cmulti **A, int LDA1, int LDA2, rmulti **B, int LDB1, int LDB2)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	cadd_cr(MAT3(C,i,j,k,LDC1,LDC2),MAT3(A,i,j,k,LDA1,LDA2),MAT3(B,i,j,k,LDB1,LDB2));
      }
    }
  }
}

/**
 @brief cmulti型の3次元配列の足し算 C=A+b.
*/
void cmat3_add_c(int m, int n, int l, cmulti **C, int LDC1, int LDC2, cmulti **A, int LDA1, int LDA2, cmulti *b)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	cadd_cc(MAT3(C,i,j,k,LDC1,LDC2),MAT3(A,i,j,k,LDA1,LDA2),b);
      }
    }
  }
}

/**
 @brief cmulti型の3次元配列の足し算 C=A+b.
*/
void cmat3_add_r(int m, int n, int l, cmulti **C, int LDC1, int LDC2, cmulti **A, int LDA1, int LDA2, rmulti *b)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	cadd_cr(MAT3(C,i,j,k,LDC1,LDC2),MAT3(A,i,j,k,LDA1,LDA2),b);
      }
    }
  }
}

/**
 @brief cmulti型の3次元配列の足し算 C=A+b.
*/
void rmat3_add_c(int m, int n, int l, cmulti **C, int LDC1, int LDC2, rmulti **A, int LDA1, int LDA2, cmulti *b)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
        cadd_cr(MAT3(C,i,j,k,LDC1,LDC2),b,MAT3(A,i,j,k,LDA1,LDA2));
      }
    }
  }
}

/**
 @brief cmulti型の3次元配列の引き算 C=A-B.
*/
void cmat3_sub(int m, int n, int l, cmulti **C, int LDC1, int LDC2, cmulti **A, int LDA1, int LDA2, cmulti **B, int LDB1, int LDB2)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	csub_cc(MAT3(C,i,j,k,LDC1,LDC2),MAT3(A,i,j,k,LDA1,LDA2),MAT3(B,i,j,k,LDB1,LDB2));
      }
    }
  }
}

/**
 @brief cmulti型の3次元配列の引き算 C=A-B.
*/
void rmat3_sub_cmat3(int m, int n, int l, cmulti **C, int LDC1, int LDC2, rmulti **A, int LDA1, int LDA2, cmulti **B, int LDB1, int LDB2)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	csub_rc(MAT3(C,i,j,k,LDC1,LDC2),MAT3(A,i,j,k,LDA1,LDA2),MAT3(B,i,j,k,LDB1,LDB2));
      }
    }
  }
}

/**
 @brief cmulti型の3次元配列の引き算 C=A-B.
*/
void cmat3_sub_rmat3(int m, int n, int l, cmulti **C, int LDC1, int LDC2, cmulti **A, int LDA1, int LDA2, rmulti **B, int LDB1, int LDB2)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	csub_cr(MAT3(C,i,j,k,LDC1,LDC2),MAT3(A,i,j,k,LDA1,LDA2),MAT3(B,i,j,k,LDB1,LDB2));
      }
    }
  }
}

/**
 @brief cmulti型の3次元配列の引き算 C=a-B.
*/
void cmat3_sub_c1(int m, int n, int l, cmulti **C, int LDC1, int LDC2, cmulti *a, cmulti **B, int LDB1, int LDB2)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	csub_cc(MAT3(C,i,j,k,LDC1,LDC2),a,MAT3(B,i,j,k,LDB1,LDB2));
      }
    }
  }
}

/**
 @brief cmulti型の3次元配列の引き算 C=A-b.
*/
void cmat3_sub_c2(int m, int n, int l, cmulti **C, int LDC1, int LDC2, cmulti **A, int LDA1, int LDA2, cmulti *b)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	csub_cc(MAT3(C,i,j,k,LDC1,LDC2),MAT3(A,i,j,k,LDA1,LDA2),b);
      }
    }
  }
}

/**
 @brief cmulti型の3次元配列の引き算 C=a-B.
*/
void cmat3_sub_r1(int m, int n, int l, cmulti **C, int LDC1, int LDC2, rmulti *a, cmulti **B, int LDB1, int LDB2)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	csub_rc(MAT3(C,i,j,k,LDC1,LDC2),a,MAT3(B,i,j,k,LDB1,LDB2));
      }
    }
  }
}

/**
 @brief cmulti型の3次元配列の引き算 C=A-b.
*/
void cmat3_sub_r2(int m, int n, int l, cmulti **C, int LDC1, int LDC2, cmulti **A, int LDA1, int LDA2, rmulti *b)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	csub_cr(MAT3(C,i,j,k,LDC1,LDC2),MAT3(A,i,j,k,LDA1,LDA2),b);
      }
    }
  }
}

/**
 @brief cmulti型の3次元配列の引き算 C=a-B.
*/
void rmat3_sub_c1(int m, int n, int l, cmulti **C, int LDC1, int LDC2, cmulti *a, rmulti **B, int LDB1, int LDB2)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	csub_cr(MAT3(C,i,j,k,LDC1,LDC2),a,MAT3(B,i,j,k,LDB1,LDB2));
      }
    }
  }
}

/**
 @brief cmulti型の3次元配列の引き算 C=A-b.
*/
void rmat3_sub_c2(int m, int n, int l, cmulti **C, int LDC1, int LDC2, rmulti **A, int LDA1, int LDA2, cmulti *b)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	csub_rc(MAT3(C,i,j,k,LDC1,LDC2),MAT3(A,i,j,k,LDA1,LDA2),b);
      }
    }
  }
}

/**
 @brief cmulti型の3次元配列の掛け算 C=A.*B.
*/
void cmat3_mul(int m, int n, int l, cmulti **C, int LDC1, int LDC2, cmulti **A, int LDA1, int LDA2, cmulti **B, int LDB1, int LDB2)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	cmul_cc(MAT3(C,i,j,k,LDC1,LDC2),MAT3(A,i,j,k,LDA1,LDA2),MAT3(B,i,j,k,LDB1,LDB2));
      }
    }
  }
}

/**
 @brief cmulti型の3次元配列の掛け算 C=A.*B.
*/
void cmat3_mul_rmat3(int m, int n, int l, cmulti **C, int LDC1, int LDC2, cmulti **A, int LDA1, int LDA2, rmulti **B, int LDB1, int LDB2)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	cmul_cr(MAT3(C,i,j,k,LDC1,LDC2),MAT3(A,i,j,k,LDA1,LDA2),MAT3(B,i,j,k,LDB1,LDB2));
      }
    }
  }
}

/**
 @brief cmulti型の3次元配列の掛け算 C=A.*b.
*/
void cmat3_mul_c(int m, int n, int l, cmulti **C, int LDC1, int LDC2, cmulti **A, int LDA1, int LDA2, cmulti *b)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	cmul_cc(MAT3(C,i,j,k,LDC1,LDC2),MAT3(A,i,j,k,LDA1,LDA2),b);
      }
    }
  }
}

/**
 @brief cmulti型の3次元配列の掛け算 C=A.*b.
*/
void cmat3_mul_r(int m, int n, int l, cmulti **C, int LDC1, int LDC2, cmulti **A, int LDA1, int LDA2, rmulti *b)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	cmul_cr(MAT3(C,i,j,k,LDC1,LDC2),MAT3(A,i,j,k,LDA1,LDA2),b);
      }
    }
  }
}

/**
 @brief cmulti型の3次元配列の掛け算 C=A.*b.
*/
void rmat3_mul_c(int m, int n, int l, cmulti **C, int LDC1, int LDC2, rmulti **A, int LDA1, int LDA2, cmulti *b)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	cmul_cr(MAT3(C,i,j,k,LDC1,LDC2),b,MAT3(A,i,j,k,LDA1,LDA2));
      }
    }
  }
}

/** @} */

///////////////////////////////////

/** @name cmulti型の3次元配列の数学関数に関する関数 */
/** @{ */

/**
 @brief cmulti型の3次元配列の絶対値 C=abs(A).
*/
void cmat3_abs(int m, int n, int l, rmulti **C, int LDC1, int LDC2, cmulti **A, int LDA1, int LDA2)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	rabs_c(MAT3(C,i,j,k,LDC1,LDC2),MAT3(A,i,j,k,LDA1,LDA2));
      }
    }
  }
}

//追加

/**
 @brief cmulti型の3次元配列の実部と虚部の絶対値 C=absc(A).
*/
void cmat3_absc(int m, int n, int l, cmulti **C, int LDC1, int LDC2, cmulti **A, int LDA1, int LDA2)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	cabsc_c(MAT3(C,i,j,k,LDC1,LDC2),MAT3(A,i,j,k,LDA1,LDA2));
      }
    }
  }
}

//ここまで

/**
 @brief cmulti型の3次元配列の偏角 B=arg(A).
*/
void cmat3_arg(int m, int n, int l, rmulti **B, int LDB1, int LDB2, cmulti **A, int LDA1, int LDA2)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	rarg_c(MAT3(B,i,j,k,LDB1,LDB2),MAT3(A,i,j,k,LDA1,LDA2));
      }
    }
  }
}

//追加

/**
 @brief cmulti型の3次元配列の和 B=sum(A).
*/
void cmat3_sum(int m, int n, int l, cmulti **B, int LDB1, int LDB2, cmulti **A, int LDA1, int LDA2)
{
  int k;
  for(k=0; k<l; k++){
    cvec_sum_cmat(m,n,&MAT3(B,0,0,k,LDB1,LDB2),&MAT3(A,0,0,k,LDA1,LDA2),LDA1);
  }
}

//ここまで

/**
 @brief cmulti型の3次元配列の最大値 B=max(A).
*/
void cmat3_max(int m, int n, int l, cmulti **B, int LDB1, int LDB2, cmulti **A, int LDA1, int LDA2)
{
  int k;
  for(k=0; k<l; k++){
    cvec_max_cmat(m,n,&MAT3(B,0,0,k,LDB1,LDB2),&MAT3(A,0,0,k,LDA1,LDA2),LDA1);
  }
}

/**
 @brief cmulti型の3次元配列の最小値 B=min(A).
*/
void cmat3_min(int m, int n, int l, cmulti **B, int LDB1, int LDB2, cmulti **A, int LDA1, int LDA2)
{
  int k;
  for(k=0; k<l; k++){
    cvec_min_cmat(m,n,&MAT3(B,0,0,k,LDB1,LDB2),&MAT3(A,0,0,k,LDA1,LDA2),LDA1);
  }
}

/**
 @brief cmulti型の3次元配列の割り算 C=A./B.
*/
void cmat3_div(int m, int n, int l, cmulti **C, int LDC1, int LDC2, cmulti **A, int LDA1, int LDA2, cmulti **B, int LDB1, int LDB2)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	cdiv_cc(MAT3(C,i,j,k,LDC1,LDC2),MAT3(A,i,j,k,LDA1,LDA2),MAT3(B,i,j,k,LDB1,LDB2));
      }
    }
  }
}

/**
 @brief cmulti型の3次元配列の割り算 C=A./B.
*/
void cmat3_div_rmat3(int m, int n, int l, cmulti **C, int LDC1, int LDC2, cmulti **A, int LDA1, int LDA2, rmulti **B, int LDB1, int LDB2)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	cdiv_cr(MAT3(C,i,j,k,LDC1,LDC2),MAT3(A,i,j,k,LDA1,LDA2),MAT3(B,i,j,k,LDB1,LDB2));
      }
    }
  }
}

/**
 @brief cmulti型の3次元配列の割り算 C=A./B.
*/
void rmat3_div_cmat3(int m, int n, int l, cmulti **C, int LDC1, int LDC2, rmulti **A, int LDA1, int LDA2, cmulti **B, int LDB1, int LDB2)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	cdiv_rc(MAT3(C,i,j,k,LDC1,LDC2),MAT3(A,i,j,k,LDA1,LDA2),MAT3(B,i,j,k,LDB1,LDB2));
      }
    }
  }
}

/**
 @brief cmulti型の3次元配列の割り算 C=A./b.
*/
void cmat3_div_c2(int m, int n, int l, cmulti **C, int LDC1, int LDC2, cmulti **A, int LDA1, int LDA2, cmulti *b)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	cdiv_cc(MAT3(C,i,j,k,LDC1,LDC2),MAT3(A,i,j,k,LDA1,LDA2),b);
      }
    }
  }
}

/**
 @brief cmulti型の3次元配列の割り算 C=A./b.
*/
void cmat3_div_r2(int m, int n, int l, cmulti **C, int LDC1, int LDC2, cmulti **A, int LDA1, int LDA2, rmulti *b)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	cdiv_cr(MAT3(C,i,j,k,LDC1,LDC2),MAT3(A,i,j,k,LDA1,LDA2),b);
      }
    }
  }
}

/**
 @brief cmulti型の3次元配列の割り算 C=A./b.
*/
void rmat3_div_c2(int m, int n, int l, cmulti **C, int LDC1, int LDC2, rmulti **A, int LDA1, int LDA2, cmulti *b)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	cdiv_rc(MAT3(C,i,j,k,LDC1,LDC2),MAT3(A,i,j,k,LDA1,LDA2),b);
      }
    }
  }
}

/**
 @brief cmulti型の3次元配列の割り算 C=a./B.
*/
void cmat3_div_c1(int m, int n, int l, cmulti **C, int LDC1, int LDC2, cmulti *a, cmulti **B, int LDB1, int LDB2)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	cdiv_cc(MAT3(C,i,j,k,LDC1,LDC2),a,MAT3(B,i,j,k,LDB1,LDB2));
      }
    }
  }
}

/**
 @brief cmulti型の3次元配列の割り算 C=a./B.
*/
void cmat3_div_r1(int m, int n, int l, cmulti **C, int LDC1, int LDC2, rmulti *a, cmulti **B, int LDB1, int LDB2)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	cdiv_rc(MAT3(C,i,j,k,LDC1,LDC2),a,MAT3(B,i,j,k,LDB1,LDB2));
      }
    }
  }
}

/**
 @brief cmulti型の3次元配列の割り算 C=a./B.
*/
void rmat3_div_c1(int m, int n, int l, cmulti **C, int LDC1, int LDC2, cmulti *a, rmulti **B, int LDB1, int LDB2)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	cdiv_cr(MAT3(C,i,j,k,LDC1,LDC2),a,MAT3(B,i,j,k,LDB1,LDB2));
      }
    }
  }
}

/**
 @brief cmulti型の3次元配列のべき乗 C=A.^B.
*/
void cmat3_pow(int m, int n, int l, cmulti **C, int LDC1, int LDC2, cmulti **A, int LDA1, int LDA2, cmulti **B, int LDB1, int LDB2)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	cpow_cc(MAT3(C,i,j,k,LDC1,LDC2),MAT3(A,i,j,k,LDA1,LDA2),MAT3(B,i,j,k,LDB1,LDB2));
      }
    }
  }
}

/**
 @brief cmulti型の3次元配列のべき乗 C=A.^B.
*/
void cmat3_pow_rmat3(int m, int n, int l, cmulti **C, int LDC1, int LDC2, cmulti **A, int LDA1, int LDA2, rmulti **B, int LDB1, int LDB2)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	cpow_cr(MAT3(C,i,j,k,LDC1,LDC2),MAT3(A,i,j,k,LDA1,LDA2),MAT3(B,i,j,k,LDB1,LDB2));
      }
    }
  }
}

/**
 @brief cmulti型の3次元配列のべき乗 C=A.^B.
*/
void rmat3_pow_cmat3(int m, int n, int l, cmulti **C, int LDC1, int LDC2, rmulti **A, int LDA1, int LDA2, cmulti **B, int LDB1, int LDB2)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	cpow_rc(MAT3(C,i,j,k,LDC1,LDC2),MAT3(A,i,j,k,LDA1,LDA2),MAT3(B,i,j,k,LDB1,LDB2));
      }
    }
  }
}

/**
 @brief cmulti型の3次元配列のべき乗 C=A.^b.
*/
void cmat3_pow_c2(int m, int n, int l, cmulti **C, int LDC1, int LDC2, cmulti **A, int LDA1, int LDA2, cmulti *b)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	cpow_cc(MAT3(C,i,j,k,LDC1,LDC2),MAT3(A,i,j,k,LDA1,LDA2),b);
      }
    }
  }
}

/**
 @brief cmulti型の3次元配列のべき乗 C=A.^b.
*/
void cmat3_pow_r2(int m, int n, int l, cmulti **C, int LDC1, int LDC2, cmulti **A, int LDA1, int LDA2, rmulti *b)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	cpow_cr(MAT3(C,i,j,k,LDC1,LDC2),MAT3(A,i,j,k,LDA1,LDA2),b);
      }
    }
  }
}

/**
 @brief cmulti型の3次元配列のべき乗 C=A.^b.
*/
void rmat3_pow_c2(int m, int n, int l, cmulti **C, int LDC1, int LDC2, rmulti **A, int LDA1, int LDA2, cmulti *b)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	cpow_rc(MAT3(C,i,j,k,LDC1,LDC2),MAT3(A,i,j,k,LDA1,LDA2),b);
      }
    }
  }
}

/**
 @brief cmulti型の3次元配列のべき乗 C=a.^B.
*/
void cmat3_pow_c1(int m, int n, int l, cmulti **C, int LDC1, int LDC2, cmulti *a, cmulti **B, int LDB1, int LDB2)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	cpow_cc(MAT3(C,i,j,k,LDC1,LDC2),a,MAT3(B,i,j,k,LDB1,LDB2));
      }
    }
  }
}

/**
 @brief cmulti型の3次元配列のべき乗 C=a.^B.
*/
void cmat3_pow_r1(int m, int n, int l, cmulti **C, int LDC1, int LDC2, rmulti *a, cmulti **B, int LDB1, int LDB2)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	cpow_rc(MAT3(C,i,j,k,LDC1,LDC2),a,MAT3(B,i,j,k,LDB1,LDB2));
      }
    }
  }
}

/**
 @brief cmulti型の3次元配列のべき乗 C=a.^B.
*/
void rmat3_pow_c1(int m, int n, int l, cmulti **C, int LDC1, int LDC2, cmulti *a, rmulti **B, int LDB1, int LDB2)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	cpow_cr(MAT3(C,i,j,k,LDC1,LDC2),a,MAT3(B,i,j,k,LDB1,LDB2));
      }
    }
  }
}

//追加

/**
 @brief cmulti型の3次元配列の指数部で評価 C=10^(floor(log10(abs(x))-y))
*/
void cmat3_exp10_floor_log10_abs_sub(int m, int n, int l, cmulti **C, int LDC1, int LDC2, cmulti **A, int LDA1, int LDA2, double b)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	cexp10_floor_log10_abs_sub(MAT3(C,i,j,k,LDC1,LDC2),MAT3(A,i,j,k,LDA1,LDA2),b);
      }
    }
  }
}

/**
 @brief cmulti型の3次元配列の指数部で評価 C=2^(floor(log2(abs(x))-y))
*/
void cmat3_exp2_floor_log2_abs_sub(int m, int n, int l, cmulti **C, int LDC1, int LDC2, cmulti **A, int LDA1, int LDA2, double b)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	cexp2_floor_log2_abs_sub(MAT3(C,i,j,k,LDC1,LDC2),MAT3(A,i,j,k,LDA1,LDA2),b);
      }
    }
  }
}

//ここまで
/** @} */

////////////////////////////////////////////////////////////////////////

/** @name cmulti型の3次元配列の要素の比較に関する関数 */
/** @{ */


/**
 @brief rmulti型の3次元配列の比較 B=isnan(A).
*/
void cmat3_isnan(int m, int n, int l, int *B, int LDB1, int LDB2, cmulti **A, int LDA1, int LDA2)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	MAT3(B,i,j,k,LDB1,LDB2)=cis_nan(MAT3(A,i,j,k,LDA1,LDA2));
      }
    }
  }
  return;
}

/**
 @brief rmulti型の3次元配列の比較 B=isinf(A).
*/
void cmat3_isinf(int m, int n, int l, int *B, int LDB1, int LDB2, cmulti **A, int LDA1, int LDA2)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	MAT3(B,i,j,k,LDB1,LDB2)=cis_inf(MAT3(A,i,j,k,LDA1,LDA2));
      }
    }
  }
  return;
}

/**
 @brief cmulti型の3次元配列の比較　C=(A==B).
*/
void cmat3_eq(int m, int n, int l, int *C, int LDC1, int LDC2, cmulti **A, int LDA1, int LDA2, cmulti **B, int LDB1, int LDB2)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	MAT3(C,i,j,k,LDC1,LDC2)=eq_cc(MAT3(A,i,j,k,LDA1,LDA2),MAT3(B,i,j,k,LDB1,LDB2));
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
	MAT3(C,i,j,k,LDC1,LDC2)=eq_cr(MAT3(A,i,j,k,LDA1,LDA2),MAT3(B,i,j,k,LDB1,LDB2));
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
	MAT3(C,i,j,k,LDC1,LDC2)=eq_cc(MAT3(A,i,j,k,LDA1,LDA2),b);
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
	MAT3(C,i,j,k,LDC1,LDC2)=eq_cr(MAT3(A,i,j,k,LDA1,LDA2),b);
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
	MAT3(C,i,j,k,LDC1,LDC2)=eq_cr(b,MAT3(A,i,j,k,LDA1,LDA2));
      }
    }
  }
  return;
}

/**
 @brief cmulti型の3次元配列の比較　C=(A!=B).
*/
void cmat3_ne(int m, int n, int l, int *C, int LDC1, int LDC2, cmulti **A, int LDA1, int LDA2, cmulti **B, int LDB1, int LDB2)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	MAT3(C,i,j,k,LDC1,LDC2)=ne_cc(MAT3(A,i,j,k,LDA1,LDA2),MAT3(B,i,j,k,LDB1,LDB2));
      }
    }
  }
  return;
}

/**
 @brief cmulti型の3次元配列の比較　C=(A!=B).
*/
void cmat3_ne_rmat3(int m, int n, int l, int *C, int LDC1, int LDC2, cmulti **A, int LDA1, int LDA2, rmulti **B, int LDB1, int LDB2)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	MAT3(C,i,j,k,LDC1,LDC2)=ne_cr(MAT3(A,i,j,k,LDA1,LDA2),MAT3(B,i,j,k,LDB1,LDB2));
      }
    }
  }
  return;
}

/**
 @brief cmulti型の3次元配列の比較　C=(A!=b).
*/
void cmat3_ne_c(int m, int n, int l, int *C, int LDC1, int LDC2, cmulti **A, int LDA1, int LDA2, cmulti *b)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	MAT3(C,i,j,k,LDC1,LDC2)=ne_cc(MAT3(A,i,j,k,LDA1,LDA2),b);
      }
    }
  }
  return;
}

/**
 @brief cmulti型の3次元配列の比較　C=(A!=b).
*/
void cmat3_ne_r(int m, int n, int l, int *C, int LDC1, int LDC2, cmulti **A, int LDA1, int LDA2, rmulti *b)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	MAT3(C,i,j,k,LDC1,LDC2)=ne_cr(MAT3(A,i,j,k,LDA1,LDA2),b);
      }
    }
  }
  return;
}

/**
 @brief cmulti型の3次元配列の比較　C=(A!=b).
*/
void rmat3_ne_c(int m, int n, int l, int *C, int LDC1, int LDC2, rmulti **A, int LDA1, int LDA2, cmulti *b)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	MAT3(C,i,j,k,LDC1,LDC2)=ne_cr(b,MAT3(A,i,j,k,LDA1,LDA2));
      }
    }
  }
  return;
}

/**
 @brief cmulti型の3次元配列の比較　C=(A>=B).
*/
void cmat3_ge(int m, int n, int l, int *C, int LDC1, int LDC2, cmulti **A, int LDA1, int LDA2, cmulti **B, int LDB1, int LDB2)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	MAT3(C,i,j,k,LDC1,LDC2)=ge_cc(MAT3(A,i,j,k,LDA1,LDA2),MAT3(B,i,j,k,LDB1,LDB2));
      }
    }
  }
  return;
}

/**
 @brief cmulti型の3次元配列の比較　C=(A>=B).
*/
void cmat3_ge_rmat3(int m, int n, int l, int *C, int LDC1, int LDC2, cmulti **A, int LDA1, int LDA2, rmulti **B, int LDB1, int LDB2)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	MAT3(C,i,j,k,LDC1,LDC2)=ge_cr(MAT3(A,i,j,k,LDA1,LDA2),MAT3(B,i,j,k,LDB1,LDB2));
      }
    }
  }
  return;
}

/**
 @brief cmulti型の3次元配列の比較　C=(A>=b).
*/
void cmat3_ge_c(int m, int n, int l, int *C, int LDC1, int LDC2, cmulti **A, int LDA1, int LDA2, cmulti *b)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	MAT3(C,i,j,k,LDC1,LDC2)=ge_cc(MAT3(A,i,j,k,LDA1,LDA2),b);
      }
    }
  }
  return;
}

/**
 @brief cmulti型の3次元配列の比較　C=(A>=b).
*/
void cmat3_ge_r(int m, int n, int l, int *C, int LDC1, int LDC2, cmulti **A, int LDA1, int LDA2, rmulti *b)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	MAT3(C,i,j,k,LDC1,LDC2)=ge_cr(MAT3(A,i,j,k,LDA1,LDA2),b);
      }
    }
  }
  return;
}

/**
 @brief cmulti型の3次元配列の比較　C=(A>=b).
*/
void rmat3_ge_c(int m, int n, int l, int *C, int LDC1, int LDC2, rmulti **A, int LDA1, int LDA2, cmulti *b)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	MAT3(C,i,j,k,LDC1,LDC2)=ge_rc(MAT3(A,i,j,k,LDA1,LDA2),b);
      }
    }
  }
  return;
}

/**
 @brief cmulti型の3次元配列の比較　C=(A>B).
*/
void cmat3_gt(int m, int n, int l, int *C, int LDC1, int LDC2, cmulti **A, int LDA1, int LDA2, cmulti **B, int LDB1, int LDB2)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	MAT3(C,i,j,k,LDC1,LDC2)=gt_cc(MAT3(A,i,j,k,LDA1,LDA2),MAT3(B,i,j,k,LDB1,LDB2));
      }
    }
  }
  return;
}

/**
 @brief cmulti型の3次元配列の比較　C=(A>B).
*/
void cmat3_gt_rmat3(int m, int n, int l, int *C, int LDC1, int LDC2, cmulti **A, int LDA1, int LDA2, rmulti **B, int LDB1, int LDB2)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	MAT3(C,i,j,k,LDC1,LDC2)=gt_cr(MAT3(A,i,j,k,LDA1,LDA2),MAT3(B,i,j,k,LDB1,LDB2));
      }
    }
  }
  return;
}

/**
 @brief cmulti型の3次元配列の比較　C=(A>b).
*/
void cmat3_gt_c(int m, int n, int l, int *C, int LDC1, int LDC2, cmulti **A, int LDA1, int LDA2, cmulti *b)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	MAT3(C,i,j,k,LDC1,LDC2)=gt_cc(MAT3(A,i,j,k,LDA1,LDA2),b);
      }
    }
  }
  return;
}

/**
 @brief cmulti型の3次元配列の比較　C=(A>b).
*/
void cmat3_gt_r(int m, int n, int l, int *C, int LDC1, int LDC2, cmulti **A, int LDA1, int LDA2, rmulti *b)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	MAT3(C,i,j,k,LDC1,LDC2)=gt_cr(MAT3(A,i,j,k,LDA1,LDA2),b);
      }
    }
  }
  return;
}

/**
 @brief cmulti型の3次元配列の比較　C=(A>b).
*/
void rmat3_gt_c(int m, int n, int l, int *C, int LDC1, int LDC2, rmulti **A, int LDA1, int LDA2, cmulti *b)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	MAT3(C,i,j,k,LDC1,LDC2)=gt_rc(MAT3(A,i,j,k,LDA1,LDA2),b);
      }
    }
  }
  return;
}

/**
 @brief cmulti型の3次元配列の比較　C=(A<=B).
*/
void cmat3_le(int m, int n, int l, int *C, int LDC1, int LDC2, cmulti **A, int LDA1, int LDA2, cmulti **B, int LDB1, int LDB2)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	MAT3(C,i,j,k,LDC1,LDC2)=le_cc(MAT3(A,i,j,k,LDA1,LDA2),MAT3(B,i,j,k,LDB1,LDB2));
      }
    }
  }
  return;
}

/**
 @brief cmulti型の3次元配列の比較　C=(A<=B).
*/
void cmat3_le_rmat3(int m, int n, int l, int *C, int LDC1, int LDC2, cmulti **A, int LDA1, int LDA2, rmulti **B, int LDB1, int LDB2)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	MAT3(C,i,j,k,LDC1,LDC2)=le_cr(MAT3(A,i,j,k,LDA1,LDA2),MAT3(B,i,j,k,LDB1,LDB2));
      }
    }
  }
  return;
}

/**
 @brief cmulti型の3次元配列の比較　C=(A<=b).
*/
void cmat3_le_c(int m, int n, int l, int *C, int LDC1, int LDC2, cmulti **A, int LDA1, int LDA2, cmulti *b)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	MAT3(C,i,j,k,LDC1,LDC2)=le_cc(MAT3(A,i,j,k,LDA1,LDA2),b);
      }
    }
  }
  return;
}

/**
 @brief cmulti型の3次元配列の比較　C=(A<=b).
*/
void cmat3_le_r(int m, int n, int l, int *C, int LDC1, int LDC2, cmulti **A, int LDA1, int LDA2, rmulti *b)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	MAT3(C,i,j,k,LDC1,LDC2)=le_cr(MAT3(A,i,j,k,LDA1,LDA2),b);
      }
    }
  }
  return;
}

/**
 @brief cmulti型の3次元配列の比較　C=(A<=b).
*/
void rmat3_le_c(int m, int n, int l, int *C, int LDC1, int LDC2, rmulti **A, int LDA1, int LDA2, cmulti *b)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	MAT3(C,i,j,k,LDC1,LDC2)=le_rc(MAT3(A,i,j,k,LDA1,LDA2),b);
      }
    }
  }
  return;
}

/**
 @brief cmulti型の3次元配列の比較　C=(A<B).
*/
void cmat3_lt(int m, int n, int l, int *C, int LDC1, int LDC2, cmulti **A, int LDA1, int LDA2, cmulti **B, int LDB1, int LDB2)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	MAT3(C,i,j,k,LDC1,LDC2)=lt_cc(MAT3(A,i,j,k,LDA1,LDA2),MAT3(B,i,j,k,LDB1,LDB2));
      }
    }
  }
  return;
}

/**
 @brief cmulti型の3次元配列の比較　C=(A<B).
*/
void cmat3_lt_rmat3(int m, int n, int l, int *C, int LDC1, int LDC2, cmulti **A, int LDA1, int LDA2, rmulti **B, int LDB1, int LDB2)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	MAT3(C,i,j,k,LDC1,LDC2)=lt_cr(MAT3(A,i,j,k,LDA1,LDA2),MAT3(B,i,j,k,LDB1,LDB2));
      }
    }
  }
  return;
}

/**
 @brief cmulti型の3次元配列の比較　C=(A < b).
*/
void cmat3_lt_c(int m, int n, int l, int *C, int LDC1, int LDC2, cmulti **A, int LDA1, int LDA2, cmulti *b)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	MAT3(C,i,j,k,LDC1,LDC2)=lt_cc(MAT3(A,i,j,k,LDA1,LDA2),b);
      }
    }
  }
  return;
}

/**
 @brief cmulti型の3次元配列の比較　C=(A < b).
*/
void cmat3_lt_r(int m, int n, int l, int *C, int LDC1, int LDC2, cmulti **A, int LDA1, int LDA2, rmulti *b)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	MAT3(C,i,j,k,LDC1,LDC2)=lt_cr(MAT3(A,i,j,k,LDA1,LDA2),b);
      }
    }
  }
  return;
}

/**
 @brief cmulti型の3次元配列の比較　C=(A < b).
*/
void rmat3_lt_c(int m, int n, int l, int *C, int LDC1, int LDC2, rmulti **A, int LDA1, int LDA2, cmulti *b)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	MAT3(C,i,j,k,LDC1,LDC2)=lt_rc(MAT3(A,i,j,k,LDA1,LDA2),b);
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
