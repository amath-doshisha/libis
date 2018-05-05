#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<stdarg.h>
#include"is_macros.h"
#include"is_strings.h"
#include"is_svec.h"
#include"is_dcomplex.h"
#include"is_rmulti.h"
#include"is_cmulti.h"
#include"is_cmat.h"
#include"is_rmat.h"
#include"is_zmat.h"
#include"is_dmat.h"
#include"is_cvec.h"
#include"is_csolve.h"
#include"is_func.h"

/**
 @file  cmat.c
 @brief 多倍長精度複素数型cmultiの行列に関する関数の定義
 @details cmulti型の関数に関する定義は@link cmulti.c@endlinkを参照のこと.
          cmulti型のベクトルの関数に関する定義は@link cvec.c@endlinkを参照のこと.
 */

/////////////////////////////////////////////////////

#define FILE_NAME_LENGTH_MAX 100000

/////////////////////////////////////////////////////

/** @name cmulti型の行列の初期化に関する関数 */
/** @{ */

/**
 @brief  cmulti型の行列の新規生成.
*/
cmulti **cmat_allocate(int LDA, int n)
{
  int i,j;
  cmulti **A=NULL;
  A=(cmulti**)malloc(sizeof(cmulti*)*LDA*n);
  for(j=0; j<n; j++){
    for(i=0; i<LDA; i++){
      MAT(A,i,j,LDA)=callocate();
    }
  }
  return A;
}

/**
 @brief cmulti型の行列の精度を指定しての新規生成.
*/
cmulti **cmat_allocate_prec(int LDA, int n, int prec)
{
  int i,j;
  cmulti **A=NULL;
  A=(cmulti**)malloc(sizeof(cmulti*)*LDA*n);
  for(j=0; j<n; j++){
    for(i=0; i<LDA; i++){
      MAT(A,i,j,LDA)=callocate_prec(prec);
    }
  }
  return A;
}

/**
 @brief cmulti型の行列の複製の生成.
*/
cmulti **cmat_allocate_clone(int LDB, int n, cmulti **B)
{
  int LDA,i,j;
  cmulti **A=NULL;
  LDA=LDB;
  A=(cmulti**)malloc(sizeof(cmulti*)*LDA*n);
  for(j=0; j<n; j++){
    for(i=0; i<LDA; i++){
      MAT(A,i,j,LDA)=callocate_clone(MAT(B,i,j,LDB));
    }
  }
  return A;
}

/**
 @brief cmulti型の行列の複製の生成.
*/
cmulti **cmat_allocate_clone_r(int LDB, int n, rmulti **B)
{
  int LDA,i,j;
  cmulti **A=NULL;
  LDA=LDB;
  A=(cmulti**)malloc(sizeof(cmulti*)*LDA*n);
  for(j=0; j<n; j++){
    for(i=0; i<LDA; i++){
      MAT(A,i,j,LDA)=callocate_clone_r(MAT(B,i,j,LDB));
    }
  }
  return A;
}

/**
 @brief cmulti型の行列の終了処理.
*/
cmulti **cmat_free(int LDA, int n, cmulti **A)
{
  int i,j;
  if(A==NULL) return NULL;
  for(j=0; j<n; j++){
    for(i=0; i<LDA; i++){
      MAT(A,i,j,LDA)=cfree(MAT(A,i,j,LDA));
    }
  }
  free(A);
  return A=NULL;
}

/**
 @brief 初期化済みのcmulti型の行列の浮動小数点数の精度(ビット数)を変更し再初期化.
*/
void cmat_round(int m, int n, cmulti **A, int LDA, int prec)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      cround(MAT(A,i,j,LDA),prec);
    }
  }
}

/**
 @brief cmulti型の行列の値を複成 B=A.
*/
void cmat_clone(int m, int n, cmulti **B, int LDB, cmulti **A, int LDA)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      cclone(MAT(B,i,j,LDB),MAT(A,i,j,LDA));
    }
  }
}

/**
 @brief cmulti型の行列の値を複成 B=A.
*/
void cmat_clone_r(int m, int n, cmulti **B, int LDB, rmulti **A, int LDA)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      cclone_r(MAT(B,i,j,LDB),MAT(A,i,j,LDA));
    }
  }
}

/**
 @brief cmulti型の行列の転置をとり値を複成 B=A^T.
*/
void cmat_clone_t(int m, int n, cmulti **B, int LDB, cmulti **A, int LDA)
{
  int i,j;
  for(i=0; i<m; i++){
    for(j=0; j<n; j++){
      cclone(MAT(B,j,i,LDB),MAT(A,i,j,LDA));
    }
  }
}

/**
 @brief cmulti型の行列の共役転置をとり値を複成 B=A'.
*/
void cmat_clone_ct(int m, int n, cmulti **B, int LDB, cmulti **A, int LDA)
{
  int i,j;
  for(i=0; i<m; i++){
    for(j=0; j<n; j++){
      cconj_clone(MAT(B,j,i,LDB),MAT(A,i,j,LDA));
    }
  }
}

/**
 @brief cmulti型の行列の値の交換 A<=>B.
*/
void cmat_swap(int m, int n, cmulti **A, int LDA, cmulti **B, int LDB)
{
  int j;
  for(j=0; j<n; j++){
    cvec_swap(m,&COL(A,j,LDA),&COL(B,j,LDB));
  }
}

/** @} */

/////////////////////////////////////////////////////

/** @name cmulti型の行列の型変換に関する関数 */
/** @{ */

/**
 @brief cmulti型の行列を文字列型に変換 y=char(x)
*/
void cmat_get_s(int m, int n, char **B, int LDB, cmulti **A, int LDA, char format, int digits)
{
  char f[1024],buf[1<<13];
  int i,j;
  sprintf(f,"%%-.%dR%c%%+.%dR%ci",digits,format,digits,format);
  if(A==NULL){ return; }
  for(i=0; i<m; i++){
    for(j=0; j<n; j++){
      mpfr_sprintf(buf,f,C_R(MAT(A,i,j,LDA)),C_I(MAT(A,i,j,LDA)));
      MAT(B,i,j,LDB)=char_renew(MAT(B,i,j,LDB),buf,NULL);
    }
  }
}


/** @} */

/////////////////////////////////////////////////////

/** @name cmulti型の行列のメンバ変数に関する関数 */
/** @{ */

/**
 @brief cmulti型の行列の浮動小数点数の精度(ビット数)を取得.
*/
void cmat_get_prec(int m, int n, int *P, int LDP, cmulti **A, int LDA)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      MAT(P,i,j,LDP)=cget_prec(MAT(A,i,j,LDA));
    }
  }
}

/**
 @brief cmulti型行列の精度(ビット数)の最大値の取得.
 */
int cmat_get_prec_max(int m, int n, cmulti **A, int LDA)
{
  int value,p,i=0,j=0;
  value=cget_prec(MAT(A,i,j,LDA));
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      p=cget_prec(MAT(A,i,j,LDA));
      if(p>=value){ value=p; }
    }
  }
  return value;
}

/**
 @brief cmulti型の行列の浮動小数点数の指数部の取得.
*/
void cmat_get_exp(int m, int n, int *P, int LDP, cmulti **A, int LDA)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      MAT(P,i,j,LDP)=cget_exp(MAT(A,i,j,LDA));
    }
  }
}

/**
 @brief cmulti型の行列の浮動小数点数の指数部の最大値の取得.
*/
int cmat_get_exp_max(int m, int n, cmulti **A, int LDA)
{
  int value,p,i=0,j=0;
  value=cget_exp(MAT(A,i,j,LDA));
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      p=cget_exp(MAT(A,i,j,LDA));
      if(p>=value){ value=p; }
    }
  }
  return value;
}

/**
 @brief cmulti型の行列が数であるかの判定.
*/
int cmat_is_number(int m, int n, cmulti **A, int LDA)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      if(!cis_number(MAT(A,i,j,LDA))){ return 0; }
    }
  }
  return 1;
}

/**
 @brief cmulti型の行列がNaNを含むかの判定.
*/
int cmat_has_nan(int m, int n, cmulti **A, int LDA)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      if(cis_nan(MAT(A,i,j,LDA))){ return 1; }
    }
  }
  return 0;
}

/**
 @brief cmulti型の行列が実数であるかの判定.
*/
int cmat_is_real(int m, int n, cmulti **A, int LDA)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      if(!cis_real(MAT(A,i,j,LDA))){ return 0; }
    }
  }
  return 1;
}


/** @} */

////////////////////////////////////////////////////////////////////////

/** @name cmulti型の行列の入出力に関する関数 */
/** @{ */

/**
 @brief cmulti型の行列の表示.
*/
void cmat_print(int m, int n, cmulti **A, int LDA, char *name, char format, int digits)
{
  char **s=NULL;
  if(A==NULL){
    if(name!=NULL){ printf("%s=NULL\n",name); }
    else          { printf("NULL\n"); }
    return;
  }
  s=svec_allocate(m*n);
  cmat_get_s(m,n,s,m,A,LDA,format,digits);
  smat_print(m,n,s,m,name);
  s=svec_free(m*n,s);
}

/**
 @brief cmulti型の行列の保存.
*/
void cmat_save(int m, int n, cmulti **A, int LDA, int digits, char* fmt, ...)
{
  int i,j;
  char fname[FILE_NAME_LENGTH_MAX+1];
  va_list argp;
  FILE *fid;
  char gmpfmt[16];
  // file name
  va_start(argp,fmt);
  vsprintf(fname,fmt,argp);
  // open file
  if((fid=fopen(fname,"w"))==0){ ERROR_AT; printf("Cant't open file: %s\n",fname); exit(0); }
  // write
  sprintf(gmpfmt,"%%d\t%%d\t%%.%dRe\t%%.%dRe\n",digits,digits);
  if(A!=NULL){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	if(mpfr_fprintf(fid,gmpfmt,i,j,C_R(MAT(A,i,j,LDA)),C_I(MAT(A,i,j,LDA)))<0){ ERROR_AT; exit(0); }
      }
    }
  }
  // close
  fclose(fid);
}

/**
 @brief cmulti型の行列の読み込み.
*/
void cmat_load(int m, int n, cmulti **A, int LDA, char* fmt, ...)
{
  int i,j;
  char fname[FILE_NAME_LENGTH_MAX+1];
  va_list argp;
  FILE *fid;
  char str1[100000],str2[100000];
  // file name
  va_start(argp,fmt);
  vsprintf(fname,fmt,argp);
  // open file
  if((fid=fopen(fname,"r"))==0){ ERROR_AT; printf("Cant't open file: %s\n",fname); exit(0); }
  // real
  for(j=0;j<n;j++){
    for(i=0;i<m;i++){
      if(fscanf(fid,"%s %s",str1,str2)==EOF){ ERROR_AT; exit(0); }
      cset_ss(MAT(A,i,j,LDA),str1,str2);
    }
  }
  // close
  fclose(fid);
}

/**
 @brief cmulti型の行列の保存.
*/
void cmat_bin_save(int m, int n, cmulti **A, int LDA, char* fmt, ...)
{
  int i,j;
  char fname[FILE_NAME_LENGTH_MAX+1];
  va_list argp;
  FILE *fid;
  if(A==NULL){ m=0; n=0; }
  // file name
  va_start(argp,fmt); vsprintf(fname,fmt,argp);
  // open
  if((fid=fopen(fname,"w"))==0){ ERROR_AT; printf("Cant't open file: %s\n",fname); exit(0); }
  // write header
  fwrite("cmat",sizeof(char),strlen("cmat"),fid);
  // write size
  fwrite(&m,sizeof(int),1,fid); // write row-size
  fwrite(&n,sizeof(int),1,fid); // write col-size
  // write data
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      cbin_save(MAT(A,i,j,LDA),fid);
    }
  }
  // close
  fclose(fid);
}


/**
 @brief cmulti型の行列の読み込み.
*/
cmulti **cmat_bin_load(int *m, int *n, char* fmt, ...)
{
  cmulti **cA=NULL;
  rmulti **rA=NULL;
  dcomplex *zA=NULL;
  double *dA=NULL;
  int l,i,j;
  size_t k;
  char fname[FILE_NAME_LENGTH_MAX],*buf=NULL;
  va_list argp;
  FILE *fid;
  // file name
  va_start(argp,fmt); vsprintf(fname,fmt,argp);
  // open
  if((fid=fopen(fname,"r"))==0){ fclose(fid); cA=NULL; (*m)=0; (*n)=0; }
  else{
    // read header
    l=strlen("cmat");
    buf=malloc(sizeof(char)*(l+1));
    k=fread(buf,sizeof(char),l,fid);
    if(k==(size_t)l && strncmp(buf,"cmat",l)==0){ /* cmat */
      // read size
      k=fread(m,sizeof(int),1,fid); if(k!=1 || (*m)<=0){ ERROR_AT; printf("Failed to load the row-size from the file '%s'.\n",fname); exit(0); }
      k=fread(n,sizeof(int),1,fid); if(k!=1 || (*n)<=0){ ERROR_AT; printf("Failed to load the col-size from the file '%s'.\n",fname); exit(0); }
      // allocate
      cA=(cmulti**)malloc(sizeof(cmulti*)*(*m)*(*n));
      // read data
      for(j=0; j<(*n); j++){
	for(i=0; i<(*m); i++){
	  MAT(cA,i,j,(*m))=cbin_load(fid);
	}
      }
      // close
      fclose(fid);
    }else if(k==(size_t)l && strncmp(buf,"rmat",l)==0){ /* rmat */
      fclose(fid);                            // close
      rA=rmat_bin_load(m,n,fname);            // load
      cA=cmat_allocate_clone_r((*m),(*n),rA); // allocate
      rA=rmat_free((*m),(*n),rA);             // free
    }else if(k==(size_t)l && strncmp(buf,"zmat",l)==0){ /* zmat */
      fclose(fid);                           // close
      zA=zmat_bin_load(m,n,fname);           // load
      cA=cmat_allocate_prec((*m),(*n),53);   // allocate
      cmat_set_z((*m),(*n),cA,(*m),zA,(*m)); // copy
      zA=zmat_free(zA);                      // free
    }else if(k==(size_t)l && strncmp(buf,"dmat",l)==0){ /* dmat */
      fclose(fid);                           // close
      dA=dmat_bin_load(m,n,fname);           // load
      cA=cmat_allocate_prec((*m),(*n),53);   // allocate
      cmat_set_d((*m),(*n),cA,(*m),dA,(*m)); // copy
      dA=dmat_free(dA);                      // free
    }else{ fclose(fid); cA=NULL; (*m)=0; (*n)=0; }
  }
  // done
  free(buf);
  return cA;
}

/** @} */

/////////////////////////////////////////////////////////

/** @name cmulti型の行列の値の設定に関する関数 */

/**
 @brief cmulti型の行列の値をNaNに設定.
*/
void cmat_set_nan(int m, int n, cmulti **A, int LDA)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      cset_nan(MAT(A,i,j,LDA));
    }
  }
}

/**
 @brief cmulti型の行列の値を倍精度複素数から設定.
*/
void cmat_set_z(int m, int n, cmulti **B, int LDB, const dcomplex *A, int LDA)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      cset_z(MAT(B,i,j,LDB),MAT(A,i,j,LDA));
    }
  }
}

/**
 @brief cmulti型の行列の値を倍精度実数から設定.
*/
void cmat_set_d(int m, int n, cmulti **B, int LDB, const double *A, int LDA)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      cset_d(MAT(B,i,j,LDB),MAT(A,i,j,LDA));
    }
  }
}

/**
 @brief cmulti型の行列の値を倍精度実数から設定.
*/
void cmat_set_dd(int m, int n, cmulti **B, int LDB, const double *Ar, int LDAr, const double *Ai, int LDAi)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      cset_dd(MAT(B,i,j,LDB),MAT(Ar,i,j,LDAr),MAT(Ai,i,j,LDAi));
    }
  }
}

/**
 @brief cmulti型の行列の値を倍精度複素数から設定.
*/
void cmat_set_all_z(int m, int n, cmulti **A, int LDA, dcomplex a)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      cset_z(MAT(A,i,j,LDA),a);
    }
  }
}

/**
 @brief cmulti型の行列の値を倍精度実数から設定.
*/
void cmat_set_all_dd(int m, int n, cmulti **A, int LDA, double a_r, double a_i)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      cset_dd(MAT(A,i,j,LDA),a_r,a_i);
    }
  }
}

/**
 @brief cmulti型の行列の値を倍精度実数から設定.
*/
void cmat_set_all_d(int m, int n, cmulti **A, int LDA, double a)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      cset_d(MAT(A,i,j,LDA),a);
    }
  }
}

/**
 @brief cmulti型の行列の値を零に設定.
*/
void cmat_set_zeros(int m, int n, cmulti **A, int LDA)
{
  return cmat_set_all_d(m,n,A,LDA,0);
}

/**
 @brief cmulti型の行列の値を1に設定.
*/
void cmat_set_ones(int m, int n, cmulti **A, int LDA)
{
  return cmat_set_all_d(m,n,A,LDA,1);
}

/**
 @brief cmulti型の行列の値を単位行列の設定.
*/
void cmat_set_eye(int m, int n, cmulti **A, int LDA)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      if(i==j){ cset_d(MAT(A,i,j,LDA),1); }
      else    { cset_d(MAT(A,i,j,LDA),0); }
    }
  }
}

/**
 @brief cmulti型の行列の値を区間(b,a+b)の疑似乱数値を設定.
*/
void cmat_set_rand(int m, int n, cmulti **A, int LDA, double a, double b)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      cset_rand(MAT(A,i,j,LDA));
      cmul_d(MAT(A,i,j,LDA),MAT(A,i,j,LDA),a);
      cadd_d(MAT(A,i,j,LDA),MAT(A,i,j,LDA),b);
    }
  }
}


/**
 @brief cmulti型の行列の対角成分にスカラーをコピー.
 @details A(k,k+offset)=a for k=0,1,..
*/
void cmat_set_diag_c(int m, int n, cmulti **A, int LDA, cmulti *a, int offset)
{
  int i,j;
  for(i=0; i<m; i++){
    j=i+offset;
    if(i>=0 && j>=0 && i<m && j<n){
      ccopy(MAT(A,i,j,LDA),a);
    }
  }
}

/**
 @brief cmulti型の行列の対角成分にスカラーをコピー.
 @details A(k,k+offset)=a for k=0,1,..
*/
void cmat_set_diag_d(int m, int n, cmulti **A, int LDA, double a, int offset)
{
  int i,j;
  for(i=0; i<m; i++){
    j=i+offset;
    if(i>=0 && j>=0 && i<m && j<n){
      cset_d(MAT(A,i,j,LDA),a);
    }
  }
}

/** @} */

////////////////////////////////////////////////////////////////////////

/** @name cmulti型の行列の型変換に関する関数 */
/** @{ */

/**
 @brief cmulti型をdcomplex型へキャストする B=double(A).
 */
void cmat_get_z(int m, int n, dcomplex *B, int LDB, cmulti **A, int LDA)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      MAT(B,i,j,LDB)=cget_z(MAT(A,i,j,LDA));
    }
  }
}

/**
 @brief cmulti型をdouble型へキャストする B=double(real(A)).
 */
void cmat_get_d(int m, int n, double *B, int LDB, cmulti **A, int LDA)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      MAT(B,i,j,LDB)=cget_d(MAT(A,i,j,LDA));
    }
  }
}


/** @} */

////////////////////////////////////////////////////////////////////////

/** @name cmulti型の行列の要素の配置に関する関数 */
/** @{ */

/**
 @brief cmulti型の行列の第k列と第l列を入れ替える.
 */
void cmat_cols_swap(int m, int n, cmulti **A, int LDA, int k, int l)
{
  if(k<n && l<n && k!=l){
    cvec_swap(m,&COL(A,k,LDA),&COL(A,l,LDA));
  }
}

/**
 @brief rmulti型の行列の列を列番号の配列Iの値に従って入れ替える.
 */
void cmat_cols_swap_index(int m, int n, cmulti **A, int LDA, const int *I)
{
  int j,LDB;
  cmulti **B=NULL;
  LDB=m; B=cmat_allocate_prec(LDB,n,53);
  for(j=0; j<n; j++){ cvec_swap(m,&COL(B,j,LDB),&COL(A,I[j],LDA)); }
  cmat_swap(m,n,A,LDA,B,LDB);
  B=cmat_free(LDB,n,B);
}


/**
 @brief cmulti型の行列の第k行と第l行を入れ替える.
 */
void cmat_swap_rows(int m, int n, cmulti **A, int LDA, int k, int l)
{
  int j;
  if(k<m && l<m){
    for(j=0; j<n; j++){
      cswap(MAT(A,k,j,LDA),MAT(A,l,j,LDA));
    }
  }
}

/**
 @brief cmulti型の行列の列を右にずらす
 */
void cmat_cols_rotate_right(int m, int n, cmulti **A, int LDA)
{
  int j;
  for(j=n-1; j>0; j--){ cvec_swap(m,&COL(A,j,LDA),&COL(A,j-1,LDA)); }
}

/**
 @brief cmulti型の行列の列を左にずらす
 */
void cmat_cols_rotate_left(int m, int n, cmulti **A, int LDA)
{
  int j;
  for(j=0; j<n-1; j++){ cvec_swap(m,&COL(A,j,LDA),&COL(A,j+1,LDA)); }
}


/** @} */

///////////////////////////////////

/** @name cmulti型の行列の数学関数に関する関数 */
/** @{ */

/**
 @brief cmulti型の行列の値のコピー B=A
*/
void cmat_copy(int m, int n, cmulti **B, int LDB, cmulti **A, int LDA)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      ccopy(MAT(B,i,j,LDB),MAT(A,i,j,LDA));
    }
  }
}

/**
 @brief cmulti型の行列の値のコピー B=A
*/
void cmat_copy_rmat(int m, int n, cmulti **B, int LDB, rmulti **A, int LDA)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      ccopy_r(MAT(B,i,j,LDB),MAT(A,i,j,LDA));
    }
  }
}

/**
 @brief rmulti型の行列を転置をとり値のコピー B=A^T.
 @param[in]  m   行列Aの行の個数.行列Bの列の個数.
 @param[in]  n   行列Aの列の個数.行列Bの行の個数.
 @param[in]  A   初期化済みのrmulti型の行列.サイズは(m,n).
 @param[in]  LDA Aの第1次元.
 @param[in]  B   初期化済みの行列.サイズは(n,m).
 @param[out] B   コピーされた結果.
 @param[in]  LDB Bの第1次元.
*/
void cmat_copy_t(int m, int n, cmulti **B, int LDB, cmulti **A, int LDA)
{
  int i,j;
  for(i=0; i<m; i++){
    for(j=0; j<n; j++){
      ccopy(MAT(B,j,i,LDB),MAT(A,i,j,LDA));
    }
  }
}

/**
 @brief rmulti型の行列を転置をとり値のコピー B=A^T.
 @param[in]  m   行列Aの行の個数.行列Bの列の個数.
 @param[in]  n   行列Aの列の個数.行列Bの行の個数.
 @param[in]  A   初期化済みのrmulti型の行列.サイズは(m,n).
 @param[in]  LDA Aの第1次元.
 @param[in]  B   初期化済みの行列.サイズは(n,m).
 @param[out] B   コピーされた結果.
 @param[in]  LDB Bの第1次元.
*/
void cmat_copy_rmat_t(int m, int n, cmulti **B, int LDB, rmulti **A, int LDA)
{
  int i,j;
  for(i=0; i<m; i++){
    for(j=0; j<n; j++){
      ccopy_r(MAT(B,j,i,LDB),MAT(A,i,j,LDA));
    }
  }
}

/**
 @brief rmulti型の行列を共役転置をとり値のコピー B=A'.
 @param[in]  m   行列Aの行の個数.行列Bの列の個数.
 @param[in]  n   行列Aの列の個数.行列Bの行の個数.
 @param[in]  A   初期化済みのrmulti型の行列.サイズは(m,n).
 @param[in]  LDA Aの第1次元.
 @param[in]  B   初期化済みの行列.サイズは(n,m).
 @param[out] B   コピーされた結果.
 @param[in]  LDB Bの第1次元.
*/
void cmat_copy_ct(int m, int n, cmulti **B, int LDB, cmulti **A, int LDA)
{
  int i,j;
  for(i=0; i<m; i++){
    for(j=0; j<n; j++){
      cconj(MAT(B,j,i,LDB),MAT(A,i,j,LDA));
    }
  }
}


/**
 @brief cmulti型の行列の列を列番号Iに従ったコピー B=A.
 @details B(:,j)=A(:,I(j)) j=0,1,2,..,n-1.
*/
void cmat_copy_col_index(int m, int n, cmulti **B, int LDB, cmulti **A, int LDA, const int *I)
{
  int j;
  for(j=0; j<n; j++){
    cvec_copy(m,&COL(B,j,LDB),&COL(A,I[j],LDA));
  }
}

/**
 @brief cmulti型の行列の値を要素を添字を指定してコピー B=A(I,J)
*/
void cmat_copy_index(int m, int n, cmulti **B, int LDB, cmulti **A, int LDA, int *I, int *J)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      if(I!=NULL && J!=NULL){ ccopy(MAT(B,i,j,LDB),MAT(A,I[i],J[j],LDA)); }
    }
  }
}

/**
 @brief cmulti型の行列の値を要素を添字を指定してコピー B(I,J)=A
*/
void cmat_index_copy(int m, int n, cmulti **B, int LDB, cmulti **A, int LDA, int *I, int *J)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      if(I!=NULL && J!=NULL){ ccopy(MAT(B,I[i],J[j],LDB),MAT(A,i,j,LDA)); }
    }
  }
}

/**
 @brief cmulti型の行列の値を要素を添字を指定してコピー B(I,J)=A
*/
void cmat_index_copy_rmat(int m, int n, cmulti **B, int LDB, rmulti **A, int LDA, int *I, int *J)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      if(I!=NULL && J!=NULL){ ccopy_r(MAT(B,I[i],J[j],LDB),MAT(A,i,j,LDA)); }
    }
  }
}

/**
 @brief cmulti型の行列の実部のコピー B=real(A)
 */
void cmat_real(int m, int n, rmulti **B, int LDB, cmulti **A, int LDA)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      rcopy(MAT(B,i,j,LDB),C_R(MAT(A,i,j,LDA)));
    }
  }
}

/**
 @brief cmulti型の行列の実部の複製 B=real(A)
 */
void cmat_real_clone(int m, int n, rmulti **B, int LDB, cmulti **A, int LDA)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      rclone(MAT(B,i,j,LDB),C_R(MAT(A,i,j,LDA)));
    }
  }
}

/**
 @brief cmulti型の行列の虚部のコピー B=imag(A)
 */
void cmat_imag(int m, int n, rmulti **B, int LDB, cmulti **A, int LDA)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      rcopy(MAT(B,i,j,LDB),C_I(MAT(A,i,j,LDA)));
    }
  }
}

/**
 @brief cmulti型の行列の虚部の複製 B=imag(A)
 */
void cmat_imag_clone(int m, int n, rmulti **B, int LDB, cmulti **A, int LDA)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      rclone(MAT(B,i,j,LDB),C_I(MAT(A,i,j,LDA)));
    }
  }
}

/**
 @brief cmulti型の行列の複素共役のコピー B=conj(A)
*/
void cmat_conj(int m, int n, cmulti **B, int LDB, cmulti **A, int LDA)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      cconj(MAT(B,i,j,LDB),MAT(A,i,j,LDA));
    }
  }
}

/**
 @brief cmulti型の行列の複素共役の複製 B=conj(A)
*/
void cmat_conj_clone(int m, int n, cmulti **B, int LDB, cmulti **A, int LDA)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      cconj_clone(MAT(B,i,j,LDB),MAT(A,i,j,LDA));
    }
  }
}

/**
 @brief cmulti型の行列の符号反転 B=-A.
*/
void cmat_neg(int m, int n, cmulti **B, int LDB, cmulti **A, int LDA)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      cneg(MAT(B,i,j,LDB),MAT(A,i,j,LDA));
    }
  }
}

/**
 @brief cmulti型の行列の要素の絶対値の平方 B=abs2(A).
*/
void cmat_abs2(int m, int n, rmulti **B, int LDB, cmulti **A, int LDA)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      cabs2(MAT(B,i,j,LDB),MAT(A,i,j,LDA));
    }
  }
}

/**
 @brief cmulti型の行列の足し算 C=A+B.
*/
void cmat_add(int m, int n, cmulti **C, int LDC, cmulti **A, int LDA, cmulti **B, int LDB)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      cadd(MAT(C,i,j,LDC),MAT(A,i,j,LDA),MAT(B,i,j,LDB));
    }
  }
}

/**
 @brief cmulti型の行列の足し算 C=A+B.
*/
void cmat_add_rmat(int m, int n, cmulti **C, int LDC, cmulti **A, int LDA, rmulti **B, int LDB)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      cadd_r(MAT(C,i,j,LDC),MAT(A,i,j,LDA),MAT(B,i,j,LDB));
    }
  }
}

/**
 @brief cmulti型の行列の足し算 C=A+b.
*/
void cmat_add_c(int m, int n, cmulti **C, int LDC, cmulti **A, int LDA, cmulti *b)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      cadd(MAT(C,i,j,LDC),MAT(A,i,j,LDA),b);
    }
  }
}

/**
 @brief cmulti型の行列の引き算 C=A-B.
*/
void cmat_sub(int m, int n, cmulti **C, int LDC, cmulti **A, int LDA, cmulti **B, int LDB)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      csub(MAT(C,i,j,LDC),MAT(A,i,j,LDA),MAT(B,i,j,LDB));
    }
  }
}

/**
 @brief cmulti型の行列の引き算 C=A-B.
*/
void cmat_sub_rmat1(int m, int n, cmulti **C, int LDC, rmulti **A, int LDA, cmulti **B, int LDB)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      csub_r1(MAT(C,i,j,LDC),MAT(A,i,j,LDA),MAT(B,i,j,LDB));
    }
  }
}

/**
 @brief cmulti型の行列の引き算 C=A-B.
*/
void cmat_sub_rmat2(int m, int n, cmulti **C, int LDC, cmulti **A, int LDA, rmulti **B, int LDB)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      csub_r2(MAT(C,i,j,LDC),MAT(A,i,j,LDA),MAT(B,i,j,LDB));
    }
  }
}


/**
 @brief cmulti型の行列の引き算 C=a-B.
*/
void cmat_sub_c1(int m, int n, cmulti **C, int LDC, cmulti *a, cmulti **B, int LDB)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      csub(MAT(C,i,j,LDC),a,MAT(B,i,j,LDB));
    }
  }
}

/**
 @brief cmulti型の行列の引き算 C=A-b.
*/
void cmat_sub_c2(int m, int n, cmulti **C, int LDC, cmulti **A, int LDA, cmulti *b)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      csub(MAT(C,i,j,LDC),MAT(A,i,j,LDA),b);
    }
  }
}


/**
 @brief cmulti型の行列の要素ごとの掛け算 C=A.*B
*/
void cmat_mul(int m, int n, cmulti **C, int LDC, cmulti **A, int LDA, cmulti **B, int LDB)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      cmul(MAT(C,i,j,LDC),MAT(A,i,j,LDA),MAT(B,i,j,LDB));
    }
  }
}

/**
 @brief cmulti型の行列の要素ごととスカラーの掛け算 C=A*b.
*/
void cmat_mul_c(int m, int n, cmulti **C, int LDC, cmulti **A, int LDA, cmulti *b)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      cmul(MAT(C,i,j,LDC),MAT(A,i,j,LDA),b);
    }
  }
}

/**
 @brief cmulti型の行列の要素ごととスカラーの掛け算 C=A*b.
*/
void cmat_mul_r(int m, int n, cmulti **C, int LDC, cmulti **A, int LDA, rmulti *b)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      cmul_r(MAT(C,i,j,LDC),MAT(A,i,j,LDA),b);
    }
  }
}

/**
 @brief cmulti型の行列の要素ごととスカラーの掛け算 C=A*b.
*/
void cmat_mul_z(int m, int n, cmulti **C, int LDC, cmulti **A, int LDA, dcomplex b)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      cmul_z(MAT(C,i,j,LDC),MAT(A,i,j,LDA),b);
    }
  }
}

/**
 @brief cmulti型の行列の要素ごととスカラーの掛け算 C=A*b.
*/
void cmat_mul_d(int m, int n, cmulti **C, int LDC, cmulti **A, int LDA, double b)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      cmul_d(MAT(C,i,j,LDC),MAT(A,i,j,LDA),b);
    }
  }
}

/**
 @brief cmulti型の行列の積 C=A*B.
 @param[in]  A 初期化済みのサイズが(l,m)の行列.
 @param[in]  B 初期化済みのサイズが(m,n)の行列.
 @param[in]  C 初期化済みのサイズが(l,n)の行列.
 @param[out] C 計算結果.
*/
void cmat_prod(int l, int m, int n, cmulti **C, int LDC, cmulti **A, int LDA, cmulti **B, int LDB)
{
  int i,j,k;
  cmulti **Z=NULL;
  Z=cmat_allocate_prec(l,n,cmat_get_prec_max(l,n,C,LDC));
  for(i=0; i<l; i++){
    for(j=0; j<n; j++){
      cset_zero(MAT(Z,i,j,l)); // Z=0
      for(k=0; k<m; k++){
	cadd_mul(MAT(Z,i,j,l),MAT(A,i,k,LDA),MAT(B,k,j,LDB)); // Z=A*B
      }
    }
  }
  cmat_copy(l,n,C,LDC,Z,l); // C=Z
  Z=cmat_free(l,n,Z);
}


/**
 @brief cmulti型の行列の積 C=A*B.
 @param[in]  A 初期化済みのサイズが(l,m)の行列.
 @param[in]  B 初期化済みのサイズが(m,n)の行列.
 @param[in]  C 初期化済みのサイズが(l,n)の行列.
 @param[out] C 計算結果.
*/
void cmat_prod_r1(int l, int m, int n, cmulti **C, int LDC, rmulti **A, int LDA, cmulti **B, int LDB)
{
  int i,j,k;
  cmulti **Z=NULL;
  Z=cmat_allocate_prec(l,n,cmat_get_prec_max(l,n,C,LDC));
  for(i=0; i<l; i++){
    for(j=0; j<n; j++){
      cset_zero(MAT(Z,i,j,l)); // Z=0
      for(k=0; k<m; k++){
	cadd_mul_r(MAT(Z,i,j,l),MAT(B,k,j,LDB),MAT(A,i,k,LDA)); // Z=A*B
      }
    }
  }
  cmat_copy(l,n,C,LDC,Z,l); // C=Z
  Z=cmat_free(l,n,Z);
}

/**
 @brief cmulti型の行列の積 C=A*B.
 @param[in]  A 初期化済みのサイズが(l,m)の行列.
 @param[in]  B 初期化済みのサイズが(m,n)の行列.
 @param[in]  C 初期化済みのサイズが(l,n)の行列.
 @param[out] C 計算結果.
*/
void cmat_prod_r2(int l, int m, int n, cmulti **C, int LDC, cmulti **A, int LDA, rmulti **B, int LDB)
{
  int i,j,k;
  cmulti **Z=NULL;
  Z=cmat_allocate_prec(l,n,cmat_get_prec_max(l,n,C,LDC));
  for(i=0; i<l; i++){
    for(j=0; j<n; j++){
      cset_zero(MAT(Z,i,j,l)); // Z=0
      for(k=0; k<m; k++){
	cadd_mul_r(MAT(Z,i,j,l),MAT(A,i,k,LDA),MAT(B,k,j,LDB)); // Z=A*B
      }
    }
  }
  cmat_copy(l,n,C,LDC,Z,l); // C=Z
  Z=cmat_free(l,n,Z);
}

/**
 @brief cmulti型の行列の積の加算 C+=A*B.
 @param[in]  A 初期化済みのサイズが(l,m)の行列.
 @param[in]  B 初期化済みのサイズが(m,n)の行列.
 @param[in]  C 初期化済みのサイズが(l,n)の行列.
 @param[out] C 計算結果.
*/
void cmat_add_prod(int l, int m, int n, cmulti **C, int LDC, cmulti **A, int LDA, cmulti **B, int LDB)
{
  int i,j,k;
  cmulti **Z=NULL;
  Z=cmat_allocate_prec(l,n,cmat_get_prec_max(l,n,C,LDC));
  for(i=0; i<l; i++){
    for(j=0; j<n; j++){
      ccopy(MAT(Z,i,j,l),MAT(C,i,j,LDC)); // Z=C
      for(k=0; k<m; k++){
	cadd_mul(MAT(Z,i,j,l),MAT(A,i,k,LDA),MAT(B,k,j,LDB)); // Z+=A*B
      }
    }
  }
  cmat_copy(l,n,C,LDC,Z,l); // C=Z
  Z=cmat_free(l,n,Z);
}

/**
 @brief cmulti型の行列の積の減算 C-=A*B.
 @param[in]  A 初期化済みのサイズが(l,m)の行列.
 @param[in]  B 初期化済みのサイズが(m,n)の行列.
 @param[in]  C 初期化済みのサイズが(l,n)の行列.
 @param[out] C 計算結果.
*/
void cmat_sub_prod(int l, int m, int n, cmulti **C, int LDC, cmulti **A, int LDA, cmulti **B, int LDB)
{
  int i,j,k;
  cmulti **Z=NULL;
  Z=cmat_allocate_prec(l,n,cmat_get_prec_max(l,n,C,LDC));
  for(i=0; i<l; i++){
    for(j=0; j<n; j++){
      ccopy(MAT(Z,i,j,l),MAT(C,i,j,LDC)); // Z=C
      for(k=0; k<m; k++){
	csub_mul(MAT(Z,i,j,l),MAT(A,i,k,LDA),MAT(B,k,j,LDB)); // Z-=A*B
      }
    }
  }
  cmat_copy(l,n,C,LDC,Z,l); // C=Z
  Z=cmat_free(l,n,Z);
}

/**
 @brief cmulti型の1ランク更新 B=A+a*x*y'
*/
void cmat_rank1op(int m, int n, cmulti **B, int LDB, cmulti **A, int LDA, cmulti *a, cmulti **x, cmulti **y)
{
  int i,j;
  cmulti *value=NULL;
  value=callocate_prec(cmat_get_prec_max(m,n,B,LDB));
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      cdot(value,y[j],x[i]);                     // vlaue=conj(y[j])*x[i]
      cmul(value,a,value);                       // value=a*conj(y[j])*x[i]
      cadd(MAT(B,i,j,LDB),MAT(A,i,j,LDA),value); // B(i,j)=A(i,j)+a*x[i]*conj(y[j])
    }
  }
  value=cfree(value);
}

/**
 @brief cmulti型の行列の対角行列 A=diag(a).
*/
void cmat_diag_copy_cvec(int m, int n, cmulti **A, int LDA, cmulti **a)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      if(i==j){ ccopy(MAT(A,i,j,LDA),a[i]); }
      else    { cset_zero(MAT(A,i,j,LDA)); }
    }
  }
}

/**
 @brief cmulti型の行列の対角行列 A=diag(a).
*/
void cmat_diag_copy_rvec(int m, int n, cmulti **A, int LDA, rmulti **a)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      if(i==j){ ccopy_r(MAT(A,i,j,LDA),a[i]); }
      else    { cset_zero(MAT(A,i,j,LDA)); }
    }
  }
}

/**
 @brief cmulti型の行列の対角行列 A=diag(a).
*/
void cmat_diag_copy_c(int m, int n, cmulti **A, int LDA, cmulti *a)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      if(i==j){ ccopy(MAT(A,i,j,LDA),a); }
      else    { cset_zero(MAT(A,i,j,LDA)); }
    }
  }
}

/**
 @brief cmulti型の行列の対角行列 A=diag(a).
*/
void cmat_diag_copy_r(int m, int n, cmulti **A, int LDA, rmulti *a)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      if(i==j){ ccopy_r(MAT(A,i,j,LDA),a); }
      else    { cset_zero(MAT(A,i,j,LDA)); }
    }
  }
}


/**
 @brief cmulti型の行列の対角行列との足し算 B=A+diag(a).
*/
void cmat_diag_add_cvec(int m, int n, cmulti **B, int LDB, cmulti **A, int LDA, cmulti **a)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      if(i==j){ cadd(MAT(B,i,j,LDB),MAT(A,i,j,LDA),a[i]); }
      else    { ccopy(MAT(B,i,j,LDB),MAT(A,i,j,LDA)); }
    }
  }
}

/**
 @brief cmulti型の行列の対角行列との足し算 B=A+diag(a).
*/
void cmat_diag_add_c(int m, int n, cmulti **B, int LDB, cmulti **A, int LDA, cmulti *a)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      if(i==j){ cadd(MAT(B,i,j,LDB),MAT(A,i,j,LDA),a); }
      else    { ccopy(MAT(B,i,j,LDB),MAT(A,i,j,LDA)); }
    }
  }
}

/**
 @brief cmulti型の行列の対角行列との足し算 B=A+diag(a).
*/
void cmat_diag_add_rvec(int m, int n, cmulti **B, int LDB, cmulti **A, int LDA, rmulti **a)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      if(i==j){ cadd_r(MAT(B,i,j,LDB),MAT(A,i,j,LDA),a[i]); }
      else    { ccopy(MAT(B,i,j,LDB),MAT(A,i,j,LDA)); }
    }
  }
}

/**
 @brief cmulti型の行列の対角行列との足し算 B=A+diag(a).
*/
void cmat_diag_add_r(int m, int n, cmulti **B, int LDB, cmulti **A, int LDA, rmulti *a)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      if(i==j){ cadd_r(MAT(B,i,j,LDB),MAT(A,i,j,LDA),a); }
      else    { ccopy(MAT(B,i,j,LDB),MAT(A,i,j,LDA)); }
    }
  }
}

/**
 @brief cmulti型の行列の対角行列との引き算 B=A-diag(a).
*/
void cmat_diag_sub_cvec(int m, int n, cmulti **B, int LDB, cmulti **A, int LDA, cmulti **a)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      if(i==j){ csub(MAT(B,i,j,LDB),MAT(A,i,j,LDA),a[i]); }
      else    { ccopy(MAT(B,i,j,LDB),MAT(A,i,j,LDA)); }
    }
  }
}

/**
 @brief cmulti型の行列の対角行列との引き算 B=A-diag(a).
*/
void cmat_diag_sub_c(int m, int n, cmulti **B, int LDB, cmulti **A, int LDA, cmulti *a)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      if(i==j){ csub(MAT(B,i,j,LDB),MAT(A,i,j,LDA),a); }
      else    { ccopy(MAT(B,i,j,LDB),MAT(A,i,j,LDA)); }
    }
  }
}

/**
 @brief cmulti型の行列の対角行列との引き算 B=A-diag(a).
*/
void cmat_diag_sub_rvec(int m, int n, cmulti **B, int LDB, cmulti **A, int LDA, rmulti **a)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      if(i==j){ csub_r2(MAT(B,i,j,LDB),MAT(A,i,j,LDA),a[i]); }
      else    { ccopy(MAT(B,i,j,LDB),MAT(A,i,j,LDA)); }
    }
  }
}

/**
 @brief cmulti型の行列の対角行列との引き算 B=A-diag(a).
*/
void cmat_diag_sub_r(int m, int n, cmulti **B, int LDB, cmulti **A, int LDA, rmulti *a)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      if(i==j){ csub_r2(MAT(B,i,j,LDB),MAT(A,i,j,LDA),a); }
      else    { ccopy(MAT(B,i,j,LDB),MAT(A,i,j,LDA)); }
    }
  }
}

/**
 @brief cmulti型の行列の要素ごとの割り算 C=A./B
*/
void cmat_div(int m, int n, cmulti **C, int LDC, cmulti **A, int LDA, cmulti **B, int LDB)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      cdiv(MAT(C,i,j,LDC),MAT(A,i,j,LDA),MAT(B,i,j,LDB));
    }
  }
}

/**
 @brief cmulti型の行列の要素ごとの割り算 C=a./B
*/
void cmat_div_c1(int m, int n, cmulti **C, int LDC, cmulti *a, cmulti **B, int LDB)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      cdiv(MAT(C,i,j,LDC),a,MAT(B,i,j,LDB));
    }
  }
}

/**
 @brief cmulti型の行列の要素ごとの割り算 C=A./b
*/
void cmat_div_c2(int m, int n, cmulti **C, int LDC, cmulti **A, int LDA, cmulti *b)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      cdiv(MAT(C,i,j,LDC),MAT(A,i,j,LDA),b);
    }
  }
}

/**
 @brief cmulti型の逆行列 B=inv(A)
 @param[in]  A 初期化済みのサイズが(n,n)の行列.
 @param[out] B 初期化済みのサイズが(n,n)の行列.
*/
void cmat_inv(int n, cmulti **B, int LDB, cmulti **A, int LDA)
{
  int info;
  cmat_set_eye(n,n,B,LDB);
  csolve(n,n,B,LDB,A,LDA,&info);
}

/**
 @brief cmulti型の行列の累乗 B=A^p.
 @param[in]  A 初期化済みのサイズが(l,l)の行列.
 @param[in]  x スカラー
 @param[out] B 初期化済みのサイズが(l,l)の行列.
*/
void cmat_power(int n, cmulti **B, int LDB, cmulti **A, int LDA, int p)
{
  int LDZ,i;
  cmulti **Z=NULL;
  LDZ=n; Z=cmat_allocate_prec(LDZ,n,cmat_get_prec_max(n,n,B,LDB));
  if(p<0){
    cmat_inv(n,Z,LDZ,A,LDA);
    cmat_power(n,B,LDB,Z,LDZ,-p);
  }else if(p==0){
    cmat_set_eye(n,n,B,LDB);
  }else if(p==1){
    cmat_copy(n,n,B,LDB,A,LDA);
  }else{
    cmat_copy(n,n,B,LDB,A,LDA);
    for(i=1; i<p; i++){
      cmat_prod(n,n,n,Z,LDZ,B,LDB,A,LDA);
      cmat_copy(n,n,B,LDB,Z,LDZ);
    }
  }
  Z=cmat_free(LDZ,n,Z);
}

/**
 @brief cmulti型の行列の要素の絶対値 B=abs(A).
*/
void cmat_abs(int m, int n, rmulti **B, int LDB, cmulti **A, int LDA)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      cabsv(MAT(B,i,j,LDB),MAT(A,i,j,LDA));
    }
  }
}

/**
 @brief cmulti型の行列の要素の偏角 B=arg(A).
*/
void cmat_arg(int m, int n, rmulti **B, int LDB, cmulti **A, int LDA)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      cargument(MAT(B,i,j,LDB),MAT(A,i,j,LDA));
    }
  }
}

//追加

/**
 @brief cmulti型の列ごとの和 B=sum(A)
*/
void cvec_sum_cmat(int m, int n, cmulti **B, cmulti **A, int LDA)
{
  int j;
  for(j=0; j<n; j++){ cvec_sum(B[j],m,&COL(A,j,LDA)); }
}

//ここまで

/**
 @brief cmulti型の列ごとの最大値 B=max(A)
*/
void cvec_max_cmat(int m, int n, cmulti **B, cmulti **A, int LDA)
{
  int j;
  for(j=0; j<n; j++){ cvec_max(B[j],m,&COL(A,j,LDA)); }
}

/**
 @brief cmulti型の列ごとの最小値 B=min(A)
*/
void cvec_min_cmat(int m, int n, cmulti **B, cmulti **A, int LDA)
{
  int j;
  for(j=0; j<n; j++){ cvec_min(B[j],m,&COL(A,j,LDA)); }
}

/**
 @brief cmulti型の行列の要素の絶対値の最大値 value=max(abs(x))
*/
void cmat_max_abs(rmulti *value, int m, int n, cmulti **A, int LDA)
{
  int i,j;
  rmulti *a=NULL;
  a=rallocate_prec(rget_prec(value));
  cabs2(value,MAT(A,0,0,LDA));  // value=abs(x[0])^2
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      cabs2(a,MAT(A,i,j,LDA));  // a=abs(x[i])^2
      if(rgt(a,value)){            // a>value
	rcopy(value,a);         // value=a
      }
    }
  }
  rsqrt(value,value); // value=sqrt(value)
  a=rfree(a);            // free
}

/**
 @brief cmulti型の行列の要素の実部，虚部の絶対値の最大値 value=max(abs(x))
*/
void cmat_max_absc(rmulti *value, int m, int n, cmulti **A, int LDA)
{
  int i,j;
  rmulti *a=NULL;
  a=rallocate_prec(rget_prec(value));
  rabs(value,C_R(MAT(A,0,0,LDA)));
  rabs(a,C_I(MAT(A,0,0,LDA)));
  if(rgt(a,value)){ rcopy(value,a); }
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      rabs(a,C_R(MAT(A,i,j,LDA)));
      if(rgt(a,value)){ rcopy(value,a); }
      rabs(a,C_I(MAT(A,i,j,LDA)));
      if(rgt(a,value)){ rcopy(value,a); }
    }
  }
  a=rfree(a);
}

/**
 @brief cmulti型の行列の要素の絶対値の最小値 value=min(abs(x))
*/
void cmat_min_abs(rmulti *value, int m, int n, cmulti **A, int LDA)
{
  int i,j;
  rmulti *a=NULL;
  a=rallocate_prec(rget_prec(value));
  cabs2(value,MAT(A,0,0,LDA));  // value=abs(x[0])^2
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      cabs2(a,MAT(A,i,j,LDA));    // a=abs(x[i])^2
      if(rlt(a,value)){ // a<value
	rcopy(value,a); // value=a
      }
    }
  }
  rsqrt(value,value); // value=sqrt(value)
  a=rfree(a);         // free
}

/**
 @brief cmulti型の行列の列ごとの正規化.
*/
void cmat_cols_normalize(int m, int n, cmulti **B, int LDB, cmulti **A, int LDA)
{
  int j;
  for(j=0; j<n; j++){
    cvec_normalize(m,&COL(B,j,LDB),&COL(A,j,LDA));
  }
}

/**
 @brief cmulti型の行列の列ごとの正規化.
*/
void cmat_cols_normalize_sgn(int m, int n, cmulti **B, int LDB, cmulti **A, int LDA)
{
  int j;
  for(j=0; j<n; j++){
    cvec_normalize_sgn(m,&COL(B,j,LDB),&COL(A,j,LDA));
  }
}

/**
 @brief cmulti型の行列の列どうしの差の絶対値の最大値 x(j)=max(abs(A(:,j)-B(:,j))).
*/
void cmat_cols_max_abs_sub(rmulti **x, int m, int n, cmulti **A, int LDA, cmulti **B, int LDB)
{
  int j;
  for(j=0; j<n; j++){
    cvec_max_abs_sub(x[j],m,&COL(A,j,LDA),&COL(B,j,LDB));
  }
}

/**
 @brief rmulti型の行列の列とベクトルとの差の絶対値の最大値 y(j)=max(abs(A(:,j)-x)).
*/
void cmat_cols_max_abs_sub_cvec(rmulti **y, int m, int n, cmulti **A, int LDA, cmulti **x)
{
  int j;
  for(j=0; j<n; j++){
    cvec_max_abs_sub(y[j],m,&COL(A,j,LDA),x);
  }
}


/** @} */

////////////////////////////////////////////////////////////////////////

/** @name cmulti型の行列の要素の比較に関する関数 */
/** @{ */

/**
 @brief cmulti型の行列の値の比較 A<==>B.
*/
int cmat_cmp(int m, int n, cmulti **A, int LDA, int k, int l, cmulti **B, int LDB)
{
  int i,j,value;
  for(j=0; j<MIN2(n,l); j++){
    for(i=0; i<MIN2(m,k); i++){
      value=ccmp(MAT(A,i,j,LDA),MAT(B,i,j,LDB));
      if(value!=0){ return value; } // X!=Y
    }
  }
  if     (m<k){ return -1; } // X<Y
  else if(m>k){ return +1; } // X>Y
  if     (n<l){ return -1; } // X<Y
  else if(n>l){ return +1; } // X>Y
  return 0; // X=Y
}

/** @} */

////////////////////////////////////////////////////////////////////

/** @name cmulti型の写像に関する関数 */
/** @{ */

/**
 @brief cmulti型の行列に関する行列写像 A=f(x).
*/
void cmat_func_list2(int m, int n, cmulti **A, int LDA, func_t *f, int l, cmulti **x)
{
  int i,j;
  for(i=0; i<m; i++){
    for(j=0; j<n; j++){
      if(func_is_list(f) && i<func_asize(f) && func_is_list(func_aget(f,i)) && j<func_asize(func_aget(f,i))){
	cvec_func(MAT(A,i,j,LDA),func_aget(func_aget(f,i),j),l,x);
      }else{ cset_nan(MAT(A,i,j,LDA)); }
    }
  }
}

/** @} */

////////////////////////////////////////////////////////////////////

/** @name 未処理 */
/** @{ */

void cmat_angle_deg_list(rmulti **angle, int m, int n, cmulti **A, int LDA)
{
  int i,j,k=0;
  for(i=0; i<n; i++){
    for(j=i+1; j<n; j++){
      cvec_angle_deg(angle[k++],m,&COL(A,i,LDA),&COL(A,j,LDA));
    }
  }
}

/** @} */

//EOF
