#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<stdarg.h>
#include"is_macros.h"
#include"is_rmulti.h"
#include"is_rmat.h"
#include"is_dmat.h"
#include"is_rvec.h"
#include"is_rsolve.h"
#include"is_func.h"

/**
 @file  rmat.c
 @brief 多倍長精度実数型rmultiの行列に関する関数の定義
 @details rmulti型の関数に関する定義は@link rmulti.c@endlinkを参照のこと.
          rmulti型のベクトルの関数に関する定義は@link rvec.c@endlinkを参照のこと.
 */

/////////////////////////////////////////////////////

#define FILE_NAME_LENGTH_MAX 100000

/////////////////////////////////////////////////////

/** @name rmulti型の行列の初期化に関する関数 */
/** @{ */

/**
 @brief  rmulti型の行列の新規生成.
*/
rmulti **rmat_allocate(int LDA, int n)
{
  int i,j;
  rmulti **A=NULL;
  A=(rmulti**)malloc(sizeof(rmulti*)*LDA*n);
  for(j=0; j<n; j++){
    for(i=0; i<LDA; i++){
      MAT(A,i,j,LDA)=rallocate();
    }
  }
  return A;
}

/**
 @brief  rmulti型の行列の精度を指定しての新規生成.
*/
rmulti **rmat_allocate_prec(int LDA, int n, int prec)
{
  int i,j;
  rmulti **A=NULL;
  A=(rmulti**)malloc(sizeof(rmulti*)*LDA*n);
  for(j=0; j<n; j++){
    for(i=0; i<LDA; i++){
      MAT(A,i,j,LDA)=rallocate_prec(prec);
    }
  }
  return A;
}

/**
 @brief  rmulti型の行列の複製の生成.
*/
rmulti **rmat_allocate_clone(int LDB, int n, rmulti **B)
{
  int LDA,i,j;
  rmulti **A=NULL;
  LDA=LDB;
  A=(rmulti**)malloc(sizeof(rmulti*)*LDA*n);
  for(j=0; j<n; j++){
    for(i=0; i<LDA; i++){
      MAT(A,i,j,LDA)=rallocate_clone(MAT(B,i,j,LDB));
    }
  }
  return A;
}

/**
 @brief   rmulti型の行列の終了処理.
*/
rmulti **rmat_free(int LDA, int n, rmulti **A)
{
  int i,j;
  if(A==NULL) return NULL;
  for(j=0; j<n; j++){
    for(i=0; i<LDA; i++){
      MAT(A,i,j,LDA)=rfree(MAT(A,i,j,LDA));
    }
  }
  free(A);
  return A=NULL;
}

/**
 @brief 初期化済みのrmulti型の行列の浮動小数点数の精度(ビット数)を変更し再初期化.
*/
void rmat_round(int m, int n, rmulti **A, int LDA, int prec)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      rround(MAT(A,i,j,LDA),prec);
    }
  }
}


/**
 @brief   rmulti型の行列の値を複成 B=A.
*/
void rmat_clone(int m, int n, rmulti **B, int LDB, rmulti **A, int LDA)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      rclone(MAT(B,i,j,LDB),MAT(A,i,j,LDA));
    }
  }
}

/**
 @brief   rmulti型の行列の転置をとり値を複成 B=A'.
 @param[in]  m   行列Aの行の個数.行列Bの列の個数.
 @param[in]  n   行列Aの列の個数.行列Bの行の個数.
 @param[in]  A   初期化済みのrmulti型の行列.サイズは(m,n).
 @param[in]  LDA Aの第1次元.
 @param[out] B   初期化済みのrmulti型の行列.サイズは(n,m).
 @param[in]  LDB Bの第1次元.
*/
void rmat_clone_t(int m, int n, rmulti **B, int LDB, rmulti **A, int LDA)
{
  int i,j;
  for(i=0; i<m; i++){
    for(j=0; j<n; j++){
      rclone(MAT(B,j,i,LDB),MAT(A,i,j,LDA));
    }
  }
}

/**
 @brief rmulti型の行列の列を列番号Iに従って複成 B(:,j)=A(:,I(j)).
 @details B(:,j)=A(:,I(j)), j=0,1,2,..,n-1.
*/
void rmat_clone_index(int m, int n, rmulti **B, int LDB, rmulti **A, int LDA, const int *I)
{
  int j;
  for(j=0; j<n; j++){
    rvec_clone(m,&COL(B,j,LDB),&COL(A,I[j],LDA));
  }
}

/**
 @brief rmulti型の行列の値の交換 A<=>B.
*/
void rmat_swap(int m, int n, rmulti **A, int LDA, rmulti **B, int LDB)
{
  int j;
  for(j=0; j<n; j++){
    rvec_swap(m,&COL(A,j,LDA),&COL(B,j,LDB));
  }
}

/** @} */


/////////////////////////////////////////////////////

/** @name rmulti型の行列のメンバ変数に関する関数 */
/** @{ */


/**
 @brief rmulti型の行列の浮動小数点数の精度(ビット数)を取得.
*/
void rmat_get_prec(int m, int n, int *P, int LDP, rmulti **A, int LDA)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      MAT(P,i,j,LDP)=rget_prec(MAT(A,i,j,LDA));
    }
  }
}

/**
 @brief rmulti型行列の精度(ビット数)の最大値の取得.
 */
int rmat_get_prec_max(int m, int n, rmulti **A, int LDA)
{
  int i,j,p,value=0;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      p=rget_prec(MAT(A,i,j,LDA));
      if(p>value){ value=p; }
    }
  }
  return value;
}


/**
 @brief rmulti型の行列の浮動小数点数の指数部の取得.
*/
void rmat_get_exp(int m, int n, int *P, int LDP, rmulti **A, int LDA)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      MAT(P,i,j,LDP)=rget_exp(MAT(A,i,j,LDA));
    }
  }
}

/**
 @brief rmulti型の行列の浮動小数点数の指数部の最大値の取得.
*/
int rmat_get_exp_max(int m, int n, rmulti **A, int LDA)
{
  int value,p,i=0,j=0;
  value=rget_exp(MAT(A,i,j,LDA));
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      p=rget_exp(MAT(A,i,j,LDA));
      if(p>=value){ value=p; }
    }
  }
  return value;
}

/**
 @brief rmulti型の行列が数であるかの判定.
*/
int rmat_is_number(int m, int n, rmulti **A, int LDA)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      if(!ris_number(MAT(A,i,j,LDA))){ return 0; }
    }
  }
  return 1;
}

/**
 @brief rmulti型の行列がNaNを含むかの判定.
*/
int rmat_has_nan(int m, int n, rmulti **A, int LDA)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      if(ris_nan(MAT(A,i,j,LDA))){ return 1; }
    }
  }
  return 0;
}


/** @} */

////////////////////////////////////////////////////////////////////////


/** @name rmulti型の行列の入出力に関する関数 */
/** @{ */

/**
 @brief rmulti型の行列の表示.
*/
void rmat_print(int m, int n, rmulti **A, int LDA, const char *name, const char *f, int digits)
{
  int i,j,k;
  char format[128];
  if(STR_EQ(f,"f") || STR_EQ(f,"F")){ k=3; }
  else if((strcmp(f,"e")==0) || (strcmp(f,"E")==0)){ k=7; }
  else{ k=3; }
  sprintf(format,"%%%d.%dR%s ",digits+k,digits,f);
  if(A==NULL){
    if(name!=NULL){ printf("%s=NULL\n",name); }
    else          { printf("NULL\n"); }
    return;
  }
  if(name!=NULL){ printf("%s\n",name); }
  if(A==NULL) return;
  for(i=0; i<m; i++){
    for(j=0; j<n; j++){
      mpfr_printf(format,MAT(A,i,j,LDA));
    }
    printf("\n");
  }
}


/**
 @brief rmulti型の行列の保存.
*/
void rmat_save(int m, int n, rmulti **A, int LDA, int digits, char* fmt, ...)
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
  sprintf(gmpfmt,"%%d\t%%d\t%%.%dRe\n",digits);
  if(A!=NULL){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	if(mpfr_fprintf(fid,gmpfmt,i,j,MAT(A,i,j,LDA))<0){ ERROR_AT; exit(0); }
      }
    }
  }
  // close
  fclose(fid);
}

/**
 @brief rmulti型の行列の読み込み.
*/
void rmat_load(int m, int n, rmulti **A, int LDA, char* fmt, ...)
{
  int i,j;
  char fname[FILE_NAME_LENGTH_MAX+1];
  va_list argp;
  FILE *fid;
  char str[100000];
  // file name
  va_start(argp,fmt);
  vsprintf(fname,fmt,argp);
  // open file
  if((fid=fopen(fname,"r"))==0){ ERROR_AT; printf("Cant't open file: %s\n",fname); exit(0); }
  // real
  for(j=0;j<n;j++){
    for(i=0;i<m;i++){
      if(fscanf(fid,"%s",str)==EOF){ ERROR_AT; }
      rset_s(MAT(A,i,j,LDA),str);
    }
  }
  // close
  fclose(fid);
}

/**
 @brief rmulti型の行列の保存.
*/
void rmat_bin_save(int m, int n, rmulti **A, int LDA, char* fmt, ...)
{
  int i,j;
  char fname[FILE_NAME_LENGTH_MAX+1];
  va_list argp;
  FILE *fid;
  // file name
  va_start(argp,fmt); vsprintf(fname,fmt,argp);
  // open
  if((fid=fopen(fname,"w"))==0){ ERROR_AT; printf("Cant't open file: %s\n",fname); exit(0); } 
  if(A!=NULL){
    // write header
    fwrite("rmat",sizeof(char),strlen("rmat"),fid);
    // write size
    fwrite(&m,sizeof(int),1,fid); // write row-size
    fwrite(&n,sizeof(int),1,fid); // write col-size
    // write data
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	rbin_save(MAT(A,i,j,LDA),fid);
      }
    }
  }
  // close
  fclose(fid);
}

/**
 @brief rmulti型の行列の読み込み.
*/
rmulti **rmat_bin_load(int *m, int *n, char* fmt, ...)
{
  rmulti **rA=NULL;
  double *dA=NULL;
  int l,i,j;
  size_t k;
  char fname[FILE_NAME_LENGTH_MAX],*buf=NULL;
  va_list argp;
  FILE *fid;
  // file name
  va_start(argp,fmt); vsprintf(fname,fmt,argp);
  // open
  if((fid=fopen(fname,"r"))==0){ fclose(fid); rA=NULL; (*m)=0; (*n)=0; }
  else{
    // read header
    l=strlen("rmat");
    buf=malloc(sizeof(char)*(l+1));
    k=fread(buf,sizeof(char),l,fid);
    if(k==(size_t)l && strncmp(buf,"rmat",l)==0){ /* rmat */
      // read size
      k=fread(m,sizeof(int),1,fid); if(k!=1 || (*m)<=0){ ERROR_AT; printf("Failed to load the row-size from the file '%s'.\n",fname); exit(0); }
      k=fread(n,sizeof(int),1,fid); if(k!=1 || (*n)<=0){ ERROR_AT; printf("Failed to load the col-size from the file '%s'.\n",fname); exit(0); }
      // allocate
      rA=(rmulti**)malloc(sizeof(rmulti*)*(*m)*(*n));
      // read data
      for(j=0; j<(*n); j++){ for(i=0; i<(*m); i++){ MAT(rA,i,j,(*m))=rbin_load(fid); } }
      // close
      fclose(fid);
    }else if(k==(size_t)l && strncmp(buf,"dmat",l)==0){ /* dmat */
      fclose(fid);                           // close
      dA=dmat_bin_load(m,n,fname);           // load
      rA=rmat_allocate_prec((*m),(*n),53);   // allocate
      rmat_set_d((*m),(*n),rA,(*m),dA,(*m)); // copy
      dA=dmat_free(dA);                      // free
    }else{ fclose(fid); rA=NULL; (*m)=0; (*n)=0; }
  }
  // done
  free(buf);
  return rA;
}


/** @} */

/////////////////////////////////////////////////////////

/** @name rmulti型の行列の値の設定に関する関数 */


/**
 @brief rmulti型の行列の値をNaNに設定.
*/
void rmat_set_nan(int m, int n, rmulti **A, int LDA)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      rset_nan(MAT(A,i,j,LDA));
    }
  }
}

/**
 @brief rmulti型の行列の値を文字列から設定.
*/
void rmat_set_s(int m, int n, rmulti **B, int LDB, char **A, int LDA)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      rset_s(MAT(B,i,j,LDB),MAT(A,i,j,LDA));
    }
  }
}

/**
 @brief rmulti型の行列の値を倍精度実数から設定.
*/
void rmat_set_d(int m, int n, rmulti **B, int LDB, const double *A, int LDA)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      rset_d(MAT(B,i,j,LDB),MAT(A,i,j,LDA));
    }
  }
}

/**
 @brief rmulti型の行列の値を倍精度複素数から設定.
*/
void rmat_set_z(int m, int n, rmulti **B, int LDB, const dcomplex *A, int LDA)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      rset_d(MAT(B,i,j,LDB),Z_R(MAT(A,i,j,LDA)));
    }
  }
}

/**
 @brief rmulti型の行列の値を倍精度浮動小数点数から設定.
*/
void rmat_set_all_d(int m, int n, rmulti **A, int LDA, double a)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      rset_d(MAT(A,i,j,LDA),a);
    }
  }
}

/**
 @brief rmulti型の行列の値を零に設定.
*/
void rmat_set_zeros(int m, int n, rmulti **A, int LDA)
{
  return rmat_set_all_d(m,n,A,LDA,0);
}

/**
 @brief rmulti型の行列の値を1に設定.
*/
void rmat_set_ones(int m, int n, rmulti **A, int LDA)
{
  return rmat_set_all_d(m,n,A,LDA,1);
}

/**
 @brief rmulti型の行列の値を単位行列の設定.
*/
void rmat_set_eye(int m, int n, rmulti **A, int LDA)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      if(i==j){ rset_d(MAT(A,i,j,LDA),1); }
      else    { rset_d(MAT(A,i,j,LDA),0); }
    }
  }
}

/**
 @brief rmulti型の行列の値を区間(b,a+b)の疑似乱数値を設定.
*/
void rmat_set_rand(int m, int n, rmulti **A, int LDA, double a, double b)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      rset_rand(MAT(A,i,j,LDA));
      rmul_d(MAT(A,i,j,LDA),MAT(A,i,j,LDA),a);
      radd_d(MAT(A,i,j,LDA),MAT(A,i,j,LDA),b);
    }
  }
}

/**
 @brief rmulti型の行列の対角成分にスカラーをコピー.
 @details A(k,k+offset)=a for k=0,1,..
*/
void rmat_set_diag_r(int m, int n, rmulti **A, int LDA, rmulti *a, int offset)
{
  int i,j;
  for(i=0; i<m; i++){
    j=i+offset;
    if(i>=0 && j>=0 && i<m && j<n){
      rcopy(MAT(A,i,j,LDA),a);
    }
  }
}

/**
 @brief rmulti型の行列の対角成分にスカラーをコピー.
 @details A(k,k+offset)=a for k=0,1,..
*/
void rmat_set_diag_d(int m, int n, rmulti **A, int LDA, double a, int offset)
{
  int i,j;
  for(i=0; i<m; i++){
    j=i+offset;
    if(i>=0 && j>=0 && i<m && j<n){
      rset_d(MAT(A,i,j,LDA),a);
    }
  }
}

/** @} */


////////////////////////////////////////////////////////////////////////


/** @name rmulti型の行列の型変換に関する関数 */
/** @{ */

/**
 @brief rmulti型の行列を倍精度型に変換.
 */
void rmat_get_d(int m, int n, double *B, int LDB, rmulti **A, int LDA)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      MAT(B,i,j,LDB)=rget_d(MAT(A,i,j,LDA));
    }
  }
}

/**
 @brief rmulti型をdcomplex型へキャストする B=double(A).
 */
void rmat_get_z(int m, int n, dcomplex *B, int LDB, rmulti **A, int LDA)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      Z_SET(MAT(B,i,j,LDB),rget_d(MAT(A,i,j,LDA)),0);
    }
  }
}


/** @} */

////////////////////////////////////////////////////////////////////////

/** @name rmulti型の行列の要素の配置に関する関数 */
/** @{ */


/**
 @brief rmulti型の行列の第k列と第l列を入れ替える.
 */
void rmat_cols_swap(int m, int n, rmulti **A, int LDA, int k, int l)
{
  if(k<n && l<n && k!=l){
    rvec_swap(m,&COL(A,k,LDA),&COL(A,l,LDA));
  }
}

/**
 @brief rmulti型の行列の列を列番号の配列Iの値に従って入れ替える.
 */
void rmat_cols_swap_index(int m, int n, rmulti **A, int LDA, const int *I)
{
  int j,LDB;
  rmulti **B=NULL;
  LDB=m; B=rmat_allocate_prec(LDB,n,53);
  for(j=0; j<n; j++){ rvec_swap(m,&COL(B,j,LDB),&COL(A,I[j],LDA)); }
  rmat_swap(m,n,A,LDA,B,LDB);
  B=rmat_free(LDB,n,B);
}

/**
 @brief rmulti型の行列の第addと第l行を入れ替える.
 */
void rmat_rows_swap(int m, int n, rmulti **A, int LDA, int k, int l)
{
  int j;
  if(k<m && l<m){
    for(j=0; j<n; j++){
      rswap(MAT(A,k,j,LDA),MAT(A,l,j,LDA));
    }
  }
}

/** @} */

////////////////////////////////////////////////////////////////////////

/** @name rmulti型の行列の数学関数に関する関数 */
/** @{ */

/**
 @brief rmulti型の行列の値のコピー B=A
*/
void rmat_copy(int m, int n, rmulti **B, int LDB, rmulti **A, int LDA)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      rcopy(MAT(B,i,j,LDB),MAT(A,i,j,LDA));
    }
  }
}

/**
 @brief rmulti型の行列を転置をとり値のコピー B=A'.
 @param[in]  m   行列Aの行の個数.行列Bの列の個数.
 @param[in]  n   行列Aの列の個数.行列Bの行の個数.
 @param[in]  A   初期化済みのrmulti型の行列.サイズは(m,n).
 @param[in]  LDA Aの第1次元.
 @param[in]  B   初期化済みの行列.サイズは(n,m).
 @param[out] B   コピーされた結果.
 @param[in]  LDB Bの第1次元.
*/
void rmat_copy_t(int m, int n, rmulti **B, int LDB, rmulti **A, int LDA)
{
  int i,j;
  for(i=0; i<m; i++){
    for(j=0; j<n; j++){
      rcopy(MAT(B,j,i,LDB),MAT(A,i,j,LDA));
    }
  }
}

/**
 @brief rmulti型の行列の列を列番号Iに従ったコピー B=A.
 @details B(:,j)=A(:,I(j)) j=0,1,2,..,n-1.
*/
void rmat_copy_col_index(int m, int n, rmulti **B, int LDB, rmulti **A, int LDA, const int *I)
{
  int j;
  for(j=0; j<n; j++){
    rvec_copy(m,&COL(B,j,LDB),&COL(A,I[j],LDA));
  }
}

/**
 @brief rmulti型の行列の値を要素を添字を指定してコピー B=A(I,J)
*/
void rmat_copy_index(int m, int n, rmulti **B, int LDB, rmulti **A, int LDA, int *I, int *J)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      if(I!=NULL && J!=NULL){ rcopy(MAT(B,i,j,LDB),MAT(A,I[i],J[j],LDA)); }
    }
  }
}

/**
 @brief rmulti型の行列の値を要素を添字を指定してコピー B(I,J)=A
*/
void rmat_index_copy(int m, int n, rmulti **B, int LDB, rmulti **A, int LDA, int *I, int *J)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      if(I!=NULL && J!=NULL){ rcopy(MAT(B,I[i],J[j],LDB),MAT(A,i,j,LDA)); }
    }
  }
}

/**
 @brief rmulti型の行列の符号反転 B=-A.
*/
void rmat_neg(int m, int n, rmulti **B, int LDB, rmulti **A, int LDA)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      rneg(MAT(B,i,j,LDB),MAT(A,i,j,LDA));
    }
  }
}

/**
 @brief rmulti型の行列の要素の絶対値 B=abs(A).
*/
void rmat_abs(int m, int n, rmulti **B, int LDB, rmulti **A, int LDA)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      rabs(MAT(B,i,j,LDB),MAT(A,i,j,LDA));
    }
  }
}


/**
 @brief rmulti型の行列の足し算 C=A+B.
*/
void rmat_add(int m, int n, rmulti **C, int LDC, rmulti **A, int LDA, rmulti **B, int LDB)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      radd(MAT(C,i,j,LDC),MAT(A,i,j,LDA),MAT(B,i,j,LDB));
    }
  }
}

/**
 @brief rmulti型の行列の足し算 C=A+b.
*/
void rmat_add_r(int m, int n, rmulti **C, int LDC, rmulti **A, int LDA, rmulti *b)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      radd(MAT(C,i,j,LDC),MAT(A,i,j,LDA),b);
    }
  }
}

/**
 @brief rmulti型の行列の引き算 C=A-B.
*/
void rmat_sub(int m, int n, rmulti **C, int LDC, rmulti **A, int LDA, rmulti **B, int LDB)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      rsub(MAT(C,i,j,LDC),MAT(A,i,j,LDA),MAT(B,i,j,LDB));
    }
  }
}

/**
 @brief rmulti型の行列の引き算 C=a-B.
*/
void rmat_sub_r1(int m, int n, rmulti **C, int LDC, rmulti *a, rmulti **B, int LDB)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      rsub(MAT(C,i,j,LDC),a,MAT(B,i,j,LDB));
    }
  }
}

/**
 @brief rmulti型の行列の引き算 C=A-b.
*/
void rmat_sub_r2(int m, int n, rmulti **C, int LDC, rmulti **A, int LDA, rmulti *b)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      rsub(MAT(C,i,j,LDC),MAT(A,i,j,LDA),b);
    }
  }
}

/**
 @brief rmulti型の行列の要素ごとの掛け算 C=A.*B.
*/
void rmat_mul(int m, int n, rmulti **C, int LDC, rmulti **A, int LDA, rmulti **B, int LDB)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      rmul(MAT(C,i,j,LDC),MAT(A,i,j,LDA),MAT(B,i,j,LDB));
    }
  }
}

/**
 @brief rmulti型の行列の要素ごととスカラーの掛け算 C=A*b.
*/
void rmat_mul_r(int m, int n, rmulti **C, int LDC, rmulti **A, int LDA, rmulti *b)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      rmul(MAT(C,i,j,LDC),MAT(A,i,j,LDA),b);
    }
  }
}

/**
 @brief rmulti型の行列の要素ごととスカラーの掛け算 C=A*b.
*/
void rmat_mul_d(int m, int n, rmulti **C, int LDC, rmulti **A, int LDA, double b)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      rmul_d(MAT(C,i,j,LDC),MAT(A,i,j,LDA),b);
    }
  }
}

/**
 @brief rmulti型の行列の積 C=A*B.
 @param[in]  A 初期化済みのサイズが(l,m)の行列.
 @param[in]  B 初期化済みのサイズが(m,n)の行列.
 @param[out] C 初期化済みのサイズが(l,n)の行列.
*/
void rmat_prod(int l, int m, int n, rmulti **C, int LDC, rmulti **A, int LDA, rmulti **B, int LDB)
{
  int i,j,k;
  rmulti **Z=NULL;
  Z=rmat_allocate_prec(l,n,rmat_get_prec_max(l,n,C,LDC));
  for(i=0; i<l; i++){
    for(j=0; j<n; j++){
      rset_zero(MAT(Z,i,j,l)); // Z=0
      for(k=0; k<m; k++){
	radd_mul(MAT(Z,i,j,l),MAT(A,i,k,LDA),MAT(B,k,j,LDB)); // Z=A*B
      }
    }
  }
  rmat_copy(l,n,C,LDC,Z,l); // C=Z
  Z=rmat_free(l,n,Z);
}

/**
 @brief rmulti型の行列の積の加算 C+=A*B.
 @param[in]  A 初期化済みのサイズが(l,m)の行列.
 @param[in]  B 初期化済みのサイズが(m,n)の行列.
 @param[out] C 初期化済みのサイズが(l,n)の行列.
*/
void rmat_add_prod(int l, int m, int n, rmulti **C, int LDC, rmulti **A, int LDA, rmulti **B, int LDB)
{
  int i,j,k;
  rmulti **Z=NULL;
  Z=rmat_allocate_prec(l,n,rmat_get_prec_max(l,n,C,LDC));
  for(i=0; i<l; i++){
    for(j=0; j<n; j++){
      rcopy(MAT(Z,i,j,l),MAT(C,i,j,LDC)); // Z=C
      for(k=0; k<m; k++){
	radd_mul(MAT(Z,i,j,l),MAT(A,i,k,LDA),MAT(B,k,j,LDB)); // Z+=A*B
      }
    }
  }
  rmat_copy(l,n,C,LDC,Z,l); // C=Z
  Z=rmat_free(l,n,Z);
}

/**
 @brief rmulti型の行列の積の減算 C-=A*B.
 @param[in]  A 初期化済みのサイズが(l,m)の行列.
 @param[in]  B 初期化済みのサイズが(m,n)の行列.
 @param[out] C 初期化済みのサイズが(l,n)の行列.
*/
void rmat_sub_prod(int l, int m, int n, rmulti **C, int LDC, rmulti **A, int LDA, rmulti **B, int LDB)
{
  int i,j,k;
  rmulti **Z=NULL;
  Z=rmat_allocate_prec(l,n,rmat_get_prec_max(l,n,C,LDC));
  for(i=0; i<l; i++){
    for(j=0; j<n; j++){
      rcopy(MAT(Z,i,j,l),MAT(C,i,j,LDC));
      for(k=0; k<m; k++){
	rsub_mul(MAT(Z,i,j,l),MAT(A,i,k,LDA),MAT(B,k,j,LDB));
      }
    }
  }
  rmat_copy(l,n,C,LDC,Z,l);
  Z=rmat_free(l,n,Z);
}

/**
 @brief rmulti型の1ランク更新 B=A+a*x*y'
*/
void rmat_rank1op(int m, int n, rmulti **B, int LDB, rmulti **A, int LDA, rmulti *a, rmulti **x, rmulti **y)
{
  int i,j;
  rmulti *value=NULL;
  value=rallocate_prec(rmat_get_prec_max(m,n,B,LDB));
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      rmul(value,x[i],y[j]);                     // value=x[i]*y[j]
      rmul(value,a,value);                       // value=a*x[i]*y[j]
      radd(MAT(B,i,j,LDB),MAT(A,i,j,LDA),value); // B(i,j)=A(i,j)+a*x[i]*y[j]
    }
  }
  value=rfree(value);
}

/**
 @brief rmulti型の行列の対角行列 A=diag(a).
*/
void rmat_diag_copy_rvec(int m, int n, rmulti **A, int LDA, rmulti **a)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      if(i==j){ rcopy(MAT(A,i,j,LDA),a[i]); }
      else    { rset_zero(MAT(A,i,j,LDA)); }
    }
  }
}

/**
 @brief rmulti型の行列の対角行列 A=diag(a).
*/
void rmat_diag_copy_r(int m, int n, rmulti **A, int LDA, rmulti *a)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      if(i==j){ rcopy(MAT(A,i,j,LDA),a); }
      else    { rset_zero(MAT(A,i,j,LDA)); }
    }
  }
}

/**
 @brief rmulti型の行列の対角行列 A=diag(a).
*/
void rmat_diag_copy_d(int m, int n, rmulti **A, int LDA, double a)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      if(i==j){ rset_d(MAT(A,i,j,LDA),a); }
      else    { rset_zero(MAT(A,i,j,LDA)); }
    }
  }
}

/**
 @brief rmulti型の行列の対角行列との足し算 B=A+diag(a).
*/
void rmat_diag_add_rvec(int m, int n, rmulti **B, int LDB, rmulti **A, int LDA, rmulti **a)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      if(i==j){ radd(MAT(B,i,j,LDB),MAT(A,i,j,LDA),a[i]); }
      else    { rcopy(MAT(B,i,j,LDB),MAT(A,i,j,LDA)); }
    }
  }
}

/**
 @brief rmulti型の行列の対角行列との足し算 B=A+diag(a).
*/
void rmat_diag_add_r(int m, int n, rmulti **B, int LDB, rmulti **A, int LDA, rmulti *a)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      if(i==j){ radd(MAT(B,i,j,LDB),MAT(A,i,j,LDA),a); }
      else    { rcopy(MAT(B,i,j,LDB),MAT(A,i,j,LDA)); }
    }
  }
}

/**
 @brief rmulti型の行列の対角行列との引き算 B=A-diag(a).
*/
void rmat_diag_sub_rvec(int m, int n, rmulti **B, int LDB, rmulti **A, int LDA, rmulti **a)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      if(i==j){ rsub(MAT(B,i,j,LDB),MAT(A,i,j,LDA),a[i]); }
      else    { rcopy(MAT(B,i,j,LDB),MAT(A,i,j,LDA)); }
    }
  }
}

/**
 @brief rmulti型の行列の対角行列との引き算 B=A-diag(a).
*/
void rmat_diag_sub_r(int m, int n, rmulti **B, int LDB, rmulti **A, int LDA, rmulti *a)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      if(i==j){ rsub(MAT(B,i,j,LDB),MAT(A,i,j,LDA),a); }
      else    { rcopy(MAT(B,i,j,LDB),MAT(A,i,j,LDA)); }
    }
  }
}

//追加

/**
 @brief rmulti型の列ごとの和 B=sum(A)
*/
void rvec_sum_rmat(int m, int n, rmulti **B, rmulti **A, int LDA)
{
  int j;
  for(j=0; j<n; j++){ rvec_sum(B[j],m,&COL(A,j,LDA)); }
}

//ここまで

/**
 @brief rmulti型の列ごとの最大値 B=max(A)
*/
void rvec_max_rmat(int m, int n, rmulti **B, rmulti **A, int LDA)
{
  int j;
  for(j=0; j<n; j++){ rvec_max(B[j],m,&COL(A,j,LDA)); }
}

/**
 @brief rmulti型の列ごとの最小値 B=min(A)
*/
void rvec_min_rmat(int m, int n, rmulti **B, rmulti **A, int LDA)
{
  int j;
  for(j=0; j<n; j++){ rvec_min(B[j],m,&COL(A,j,LDA)); }
}

/**
 @brief rmulti型の行列の要素の絶対値の最大値 value=max(abs(x))
*/
void rmat_max_abs(rmulti *value, int m, int n, rmulti **A, int LDA)
{
  int i,j;
  rmulti *a=NULL;
  a=rallocate_prec(rget_prec(value));
  rabs(value,MAT(A,0,0,LDA));  // value=abs(x[0])
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      rabs(a,MAT(A,i,j,LDA));  // a=abs(x[i])
      if(rgt(a,value)){            // a>value
	rcopy(value,a);         // value=a
      }
    }
  }
  a=rfree(a);
}

/**
 @brief rmulti型の行列の要素の絶対値の最小値 value=min(abs(x))
*/
void rmat_min_abs(rmulti *value, int m, int n, rmulti **A, int LDA)
{
  int i,j;
  rmulti *a=NULL;
  a=rallocate_prec(rget_prec(value));
  rabs(value,MAT(A,0,0,LDA));  // value=abs(x[0])
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      rabs(a,MAT(A,i,j,LDA));    // a=abs(x[i])
      if(rlt(a,value)){ // a<value
	rcopy(value,a); // value=a
      }
    }
  }
  a=rfree(a);
}


/**
 @brief rmulti型の行列の列どうしの差の絶対値の最大値 x(j)=max(abs(A(:,j)-B(:,j))).
*/
void rmat_cols_max_abs_sub(rmulti **x, int m, int n, rmulti **A, int LDA, rmulti **B, int LDB)
{
  int j;
  for(j=0; j<n; j++){
    rvec_max_abs_sub(x[j],m,&COL(A,j,LDA),&COL(B,j,LDB));
  }
}

/**
 @brief rmulti型の行列の列とベクトルとの差の絶対値の最大値 y(j)=max(abs(A(:,j)-x)).
*/
void rmat_cols_max_abs_sub_rvec(rmulti **y, int m, int n, rmulti **A, int LDA, rmulti **x)
{
  int j;
  for(j=0; j<n; j++){
    rvec_max_abs_sub(y[j],m,&COL(A,j,LDA),x);
  }
}

/**
 @brief rmulti型の行列の要素ごとの割り算 C=A./B.
*/
void rmat_div(int m, int n, rmulti **C, int LDC, rmulti **A, int LDA, rmulti **B, int LDB)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      rdiv(MAT(C,i,j,LDC),MAT(A,i,j,LDA),MAT(B,i,j,LDB));
    }
  }
}

/**
 @brief rmulti型の行列の要素ごととスカラーの割り算 C=A./b.
*/
void rmat_div_r2(int m, int n, rmulti **C, int LDC, rmulti **A, int LDA, rmulti *b)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      rdiv(MAT(C,i,j,LDC),MAT(A,i,j,LDA),b);
    }
  }
}

/**
 @brief rmulti型の行列の要素ごととスカラーの割り算 C=a./B.
*/
void rmat_div_r1(int m, int n, rmulti **C, int LDC, rmulti *a, rmulti **B, int LDB)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      rdiv(MAT(C,i,j,LDC),a,MAT(B,i,j,LDB));
    }
  }
}

/**
 @brief rmulti型の行列の要素ごととスカラーの割り算 C=A./b.
*/
void rmat_div_d2(int m, int n, rmulti **C, int LDC, rmulti **A, int LDA, double b)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      rdiv_d2(MAT(C,i,j,LDC),MAT(A,i,j,LDA),b);
    }
  }
}

/**
 @brief rmulti型の逆行列 B=inv(A)
 @param[in]  A 初期化済みのサイズが(n,n)の行列.
 @param[out] B 初期化済みのサイズが(n,n)の行列.
*/
void rmat_inv(int n, rmulti **B, int LDB, rmulti **A, int LDA)
{
  int info;
  rmat_set_eye(n,n,B,LDB);
  rsolve(n,n,B,LDB,A,LDA,&info);
}

/**
 @brief rmulti型の行列の累乗 B=A^p.
 @param[in]  A 初期化済みのサイズが(n,n)の行列.
 @param[in]  x スカラー
 @param[out] B 初期化済みのサイズが(n,n)の行列.
*/
void rmat_power(int n, rmulti **B, int LDB, rmulti **A, int LDA, int p)
{
  int LDZ,i;
  rmulti **Z=NULL;
  LDZ=n; Z=rmat_allocate_prec(LDZ,n,rmat_get_prec_max(n,n,B,LDB));
  if(p<0){
    rmat_inv(n,Z,LDZ,A,LDA);
    rmat_power(n,B,LDB,Z,LDZ,-p);
  }else if(p==0){
    rmat_set_eye(n,n,B,LDB);
  }else if(p==1){
    rmat_copy(n,n,B,LDB,A,LDA);
  }else{
    rmat_copy(n,n,B,LDB,A,LDA);
    for(i=1; i<p; i++){
      rmat_prod(n,n,n,Z,LDZ,B,LDB,A,LDA);
      rmat_copy(n,n,B,LDB,Z,LDZ);
    }
  }
  Z=rmat_free(LDZ,n,Z);
}

/**
 @brief rmulti型の行列の列ごとの正規化.
*/
void rmat_cols_normalize(int m, int n, rmulti **B, int LDB, rmulti **A, int LDA)
{
  int j;
  for(j=0; j<n; j++){
    rvec_normalize(m,&COL(B,j,LDB),&COL(A,j,LDA));
  }
}

/**
 @brief rmulti型の行列の列ごとの正規化.
*/
void rmat_cols_normalize_sgn(int m, int n, rmulti **B, int LDB, rmulti **A, int LDA)
{
  int j;
  for(j=0; j<n; j++){
    rvec_normalize_sgn(m,&COL(B,j,LDB),&COL(A,j,LDA));
  }
}

/** @} */


////////////////////////////////////////////////////////////////////////

/** @name rmulti型の行列の要素の比較に関する関数 */
/** @{ */

/**
 @brief rmulti型の行列の値の比較 A<==>B.
*/
int rmat_cmp(int m, int n, rmulti **A, int LDA, int k, int l, rmulti **B, int LDB)
{
  int i,j,value;
  for(j=0; j<MIN2(n,l); j++){
    for(i=0; i<MIN2(m,k); i++){
      value=rcmp(MAT(A,i,j,LDA),MAT(B,i,j,LDB));
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

/** @name rmulti型の写像に関する関数 */
/** @{ */

/**
 @brief rmulti型の行列に関する行列写像 A=f(x).
*/
void rmat_func_list2(int m, int n, rmulti **A, int LDA, func_t *f, int l, rmulti **x)
{
  int i,j;
  for(i=0; i<m; i++){
    for(j=0; j<n; j++){
      if(func_is_list(f) && i<func_asize(f) && func_is_list(func_aget(f,i)) && j<func_asize(func_aget(f,i))){
	rvec_func(MAT(A,i,j,LDA),func_aget(func_aget(f,i),j),l,x);
      }else{ rset_nan(MAT(A,i,j,LDA)); }
    }
  }
}

/** @} */



//EOF
