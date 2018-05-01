#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<stdarg.h>
#include"is_macros.h"
#include"is_rmulti.h"
#include"is_ivec.h"
#include"is_rvec.h"
#include"is_zvec.h"
#include"is_dvec.h"
#include"is_func.h"


/**
 @file  rvec.c
 @brief 多倍長精度実数型rmultiのベクトルに関する関数の定義
 @details rmulti型の関数に関する定義は@link rmulti.c@endlinkを参照のこと.
          rmulti型の行列の関数に関する定義は@link rmat.c@endlinkを参照のこと.
 */


/////////////////////////////////////////////////////

#define FR(F) func_retain(F)
#define RAp(X,Y)   ((X)=rallocate_prec(rget_prec(Y)))
#define RAP(X,Y,N) ((X)=rallocate_prec(rvec_get_prec_max(N,Y)))
#define RF(X)      ((X)=rfree(X))
#define RVAp(X,Y,N){ X=rvec_allocate_prec(N,rvec_get_prec_max(N,Y)); }
#define RVF(X,N)   { X=rvec_free(N,X); }

#define FILE_NAME_LENGTH_MAX 100000

/////////////////////////////////////////////////////

/** @name rmulti型のベクトルの初期化に関する関数 */
/** @{ */

/**
 @brief  rmulti型のベクトルの新規生成.
*/
rmulti **rvec_allocate(int n)
{
  int i;
  rmulti **x;
  x=(rmulti**)malloc(sizeof(rmulti*)*n);
  for(i=0; i<n; i++){ x[i]=rallocate(); }
  return x;
}

/**
 @brief  rmulti型のベクトルの精度を指定しての新規生成.
*/
rmulti **rvec_allocate_prec(int n, int prec)
{
  int i;
  rmulti **x;
  x=(rmulti**)malloc(sizeof(rmulti*)*n);
  for(i=0; i<n; i++){ x[i]=rallocate_prec(prec); }
  return x;
}

/**
 @brief  rmulti型のベクトルの複製の生成.
*/
rmulti **rvec_allocate_clone(int n, rmulti **y)
{
  int i;
  rmulti **x;
  x=(rmulti**)malloc(sizeof(rmulti*)*n);
  for(i=0; i<n; i++){ x[i]=rallocate_clone(y[i]); }
  return x;
}

/**
 @brief   rmulti型のベクトルの終了処理.
*/
rmulti **rvec_free(int n, rmulti **x)
{
  int i;
  if(x==NULL) return NULL;
  for(i=0; i<n; i++){ x[i]=rfree(x[i]); }
  free(x);
  x=NULL;
  return x;
}

/**
 @brief 初期化済みのrmulti型のベクトルの浮動小数点数の精度(ビット数)を変更し再初期化.
*/
void rvec_round(int n, rmulti **x, int prec)
{
  int i;
  for(i=0; i<n; i++){ rround(x[i],prec); }
}

/**
 @brief rmulti型のベクトルの値を複成 y=x.
*/
void rvec_clone(int n, rmulti **y, rmulti **x)
{
  int i;
  for(i=0; i<n; i++){ rclone(y[i],x[i]); }
}

/**
 @brief rmulti型のベクトルの値を添字を指定して複成 y=x(I).
*/
void rvec_clone_index(int n, rmulti **y, rmulti **x, int *I)
{
  int i;
  for(i=0; i<n; i++){ rclone(y[i],x[I[i]]); }
}

/**
 @brief rmulti型のベクトルの値を添字を指定して複成 y(I)=x.
*/
void rvec_index_clone(int n, rmulti **y, rmulti **x, int *I)
{
  int i;
  for(i=0; i<n; i++){ rclone(y[I[i]],x[i]); }
}

/**
 @brief rmulti型のベクトルの値の交換 x<==>y.
*/
void rvec_swap(int n, rmulti **x, rmulti **y)
{
  int i;
  for(i=0; i<n; i++){ rswap(x[i],y[i]); }
}

/** @} */

/////////////////////////////////////////////////////

/** @name rmulti型のベクトルのメンバ変数に関する関数 */
/** @{ */

/**
 @brief rmulti型のベクトルの浮動小数点数の精度(ビット数)を取得.
*/
void rvec_get_prec(int n, int *p, rmulti **x)
{
  int i;
  for(i=0; i<n; i++){ p[i]=rget_prec(x[i]); }
}

/**
 @brief rmulti型のベクトルの浮動小数点数の精度(ビット数)の最大値を取得.
*/
int rvec_get_prec_max(int n, rmulti **x)
{
  int value,p,i=0;
  value=rget_prec(x[i]);
  for(i=1; i<n; i++){
    p=rget_prec(x[i]);
    if(p>=value){ value=p; }
  }
  return value;
}

/**
 @brief rmulti型のベクトルの浮動小数点数の指数部の取得.
*/
void rvec_get_exp(int n, int *p, rmulti **x)
{
  int i;
  for(i=0; i<n; i++){ p[i]=rget_exp(x[i]); }
}

/**
 @brief rmulti型のベクトルの浮動小数点数の指数部の最大値の取得.
*/
int rvec_get_exp_max(int n, rmulti **x)
{
  int value,p,i=0;
  value=rget_exp(x[i]);
  for(i=1; i<n; i++){
    p=rget_exp(x[i]);
    if(p>=value){ value=p; }
  }
  return value;
}

/**
 @brief rmulti型のベクトルが数であるかの判定.
*/
int rvec_is_number(int n, rmulti **x)
{
  int i;
  for(i=0; i<n; i++){ if(!ris_number(x[i])){ return 0; } }
  return 1;
}

/**
 @brief rmulti型のベクトルがNaNを含むかの判定.
*/
int rvec_has_nan(int n, rmulti **x)
{
  int i;
  for(i=0; i<n; i++){ if(ris_nan(x[i])){ return 1; } }
  return 0;
}

/** @} */

////////////////////////////////////////////////////////////////////////


/** @name rmulti型のベクトルの入出力に関する関数 */
/** @{ */

/**
 @brief rmulti型のベクトルの表示.
*/
void rvec_print(int n, rmulti **x, const char *name, const char *f, int digits)
{
  int i,k;
  char format[128];
  if(STR_EQ(f,"f") || STR_EQ(f,"F")){ k=3; }
  else if((strcmp(f,"e")==0) || (strcmp(f,"E")==0)){ k=8; }
  else if((strcmp(f,"g")==0) || (strcmp(f,"G")==0)){ k=7; }
  else{ k=3; }
  sprintf(format,"%%%d.%dR%s\n",digits+k,digits,f);
  if(name!=NULL){ printf("%s\n",name); }
  if(x==NULL){ printf("NULL\n"); return; }
  for(i=0; i<n; i++){ mpfr_printf(format,x[i]); }
}

/**
 @brief rmulti型のベクトルの表示.
*/
void rvec_print_bin(int n, rmulti **x, const char *name, int digits)
{
  int i;
  char format[128];
  rmulti *a;
  a=rallocate();
  sprintf(format,"%%%d.%dRfb%%d\n",digits+4,digits);
  if(x==NULL){
    printf("%s=NULL\n",name);
    return;
  }
  printf("%s\n",name);
  for(i=0; i<n; i++){
    rdiv_2exp(a,x[i],rget_exp(x[i])-1);
    mpfr_printf(format,a,rget_exp(x[i])-1);
  }
  a=rfree(a);
}

/**
 @brief rmulti型のベクトルの表示.
*/
void rvec_print_prec(int n, rmulti **x, const char *name, const char *f, int digits)
{
  int i,k;
  char format[128];
  if(STR_EQ(f,"f") || STR_EQ(f,"F")){ k=3; }
  else if((strcmp(f,"e")==0) || (strcmp(f,"E")==0)){ k=8; }
  else{ k=3; }
  sprintf(format,"%%%d.%dR%s [%%d]\n",digits+k,digits,f);
  if(x==NULL){
    printf("%s=NULL\n",name);
    return;
  }
  printf("%s\n",name);
  for(i=0; i<n; i++){
    mpfr_printf(format,x[i],rget_prec(x[i]));
  }
}

/**
 @brief rmulti型のベクトルの表示.
*/
void rvec_print_exp(int n, rmulti **x, const char *name)
{
  int i;
  if(x==NULL){
    printf("%s=NULL\n",name);
    return;
  }
  printf("%s\n",name);
  for(i=0; i<n; i++){
    if(ris_zero(x[i])){
      mpfr_printf("-Inf\n");
    }else{
      mpfr_printf("%4d\n",rget_exp(x[i]));
    }
  }
}

/**
 @brief rmulti型のベクトルの保存.
*/
void rvec_save(int n, rmulti **x, int offset, int digits, char* fmt, ...)
{
  int i;
  char fname[FILE_NAME_LENGTH_MAX];
  va_list argp;
  FILE *fid;
  char gmpfmt[16];
  // file name
  va_start(argp,fmt);
  vsprintf(fname,fmt,argp);
  // open
  if((fid=fopen(fname,"w"))==0){ ERROR_AT; printf("Cant't open file: %s\n",fname); exit(0); }
  // write
  sprintf(gmpfmt,"%%d\t%%+.%dRe\n",digits);
  for(i=0;i<n;i++){ if(mpfr_fprintf(fid,gmpfmt,i+offset,x[i])<0){ ERROR_AT; exit(0); } }
  // close
  fclose(fid);
}

/**
 @brief rmulti型のベクトルの保存.
*/
void rvec_save_log2_abs(int n, rmulti **x,int offset, int digits, char* fmt, ...)
{
  int i;
  char fname[FILE_NAME_LENGTH_MAX];
  va_list argp;
  FILE *fid;
  char gmpfmt[16];
  rmulti *value=NULL;  
  // allocation
  value=rallocate();
  // file name
  va_start(argp,fmt);
  vsprintf(fname,fmt,argp);
  // open file
  if((fid=fopen(fname,"w"))==0){ ERROR_AT; printf("Cant't open file: %s\n",fname); exit(0); } 
  // write
  sprintf(gmpfmt,"%%d\t%%.%dRe\n",digits);
  for(i=0;i<n;i++){
    rabs(value,x[i]);
    rlog2(value,value);
    if(mpfr_fprintf(fid,gmpfmt,i+offset,value)<0){ ERROR_AT; exit(0); }
  }
  // close
  fclose(fid);
  value=rfree(value);
}

/**
 @brief rmulti型のベクトルの保存.
*/
void rvec_save_itrmap(int n, rmulti **x, int digits, char* fmt, ...)
{
  int i;
  char fname[FILE_NAME_LENGTH_MAX];
  va_list argp;
  FILE *fid;
  char format[16];
  // file name
  va_start(argp,fmt);
  vsprintf(fname,fmt,argp);
  // open file
  if((fid=fopen(fname,"w"))==0){ ERROR_AT; printf("Cant't open file: %s\n",fname); exit(0); }
  // loop
  sprintf(format,"%%+.%dRe\t%%+.%dRe\n",digits,digits);
  for(i=1;i<n;i++){
    if(mpfr_fprintf(fid,format,x[i-1],x[i])<0){ ERROR_AT; exit(0); }
  }
  // close
  fclose(fid);
}

/**
 @brief rmulti型のベクトルの読み込み.
*/
void rvec_load(int n, rmulti **x, char* fmt, ...)
{
  int i;
  char fname[FILE_NAME_LENGTH_MAX];
  va_list argp;
  FILE *fid;
  char str[100000];
  // file name
  va_start(argp,fmt);
  vsprintf(fname,fmt,argp);
  // open file
  if((fid=fopen(fname,"r"))==0){ ERROR_AT; printf("Cant't open file: %s\n",fname); exit(0); }
  // real
  for(i=0;i<n;i++){
    if(fscanf(fid,"%s",str)==EOF){ ERROR_EXIT("Error in rvec_read(), fscanf(), fname=%s\n",fname); }
    rset_s(x[i],str);
  }
  // close
  fclose(fid);
}

/**
 @brief rmulti型のベクトルの保存.
*/
void rvec_bin_save(int n, rmulti **x, char* fmt, ...)
{
  int i;
  char fname[FILE_NAME_LENGTH_MAX+1];
  va_list argp;
  FILE *fid;
  // file name
  va_start(argp,fmt); vsprintf(fname,fmt,argp);
  // open
  if((fid=fopen(fname,"w"))==0){ ERROR_AT; printf("Cant't open file: %s\n",fname); exit(0); }
  if(x!=NULL){
    // write header
    fwrite("rvec",sizeof(char),strlen("rvec"),fid);
    // write size
    fwrite(&n,sizeof(int),1,fid);
    // write data
    for(i=0; i<n; i++){ rbin_save(x[i],fid); }  
  }
  // close
  fclose(fid);
}

/**
 @brief rmulti型のベクトルの読み込み.
*/
rmulti **rvec_bin_load(int *n, char* fmt, ...)
{
  int l,i;
  size_t k;
  rmulti **rx=NULL;
  dcomplex *zx=NULL;
  double *dx=NULL;
  char fname[FILE_NAME_LENGTH_MAX+1],*buf=NULL;
  va_list argp;
  FILE *fid;
  // file name
  va_start(argp,fmt); vsprintf(fname,fmt,argp);
  // open
  if((fid=fopen(fname,"r"))==0){ fclose(fid); rx=NULL; (*n)=0; }
  else{
    // read header
    l=strlen("rvec");
    buf=(char*)malloc(sizeof(char)*(l+1));
    k=fread(buf,sizeof(char),l,fid);
    if(k==(size_t)l && strncmp(buf,"rvec",l)==0){ /* rvec */
      k=fread(n,sizeof(int),1,fid);                  // read size
      if(k!=1 || (*n)<=0){ ERROR_AT; printf("Failed to load the size from the file '%s'.\n",fname); exit(0); }
      rx=(rmulti**)malloc(sizeof(rmulti*)*(*n));     // allocate
      for(i=0; i<(*n); i++){ rx[i]=rbin_load(fid); } // read data
      fclose(fid);                                   // close
    }else if(k==(size_t)l && strncmp(buf,"zvec",l)==0){ /* zvec */
      fclose(fid);                    // close
      zx=zvec_bin_load(n,fname);      // load
      rx=rvec_allocate_prec((*n),53); // allocate
      rvec_set_z((*n),rx,zx);         // copy
      zx=zvec_free(zx);               // free
    }else if(k==(size_t)l && strncmp(buf,"dvec",l)==0){ /* dvec */
      fclose(fid);                    // close
      dx=dvec_bin_load(n,fname);      // load
      rx=rvec_allocate_prec((*n),53); // allocate
      rvec_set_d((*n),rx,dx);         // copy
      dx=dvec_free(dx);               // free
    }else{ fclose(fid); rx=NULL; (*n)=0; }
  }
  // done
  free(buf);
  return rx;
}

/** @} */

/////////////////////////////////////////////////////////

/** @name rmulti型のベクトルの値の設定に関する関数 */
/** @{ */

/**
 @brief rmulti型のベクトルの値をNaNに設定.
*/
void rvec_set_nan(int n, rmulti **x)
{
  int i;
  for(i=0; i<n; i++){ rset_nan(x[i]); }
}

/**
 @brief rmulti型のベクトルの値を文字列から設定.
*/
void rvec_set_s(int n, rmulti **x, const char **str)
{
  int i;
  for(i=0; i<n; i++){ rset_s(x[i],str[i]); }
}

/**
 @brief rmulti型のベクトルの値を倍精度実数から設定.
*/
void rvec_set_d(int n, rmulti **y, double *x)
{
  int i;
  for(i=0; i<n; i++){ rset_d(y[i],x[i]); }
}

/**
 @brief rmulti型のベクトルの値を倍精度複素数から設定.
*/
void rvec_set_z(int n, rmulti **y, dcomplex *x)
{
  int i;
  for(i=0; i<n; i++){ rset_d(y[i],Z_R(x[i])); }
}

/**
 @brief rmulti型のベクトルの全ての値を倍精度実数から設定.
*/
void rvec_set_all_d(int n, rmulti **x, double a)
{
  int i;
  for(i=0; i<n; i++){ rset_d(x[i],a); }
}

/**
 @brief rmulti型のベクトルの値を零に設定.
*/
void rvec_set_zeros(int n, rmulti **x)
{
  return rvec_set_all_d(n,x,0);
}

/**
 @brief rmulti型のベクトルの値を1に設定.
*/
void rvec_set_ones(int n, rmulti **x)
{
  return rvec_set_all_d(n,x,1);
}

/**
 @brief rmulti型のベクトルの値を第k要素が1の基本ベクトルに設定.
*/
void rvec_set_unit(int n, rmulti **x, int k)
{
  rvec_set_zeros(n,x);
  rset_d(x[k],1);
}

/**
 @brief rmulti型のベクトルの値を0,1,2,...に設定.
*/
void rvec_set_grid(int n, rmulti **x)
{
  int i;
  for(i=0; i<n; i++){ rset_d(x[i],i); }
}

/**
 @brief rmulti型のベクトルの値を区間(b,a+b)の疑似乱数値を設定.
*/
void rvec_set_rand(int n, rmulti **x, double a, double b)
{
  int i;
  for(i=0; i<n; i++){
    rset_rand(x[i]);
    rmul_d(x[i],x[i],a);
    radd_d(x[i],x[i],b);
  }
}

/** @} */


////////////////////////////////////////////////////////////////////////


/** @name rmulti型のベクトルの型変換に関する関数 */
/** @{ */

/**
 @brief rmulti型のベクトルを倍精度実数型に変換.
 */
void rvec_get_d(int n, double *y, rmulti **x)
{
  int i;
  for(i=0; i<n; i++){ y[i]=rget_d(x[i]); }
}

/**
 @brief rmulti型のベクトルを倍精度複素数型に変換.
 */
void rvec_get_z(int n, dcomplex *y, rmulti **x)
{
  int i;
  for(i=0; i<n; i++){ Z_SET(y[i],rget_d(x[i]),0); }
}

/** @} */

////////////////////////////////////////////////////////////////////////

/** @name rmulti型のベクトルの要素の配置に関する関数 */
/** @{ */

/**
 @brief rmulti型のベクトルの要素の順を逆に入れ替える.
 */
void rvec_reverse(int n, rmulti **x)
{
  int i;
  for(i=0; i<n/2; i++){
    rswap(x[i],x[n-i-1]);
  }
}

/**
 @brief rmulti型のベクトルの第i要素と第j要素を入れ替える.
 */
void rvec_swap_at(rmulti **x, int i, int j)
{
  rswap(x[i],x[j]);
}

/**
 @brief rmulti型のベクトルの要素を要素番号の配列Iの値に従って入れ替える.
 */
void rvec_swap_index(int n, rmulti **x, const int *I)
{
  int i;
  rmulti **y=NULL;
  y=rvec_allocate_prec(n,53);
  for(i=0; i<n; i++){ rswap(y[i],x[I[i]]); }
  rvec_swap(n,x,y);
  y=rvec_free(n,y);
}

/**
 @brief rmulti型のベクトルの要素を小さい順に入れ替える.
 @details 要素x[i]の値の小さい順に並べ替え，要素番号を配列Iに保存する.
          ただし，I==NULLの場合は保存は行われない.
 @param[in]  n ベクトルのサイズ.
 @param[in]  x 初期化済みのrmulti型のベクトル.
 @param[out] x ソートされたrmulti型のベクトル.
 @param[in]  I 入れ替えられた要素の要素番号の保存先.
 @param[out] I 保存された要素番号.
 */
void rvec_sort(int n, rmulti **x, int *I)
{
  if(I!=NULL){ ivec_grid(n,I); }
  rvec_quick_sort(n,x,I,0,n-1);
}

/**
 @brief rmulti型のベクトルの要素に関するクイックソートアルゴリズム(ユーザは直接使用してはならない).
*/
void rvec_quick_sort(int n, rmulti **x, int *I, int left, int right)
{
  int i,last;
  if(left>=right) return;
  rswap(x[left],x[(left+right)/2]);
  if(I!=NULL) ivec_swap_at(I,left,(left+right)/2);
  last=left;
  for(i=left+1; i<=right; i++){
    if(rlt(x[i],x[left])){
      ++last;
      rswap(x[last],x[i]);
      if(I!=NULL) ivec_swap_at(I,last,i);      
    }
  }
  rswap(x[left],x[last]);
  if(I!=NULL) ivec_swap_at(I,left,last);
  rvec_quick_sort(n,x,I,left,last-1);
  rvec_quick_sort(n,x,I,last+1,right);
}

/**
 @brief rmulti型のベクトルの要素を小さい順に入れ替えたときの要素番号の配列の取得.
*/
void rvec_sort_index(int *I, int n, rmulti **X)
{
  rmulti **Y=NULL;
  Y=rvec_allocate_clone(n,X);
  rvec_sort(n,Y,I);
  Y=rvec_free(n,Y);
}


/** @} */

////////////////////////////////////////////////////////////////////////

/** @name rmulti型のベクトルの数学関数に関する関数 */
/** @{ */

/**
 @brief rmulti型のベクトルの値のコピー y=x.
*/
void rvec_copy(int n, rmulti **y, rmulti **x)
{
  int i;
  if(y==x){ return; }
  for(i=0; i<n; i++){ rcopy(y[i],x[i]); }
}

/**
 @brief rmulti型のベクトルの値の要素番号Iに従ったコピー y=x.
 @details y[i]=x[I[i]], i=0,1,2,..,n-1.
*/
void rvec_copy_index(int n, rmulti **y, rmulti **x, const int *I)
{
  int i;
  if(y==x){ rvec_swap_index(n,x,I); return; }
  for(i=0; i<n; i++){ rcopy(y[i],x[I[i]]); }
}

/**
 @brief rmulti型のベクトルの値を添字を指定してコピー y(I)=x.
*/
void rvec_index_copy(int n, rmulti **y, rmulti **x, int *I)
{
  int i;
  for(i=0; i<n; i++){ rcopy(y[I[i]],x[i]); }
}

/**
 @brief rmulti型のベクトルの指数部の足し算 y=x*2^p
*/
void rvec_mul_2exp(int n, rmulti **y, rmulti **x, int p)
{
  int i;
  for(i=0; i<n; i++){ rmul_2exp(y[i],x[i],p); }
}

/**
 @brief rmulti型のベクトルの指数部の引き算 y=x/2^p
*/
void rvec_div_2exp(int n, rmulti **y, rmulti **x, int p)
{
  int i;
  for(i=0; i<n; i++){ rdiv_2exp(y[i],x[i],p); }
}

/**
 @brief rmulti型のベクトルの符号反転 y=-x
*/
void rvec_neg(int n, rmulti **y, rmulti **x)
{
  int i;
  for(i=0; i<n; i++){ rneg(y[i],x[i]); }
}

/**
 @brief rmulti型のベクトルの要素の絶対値 y=abs(x)
*/
void rvec_abs(int n, rmulti **y, rmulti **x)
{
  int i;
  for(i=0; i<n; i++){ rabs(y[i],x[i]); }
}

/**
 @brief rmulti型のベクトルの足し算 z=x+y
*/
void rvec_add(int n, rmulti **z, rmulti **x, rmulti **y)
{
  int i;
  for(i=0; i<n; i++){ radd(z[i],x[i],y[i]); }
}

/**
 @brief rmulti型のベクトルの足し算 z=x+y
*/
void rvec_add_r(int n, rmulti **z, rmulti **x, rmulti *y)
{
  int i;
  for(i=0; i<n; i++){ radd(z[i],x[i],y); }
}

/**
 @brief rmulti型のベクトルの足し算 z=x+y
*/
void rvec_add_d(int n, rmulti **z, rmulti **x, double y)
{
  int i;
  for(i=0; i<n; i++){ radd_d(z[i],x[i],y); }
}

/**
 @brief rmulti型のベクトルの引き算 z=x-y
*/
void rvec_sub(int n, rmulti **z, rmulti **x, rmulti **y)
{
  int i;
  for(i=0; i<n; i++){ rsub(z[i],x[i],y[i]); }
}

/**
 @brief rmulti型のベクトルの引き算 z=x-y
*/
void rvec_sub_r1(int n, rmulti **z, rmulti *x, rmulti **y)
{
  int i;
  for(i=0; i<n; i++){ rsub(z[i],x,y[i]); }
}

/**
 @brief rmulti型のベクトルの引き算 z=x-y
*/
void rvec_sub_r2(int n, rmulti **z, rmulti **x, rmulti *y)
{
  int i;
  for(i=0; i<n; i++){ rsub(z[i],x[i],y); }
}


/**
 @brief rmulti型のベクトルの要素ごとの掛け算 z=x.*y
*/
void rvec_mul(int n, rmulti **z, rmulti **x, rmulti **y)
{
  int i;
  for(i=0; i<n; i++){ rmul(z[i],x[i],y[i]); }
}

/**
 @brief rmulti型のベクトルxとスカラーyの掛け算 z=x*y
*/
void rvec_mul_r(int n, rmulti **z, rmulti **x, rmulti *y)
{
  int i;
  for(i=0; i<n; i++){ rmul(z[i],x[i],y); }
}

/**
 @brief rmulti型のベクトルxとスカラーyの掛け算 z=x*y
*/
void rvec_mul_d(int n, rmulti **z, rmulti **x, double y)
{
  int i;
  for(i=0; i<n; i++){ rmul_d(z[i],x[i],y); }
}

/**
 @brief rmulti型のベクトルの要素ごとの掛け算の加算 z+=x.*y
*/
void rvec_add_mul(int n, rmulti **z, rmulti **x, rmulti **y)
{
  int i;
  for(i=0; i<n; i++){ radd_mul(z[i],x[i],y[i]); }
}

/**
 @brief rmulti型のベクトルxとスカラーyの掛け算の加算 z+=x*y
*/
void rvec_add_mul_r(int n, rmulti **z, rmulti **x, rmulti *y)
{
  int i;
  for(i=0; i<n; i++){ radd_mul(z[i],x[i],y); }
}

/**
 @brief rmulti型のベクトルxとスカラーyの掛け算の加算 z+=x*y
*/
void rvec_add_mul_d(int n, rmulti **z, rmulti **x, double y)
{
  int i;
  for(i=0; i<n; i++){ radd_mul_d(z[i],x[i],y); }
}

/**
 @brief rmulti型のベクトルの要素ごとの掛け算の減算 z-=x.*y
*/
void rvec_sub_mul(int n, rmulti **z, rmulti **x, rmulti **y)
{
  int i;
  for(i=0; i<n; i++){ rsub_mul(z[i],x[i],y[i]); }
}

/**
 @brief rmulti型のベクトルxとスカラーyの掛け算の減算 z-=x*y
*/
void rvec_sub_mul_r(int n, rmulti **z, rmulti **x, rmulti *y)
{
  int i;
  for(i=0; i<n; i++){ rsub_mul(z[i],x[i],y); }
}

/**
 @brief rmulti型のベクトルxとスカラーyの掛け算の減算 z-=x*y
*/
void rvec_sub_mul_d(int n, rmulti **z, rmulti **x, double y)
{
  int i;
  for(i=0; i<n; i++){ rsub_mul_d(z[i],x[i],y); }
}

/**
 @brief rmulti型のベクトルの要素の差の絶対値 z=abs(x-y)
*/
void rvec_abs_sub(int n, rmulti **z, rmulti **x, rmulti **y)
{
  int i;
  for(i=0; i<n; i++){
    rabs_sub(z[i],x[i],y[i]);
  }
}

/**
 @brief rmulti型のベクトルの要素とスカラーの差の絶対値 z=abs(x-y)
*/
void rvec_abs_sub_r(int n, rmulti **z, rmulti **x, rmulti *y)
{
  int i;
  for(i=0; i<n; i++){
    rabs_sub(z[i],x[i],y);
  }
}

/**
 @brief rmulti型のべき乗 z=x^y
*/
void rvec_pow_ui(int n, rmulti **z, rmulti **x, ulong y)
{
  int i;
  for(i=0; i<n; i++){ rpow_ui(z[i],x[i],y); }
}

/**
 @brief rmulti型のベクトルxに関する行列Aによる1次変換 y=A*x
 @param[in]  m   ベクトルyのサイズ.行列Aの行のサイズ.
 @param[in]  n   ベクトルxのサイズ.行列Aの列のサイズ.
 @param[int] y   初期化済みのサイズmのベクトル.
 @param[out] y   y=A*xの計算結果.
 @param[in]  A   サイズ(m,n)の行列.
 @param[in]  LDA 行列Aの第1次元.
 @param[in]  x   サイズnのベクトル.
 @return         丸めの発生回数.
*/
void rvec_lintr(int m, int n, rmulti **y, rmulti **A, int LDA, rmulti **x)
{
  int i,j;
  rmulti **z=NULL;
  RVAp(z,y,m);
  for(i=0; i<m; i++){
    rset_zero(z[i]); // z=0
    for(j=0; j<n; j++){
      radd_mul(z[i],MAT(A,i,j,LDA),x[j]); // z+=A*x
    }
  }
  rvec_copy(m,y,z); // y=z
  RVF(z,m);
}

/**
 @brief rmulti型のベクトルxに関する行列Aによる1次変換の加算 y+=A*x
 @param[in]  m   ベクトルyのサイズ.行列Aの行のサイズ.
 @param[in]  n   ベクトルxのサイズ.行列Aの列のサイズ.
 @param[int] y   初期化済みのサイズmのベクトル.
 @param[out] y   y+=A*xの計算結果.
 @param[in]  A   サイズ(m,n)の行列.
 @param[in]  LDA 行列Aの第1次元.
 @param[in]  x   サイズnのベクトル.
 @return         丸めの発生回数.
*/
void rvec_add_lintr(int m, int n, rmulti **y, rmulti **A, int LDA, rmulti **x)
{
  int i,j;
  rmulti **z=NULL;
  RVAp(z,y,m);
  for(i=0; i<m; i++){
    rcopy(z[i],y[i]); // z=y
    for(j=0; j<n; j++){
      radd_mul(z[i],MAT(A,i,j,LDA),x[j]); // z+=A*x
    }
  }
  rvec_copy(m,y,z); // y=z
  RVF(z,m);
}

/**
 @brief rmulti型のベクトルxに関する行列Aによる1次変換の減算 y-=A*x
 @param[in]  m   ベクトルyのサイズ.行列Aの行のサイズ.
 @param[in]  n   ベクトルxのサイズ.行列Aの列のサイズ.
 @param[int] y   初期化済みのサイズmのベクトル.
 @param[out] y   y-=A*xの計算結果.
 @param[in]  A   サイズ(m,n)の行列.
 @param[in]  LDA 行列Aの第1次元.
 @param[in]  x   サイズnのベクトル.
 @return         丸めの発生回数.
*/
void rvec_sub_lintr(int m, int n, rmulti **y, rmulti **A, int LDA, rmulti **x)
{
  int i,j;
  rmulti **z=NULL;
  RVAp(z,y,m);
  for(i=0; i<m; i++){
    rcopy(z[i],y[i]); // z=y
    for(j=0; j<n; j++){
      rsub_mul(z[i],MAT(A,i,j,LDA),x[j]); // z-=A*x
    }
  }
  rvec_copy(m,y,z); // y=z
  RVF(z,m);
}

/**
 @brief rmulti型のベクトルxに関する転置行列Aによる1次変換 y=A'*x
 @param[in]  m   ベクトルxのサイズ.行列Aの行のサイズ.
 @param[in]  n   ベクトルyのサイズ.行列Aの列のサイズ.
 @param[int] y   初期化済みのサイズnのベクトル.
 @param[out] y   y=A'*xの計算結果.
 @param[in]  A   サイズ(m,n)の行列.
 @param[in]  LDA 行列Aの第1次元.
 @param[in]  x   サイズmのベクトル.
 @return         丸めの発生回数.
*/
void rvec_lintr_t(int m, int n, rmulti **y, rmulti **A, int LDA, rmulti **x)
{
  int i,j;
  rmulti **z=NULL;
  RVAp(z,y,n);
  for(j=0; j<n; j++){
    rset_zero(z[j]); // z=0
    for(i=0; i<m; i++){
      radd_mul(z[j],MAT(A,i,j,LDA),x[i]); // z=A'*x
    }
  }
  rvec_copy(n,y,z); // y=z
  RVF(z,n);
}

/**
 @brief rmulti型のベクトルxに関する転置行列Aによる1次変換の加算 y+=A'*x
 @param[in]  m   ベクトルxのサイズ.行列Aの行のサイズ.
 @param[in]  n   ベクトルyのサイズ.行列Aの列のサイズ.
 @param[int] y   初期化済みのサイズnのベクトル.
 @param[out] y   y+=A'*xの計算結果.
 @param[in]  A   サイズ(m,n)の行列.
 @param[in]  LDA 行列Aの第1次元.
 @param[in]  x   サイズmのベクトル.
 @return         丸めの発生回数.
*/
void rvec_add_lintr_t(int m, int n, rmulti **y, rmulti **A, int LDA, rmulti **x)
{
  int i,j;
  rmulti **z=NULL;
  RVAp(z,y,n);
  for(j=0; j<n; j++){
    rcopy(z[j],y[j]); // z=y
    for(i=0; i<m; i++){
      radd_mul(z[j],MAT(A,i,j,LDA),x[i]); // z+=A'*x
    }
  }
  rvec_copy(n,y,z); // y=z
  RVF(z,n);
}

/**
 @brief rmulti型のベクトルxに関する転置行列Aによる1次変換の減算 y-=A'*x
 @param[in]  m   ベクトルxのサイズ.行列Aの行のサイズ.
 @param[in]  n   ベクトルyのサイズ.行列Aの列のサイズ.
 @param[int] y   初期化済みのサイズnのベクトル.
 @param[out] y   y-=A'*xの計算結果.
 @param[in]  A   サイズ(m,n)の行列.
 @param[in]  LDA 行列Aの第1次元.
 @param[in]  x   サイズmのベクトル.
 @return         丸めの発生回数.
*/
void rvec_sub_lintr_t(int m, int n, rmulti **y, rmulti **A, int LDA, rmulti **x)
{
  int i,j;
  rmulti **z=NULL;
  RVAp(z,y,n);
  for(j=0; j<n; j++){
    rcopy(z[j],y[j]); // z=y
    for(i=0; i<m; i++){
      rsub_mul(z[j],MAT(A,i,j,LDA),x[i]); // z-=A'*x
    }
  }
  rvec_copy(n,y,z); // y=z
  RVF(z,n);
}

/**
 @brief rmulti型のベクトルの要素の総和 value=sum(x)
 */
void rvec_sum(rmulti *value, int n, rmulti **x)
{
  int i;
  rset_zero(value);
  for(i=0; i<n; i++){ radd(value,value,x[i]); }
}

/**
 @brief rmulti型のベクトルの要素の絶対値の総和 value=sum(abs(x))
 */
void rvec_sum_abs(rmulti *value, int n, rmulti **x)
{
  int i;
  rset_zero(value); // value=0
  for(i=0; i<n; i++){ radd_abs(value,x[i]); } // value+=abs(x[i])
}

/**
 @brief rmulti型のベクトルの要素の平方の総和 value=sum(x.^2)
*/
void rvec_sum_pow2(rmulti *value, int n, rmulti **x)
{
  int i;
  rset_zero(value);
  for(i=0; i<n; i++){ radd_mul(value,x[i],x[i]); }
}

/**
 @brief rmulti型のベクトルの内積 value=sum(x.*y)
 */
void rvec_sum_mul(rmulti *value, int n, rmulti **x, rmulti **y)
{
  int i;
  rset_zero(value);
  for(i=0; i<n; i++){ radd_mul(value,x[i],y[i]); }
}

/**
 @brief rmulti型のベクトルの要素の差の絶対値の総和 value=sum(abs(x-y))
*/
void rvec_sum_abs_sub(rmulti *value, int n, rmulti **x, rmulti **y)
{
  int i;
  rmulti *a=NULL;
  RAp(a,value);
  rset_zero(value); // value=0
  for(i=0; i<n; i++){
    rsub(a,x[i],y[i]); // a=x[i]-y[i]
    radd_abs(value,a); // value+=abs(x[i]-y[i])
  }
  RF(a);
}


/**
 @brief rmulti型のベクトルの要素の最大値 value=max(x)
*/
void rvec_max(rmulti *value, int n, rmulti **x)
{
  int i;
  rcopy(value,x[0]);      // value=x[0]
  for(i=1; i<n; i++){
    if(rgt(x[i],value)){     // x[i]>value
      rcopy(value,x[i]);  // value=x[i]
    }
  }
}

/**
 @brief rmulti型のベクトルの要素の絶対値の最大値 value=max(abs(x))
*/
void rvec_max_abs(rmulti *value, int n, rmulti **x)
{
  int i;
  rmulti *a=NULL;
  RAp(a,value);
  rabs(value,x[0]);   // value=abs(x[0])
  for(i=1; i<n; i++){
    rabs(a,x[i]);     // a=abs(x[i])
    if(rgt(a,value)){    // a>value
      rcopy(value,a); // value=a
    }
  }
  RF(a);
}

/**
 @brief rmulti型のベクトルの要素の絶対値の最大値 value=max(abs(x)) とその添え字 I
*/
void rvec_max_abs_index(rmulti *value, int n, rmulti **x, int *I)
{
  int i;
  rmulti *a=NULL;
  RAp(a,value);
  rabs(value,x[0]);       // value=abs(x[0])
  if(I!=NULL){ (*I)=0; }     // I=i
  for(i=1; i<n; i++){
    rabs(a,x[i]);         // value=abs(x[i])
    if(rgt(a,value)){        // a>value
      rcopy(value,a);        // value=a
      if(I!=NULL){ (*I)=i; } // I=i
    }
  }
  RF(a);
}

/**
 @brief rmulti型のベクトルの要素の差の絶対値の最大値 value=max(abs(x-y))
*/
void rvec_max_abs_sub(rmulti *value, int n, rmulti **x, rmulti **y)
{
  int i;
  rmulti *a=NULL;
  RAp(a,value);
  rsub(a,x[0],y[0]);   // a=x[0]-y[0]
  rabs(value,a);       // value=abs(x[0]-y[0])
  for(i=1; i<n; i++){
    rsub(a,x[i],y[i]); // a=x[i]-y[i]
    rabs(a,a);         // a=abs(x[i]-y[i])
    if(rgt(a,value)){
      rcopy(value,a);  // value=a
    }
  }
  RF(a);
}

/**
 @brief rmulti型のベクトルの要素の最小値 value=min(x)
*/
void rvec_min(rmulti *value, int n, rmulti **x)
{
  int i;
  rcopy(value,x[0]);      // value=x[0]
  for(i=1; i<n; i++){
    if(rlt(x[i],value)){     // x[i]<value
      rcopy(value,x[i]);  // value=x[i]
    }
  }
}

/**
 @brief rmulti型のベクトルの要素の絶対値の最小値 value=min(abs(x))
*/
void rvec_min_abs(rmulti *value, int n, rmulti **x)
{
  int i;
  rmulti *a=NULL;
  RAp(a,value);
  rabs(value,x[0]);       // value=abs(x[0])
  for(i=1; i<n; i++){
    rabs(a,x[i]);         // value=abs(x[i])
    if(rlt(a,value)){     // a<value
      rcopy(value,a);     // value=a
    }
  }
  RF(a);
}

/**
 @brief rmulti型のベクトルの要素の絶対値の最小値 value=min(abs(x)) とその添え字 I
*/
void rvec_min_abs_index(rmulti *value, int n, rmulti **x, int *I)
{
  int i;
  rmulti *a=NULL;
  RAp(a,value);
  rabs(value,x[0]);       // value=abs(x[0])
  if(I!=NULL){ (*I)=0; }     // I=i
  for(i=1; i<n; i++){
    rabs(a,x[i]);         // value=abs(x[i])
    if(rlt(a,value)){        // a<value
      rcopy(value,a);     // value=a
      if(I!=NULL){ (*I)=i; } // I=i
    }
  }
  RF(a);
}

/**
 @brief rmulti型のベクトルの要素ごとの割り算 z=x./y
*/
void rvec_div(int n, rmulti **z, rmulti **x, rmulti **y)
{
  int i;
  for(i=0; i<n; i++){ rdiv(z[i],x[i],y[i]); }
}

/**
 @brief rmulti型のスカラーとベクトルの要素との割り算 z=x./y
*/
void rvec_div_r1(int n, rmulti **z, rmulti *x, rmulti **y)
{
  int i;
  for(i=0; i<n; i++){ rdiv(z[i],x,y[i]); }
}

/**
 @brief rmulti型のベクトルの要素とスカラーの割り算 z=x/y
*/
void rvec_div_r2(int n, rmulti **z, rmulti **x, rmulti *y)
{
  int i;
  for(i=0; i<n; i++){ rdiv(z[i],x[i],y); }
}

/**
 @brief rmulti型のスカラーとベクトルの要素との割り算 z=x/y
*/
void rvec_div_d1(int n, rmulti **z, double x, rmulti **y)
{
  int i;
  for(i=0; i<n; i++){ rdiv_d1(z[i],x,y[i]); }
}

/**
 @brief rmulti型のベクトルの要素とスカラーの割り算 z=x/y
*/
void rvec_div_d2(int n, rmulti **z, rmulti **x, double y)
{
  int i;
  for(i=0; i<n; i++){ rdiv_d2(z[i],x[i],y); }
}

/**
 @brief rmulti型のべき乗 z=x.^y
*/
void rvec_pow(int n, rmulti **z, rmulti **x, rmulti **y) 
{
  int i;
  for(i=0; i<n; i++){ rpow(z[i],x[i],y[i]); }
}

/**
 @brief rmulti型のべき乗 z=x^y
*/
void rvec_pow_si(int n, rmulti **z, rmulti **x, long y)
{
  int i;
  for(i=0; i<n; i++){ rpow_si(z[i],x[i],y); }
}

/**
 @brief rmulti型のべき乗 z=x^y
*/
void rvec_pow_r(int n, rmulti **z, rmulti **x, rmulti *y)
{
  int i;
  for(i=0; i<n; i++){ rpow(z[i],x[i],y); }
}

/**
 @brief rmulti型のベクトルの要素ごとの絶対値の対数 y=log2(abs(x))
*/
void rvec_log2_abs(int n, rmulti **y, rmulti **x)
{
  int i;
  rmulti *a=NULL;
  RAP(a,y,n);
  for(i=0; i<n; i++){
    rabs(a,x[i]);
    rlog2(y[i],a);
  }
  RF(a);
}

/**
 @brief rmulti型のベクトルの規格化 y=x/sqrt(x'*x)
*/
void rvec_normalize(int n, rmulti **y, rmulti **x)
{
  rmulti *a=NULL;
  RAP(a,y,n);
  rvec_norm2(a,n,x);   // a=sqrt(x'*x)
  rinv(a,a);           // a=1/sqrt(x'*x)
  rvec_mul_r(n,y,x,a); // y=x/sqrt(x'*x)
  RF(a);
}

/**
 @brief rmulti型のベクトルの規格化 y=x/sqrt(x'*x)
*/
void rvec_normalize_sgn(int n, rmulti **y, rmulti **x)
{
  int k=-1;
  rmulti *a=NULL;
  RAP(a,y,n);
  rvec_max_abs_index(a,n,x,&k);        // a=abs(x[k])
  rvec_norm2(a,n,x);                   // a=sqrt(x'*x)
  rinv(a,a);                           // a=1/a
  if(ris_negative(x[k])){ rneg(a,a); } // a=-a
  rvec_mul_r(n,y,x,a);                 // y=x/sqrt(x'*x)
  RF(a);
}

/**
 @brief rmulti型のベクトルの直交化 y-=(x'*y)*x where x'*x=1
*/
void rvec_orthogonalize(int n, rmulti **y, rmulti **x)
{
  rmulti *a=NULL;
  RAP(a,y,n);
  rvec_sum_mul(a,n,x,y);    // a=x'*y
  rvec_sub_mul_r(n,y,x,a);  // y-=x*a
  RF(a);
}

/**
 @brief rmulti型のベクトルの2乗ノルム value=sqrt(sum(x.^2))
*/
void rvec_norm2(rmulti *value, int n, rmulti **x)
{
  rvec_sum_pow2(value,n,x); // value=sum(x.^2)
  rsqrt(value,value);       // value=sqrt(value)
}

/**
 @brief rmulti型のベクトルの要素の平均 value=sum(x)/n
*/
void rvec_average(rmulti *value, int n, rmulti **x)
{
  rvec_sum(value,n,x);
  rdiv_si2(value,value,n);
}

/**
 @brief rmulti型のベクトルの要素の絶対値の比の最大値 value=max(abs(x)./abs(y))
*/
void rvec_max_div_abs(rmulti *value, int n, rmulti **x, rmulti **y)
{
  int i=0;
  rmulti *a=NULL;
  RAp(a,value);
  rdiv(a,x[i],y[i]);   // value=x[0]/y[0]
  rabs(value,a);       // value=abs(x[0]/y[0])
  for(i=1; i<n; i++){
    rdiv(a,x[i],y[i]); // a=x[i]/y[i]
    rabs(a,a);         // a=abs(x[i]/y[i])
    if(rgt(a,value)){
      rcopy(value,a);     // value=a
    }
  }
  RF(a);
}

/**
 @brief rmulti型のベクトルの方向余弦 (x'*y)/sqrt(x'*x)/sqrt(y'*y)
*/
void rvec_dcos(rmulti *value, int n, rmulti **x, rmulti **y)
{
  rmulti *a=NULL,*b=NULL;
  RAp(a,value); RAp(b,value);
  rvec_norm2(a,n,x); // a=sqrt(x'*x)
  rvec_norm2(b,n,y); // b=sqrt(y'*y)
  rmul(a,a,b);       // a=sqrt(x'*x)*sqrt(y'*y)
  rvec_sum_mul(b,n,x,y); // b=x'*y
  rdiv(value,b,a);   // valu=(x'*y)/sqrt(x'*x)/sqrt(y'*y)
  RF(a); RF(b);
}

/**
 @brief rmulti型のベクトルの方向余弦の絶対値 abs(x'*y)/sqrt(x'*x)/sqrt(y'*y)
*/
void rvec_abs_dcos(rmulti *value, int n, rmulti **x, rmulti **y)
{
  rvec_dcos(value,n,x,y);
  rabs(value,value);
}

/**
 @brief rmulti型のベクトルの角度 theta=acos(abs(x'*y)/sqrt(x'*x)/sqrt(y'*y))
*/
void rvec_angle(rmulti *theta, int n, rmulti **x, rmulti **y)
{
  rmulti *a=NULL;
  RAp(a,theta);
  rvec_abs_dcos(a,n,x,y); // a=abs(x'*y)/sqrt(x'*x)/sqrt(y'*y)
  if(rgt_d2(a,1)){ rset_d(a,1); }
  racos(theta,a);         // theta=acos(dcos)
  RF(a);
}

/**
 @brief rmulti型のベクトルの角度[deg] theta=acos(abs(x'*y)/sqrt(x'*x)/sqrt(y'*y))
*/
void rvec_angle_deg(rmulti *theta, int n, rmulti **x, rmulti **y)
{
  rvec_angle(theta,n,x,y);
  rmul_d(theta,theta,(180.0/M_PI));
}

/** @} */


////////////////////////////////////////////////////////////////////////

/** @name rmulti型の値の比較に関する関数 */
/** @{ */

/**
 @brief rmulti型のベクトルの値の比較 X<=>Y
*/
int rvec_cmp(int n, rmulti **x, int m, rmulti **y)
{
  int i,value;
  for(i=0; i<MIN2(n,m); i++){
    value=rcmp(x[i],y[i]);
    if(value!=0){ return value; } // X!=Y
  }
  if     (n<m){ return -1; } // X<Y
  else if(n>m){ return +1; } // X>Y
  else        { return  0; } // X=Y
}

/**
 @brief rmulti型のベクトルの値の比較 X==Y
*/
int rvec_eq(int n, rmulti **x, rmulti **y)
{
  int i;
  for(i=0; i<n; i++){ if(!req(x[i],y[i])){ return 0; } }
  return 1;
}

/**
 @brief rmulti型のベクトルの値の比較 X>Y
*/
int rvec_gt(int n, rmulti **x, rmulti **y)
{
  int i;
  for(i=0; i<n; i++){ if(!rgt(x[i],y[i])){ return 0; } }
  return 1;
}

/**
 @brief rmulti型のベクトルの値の比較 X>=Y
*/
int rvec_ge(int n, rmulti **x, rmulti **y)
{
  int i;
  for(i=0; i<n; i++){ if(!rge(x[i],y[i])){ return 0; } }
  return 1;
}

/**
 @brief rmulti型のベクトルの値の比較 X<Y
*/
int rvec_lt(int n, rmulti **x, rmulti **y)
{
  int i;
  for(i=0; i<n; i++){ if(!rlt(x[i],y[i])){ return 0; } }
  return 1;
}

/**
 @brief rmulti型のベクトルの値の比較 X<=Y
*/
int rvec_le(int n, rmulti **x, rmulti **y)
{
  int i;
  for(i=0; i<n; i++){ if(!rle(x[i],y[i])){ return 0; } }
  return 1;
}

/**
 @brief rmulti型のベクトルの値の比較 X==Y
*/
int rvec_eq_d(int n, rmulti **x, double y)
{
  int i;
  for(i=0; i<n; i++){ if(!req_d(x[i],y)){ return 0; } }
  return 1;
}

/**
 @brief rmulti型のベクトルの値の比較 X>Y
*/
int rvec_gt_d2(int n, rmulti **x, double y)
{
  int i;
  for(i=0; i<n; i++){ if(!rgt_d2(x[i],y)){ return 0; } }
  return 1;
}

/**
 @brief rmulti型のベクトルの値の比較 X>=Y
*/
int rvec_ge_d2(int n, rmulti **x, double y)
{
  int i;
  for(i=0; i<n; i++){ if(!rge_d2(x[i],y)){ return 0; } }
  return 1;
}

/**
 @brief rmulti型のベクトルの値の比較 X<Y
*/
int rvec_lt_d2(int n, rmulti **x, double y)
{
  int i;
  for(i=0; i<n; i++){ if(!rlt_d2(x[i],y)){ return 0; } }
  return 1;
}

/**
 @brief rmulti型のベクトルの値の比較 X<=Y
*/
int rvec_le_d2(int n, rmulti **x, double y)
{
  int i;
  for(i=0; i<n; i++){ if(!rle_d2(x[i],y)){ return 0; } }
  return 1;
}

/** @} */


///////////////////////////////////////////////////////////////////////////////////////////////////////////

/** @name rmulti型の写像に関する関数 */
/** @{ */

/**
 @brief rmulti型のベクトルに関する写像 y=f(x)
*/
void rvec_func(rmulti *y, func_t *f, int n, rmulti **x)
{
  int i;
  rmulti *a=NULL,*b=NULL,*z=NULL;
  RAp(a,y); RAp(b,y); RAp(z,y);
  if(f==NULL)                { FUNC_ERROR_ARG1("rvec_func",f); }
  else if(func_is(f,"nan"))  { rset_nan(z); }
  else if(func_is(f,"inf"))  { rset_inf(z,1); }
  else if(func_is_zero(f))   { rset_si(z,0); }
  else if(func_is_one(f))    { rset_si(z,1); }
  else if(func_is_bigint(f)) { bigint_get_rmulti(z,func_bigint_p(f)); }
  else if(func_is_real(f))   { rcopy(z,func_real_p(f)); }
  else if(func_is_complex(f)){ rset_si(a,0); rset_si(b,0); rdiv(z,a,b); }
  else if(func_is_var(f))    {
    rset_d(z,1);
    for(i=0; i<func_var_size(f); i++){
      if(func_var_pow(f,i)!=0){
	if(0<=func_var_num(f,i) && func_var_num(f,i)<n){
	  rpow_si(a,x[func_var_num(f,i)],func_var_pow(f,i));
	}else{ rset_nan(a); }
	rmul(z,z,a);
      }
    }
  }
  else if(func_is_add(f))    { rset_d(z,0); for(i=0; i<func_asize(f); i++){ rvec_func(a,func_aget(f,i),n,x); radd(z,z,a); } }
  else if(func_is_mul(f))    { rset_d(z,1); for(i=0; i<func_asize(f); i++){ rvec_func(a,func_aget(f,i),n,x); rmul(z,z,a); } }
  else if(func_is(f,"sqrt")) { rvec_func(a,func_aget(f,0),n,x); rsqrt(z,a); }
  else if(func_is(f,"exp"))  { rvec_func(a,func_aget(f,0),n,x); rexp(z,a); }
  else if(func_is(f,"log"))  { rvec_func(a,func_aget(f,0),n,x); rlog(z,a); }
  else if(func_is(f,"sin"))  { rvec_func(a,func_aget(f,0),n,x); rsin(z,a); }
  else if(func_is(f,"cos"))  { rvec_func(a,func_aget(f,0),n,x); rcos(z,a); }
  else if(func_is(f,"tan"))  { rvec_func(a,func_aget(f,0),n,x); rtan(z,a); }
  else if(func_is(f,"asin")) { rvec_func(a,func_aget(f,0),n,x); rasin(z,a); }
  else if(func_is(f,"acos")) { rvec_func(a,func_aget(f,0),n,x); racos(z,a); }
  else if(func_is(f,"atan")) { rvec_func(a,func_aget(f,0),n,x); ratan(z,a); }
  else if(func_is(f,"sinh")) { rvec_func(a,func_aget(f,0),n,x); rsinh(z,a); }
  else if(func_is(f,"cosh")) { rvec_func(a,func_aget(f,0),n,x); rcosh(z,a); }
  else if(func_is(f,"tanh")) { rvec_func(a,func_aget(f,0),n,x); rtanh(z,a); }
  else if(func_is(f,"asinh")){ rvec_func(a,func_aget(f,0),n,x); rasinh(z,a); }
  else if(func_is(f,"acosh")){ rvec_func(a,func_aget(f,0),n,x); racosh(z,a); }
  else if(func_is(f,"atanh")){ rvec_func(a,func_aget(f,0),n,x); ratanh(z,a); }
  else if(func_is(f,"pow"))  { rvec_func(a,func_aget(f,0),n,x); rvec_func(b,func_aget(f,1),n,x); rpow(z,a,b); }
  else if(func_is_list(f))   { rset_nan(a); }
  else if(func_is_rvec(f))   { rset_nan(a); }
  else if(func_is_cvec(f))   { rset_nan(a); }
  else if(func_is_rmat(f))   { rset_nan(a); }
  else if(func_is_cmat(f))   { rset_nan(a); }
  else                       { rset_nan(a); }
  if(func_has_power(f))      { rpow_si(z,z,func_power(f)); }
  rcopy(y,z);
  RF(a); RF(b); RF(z);
}

/**
 @brief rmulti型のベクトルに関するベクトル写像 y=f(x) 
*/
void rvec_func_list(int m, rmulti **y, func_t *f, int n, rmulti **x)
{
  int i;
  for(i=0; i<m; i++){
    if(func_is_list(f) && i<func_asize(f)){
      rvec_func(y[i],func_aget(f,i),n,x);
    }else{
      rset_nan(y[i]);
    }
  }
}

/**
 @brief rmulti型に関する写像 y=f(x0)
*/
void r1_func(rmulti *y, func_t *f, rmulti *x0)
{
  int n=1;
  rmulti **x=NULL;
  x=rvec_allocate_prec(n,rget_prec(y));
  rclone(x[0],x0);
  rvec_func(y,f,n,x);
  x=rvec_free(n,x);
}

/**
 @brief rmulti型に関する写像 y=f(x0,x1)
*/
void r2_func(rmulti *y, func_t *f, rmulti *x0, rmulti *x1)
{
  int n=2;
  rmulti **x=NULL;
  x=rvec_allocate_prec(n,rget_prec(y));
  rclone(x[0],x0);
  rclone(x[1],x1);
  rvec_func(y,f,n,x);
  x=rvec_free(n,x);
}

/**
 @brief rmulti型に関する写像 y=f(x0,x1,x2)
*/
void r3_func(rmulti *y, func_t *f, rmulti *x0, rmulti *x1, rmulti *x2)
{
  int n=3;
  rmulti **x=NULL;
  x=rvec_allocate_prec(n,rget_prec(y));
  rclone(x[0],x0);
  rclone(x[1],x1);
  rclone(x[2],x2);
  rvec_func(y,f,n,x);
  x=rvec_free(n,x);
}

/** @} */



//EOF
