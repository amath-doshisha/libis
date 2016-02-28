#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<stdarg.h>
#include"is_macros.h"
#include"is_rmulti.h"
#include"is_cmulti.h"
#include"is_cvec.h"
#include"is_rvec.h"
#include"is_zvec.h"
#include"is_dvec.h"
#include"is_ivec.h"
#include"is_func.h"


/**
 @file  cvec.c
 @brief 多倍長精度複素数型cmultiのベクトルに関する関数の定義.
 @details cmulti型の関数に関する定義は@link cmulti.c@endlinkを参照のこと.
          cmulti型の行列の関数に関する定義は@link cmat.c@endlinkを参照のこと.
 */

/////////////////////////////////////////////////////

#define FR(F)     (func_retain(F))

/////////////////////////////////////////////////////

#define FILE_NAME_LENGTH_MAX 100000

/////////////////////////////////////////////////////

/** @name cmulti型のベクトルの初期化に関する関数 */
/** @{ */

/**
 @brief cmulti型のベクトルの新規生成.
*/
cmulti **cvec_allocate(int n)
{
  int i;
  cmulti **x=NULL;
  x=(cmulti**)malloc(sizeof(cmulti*)*n);
  for(i=0; i<n; i++){ x[i]=callocate(); }
  return x;
}

/**
 @brief cmulti型のベクトルの精度を指定しての新規生成.
*/
cmulti **cvec_allocate_prec(int n, int prec)
{
  int i;
  cmulti **x=NULL;
  x=(cmulti**)malloc(sizeof(cmulti*)*n);
  for(i=0; i<n; i++){ x[i]=callocate_prec(prec); }
  return x;
}

/**
 @brief cmulti型のベクトルの複製の生成.
*/
cmulti **cvec_allocate_clone(int n, cmulti **y)
{
  int i;
  cmulti **x;
  x=(cmulti**)malloc(sizeof(cmulti*)*n);
  for(i=0; i<n; i++){ x[i]=callocate_clone(y[i]); }
  return x;
}

/**
 @brief cmulti型のベクトルのrmulti型ベクトルからの複製の生成.
*/
cmulti **cvec_allocate_clone_r(int n, rmulti **y)
{
  int i;
  cmulti **x;
  x=(cmulti**)malloc(sizeof(cmulti*)*n);
  for(i=0; i<n; i++){ x[i]=callocate_clone_r(y[i]); }
  return x;
}

/**
 @brief cmulti型のベクトルの終了処理.
*/
cmulti **cvec_free(int n, cmulti **x)
{
  int i;
  if(x==NULL) return NULL;
  for(i=0; i<n; i++){ x[i]=cfree(x[i]); }
  free(x);
  x=NULL;
  return x;
}

/**
 @brief 初期化済みのcmulti型のベクトルの浮動小数点数の精度(ビット数)を変更し再初期化.
*/
int cvec_round(int n, cmulti **x, int prec)
{
  int i,e=0;
  for(i=0; i<n; i++){ e+=cround(x[i],prec); }
  return e;
}

/**
 @brief cmulti型のベクトルの値を複成 y=x.
*/
int cvec_clone(int n, cmulti **y, cmulti **x)
{
  int i,e=0;
  for(i=0; i<n; i++){ e+=cclone(y[i],x[i]); }
  return e;
}

/**
 @brief cmulti型のベクトルの値を添字を指定して複成 y=x(I).
*/
int cvec_clone_index(int n, cmulti **y, cmulti **x, int *I)
{
  int i,e=0;
  for(i=0; i<n; i++){ e+=cclone(y[i],x[I[i]]); }
  return e;
}

/**
 @brief cmulti型のベクトルの値を添字を指定して複成 y(I)=x.
*/
int cvec_index_clone(int n, cmulti **y, cmulti **x, int *I)
{
  int i,e=0;
  for(i=0; i<n; i++){ e+=cclone(y[I[i]],x[i]); }
  return e;
}

/**
 @brief cmulti型のベクトルの値をrmulti型のベクトルから複成 y=x.
*/
int cvec_clone_r(int n, cmulti **y, rmulti **x)
{
  int i,e=0;
  for(i=0; i<n; i++){ e+=cclone_r(y[i],x[i]); }
  return e;
}

/**
 @brief cmulti型のベクトルの値をrmulti型の2つのベクトルから複成 y=x_r+i*x_i.
*/
int cvec_clone_rr(int n, cmulti **y, rmulti **x_r, rmulti **x_i)
{
  int i,e=0;
  for(i=0; i<n; i++){ e+=cclone_rr(y[i],x_r[i],x_i[i]); }
  return e;
}

/**
 @brief cmulti型のベクトルの値の交換 x <=> y.
*/
void cvec_swap(int n, cmulti **x, cmulti **y)
{
  int i;
  for(i=0; i<n; i++){ cswap(x[i],y[i]); }
}

/** @} */

////////////////////////////////////////////////////////////////////////

/** @name cmulti型のベクトルのメンバ変数に関する関数 */
/** @{ */

/**
 @brief cmulti型のベクトルの浮動小数点数の精度(ビット数)を取得.
*/
void cvec_get_prec(int n, int *p, cmulti **x)
{
  int i;
  for(i=0; i<n; i++){ p[i]=cget_prec(x[i]); }
}

/**
 @brief cmulti型のベクトルの浮動小数点数の精度(ビット数)の最大値を取得.
*/
int cvec_get_prec_max(int n, cmulti **x)
{
  int value,p,i=0;
  value=cget_prec(x[i]);
  for(i=1; i<n; i++){
    p=cget_prec(x[i]);
    if(p>=value){ value=p; }
  }
  return value;
}

/**
 @brief cmulti型のベクトルの浮動小数点数の指数部の取得.
*/
void cvec_get_exp(int n, int *p, cmulti **x)
{
  int i;
  for(i=0; i<n; i++){ p[i]=cget_exp(x[i]); }
}

/**
 @brief cmulti型のベクトルの浮動小数点数の指数部の最大値の取得.
*/
int cvec_get_exp_max(int n, cmulti **x)
{
  int value,p,i=0;
  value=cget_exp(x[i]);
  for(i=1; i<n; i++){
    p=cget_exp(x[i]);
    if(p>=value){ value=p; }
  }
  return value;
}

/**
 @brief cmulti型のベクトルが数であるかの判定.
*/
int cvec_is_number(int n, cmulti **x)
{
  int i;
  for(i=0; i<n; i++){ if(!cis_number(x[i])){ return 0; } }
  return 1;
}

/**
 @brief cmulti型のベクトルがNaNを含むかの判定.
*/
int cvec_has_nan(int n, cmulti **x)
{
  int i;
  for(i=0; i<n; i++){ if(cis_nan(x[i])){ return 1; } }
  return 0;
}

/** @} */

/////////////////////////////////////////////////////

/** @name cmulti型のベクトルの入出力に関する関数 */
/** @{ */

/**
 @brief cmulti型のベクトルの表示.
*/
void cvec_print(int n, cmulti **x, const char *name, const char *f, int digits)
{
  int i,k;
  char format[128];
  if(STR_EQ(f,"f") || STR_EQ(f,"F")){ k=3; }
  else if((strcmp(f,"e")==0) || (strcmp(f,"E")==0)){ k=9; }
  else{ k=3; }
  sprintf(format,"%%%d.%dR%s %%%d.%dR%s\n",digits+k,digits,f,digits+k,digits,f);
  if(x==NULL){
    if(name!=NULL){ printf("%s=NULL\n",name); }
    else          { printf("NULL\n"); }
    return;
  }
  if(name!=NULL){ printf("%s\n",name); }
  if(x==NULL) return;
  for(i=0; i<n; i++){ mpfr_printf(format,C_R(x[i]),C_I(x[i])); }
}

/**
 @brief cmulti型のベクトルの表示.
*/
void cvec_print_prec(int n, cmulti **x, const char *name, const char *f, int digits)
{
  int i,k;
  char format[128];
  if(STR_EQ(f,"f") || STR_EQ(f,"F")){ k=3; }
  else if((strcmp(f,"e")==0) || (strcmp(f,"E")==0)){ k=8; }
  else{ k=3; }
  sprintf(format,"%%%d.%dR%s %%%d.%dR%s [%%d] [%%d]\n",digits+k,digits,f,digits+k,digits,f);
  if(x==NULL){
    if(name!=NULL){ printf("%s=NULL\n",name); }
    else          { printf("NULL\n"); }
    return;
  }
  if(name!=NULL){ printf("%s\n",name); }
  if(x==NULL) return;
  for(i=0; i<n; i++){
    mpfr_printf(format,C_R(x[i]),C_I(x[i]),rget_prec(C_R(x[i])),rget_prec(C_I(x[i])));
  }
}

/**
 @brief cmulti型のベクトルの表示.
*/
void cvec_print_exp(int n, cmulti **x, const char *name)
{
  int i;
  if(x==NULL){
    printf("%s=NULL\n",name);
    return;
  }
  printf("%s\n",name);
  for(i=0; i<n; i++){
    if(ris_zero(C_R(x[i]))){
      mpfr_printf("-Inf ");
    }else{
      mpfr_printf("%4d ",rget_exp(C_R(x[i])));
    }
    if(ris_zero(C_I(x[i]))){
      mpfr_printf("-Inf\n");
    }else{
      mpfr_printf("%4d\n",rget_exp(C_I(x[i])));
    }
  }
}


/**
 @brief cmulti型のベクトルの保存.
*/
void cvec_save(int n, cmulti **x, int offset, int digits, char* fmt, ...)
{
  int i;
  char fname[FILE_NAME_LENGTH_MAX];
  va_list argp;
  FILE *fid;
  char format[128];
  // file name
  va_start(argp,fmt); vsprintf(fname,fmt,argp);
  // open file
  if((fid=fopen(fname,"w"))==0){ ERROR_AT; printf("Cant't open file: %s\n",fname); exit(0); }
  // write
  sprintf(format,"%%d %%+.%dRe %%+.%dRe\n",digits,digits);
  for(i=0;i<n;i++){ if(mpfr_fprintf(fid,format,i+offset,C_R(x[i]),C_I(x[i]))<0){ ERROR_AT; exit(0); } }
  // close
  fclose(fid);
}

/**
 @brief cmulti型のベクトルの保存.
*/
void cvec_save_cplane(int n, cmulti **x, int digits, char* fmt, ...)
{
  int i;
  char fname[FILE_NAME_LENGTH_MAX];
  va_list argp;
  FILE *fid;
  char format[128];
  // file name
  va_start(argp,fmt); vsprintf(fname,fmt,argp);
  // open file
  if((fid=fopen(fname,"w"))==0){ ERROR_AT; printf("Cant't open file: %s\n",fname); exit(0); }
  // write
  sprintf(format,"%%+.%dRe %%+.%dRe\n",digits,digits);
  for(i=0;i<n;i++){ if(mpfr_fprintf(fid,format,C_R(x[i]),C_I(x[i]))<0){ ERROR_AT; exit(0); } }
  // close
  fclose(fid);
}

/**
 @brief cmulti型のベクトルの読み込み.
*/
void cvec_load(int n, cmulti **x, char* fmt, ...)
{
  int i;
  char fname[FILE_NAME_LENGTH_MAX];
  va_list argp;
  FILE *fid;
  char str1[100000],str2[100000];
  // file name
  va_start(argp,fmt);
  vsprintf(fname,fmt,argp);
  // open file
  if((fid=fopen(fname,"r"))==0){ ERROR_EXIT("Cant't open file: %s\n",fname); }
  // real
  for(i=0;i<n;i++){
    if(fscanf(fid,"%s %s",str1,str2)==EOF){ ERROR_EXIT("Error in cvec_read(), fscanf(), fname=%s\n",fname); }
    cset_ss(x[i],str1,str2);
  }
  // close
  fclose(fid);
  fprintf(stderr,"loaded: %s\n",fname);
}

/**
 @brief cmulti型のベクトルの保存.
*/
void cvec_bin_save(int n, cmulti **x, char* fmt, ...)
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
    fwrite("cvec",sizeof(char),strlen("cvec"),fid);
    // write size
    fwrite(&n,sizeof(int),1,fid);
    // write data
    for(i=0; i<n; i++){ cbin_save(x[i],fid); }
  }
  // close
  fclose(fid);
}

/**
 @brief cmulti型のベクトルの読み込み.
*/
cmulti **cvec_bin_load(int *n, char* fmt, ...)
{
  int l,i;
  size_t k;
  cmulti **cx=NULL;
  rmulti **rx=NULL;
  dcomplex *zx=NULL;
  double *dx=NULL;
  char fname[FILE_NAME_LENGTH_MAX+1],*buf=NULL;
  va_list argp;
  FILE *fid;
  // file name
  va_start(argp,fmt); vsprintf(fname,fmt,argp);
  // open
  if((fid=fopen(fname,"r"))==0){ fclose(fid); cx=NULL; (*n)=0; }
  else{
    // read header
    l=strlen("cvec");
    buf=(char*)malloc(sizeof(char)*(l+1));
    k=fread(buf,sizeof(char),l,fid);
    if(k==(size_t)l && strncmp(buf,"cvec",l)==0){ /* cvec */
      k=fread(n,sizeof(int),1,fid);                  // read size
      if(k!=1 || (*n)<=0){ ERROR_AT; printf("Failed to load the size from the file '%s'.\n",fname); exit(0); }
      cx=(cmulti**)malloc(sizeof(cmulti*)*(*n));     // allocate
      for(i=0; i<(*n); i++){ cx[i]=cbin_load(fid); } // read data
      fclose(fid);                                   // close
    }else if(k==(size_t)l && strncmp(buf,"rvec",l)==0){ /* rvec */
      fclose(fid);                       // close
      rx=rvec_bin_load(n,fname);         // load
      cx=cvec_allocate_clone_r((*n),rx); // allocate
      rx=rvec_free((*n),rx);             // free
    }else if(k==(size_t)l && strncmp(buf,"zvec",l)==0){ /* zvec */
      fclose(fid);                    // close
      zx=zvec_bin_load(n,fname);      // load
      cx=cvec_allocate_prec((*n),53); // allocate
      cvec_set_z((*n),cx,zx);         // copy
      zx=zvec_free(zx);               // free
    }else if(k==(size_t)l && strncmp(buf,"dvec",l)==0){ /* dvec */
      fclose(fid);                    // close
      dx=dvec_bin_load(n,fname);      // load
      cx=cvec_allocate_prec((*n),53); // allocate
      cvec_set_d((*n),cx,dx);         // copy
      dx=dvec_free(dx);               // free
    }else{ fclose(fid); cx=NULL; (*n)=0; }
  }
  // done
  free(buf);
  return cx;
}

/** @} */

/////////////////////////////////////////////////////

/** @name cmulti型のベクトルの値の設定に関する関数 */
/** @{ */

/**
 @brief cmulti型のベクトルの値をNaNに設定.
*/
void cvec_set_nan(int n, cmulti **x)
{
  int i;
  for(i=0; i<n; i++){ cset_nan(x[i]); }
}

/**
 @brief cmulti型のベクトルの値を文字列から設定.
*/
void cvec_set_s(int n, cmulti **x, const char **str)
{
  int i;
  for(i=0; i<n; i++){ cset_s(x[i],str[i]); }
}

/**
 @brief cmulti型のベクトルの値を文字列から設定.
*/
void cvec_set_ss(int n, cmulti **x, const char **str)
{
  int i;
  for(i=0; i<n; i++){ cset_ss(x[i],str[2*i],str[2*i+1]); }
}

/**
 @brief cmulti型のベクトルの値を倍精度複素数から設定.
*/
int cvec_set_z(int n, cmulti **y, dcomplex *x)
{
  int i,e=0;
  for(i=0; i<n; i++){ e+=cset_z(y[i],x[i]); }
  return e;
}

/**
 @brief cmulti型のベクトルの値を倍精度実数から設定.
*/
int cvec_set_d(int n, cmulti **y, double *x)
{
  int i,e=0;
  for(i=0; i<n; i++){ e+=cset_d(y[i],x[i]); }
  return e;
}

/**
 @brief cmulti型のベクトルの値を倍精度浮実数から設定.
*/
int cvec_set_dd(int n, cmulti **y, double *x_r, double *x_i)
{
  int i,e=0;
  for(i=0; i<n; i++){ e+=cset_dd(y[i],x_r[i],x_i[i]); }
  return e;
}

/**
 @brief cmulti型のベクトルの全ての値を倍精度複素数から設定.
*/
int cvec_set_all_z(int n, cmulti **x, dcomplex a)
{
  int i,e=0;
  for(i=0; i<n; i++){ e+=cset_z(x[i],a); }
  return e;
}

/**
 @brief cmulti型のベクトルの全ての値を倍精度実数から設定.
*/
int cvec_set_all_d(int n, cmulti **x, double a)
{
  int i,e=0;
  for(i=0; i<n; i++){ e+=cset_d(x[i],a); }
  return e;
}

/**
 @brief cmulti型のベクトルの全ての値を倍精度実数から設定.
*/
int cvec_set_all_dd(int n, cmulti **x, double a_r, double a_i)
{
  int i,e=0;
  for(i=0; i<n; i++){ e+=cset_dd(x[i],a_r,a_i); }
  return e;
}

/**
 @brief cmulti型のベクトルの値を零に設定.
*/
int cvec_set_zeros(int n, cmulti **x)
{
  return cvec_set_all_d(n,x,0);
}

/**
 @brief cmulti型のベクトルの値を1に設定.
*/
int cvec_set_ones(int n, cmulti **x)
{
  return cvec_set_all_d(n,x,1);
}

/**
 @brief cmulti型のベクトルの値を第k要素が1の基本ベクトルに設定.
*/
int cvec_set_unit(int n, cmulti **x, int k)
{
  int e=0;
  e+=cvec_set_zeros(n,x);
  e+=cset_d(x[k],1);
  return e;
}

/**
 @brief cmulti型のベクトルの値を0,1,2,...に設定.
*/
int cvec_set_grid(int n, cmulti **x)
{
  int i,e=0;
  for(i=0; i<n; i++){ e+=cset_d(x[i],i); }
  return e;
}

/**
 @brief cmulti型のベクトルの値を区間(b,a+b)の疑似乱数値を設定.
*/
void cvec_set_rand(int n, cmulti **x, double a, double b)
{
  int i;
  for(i=0; i<n; i++){
    cset_rand(x[i]);
    cmul_d(x[i],x[i],a);
    cadd_d(x[i],x[i],b);
  }
}

/** @} */

////////////////////////////////////////////////////////////////////////

/** @name cmulti型のベクトルの型変換に関する関数 */
/** @{ */

/**
 @brief cmulti型のベクトルを倍精度複素数型に変換.
 */
void cvec_get_z(int n, dcomplex *y, cmulti **x)
{
  int i;
  for(i=0; i<n; i++){ y[i]=cget_z(x[i]); }
}

/**
 @brief cmulti型のベクトルを倍精度実数型に変換.
 */
void cvec_get_d(int n, double *y, cmulti **x)
{
  int i;
  for(i=0; i<n; i++){ y[i]=cget_d(x[i]); }
}

/** @} */

////////////////////////////////////////////////////////////////////

/** @name cmulti型のベクトルの要素の配置に関する関数 */
/** @{ */

/**
 @brief cmulti型のベクトルの要素の順を逆に入れ替える.
 */
void cvec_reverse(int n, cmulti **x)
{
  int i;
  for(i=0; i<n/2; i++){
    cswap(x[i],x[n-i-1]);
  }
}

/**
 @brief cmulti型のベクトルの第i要素と第j要素を入れ替える.
 */
void cvec_swap_at(cmulti **x, int i, int j)
{
  cswap(x[i],x[j]);
}

/**
 @brief cmulti型のベクトルの要素を要素番号の配列Iの値に従って入れ替える.
 */
void cvec_swap_index(int n, cmulti **x, const int *I)
{
  int i;
  cmulti **y=NULL;
  y=cvec_allocate_prec(n,53);
  for(i=0; i<n; i++){
    cswap(y[i],x[I[i]]);
  }
  cvec_swap(n,x,y);
  y=cvec_free(n,y);
}

/**
 @brief cmulti型のベクトルの要素を小さい順に入れ替える.
 @details 要素x[i]の値の小さい順に並べ替え，要素番号を配列Iに保存する.
          ただし，I==NULLの場合は保存は行われない.
 @param[in]  n ベクトルのサイズ.
 @param[in]  x 初期化済みのrmulti型のベクトル.
 @param[out] x ソートされたrmulti型のベクトル.
 @param[in]  I 入れ替えられた要素の要素番号の保存先.
 @param[out] I 保存された要素番号.
*/
void cvec_sort(int n, cmulti **x, int *I)
{
  if(I!=NULL){ ivec_grid(n,I); }
  cvec_quick_sort(n,x,I,0,n-1);
}

/**
 @brief cmulti型のベクトルの要素に関するクイックソートアルゴリズム(ユーザは直接使用してはならない).
*/
void cvec_quick_sort(int n, cmulti **x, int *I, int left, int right)
{
  int i,last;
  if(left>=right) return;
  cswap(x[left],x[(left+right)/2]);
  if(I!=NULL) ivec_swap_at(I,left,(left+right)/2);
  last=left;
  for(i=left+1; i<=right; i++){
    if(clt(x[i],x[left])){
      ++last;
      cswap(x[last],x[i]);
      if(I!=NULL) ivec_swap_at(I,last,i);      
    }
  }
  cswap(x[left],x[last]);
  if(I!=NULL) ivec_swap_at(I,left,last);
  cvec_quick_sort(n,x,I,left,last-1);
  cvec_quick_sort(n,x,I,last+1,right);
}

/**
 @brief cmulti型のベクトルの要素を小さい順に入れ替えたときの要素番号の配列の取得.
*/
void cvec_sort_index(int *I, int n, cmulti **X)
{
  cmulti **Y=NULL;
  Y=cvec_allocate_clone(n,X);
  cvec_sort(n,Y,I);
  Y=cvec_free(n,Y);
}


/** @} */

////////////////////////////////////////////////////////////////////

/** @name cmulti型のベクトルの自動精度調整モードが機能する関数 */
/** @{ */

/**
 @brief cmulti型のベクトルの値のコピー y=x.
*/
int cvec_copy(int n, cmulti **y, cmulti **x)
{
  int i,e=0;
  if(y==x){ return 0; }
  for(i=0; i<n; i++){ e+=ccopy(y[i],x[i]); }
  return e;
}

/**
 @brief cmulti型のベクトルの値のrmulti型からのコピー y=x.
*/
int cvec_copy_r(int n, cmulti **y, rmulti **x)
{
  int i,e=0;
  for(i=0; i<n; i++){ e+=ccopy_r(y[i],x[i]); }
  return e;
}

/**
 @brief cmulti型のベクトルの値のrmulti型からのコピー y=x.
*/
int cvec_copy_rr(int n, cmulti **y, rmulti **x_r, rmulti **x_i)
{
  int i,e=0;
  for(i=0; i<n; i++){ e+=ccopy_rr(y[i],x_r[i],x_i[i]); }
  return e;
}

/**
 @brief cmulti型のベクトルの値の要素番号Iに従ったコピー y=x.
 @details y[i]=x[I[i]], i=0,1,2,..,n-1.
*/
int cvec_copy_index(int n, cmulti **y, cmulti **x, const int *I)
{
  int i,e=0;
  if(y==x){ cvec_swap_index(n,x,I); return 0; }
  for(i=0; i<n; i++){ e+=ccopy(y[i],x[I[i]]); }
  return e;
}

/**
 @brief cmulti型のベクトルの値を添字を指定してコピー y(I)=x.
*/
int cvec_index_copy(int n, cmulti **y, cmulti **x, int *I)
{
  int i,e=0;
  for(i=0; i<n; i++){ e+=ccopy(y[I[i]],x[i]); }
  return e;
}

/**
 @brief cmulti型の実部のコピー y=real(x)
 */
int cvec_real(int n, rmulti **y, cmulti **x)
{
  int i,e=0;
  for(i=0; i<n; i++){ e+=rcopy(y[i],C_R(x[i])); }
  return e;
}

/**
 @brief cmulti型の実部の複製 y=real(x)
 */
int cvec_real_clone(int n, rmulti **y, cmulti **x)
{
  int i,e=0;
  for(i=0; i<n; i++){ e+=rclone(y[i],C_R(x[i])); }
  return e;
}

/**
 @brief cmulti型の虚部のコピー y=imag(x)
 */
int cvec_imag(int n, rmulti **y, cmulti **x)
{
  int i,e=0;
  for(i=0; i<n; i++){ e+=rcopy(y[i],C_I(x[i])); }
  return e;
}

/**
 @brief cmulti型の虚部の複製 y=imag(x)
 */
int cvec_imag_clone(int n, rmulti **y, cmulti **x)
{
  int i,e=0;
  for(i=0; i<n; i++){ e+=rclone(y[i],C_I(x[i])); }
  return e;
}

/**
 @brief cmulti型の複素共役のコピー y=imag(x)
 */
int cvec_conj(int n, cmulti **y, cmulti **x)
{
  int e=0,i;
  for(i=0; i<n; i++){ e+=cconj(y[i],x[i]); }
  return e;
}

/**
 @brief cmulti型の複素共役の複製 y=imag(x)
 */
int cvec_conj_clone(int n, cmulti **y, cmulti **x)
{
  int e=0,i;
  for(i=0; i<n; i++){ e+=cconj_clone(y[i],x[i]); }
  return e;
}


/**
 @brief cmulti型のベクトルの指数部の足し算 y=x*2^p
*/
int cvec_mul_2exp(int n, cmulti **y, cmulti **x, int pr, int pi)
{
  int i,e=0;
  for(i=0; i<n; i++){ e+=cmul_2exp(y[i],x[i],pr,pi); }
  return e;
}

/**
 @brief cmulti型のベクトルの指数部の引き算 y=x/2^p
*/
int cvec_div_2exp(int n, cmulti **y, cmulti **x, int pr, int pi)
{
  int i,e=0;
  for(i=0; i<n; i++){ e+=cdiv_2exp(y[i],x[i],pr,pi); }
  return e;
}

/**
 @brief cmulti型のベクトルの符号反転 y=-x
*/
int cvec_neg(int n, cmulti **y, cmulti **x)
{
  int i,e=0;
  for(i=0; i<n; i++){ e+=cneg(y[i],x[i]); }
  return e;
}

/**
 @brief cmulti型のベクトルの要素の実部，虚部の絶対値 y=abs(real(x)+i*imag(x))
*/
int cvec_absc(int n, cmulti **y, cmulti **x)
{
  int i,e=0;
  for(i=0; i<n; i++){ e+=cabsc(y[i],x[i]); }
  return e;
}

/**
 @brief cmulti型のベクトルの足し算 z=x+y
*/
int cvec_add(int n, cmulti **z, cmulti **x, cmulti **y)
{
  int i,e=0;
  for(i=0; i<n; i++){ e+=cadd(z[i],x[i],y[i]); }
  return e;
}

/**
 @brief cmulti型のベクトルの足し算 z=x+y
*/
int cvec_add_rvec(int n, cmulti **z, cmulti **x, rmulti **y)
{
  int i,e=0;
  for(i=0; i<n; i++){ e+=cadd_r(z[i],x[i],y[i]); }
  return e;
}

/**
 @brief cmulti型のベクトルの足し算 z=x+y
*/
int cvec_add_r(int n, cmulti **z, cmulti **x, rmulti *y)
{
  int i,e=0;
  for(i=0; i<n; i++){ e+=cadd_r(z[i],x[i],y); }
  return e;
}

/**
 @brief cmulti型のベクトルの足し算 z=x+y
*/
int cvec_add_z(int n, cmulti **z, cmulti **x, dcomplex y)
{
  int i,e=0;
  for(i=0; i<n; i++){ e+=cadd_z(z[i],x[i],y); }
  return e;
}

/**
 @brief cmulti型のベクトルの足し算 z=x+y
*/
int cvec_add_d(int n, cmulti **z, cmulti **x, double y)
{
  int i,e=0;
  for(i=0; i<n; i++){ e+=cadd_d(z[i],x[i],y); }
  return e;
}

/**
 @brief cmulti型のベクトルの引き算 z=x-y
*/
int cvec_sub(int n, cmulti **z, cmulti **x, cmulti **y)
{
  int i,e=0;
  for(i=0; i<n; i++){ e+=csub(z[i],x[i],y[i]); }
  return e;
}

/**
 @brief cmulti型のベクトルの引き算 z=x-y
*/
int cvec_sub_rvec1(int n, cmulti **z, rmulti **x, cmulti **y)
{
  int i,e=0;
  for(i=0; i<n; i++){ e+=csub_r1(z[i],x[i],y[i]); }
  return e;
}

/**
 @brief cmulti型のベクトルの引き算 z=x-y
*/
int cvec_sub_rvec2(int n, cmulti **z, cmulti **x, rmulti **y)
{
  int i,e=0;
  for(i=0; i<n; i++){ e+=csub_r2(z[i],x[i],y[i]); }
  return e;
}

/**
 @brief cmulti型のベクトルの引き算 z=x-y
*/
int cvec_sub_c1(int n, cmulti **z, cmulti *x, cmulti **y)
{
  int i,e=0;
  for(i=0; i<n; i++){ e+=csub(z[i],x,y[i]); }
  return e;
}

/**
 @brief cmulti型のベクトルの引き算 z=x-y
*/
int cvec_sub_c2(int n, cmulti **z, cmulti **x, cmulti *y)
{
  int i,e=0;
  for(i=0; i<n; i++){ e+=csub(z[i],x[i],y); }
  return e;
}

/**
 @brief cmulti型のベクトルの引き算 z=x-y
*/
int cvec_sub_r1(int n, cmulti **z, rmulti *x, cmulti **y)
{
  int i,e=0;
  for(i=0; i<n; i++){ e+=csub_r1(z[i],x,y[i]); }
  return e;
}

/**
 @brief cmulti型のベクトルの引き算 z=x-y
*/
int cvec_sub_r2(int n, cmulti **z, cmulti **x, rmulti *y)
{
  int i,e=0;
  for(i=0; i<n; i++){ e+=csub_r2(z[i],x[i],y); }
  return e;
}

/**
 @brief cmulti型のベクトルの引き算 z=x-y
*/
int cvec_sub_z1(int n, cmulti **z, dcomplex x, cmulti **y)
{
  int i,e=0;
  for(i=0; i<n; i++){ e+=csub_z1(z[i],x,y[i]); }
  return e;
}

/**
 @brief cmulti型のベクトルの引き算 z=x-y
*/
int cvec_sub_z2(int n, cmulti **z, cmulti **x, dcomplex y)
{
  int i,e=0;
  for(i=0; i<n; i++){ e+=csub_z2(z[i],x[i],y); }
  return e;
}

/**
 @brief cmulti型のベクトルの引き算 z=x-y
*/
int cvec_sub_d1(int n, cmulti **z, double x, cmulti **y)
{
  int i,e=0;
  for(i=0; i<n; i++){ e+=csub_d1(z[i],x,y[i]); }
  return e;
}

/**
 @brief cmulti型のベクトルの引き算 z=x-y
*/
int cvec_sub_d2(int n, cmulti **z, cmulti **x, double y)
{
  int i,e=0;
  for(i=0; i<n; i++){ e+=csub_d2(z[i],x[i],y); }
  return e;
}

/**
 @brief cmulti型のベクトルの要素ごとの掛け算 z=x.*y
*/
int cvec_mul(int n, cmulti **z, cmulti **x, cmulti **y)
{
  int i,e=0;
  for(i=0; i<n; i++){ e+=cmul(z[i],x[i],y[i]); }
  return e;
}

/**
 @brief cmulti型のベクトルxとスカラーyの掛け算 z=x*y
*/
int cvec_mul_c(int n, cmulti **z, cmulti **x, cmulti *y)
{
  int i,e=0;
  for(i=0; i<n; i++){ e+=cmul(z[i],x[i],y); }
  return e;
}

/**
 @brief cmulti型のベクトルxとスカラーyの掛け算 z=x*y
*/
int cvec_mul_z(int n, cmulti **z, cmulti **x, dcomplex y)
{
  int i,e=0;
  for(i=0; i<n; i++){ e+=cmul_z(z[i],x[i],y); }
  return e;
}

/**
 @brief cmulti型のベクトルxとスカラーyの掛け算 z=x*y
*/
int cvec_mul_r(int n, cmulti **z, cmulti **x, rmulti *y)
{
  int i,e=0;
  for(i=0; i<n; i++){ e+=cmul_r(z[i],x[i],y); }
  return e;
}

/**
 @brief cmulti型のベクトルxとスカラーyの掛け算 z=x*y
*/
int cvec_mul_d(int n, cmulti **z, cmulti **x, double y)
{
  int i,e=0;
  for(i=0; i<n; i++){ e+=cmul_d(z[i],x[i],y); }
  return e;
}

/**
 @brief cmulti型のベクトルの要素ごとの掛け算の加算 z+=x.*y
*/
int cvec_add_mul(int n, cmulti **z, cmulti **x, cmulti **y)
{
  int i,e=0;
  for(i=0; i<n; i++){ e+=cadd_mul(z[i],x[i],y[i]); }
  return e;
}

/**
 @brief cmulti型のベクトルxとスカラーyの掛け算の加算 z+=x*y
*/
int cvec_add_mul_c(int n, cmulti **z, cmulti **x, cmulti *y)
{
  int i,e=0;
  for(i=0; i<n; i++){ e+=cadd_mul(z[i],x[i],y); }
  return e;
}

/**
 @brief cmulti型のベクトルxとスカラーyの掛け算の加算 z+=x*y
*/
int cvec_add_mul_r(int n, cmulti **z, cmulti **x, rmulti *y)
{
  int i,e=0;
  for(i=0; i<n; i++){ e+=cadd_mul_r(z[i],x[i],y); }
  return e;
}

/**
 @brief cmulti型のベクトルxとスカラーyの掛け算の加算 z+=x*y
*/
int cvec_add_mul_z(int n, cmulti **z, cmulti **x, dcomplex y)
{
  int i,e=0;
  for(i=0; i<n; i++){ e+=cadd_mul_z(z[i],x[i],y); }
  return e;
}

/**
 @brief cmulti型のベクトルxとスカラーyの掛け算の加算 z+=x*y
*/
int cvec_add_mul_d(int n, cmulti **z, cmulti **x, double y)
{
  int i,e=0;
  for(i=0; i<n; i++){ e+=cadd_mul_d(z[i],x[i],y); }
  return e;
}

/**
 @brief cmulti型のベクトルの要素ごとの掛け算の減算 z-=x.*y
*/
int cvec_sub_mul(int n, cmulti **z, cmulti **x, cmulti **y)
{
  int i,e=0;
  for(i=0; i<n; i++){ e+=csub_mul(z[i],x[i],y[i]); }
  return e;
}

/**
 @brief cmulti型のベクトルxとスカラーyの掛け算の減算 z-=x*y
*/
int cvec_sub_mul_c(int n, cmulti **z, cmulti **x, cmulti *y)
{
  int i,e=0;
  for(i=0; i<n; i++){ e+=csub_mul(z[i],x[i],y); }
  return e;
}

/**
 @brief cmulti型のベクトルxとスカラーyの掛け算の減算 z-=x*y
*/
int cvec_sub_mul_r(int n, cmulti **z, cmulti **x, rmulti *y)
{
  int i,e=0;
  for(i=0; i<n; i++){ e+=csub_mul_r(z[i],x[i],y); }
  return e;
}

/**
 @brief cmulti型のベクトルxとスカラーyの掛け算の減算 z-=x*y
*/
int cvec_sub_mul_z(int n, cmulti **z, cmulti **x, dcomplex y)
{
  int i,e=0;
  for(i=0; i<n; i++){ e+=csub_mul_z(z[i],x[i],y); }
  return e;
}

/**
 @brief cmulti型のベクトルxとスカラーyの掛け算の減算 z-=x*y
*/
int cvec_sub_mul_d(int n, cmulti **z, cmulti **x, double y)
{
  int i,e=0;
  for(i=0; i<n; i++){ e+=csub_mul_d(z[i],x[i],y); }
  return e;
}

/**
 @brief cmulti型のベクトルの要素ごとの掛け算 z=conj(x).*y
*/
int cvec_dot(int n, cmulti **z, cmulti **x, cmulti **y)
{
  int i,e=0;
  for(i=0; i<n; i++){ e+=cdot(z[i],x[i],y[i]); }
  return e;
}

/**
 @brief cmulti型のベクトルyとスカラーxの掛け算 z=conj(x)*y
*/
int cvec_dot_c1(int n, cmulti **z, cmulti *x, cmulti **y)
{
  int i,e=0;
  for(i=0; i<n; i++){ e+=cdot(z[i],x,y[i]); }
  return e;
}

/**
 @brief cmulti型のベクトルxとスカラーyの掛け算 z=conj(x)*y
*/
int cvec_dot_c2(int n, cmulti **z, cmulti **x, cmulti *y)
{
  int i,e=0;
  for(i=0; i<n; i++){ e+=cdot(z[i],x[i],y); }
  return e;
}

/**
 @brief cmulti型のベクトルyとスカラーxの掛け算 z=conj(x)*y
*/
int cvec_dot_z1(int n, cmulti **z, dcomplex x, cmulti **y)
{
  int i,e=0;
  for(i=0; i<n; i++){ e+=cdot_z1(z[i],x,y[i]); }
  return e;
}

/**
 @brief cmulti型のベクトルxとスカラーyの掛け算 z=conj(x)*y
*/
int cvec_dot_z2(int n, cmulti **z, cmulti **x, dcomplex y)
{
  int i,e=0;
  for(i=0; i<n; i++){ e+=cdot_z2(z[i],x[i],y); }
  return e;
}

/**
 @brief cmulti型のベクトルの要素ごとの掛け算の加算 z+=conj(x).*y
*/
int cvec_add_dot(int n, cmulti **z, cmulti **x, cmulti **y)
{
  int i,e=0;
  for(i=0; i<n; i++){ e+=cadd_dot(z[i],x[i],y[i]); }
  return e;
}

/**
 @brief cmulti型のベクトルyとスカラーxの掛け算の加算 z+=conj(x)*y
*/
int cvec_add_dot_c1(int n, cmulti **z, cmulti *x, cmulti **y)
{
  int i,e=0;
  for(i=0; i<n; i++){ e+=cadd_dot(z[i],x,y[i]); }
  return e;
}

/**
 @brief cmulti型のベクトルxとスカラーyの掛け算の加算 z+=conj(x)*y
*/
int cvec_add_dot_c2(int n, cmulti **z, cmulti **x, cmulti *y)
{
  int i,e=0;
  for(i=0; i<n; i++){ e+=cadd_dot(z[i],x[i],y); }
  return e;
}

/**
 @brief cmulti型のベクトルyとスカラーxの掛け算の加算 z+=conj(x)*y
*/
int cvec_add_dot_z1(int n, cmulti **z, dcomplex x, cmulti **y)
{
  int i,e=0;
  for(i=0; i<n; i++){ e+=cadd_dot_z1(z[i],x,y[i]); }
  return e;
}

/**
 @brief cmulti型のベクトルxとスカラーyの掛け算の加算 z+=conj(x)*y
*/
int cvec_add_dot_z2(int n, cmulti **z, cmulti **x, dcomplex y)
{
  int i,e=0;
  for(i=0; i<n; i++){ e+=cadd_dot_z2(z[i],x[i],y); }
  return e;
}

/**
 @brief cmulti型のベクトルの要素ごとの掛け算の減算 z-=conj(x).*y
*/
int cvec_sub_dot(int n, cmulti **z, cmulti **x, cmulti **y)
{
  int i,e=0;
  for(i=0; i<n; i++){ e+=csub_dot(z[i],x[i],y[i]); }
  return e;
}

/**
 @brief cmulti型のベクトルyとスカラーxの掛け算の減算 z-=conj(x)*y
*/
int cvec_sub_dot_c1(int n, cmulti **z, cmulti *x, cmulti **y)
{
  int i,e=0;
  for(i=0; i<n; i++){ e+=csub_dot(z[i],x,y[i]); }
  return e;
}

/**
 @brief cmulti型のベクトルxとスカラーyの掛け算の減算 z-=conj(x)*y
*/
int cvec_sub_dot_c2(int n, cmulti **z, cmulti **x, cmulti *y)
{
  int i,e=0;
  for(i=0; i<n; i++){ e+=csub_dot(z[i],x[i],y); }
  return e;
}

/**
 @brief cmulti型のベクトルyとスカラーxの掛け算の減算 z-=conj(x)*y
*/
int cvec_sub_dot_z1(int n, cmulti **z, dcomplex x, cmulti **y)
{
  int i,e=0;
  for(i=0; i<n; i++){ e+=csub_dot_z1(z[i],x,y[i]); }
  return e;
}

/**
 @brief cmulti型のベクトルxとスカラーyの掛け算の減算 z-=conj(x)*y
*/
int cvec_sub_dot_z2(int n, cmulti **z, cmulti **x, dcomplex y)
{
  int i,e=0;
  for(i=0; i<n; i++){ e+=csub_dot_z2(z[i],x[i],y); }
  return e;
}

/**
 @brief cmulti型のべき乗 z=x^y
*/
int cvec_pow_ui(int n, cmulti **z, cmulti **x, ulong y)
{
  int i,e=0;
  for(i=0; i<n; i++){ e+=cpow_ui(z[i],x[i],y); }
  return e;
}

/**
 @brief cmulti型のベクトルxに関する行列Aによる1次変換 y=A*x
 @param[in]  m   ベクトルyのサイズ.行列Aの行のサイズ.
 @param[in]  n   ベクトルxのサイズ.行列Aの列のサイズ.
 @param[int] y   初期化済みのサイズmのベクトル.
 @param[out] y   y=A*xの計算結果.
 @param[in]  A   サイズ(m,n)の行列.
 @param[in]  LDA 行列Aの第1次元.
 @param[in]  x   サイズnのベクトル.
 @return         丸めの発生回数.
*/
int cvec_lintr(int m, int n, cmulti **y, cmulti **A, int LDA, cmulti **x)
{
  int i,j,e=0;
  cmulti **z=NULL;
  z=cvec_allocate_prec(m,cvec_get_prec_max(m,y));
  for(i=0; i<m; i++){
    e+=cset_zero(z[i]); // z=0
    for(j=0; j<n; j++){
      e+=cadd_mul(z[i],MAT(A,i,j,LDA),x[j]); // z+=A*x
    }
  }
  e+=cvec_copy(m,y,z); // y=z
  z=cvec_free(m,z);
  return e;
}

/**
 @brief cmulti型のベクトルxに関する行列Aによる1次変換の加算 y+=A*x
 @param[in]  m   ベクトルyのサイズ.行列Aの行のサイズ.
 @param[in]  n   ベクトルxのサイズ.行列Aの列のサイズ.
 @param[int] y   初期化済みのサイズmのベクトル.
 @param[out] y   y+=A*xの計算結果.
 @param[in]  A   サイズ(m,n)の行列.
 @param[in]  LDA 行列Aの第1次元.
 @param[in]  x   サイズnのベクトル.
 @return         丸めの発生回数.
*/
int cvec_add_lintr(int m, int n, cmulti **y, cmulti **A, int LDA, cmulti **x)
{
  int i,j,e=0;
  cmulti **z=NULL;
  z=cvec_allocate_prec(m,cvec_get_prec_max(m,y));
  for(i=0; i<m; i++){
    e+=ccopy(z[i],y[i]); // z=y
    for(j=0; j<n; j++){
      e+=cadd_mul(z[i],MAT(A,i,j,LDA),x[j]); // z+=A*x
    }
  }
  e+=cvec_copy(m,y,z); // y=z
  z=cvec_free(m,z);
  return e;
}

/**
 @brief cmulti型のベクトルxに関する行列Aによる1次変換の減算 y-=A*x
 @param[in]  m   ベクトルyのサイズ.行列Aの行のサイズ.
 @param[in]  n   ベクトルxのサイズ.行列Aの列のサイズ.
 @param[int] y   初期化済みのサイズmのベクトル.
 @param[out] y   y-=A*xの計算結果.
 @param[in]  A   サイズ(m,n)の行列.
 @param[in]  LDA 行列Aの第1次元.
 @param[in]  x   サイズnのベクトル.
 @return         丸めの発生回数.
*/
int cvec_sub_lintr(int m, int n, cmulti **y, cmulti **A, int LDA, cmulti **x)
{
  int i,j,e=0;
  cmulti **z=NULL;
  z=cvec_allocate_prec(m,cvec_get_prec_max(m,y));
  for(i=0; i<m; i++){
    e+=ccopy(z[i],y[i]); // z=y
    for(j=0; j<n; j++){
      e+=csub_mul(z[i],MAT(A,i,j,LDA),x[j]); // z-=A*x
    }
  }
  e+=cvec_copy(m,y,z); // y=z
  z=cvec_free(m,z);
  return e;
}

/**
 @brief cmulti型のベクトルxに関する転置行列Aによる1次変換 y=A^T*x
 @param[in]  m   ベクトルxのサイズ.行列Aの行のサイズ.
 @param[in]  n   ベクトルyのサイズ.行列Aの列のサイズ.
 @param[int] y   初期化済みのサイズnのベクトル.
 @param[out] y   y=A^T*xの計算結果.
 @param[in]  A   サイズ(m,n)の行列.
 @param[in]  LDA 行列Aの第1次元.
 @param[in]  x   サイズmのベクトル.
 @return         丸めの発生回数.
*/
int cvec_lintr_t(int m, int n, cmulti **y, cmulti **A, int LDA, cmulti **x)
{
  int i,j,e=0;
  cmulti **z=NULL;
  z=cvec_allocate_prec(n,cvec_get_prec_max(n,y));
  for(j=0; j<n; j++){
    e+=cset_zero(z[j]); // z=0
    for(i=0; i<m; i++){
      e+=cadd_mul(z[j],MAT(A,i,j,LDA),x[i]); // z=A^T*x
    }
  }
  e+=cvec_copy(n,y,z); // y=z
  z=cvec_free(n,z);
  return e;
}

/**
 @brief cmulti型のベクトルxに関する転置行列Aによる1次変換の加算 y+=A^T*x
 @param[in]  m   ベクトルxのサイズ.行列Aの行のサイズ.
 @param[in]  n   ベクトルyのサイズ.行列Aの列のサイズ.
 @param[int] y   初期化済みのサイズnのベクトル.
 @param[out] y   y+=A^T*xの計算結果.
 @param[in]  A   サイズ(m,n)の行列.
 @param[in]  LDA 行列Aの第1次元.
 @param[in]  x   サイズmのベクトル.
 @return         丸めの発生回数.
*/
int cvec_add_lintr_t(int m, int n, cmulti **y, cmulti **A, int LDA, cmulti **x)
{
  int i,j,e=0;
  cmulti **z=NULL;
  z=cvec_allocate_prec(n,cvec_get_prec_max(n,y));
  for(j=0; j<n; j++){
    e+=ccopy(z[j],y[j]); // z=y
    for(i=0; i<m; i++){
      e+=cadd_mul(z[j],MAT(A,i,j,LDA),x[i]); // z+=A^T*x
    }
  }
  e+=cvec_copy(n,y,z); // y=z
  z=cvec_free(n,z);
  return e;
}

/**
 @brief cmulti型のベクトルxに関する転置行列Aによる1次変換の減算 y-=A^T*x
 @param[in]  m   ベクトルxのサイズ.行列Aの行のサイズ.
 @param[in]  n   ベクトルyのサイズ.行列Aの列のサイズ.
 @param[int] y   初期化済みのサイズnのベクトル.
 @param[out] y   y-=A^T*xの計算結果.
 @param[in]  A   サイズ(m,n)の行列.
 @param[in]  LDA 行列Aの第1次元.
 @param[in]  x   サイズmのベクトル.
 @return         丸めの発生回数.
*/
int cvec_sub_lintr_t(int m, int n, cmulti **y, cmulti **A, int LDA, cmulti **x)
{
  int i,j,e=0;
  cmulti **z=NULL;
  z=cvec_allocate_prec(n,cvec_get_prec_max(n,y));
  for(j=0; j<n; j++){
    e+=ccopy(z[j],y[j]); // z=y
    for(i=0; i<m; i++){
      e+=csub_mul(z[j],MAT(A,i,j,LDA),x[i]); // z-=A^T*x
    }
  }
  e+=cvec_copy(n,y,z); // y=z
  z=cvec_free(n,z);
  return e;
}

/**
 @brief cmulti型のベクトルxに関する共役転置行列Aによる1次変換 y=A'*x
 @param[in]  m   ベクトルxのサイズ.行列Aの行のサイズ.
 @param[in]  n   ベクトルyのサイズ.行列Aの列のサイズ.
 @param[int] y   初期化済みのサイズnのベクトル.
 @param[out] y   y=A'*xの計算結果.
 @param[in]  A   サイズ(m,n)の行列.
 @param[in]  LDA 行列Aの第1次元.
 @param[in]  x   サイズmのベクトル.
 @return         丸めの発生回数.
*/
int cvec_lintr_ct(int m, int n, cmulti **y, cmulti **A, int LDA, cmulti **x)
{
  int i,j,e=0;
  cmulti **z=NULL;
  z=cvec_allocate_prec(n,cvec_get_prec_max(n,y));
  for(j=0; j<n; j++){
    e+=cset_zero(z[j]); // z=0
    for(i=0; i<m; i++){
      e+=cadd_dot(z[j],MAT(A,i,j,LDA),x[i]); // z=A'*x
    }
  }
  e+=cvec_copy(n,y,z); // y=z
  z=cvec_free(n,z);
  return e;
}

/**
 @brief cmulti型のベクトルxに関する共役転置行列Aによる1次変換の加算 y+=A'*x
 @param[in]  m   ベクトルxのサイズ.行列Aの行のサイズ.
 @param[in]  n   ベクトルyのサイズ.行列Aの列のサイズ.
 @param[int] y   初期化済みのサイズnのベクトル.
 @param[out] y   y+=A'*xの計算結果.
 @param[in]  A   サイズ(m,n)の行列.
 @param[in]  LDA 行列Aの第1次元.
 @param[in]  x   サイズmのベクトル.
 @return         丸めの発生回数.
*/
int cvec_add_lintr_ct(int m, int n, cmulti **y, cmulti **A, int LDA, cmulti **x)
{
  int i,j,e=0;
  cmulti **z=NULL;
  z=cvec_allocate_prec(n,cvec_get_prec_max(n,y));
  for(j=0; j<n; j++){
    e+=ccopy(z[j],y[j]); // z=y
    for(i=0; i<m; i++){
      e+=cadd_dot(z[j],MAT(A,i,j,LDA),x[i]); // z+=A'*x
    }
  }
  e+=cvec_copy(n,y,z); // y=z
  z=cvec_free(n,z);
  return e;
}

/**
 @brief cmulti型のベクトルxに関する共役転置行列Aによる1次変換の減算 y-=A'*x
 @param[in]  m   ベクトルxのサイズ.行列Aの行のサイズ.
 @param[in]  n   ベクトルyのサイズ.行列Aの列のサイズ.
 @param[int] y   初期化済みのサイズnのベクトル.
 @param[out] y   y-=A'*xの計算結果.
 @param[in]  A   サイズ(m,n)の行列.
 @param[in]  LDA 行列Aの第1次元.
 @param[in]  x   サイズmのベクトル.
 @return         丸めの発生回数.
*/
int cvec_sub_lintr_ct(int m, int n, cmulti **y, cmulti **A, int LDA, cmulti **x)
{
  int i,j,e=0;
  cmulti **z=NULL;
  z=cvec_allocate_prec(n,cvec_get_prec_max(n,y));
  for(j=0; j<n; j++){
    e+=ccopy(z[j],y[j]); // z=y
    for(i=0; i<m; i++){
      e+=csub_dot(z[j],MAT(A,i,j,LDA),x[i]); // z-=A'*x
    }
  }
  e+=cvec_copy(n,y,z); // y=z
  z=cvec_free(n,z);
  return e;
}

/**
 @brief cmulti型のベクトルの要素の総和 value=sum(x)
 */
int cvec_sum(cmulti *value, int n, cmulti **x)
{
  int i,e=0;
  cset_zero(value);
  for(i=0; i<n; i++){ e+=cadd(value,value,x[i]); }
  return e;
}

/**
 @brief cmulti型のベクトルの要素の絶対値の平方の総和 value=sum(abs(x).^2)
*/
int cvec_sum_abs2(rmulti *value, int n, cmulti **x)
{
  int i,e=0;
  e+=rset_zero(value);
  for(i=0; i<n; i++){ e+=cadd_abs2(value,x[i]); }
  return e;
}

/**
 @brief cmulti型のベクトルの内積 value=sum(conj(x).*y)
 */
int cvec_sum_dot(cmulti *value, int n, cmulti **x, cmulti **y)
{
  int i,e=0;
  e+=cset_zero(value);
  for(i=0; i<n; i++){ e+=cadd_dot(value,x[i],y[i]); }
  return e;
}

/**
 @brief cmulti型のベクトルの要素の最大値 value=max(x)
*/
int cvec_max(cmulti *value, int n, cmulti **x)
{
  int i,e=0;
  e+=ccopy(value,x[0]);     // value=x[0]
  for(i=1; i<n; i++){
    if(cgt(x[i],value)){ // x[i]<value
      e+=ccopy(value,x[i]); // value=x[i]
    }
  }
  return e;
}

/**
 @brief rmulti型のベクトルの要素の最小値 value=min(x)
*/
int cvec_min(cmulti *value, int n, cmulti **x)
{
  int i,e=0;
  e+=ccopy(value,x[0]);     // value=x[0]
  for(i=1; i<n; i++){
    if(clt(x[i],value)){ // x[i]<value
      e+=ccopy(value,x[i]); // value=x[i]
    }
  }
  return e;
}

/** @} */

///////////////////////////////////////////////////////////////////////

/** @name cmulti型のベクトルの数学関数に関する関数 */
/** @{ */

/**
 @brief cmulti型のベクトルの要素ごとの割り算 z=x./y
*/
int cvec_div(int n, cmulti **z, cmulti **x, cmulti **y)
{
  int i,e=0;
  for(i=0; i<n; i++){ e+=cdiv(z[i],x[i],y[i]); }
  return e;
}

/**
 @brief cmulti型のスカラーとベクトルの要素との割り算 z=x./y
*/
int cvec_div_c1(int n, cmulti **z, cmulti *x, cmulti **y)
{
  int i,e=0;
  for(i=0; i<n; i++){ e+=cdiv(z[i],x,y[i]); }
  return e;
}

/**
 @brief cmulti型のスカラーとベクトルの要素との割り算 z=x/y
*/
int cvec_div_c2(int n, cmulti **z, cmulti **x, cmulti *y)
{
  int i,e=0;
  for(i=0; i<n; i++){ e+=cdiv(z[i],x[i],y); }
  return e;
}

/**
 @brief cmulti型のスカラーとベクトルの要素との割り算 z=x./y
*/
int cvec_div_r1(int n, cmulti **z, rmulti *x, cmulti **y)
{
  int i,e=0;
  for(i=0; i<n; i++){ e+=cdiv_r1(z[i],x,y[i]); }
  return e;
}

/**
 @brief cmulti型のスカラーとベクトルの要素との割り算 z=x/y
*/
int cvec_div_r2(int n, cmulti **z, cmulti **x, rmulti *y)
{
  int i,e=0;
  for(i=0; i<n; i++){ e+=cdiv_r2(z[i],x[i],y); }
  return e;
}

/**
 @brief cmulti型のスカラーとベクトルの要素との割り算 z=x./y
*/
int cvec_div_z1(int n, cmulti **z, dcomplex x, cmulti **y)
{
  int i,e=0;
  for(i=0; i<n; i++){ e+=cdiv_z1(z[i],x,y[i]); }
  return e;
}

/**
 @brief cmulti型のスカラーとベクトルの要素との割り算 z=x/y
*/
int cvec_div_z2(int n, cmulti **z, cmulti **x, dcomplex y)
{
  int i,e=0;
  for(i=0; i<n; i++){ e+=cdiv_z2(z[i],x[i],y); }
  return e;
}

/**
 @brief cmulti型のスカラーとベクトルの要素との割り算 z=x./y
*/
int cvec_div_d1(int n, cmulti **z, double x, cmulti **y)
{
  int i,e=0;
  for(i=0; i<n; i++){ e+=cdiv_d1(z[i],x,y[i]); }
  return e;
}

/**
 @brief cmulti型のスカラーとベクトルの要素との割り算 z=x/y
*/
int cvec_div_d2(int n, cmulti **z, cmulti **x, double y)
{
  int i,e=0;
  for(i=0; i<n; i++){ e+=cdiv_d2(z[i],x[i],y); }
  return e;
}

/**
 @brief cmulti型のべき乗 z=x.^y
*/
int cvec_pow(int n, cmulti **z, cmulti **x, cmulti **y) 
{
  int i,e=0;
  for(i=0; i<n; i++){ e+=cpow_c(z[i],x[i],y[i]); }
  return e;
}

/**
 @brief cmulti型のべき乗 z=x^y
*/
int cvec_pow_si(int n, cmulti **z, cmulti **x, long y)
{
  int i,e=0;
  for(i=0; i<n; i++){ e+=cpow_si(z[i],x[i],y); }
  return e;
}

/**
 @brief cmulti型のべき乗 z=x^y
*/
int cvec_pow_r(int n, cmulti **z, cmulti **x, rmulti *y)
{
  int i,e=0;
  for(i=0; i<n; i++){ e+=cpow_r2(z[i],x[i],y); }
  return e;
}

/**
 @brief cmulti型のベクトルの要素ごとの絶対値の対数 y=log2(abs(x))
*/
int cvec_log2_abs(int n, rmulti **y, cmulti **x)
{
  int i,e=0;
  for(i=0; i<n; i++){
    e+=cabsv(y[i],x[i]);
    e+=rlog2(y[i],y[i]);
  }
  return e;
}

/**
 @brief cmulti型のベクトルの規格化 y=x/sqrt(x'*x)
*/
int cvec_normalize(int n, cmulti **y, cmulti **x)
{
  int e=0;
  rmulti *a=NULL;
  a=rallocate_prec(cvec_get_prec_max(n,y));
  e+=cvec_norm2(a,n,x);   // a=sqrt(x'*x)
  e+=rinv(a,a);           // a=1/sqrt(x'*x)
  e+=cvec_mul_r(n,y,x,a); // y=x/sqrt(x'*x)
  a=rfree(a);
  return e;
}

/**
 @brief cmulti型のベクトルの規格化 y=x/sqrt(x'*x)
*/
int cvec_normalize_sgn(int n, cmulti **y, cmulti **x)
{
  int prec,k=-1,e=0;
  rmulti *a=NULL;
  cmulti *b=NULL;
  // allocate
  prec=cvec_get_prec_max(n,y);
  a=rallocate_prec(prec);
  b=callocate_prec(prec);
  // compute
  e+=cvec_max_abs_index(a,n,x,&k); // a=abs(x[k])
  e+=cnormalize(b,x[k]);     // b=sgn(x[k])
  e+=cconj(b,b);             // b=conj(sgn(x[k]))
  e+=cvec_norm2(a,n,x);      // a=sqrt(x'*x)
  e+=cdiv_r2(b,b,a);         // b=conj(sgn(x[k]))/sqrt(x'*x)
  e+=cvec_mul_c(n,y,x,b);    // y=x*(conj(sgn(x[k]))/sqrt(x'*x))
  // done
  a=rfree(a);
  b=cfree(b);
  return e;
}

/**
 @brief cmulti型のベクトルの直交化 y-=(x'*y)*x where x'*x=1
*/
int cvec_orthogonalize(int n, cmulti **y, cmulti **x)
{
  int e=0;
  cmulti *a=NULL;
  a=callocate_prec(cvec_get_prec_max(n,y));
  e+=cvec_sum_dot(a,n,x,y);   // a=x'*y
  e+=cvec_sub_mul_c(n,y,x,a); // y-=x*a
  a=cfree(a);
  return e;
}

/**
 @brief cmulti型のベクトルの2乗ノルム value=sqrt(sum(abs(x).^2))
*/
int cvec_norm2(rmulti *value, int n, cmulti **x)
{
  int e=0;
  e+=cvec_sum_abs2(value,n,x); // value=sum(abs(x).^2)
  e+=rsqrt(value,value);       // value=sqrt(value)
  return e;
}

/**
 @brief cmulti型のベクトルの要素の平均 value=sum(x)/n
*/
int cvec_average(cmulti *value, int n, cmulti **x)
{
  int e=0;
  e+=cvec_sum(value,n,x);
  e+=cdiv_si2(value,value,n);
  return e;
}

/**
 @brief cmulti型のベクトルの要素の差の絶対値の平方の最大値 value=max(abs(x-y).^2)
*/
int cvec_max_abs2_sub(rmulti *value, int n, cmulti **x, cmulti **y)
{
  int i,e=0;
  cmulti *a=NULL;
  rmulti *b=NULL;
  a=callocate_prec(rget_prec(value));
  b=rallocate_prec(rget_prec(value));
  e+=csub(a,x[0],y[0]);  // value=x[0]-y[0]
  e+=cabs2(value,a);     // value=(abs(x[0]-y[0]))^2
  for(i=1; i<n; i++){
    e+=csub(a,x[i],y[i]); // a=x[i]-y[i]
    e+=cabs2(b,a);        // b=(abs(x[i]-y[i]))^2
    if(rgt(b,value)){
      rcopy(value,b);     // value=b
    }
  }
  a=cfree(a);
  b=rfree(b);
  return e;
}

/**
 @brief cmulti型のベクトルの要素の絶対値 y=abs(x)
*/
int cvec_abs(int n, rmulti **y, cmulti **x)
{
  int i,e=0;
  for(i=0; i<n; i++){ e+=cabsv(y[i],x[i]); }
  return e;
}

/**
 @brief cmulti型のベクトルの要素の差の絶対値 z=abs(x-y)
*/
int cvec_abs_sub(int n, rmulti **z, cmulti **x, cmulti **y)
{
  int i,e=0;
  for(i=0; i<n; i++){
    e+=cabs_sub(z[i],x[i],y[i]);
  }
  return e;
}

/**
 @brief cmulti型のベクトルの要素とスカラーの差の絶対値 z=abs(x-y)
*/
int cvec_abs_sub_c(int n, rmulti **z, cmulti **x, cmulti *y)
{
  int i,e=0;
  for(i=0; i<n; i++){
    e+=cabs_sub(z[i],x[i],y);
  }
  return e;
}

/**
 @brief cmulti型のベクトルの要素とスカラーの差の絶対値 z=abs(x-y)
*/
int cvec_abs_sub_r(int n, rmulti **z, cmulti **x, rmulti *y)
{
  int i,e=0;
  for(i=0; i<n; i++){
    e+=cabs_sub_r(z[i],x[i],y);
  }
  return e;
}

/**
 @brief rmulti型のベクトルの要素の絶対値の総和 value=sum(abs(x))
 */
int cvec_sum_abs(rmulti *value, int n, cmulti **x)
{
  int i,e=0;
  e+=rset_zero(value); // value=0
  for(i=0; i<n; i++){ e+=cadd_abs(value,x[i]); } // value+=abs(x[i])
  return e;
}

/**
 @brief cmulti型のベクトルの要素の差の絶対値の総和 value=sum(abs(x-y))
*/
int cvec_sum_abs_sub(rmulti *value, int n, cmulti **x, cmulti **y)
{
  int i,e=0;
  cmulti *a=NULL;
  a=callocate_prec(rget_prec(value));
  e+=rset_zero(value); // value=0
  for(i=0; i<n; i++){
    e+=csub(a,x[i],y[i]); // a=x[i]-y[i]
    e+=cadd_abs(value,a); // value+=abs(x[i]-y[i])
  }
  a=cfree(a);
  return e;
}

/**
 @brief cmulti型のベクトルの要素の絶対値の最大値 value=max(abs(x))
*/
int cvec_max_abs(rmulti *value, int n, cmulti **x)
{
  int i,e=0;
  rmulti *a=NULL;
  a=rallocate_prec(rget_prec(value));
  e+=cabs2(value,x[0]);  // value=abs(x[0])^2
  for(i=1; i<n; i++){
    e+=cabs2(a,x[i]);    // a=abs(x[i])^2
    if(rgt(a,value)){    // a>value
      e+=rcopy(value,a); // value=a
    }
  }
  e+=rsqrt(value,value); // value=sqrt(value)
  a=rfree(a);            // free
  return e;
}

/**
 @brief cmulti型のベクトルの要素の実部，虚部の絶対値の最大値 value=max(abs(x))
*/
int cvec_max_absc(rmulti *value, int n, cmulti **x)
{
  int i,e=0;
  rmulti *a=NULL;
  a=rallocate_prec(rget_prec(value));
  e+=rabs(value,C_R(x[0]));
  e+=rabs(a,C_I(x[0]));
  if(rgt(a,value)){ e+=rcopy(value,a); }
  for(i=1; i<n; i++){
    e+=rabs(a,C_R(x[i]));
    if(rgt(a,value)){ e+=rcopy(value,a); }
    e+=rabs(a,C_I(x[i]));
    if(rgt(a,value)){ e+=rcopy(value,a); }
  }
  a=rfree(a);
  return e;
}

/**
 @brief cmulti型のベクトルの要素の絶対値の最大値 value=max(abs(x)) とその添え字 I
*/
int cvec_max_abs_index(rmulti *value, int n, cmulti **x, int *I)
{
  int i,e=0;
  rmulti *a=NULL;
  a=rallocate_prec(rget_prec(value));
  e+=cabs2(value,x[0]);   // value=abs(x[0])^2
  if(I!=NULL){ (*I)=0; }  // I=i
  for(i=1; i<n; i++){
    e+=cabs2(a,x[i]);     // a=abs(x[i])^2
    if(rgt(a,value)){     // a>value
      e+=rcopy(value,a);  // value=a
      if(I!=NULL){ (*I)=i; } // I=i
    }
  }
  e+=rsqrt(value,value); // value=sqrt(value)
  a=rfree(a);
  return e;
}


/**
 @brief cmulti型のベクトルの要素の差の絶対値の最大値 value=max(abs(x-y))
*/
int cvec_max_abs_sub(rmulti *value, int n, cmulti **x, cmulti **y)
{
  int i,e=0;
  cmulti *a=NULL;
  rmulti *b=NULL;
  a=callocate_prec(rget_prec(value));
  b=rallocate_prec(rget_prec(value));
  e+=csub(a,x[0],y[0]);  // a=x[0]-y[0]
  e+=cabs2(value,a);     // value=abs(x[0]-y[0])^2
  for(i=1; i<n; i++){
    e+=csub(a,x[i],y[i]); // a=x[i]-y[i]
    e+=cabs2(b,a);        // b=abs(x[i]-y[i])^2
    if(rgt(b,value)){
      e+=rcopy(value,b);  // value=b
    }
  }
  e+=rsqrt(value,value); // value=sqrt(value)
  a=cfree(a);
  b=rfree(b);
  return e;
}

/**
 @brief cmulti型のベクトルの要素の絶対値の最小値 value=min(abs(x))
*/
int cvec_min_abs(rmulti *value, int n, cmulti **x)
{
  int i,e=0;
  rmulti *a=NULL;
  a=rallocate_prec(rget_prec(value));
  e+=cabs2(value,x[0]);  // value=abs(x[0])^2
  for(i=1; i<n; i++){
    e+=cabs2(a,x[i]);    // a=abs(x[i])^2
    if(rlt(a,value)){    // a<value
      e+=rcopy(value,a); // value=a
    }
  }
  e+=rsqrt(value,value); // value=sqrt(value)
  a=rfree(a);
  return e;
}

/**
 @brief cmulti型のベクトルの要素の絶対値の最小値 value=min(abs(x)) とその添え字 I
*/
int cvec_min_abs_index(rmulti *value, int n, cmulti **x, int *I)
{
  int i,e=0;
  rmulti *a=NULL;
  a=rallocate_prec(rget_prec(value));
  e+=cabs2(value,x[0]);   // value=abs(x[0])^2
  if(I!=NULL){ (*I)=0; }  // I=i
  for(i=1; i<n; i++){
    e+=cabs2(a,x[i]);     // a=abs(x[i])^2
    if(rlt(a,value)){     // a<value
      e+=rcopy(value,a);  // value=a
      if(I!=NULL){ (*I)=i; } // I=i
    }
  }
  e+=rsqrt(value,value); // value=sqrt(value)
  a=rfree(a);
  return e;
}

/**
 @brief cmulti型のベクトルの方向余弦 value=x'*y/sqrt(x'*x)/sqrt(y'*y)
*/
int cvec_dcos(cmulti *value, int n, cmulti **x, cmulti **y)
{
  int e=0;
  rmulti *a=NULL;
  a=rallocate_prec(cget_prec(value));
  e+=cvec_sum_dot(value,n,x,y); // value=x'*y
  e+=cvec_norm2(a,n,x);         // a=sqrt(x'*x)
  e+=cdiv_r2(value,value,a);    // value=x'*y/sqrt(x'*x)
  e+=cvec_norm2(a,n,y);         // a=sqrt(y'*y)
  e+=cdiv_r2(value,value,a);    // value=x'*y/sqrt(x'*x)/sqrt(y'*y)
  a=rfree(a);
  return e;
}

/**
 @brief cmulti型のベクトルの方向余弦の絶対値 value=abs(x'*y)/sqrt(x'*x)/sqrt(y'*y)
*/
int cvec_abs_dcos(rmulti *value, int n, cmulti **x, cmulti **y)
{
  int e=0;
  cmulti *a=NULL;
  a=callocate_prec(rget_prec(value));
  e+=cvec_dcos(a,n,x,y);
  e+=cabsv(value,a);
  a=cfree(a);
  return e;
}

/**
 @brief cmulti型のベクトルの鋭角の角度 theta=acos(abs(x'*y)/sqrt(x'*x)/sqrt(y'*y))
*/
int cvec_angle(rmulti *theta, int n, cmulti **x, cmulti **y)
{
  int e=0;
  rmulti *a=NULL;
  a=rallocate_prec(rget_prec(theta));
  e+=cvec_abs_dcos(a,n,x,y); // a=abs(x'*y)/sqrt(x'*x)/sqrt(y'*y)
  if(rgt_d2(a,1)){ e+=rset_d(a,1); }
  e+=racos(theta,a); // theta=acos(dcos)
  a=rfree(a);
  return e;
}

/**
 @brief cmulti型のベクトルの角度[deg] theta=acos(abs(x'*y)/sqrt(x'*x)/sqrt(y'*y))
*/
int cvec_angle_deg(rmulti *theta, int n, cmulti **x, cmulti **y)
{
  int e=0;
  e+=cvec_angle(theta,n,x,y);
  e+=rmul_d(theta,theta,(180.0/M_PI));
  return e;
}

/** @} */


////////////////////////////////////////////////////////////////////////

/** @name cmulti型の値の比較に関する関数 */
/** @{ */

/**
 @brief cmulti型のベクトルの値の比較 X<=>Y
*/
int cvec_cmp(int n, cmulti **x, int m, cmulti **y)
{
  int i,value;
  for(i=0; i<MIN2(n,m); i++){
    value=ccmp(x[i],y[i]);
    if(value!=0){ return value; } // X!=Y
  }
  if     (n<m){ return -1; } // X<Y
  else if(n>m){ return +1; } // X>Y
  else        { return  0; } // X=Y
}

/**
 @brief cmulti型のベクトルの値の比較 X==Y
*/
int cvec_eq(int n, cmulti **x, cmulti **y)
{
  int i;
  for(i=0; i<n; i++){ if(!ceq(x[i],y[i])){ return 0; } }
  return 1;
}

/**
 @brief cmulti型のベクトルの値の比較 X>Y
*/
int cvec_gt(int n, cmulti **x, cmulti **y)
{
  int i;
  for(i=0; i<n; i++){ if(!cgt(x[i],y[i])){ return 0; } }
  return 1;
}

/**
 @brief cmulti型のベクトルの値の比較 X>=Y
*/
int cvec_ge(int n, cmulti **x, cmulti **y)
{
  int i;
  for(i=0; i<n; i++){ if(!cge(x[i],y[i])){ return 0; } }
  return 1;
}

/**
 @brief cmulti型のベクトルの値の比較 X<Y
*/
int cvec_lt(int n, cmulti **x, cmulti **y)
{
  int i;
  for(i=0; i<n; i++){ if(!clt(x[i],y[i])){ return 0; } }
  return 1;
}

/**
 @brief cmulti型のベクトルの値の比較 X<=Y
*/
int cvec_le(int n, cmulti **x, cmulti **y)
{
  int i;
  for(i=0; i<n; i++){ if(!cle(x[i],y[i])){ return 0; } }
  return 1;
}

/**
 @brief cmulti型のベクトルの値の比較 X==Y
*/
int cvec_eqc(int n, cmulti **x, cmulti **y)
{
  int i;
  for(i=0; i<n; i++){ if(!ceqc(x[i],y[i])){ return 0; } }
  return 1;
}

/**
 @brief cmulti型のベクトルの値の比較 X>Y
*/
int cvec_gtc(int n, cmulti **x, cmulti **y)
{
  int i;
  for(i=0; i<n; i++){ if(!cgtc(x[i],y[i])){ return 0; } }
  return 1;
}

/**
 @brief cmulti型のベクトルの値の比較 X>=Y
*/
int cvec_gec(int n, cmulti **x, cmulti **y)
{
  int i;
  for(i=0; i<n; i++){ if(!cgec(x[i],y[i])){ return 0; } }
  return 1;
}

/**
 @brief cmulti型のベクトルの値の比較 X<Y
*/
int cvec_ltc(int n, cmulti **x, cmulti **y)
{
  int i;
  for(i=0; i<n; i++){ if(!cltc(x[i],y[i])){ return 0; } }
  return 1;
}

/**
 @brief cmulti型のベクトルの値の比較 X<=Y
*/
int cvec_lec(int n, cmulti **x, cmulti **y)
{
  int i;
  for(i=0; i<n; i++){ if(!clec(x[i],y[i])){ return 0; } }
  return 1;
}

/** @} */


///////////////////////////////////////////////////////////////////////////////////////////////////////////

/** @name rmulti型の写像に関する関数 */
/** @{ */

/**
 @brief cmulti型のベクトルに関する写像 y=f(x)
*/
int cvec_func(cmulti *y, func_t *f, int n, cmulti **x)
{
  int i,e=0;
  cmulti *a=NULL,*b=NULL,*z=NULL;
  a=callocate_prec(cget_prec(y));
  b=callocate_prec(cget_prec(y));
  z=callocate_prec(cget_prec(y));
  if(f==NULL)                { FUNC_ERROR_ARG1("cvec_fun",f); }
  else if(func_is(f,"nan"))  { cset_nan(z); }
  else if(func_is(f,"inf"))  { cset_inf(z,1,1); }
  else if(func_is_zero(f))   { e+=cset_si(z,0); }
  else if(func_is_one(f))    { e+=cset_si(z,1); }
  else if(func_is_bigint(f)) { e+=bigint_get_cmulti(z,func_bigint_p(f)); }
  else if(func_is_real(f))   { e+=ccopy_r(z,func_real_p(f)); }
  else if(func_is_complex(f)){ e+=ccopy(z,func_complex_p(f)); }
  else if(func_is_var(f))    {
    cset_d(z,1);
    for(i=0; i<func_var_size(f); i++){
      if(func_var_pow(f,i)!=0){
	if(0<=func_var_num(f,i) && func_var_num(f,i)<n){
	  e+=cpow_si(a,x[func_var_num(f,i)],func_var_pow(f,i));
	}else{ cset_nan(a); }
	e+=cmul(z,z,a);
      }
    }
  }
  else if(func_is_add(f))    { cset_d(z,0); for(i=0; i<func_asize(f); i++){ e+=cvec_func(a,func_aget(f,i),n,x); e+=cadd(z,z,a); } }
  else if(func_is_mul(f))    { cset_d(z,1); for(i=0; i<func_asize(f); i++){ e+=cvec_func(a,func_aget(f,i),n,x); e+=cmul(z,z,a); } }
  else if(func_is(f,"sqrt")) { e+=cvec_func(a,func_aget(f,0),n,x); e+=csqrt_c(z,a); }
  else if(func_is(f,"exp"))  { e+=cvec_func(a,func_aget(f,0),n,x); e+=cexp_c(z,a); }
  else if(func_is(f,"log"))  { e+=cvec_func(a,func_aget(f,0),n,x); e+=clog_c(z,a); }
  else if(func_is(f,"sin"))  { e+=cvec_func(a,func_aget(f,0),n,x); e+=csin_c(z,a); }
  else if(func_is(f,"cos"))  { e+=cvec_func(a,func_aget(f,0),n,x); e+=ccos_c(z,a); }
  else if(func_is(f,"tan"))  { e+=cvec_func(a,func_aget(f,0),n,x); e+=ctan_c(z,a); }
  else if(func_is(f,"asin")) { e+=cvec_func(a,func_aget(f,0),n,x); e+=casin_c(z,a); }
  else if(func_is(f,"acos")) { e+=cvec_func(a,func_aget(f,0),n,x); e+=cacos_c(z,a); }
  else if(func_is(f,"atan")) { e+=cvec_func(a,func_aget(f,0),n,x); e+=catan_c(z,a); }
  else if(func_is(f,"sinh")) { e+=cvec_func(a,func_aget(f,0),n,x); e+=csinh_c(z,a); }
  else if(func_is(f,"cosh")) { e+=cvec_func(a,func_aget(f,0),n,x); e+=ccosh_c(z,a); }
  else if(func_is(f,"tanh")) { e+=cvec_func(a,func_aget(f,0),n,x); e+=ctanh_c(z,a); }
  else if(func_is(f,"asinh")){ e+=cvec_func(a,func_aget(f,0),n,x); e+=casinh_c(z,a); }
  else if(func_is(f,"acosh")){ e+=cvec_func(a,func_aget(f,0),n,x); e+=cacosh_c(z,a); }
  else if(func_is(f,"atanh")){ e+=cvec_func(a,func_aget(f,0),n,x); e+=catanh_c(z,a); }
  else if(func_is(f,"pow"))  { e+=cvec_func(a,func_aget(f,0),n,x); e+=cvec_func(b,func_aget(f,1),n,x); e+=cpow_c(z,a,b); }
  else if(func_is_list(f))   { cset_nan(a); }
  else if(func_is_rvec(f))   { cset_nan(a); }
  else if(func_is_cvec(f))   { cset_nan(a); }
  else if(func_is_rmat(f))   { cset_nan(a); }
  else if(func_is_cmat(f))   { cset_nan(a); }
  else                       { cset_nan(a); }
  if(func_has_power(f))      { e+=cpow_si(z,z,func_power(f)); }
  e+=ccopy(y,z);
  a=cfree(a);
  b=cfree(b);
  z=cfree(z);
  return e;
}

/**
 @brief cmulti型のベクトルに関するベクトル写像 y=f(x) 
*/
int cvec_func_list(int m, cmulti **y, func_t *f, int n, cmulti **x)
{
  int i,e=0;
  if(func_is_list(f)){
    for(i=0; i<m && i<func_asize(f); i++){
      e+=cvec_func(y[i],func_aget(f,i),n,x);
    }
  }
  return e;
}

/**
 @brief cmulti型に関する写像 y=f(x0)
*/
int c1_func(cmulti *y, func_t *f, cmulti *x0)
{
  int e=0,n=1;
  cmulti **x=NULL;
  x=cvec_allocate_prec(n,cget_prec(y));
  e+=ccopy(x[0],x0);
  e+=cvec_func(y,f,n,x);
  x=cvec_free(n,x);
  return e;
}

/**
 @brief cmulti型に関する写像 y=f(x0,x1)
*/
int c2_func(cmulti *y, func_t *f, cmulti *x0, cmulti *x1)
{
  int e=0,n=2;
  cmulti **x=NULL;
  x=cvec_allocate_prec(n,cget_prec(y));
  e+=ccopy(x[0],x0);
  e+=ccopy(x[1],x1);
  e+=cvec_func(y,f,n,x);
  x=cvec_free(n,x);
  return e;
}

/**
 @brief cmulti型に関する写像 y=f(x0,x1,x2)
*/
int c3_func(cmulti *y, func_t *f, cmulti *x0, cmulti *x1, cmulti *x2)
{
  int e=0,n=3;
  cmulti **x=NULL;
  x=cvec_allocate_prec(n,cget_prec(y));
  e+=ccopy(x[0],x0);
  e+=ccopy(x[1],x1);
  e+=ccopy(x[2],x2);
  e+=cvec_func(y,f,n,x);
  x=cvec_free(n,x);
  return e;
}

/* @} */

//EOF
