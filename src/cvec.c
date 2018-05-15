#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<stdarg.h>
#include"is_macros.h"
#include"is_svec.h"
#include"is_rmulti.h"
#include"is_cmulti.h"
#include"is_irmulti.h"
#include"is_icmulti.h"
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
cmulti **cvec_allocate_clone_cvec(int n, cmulti **y)
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
cmulti **cvec_allocate_clone_rvec(int n, rmulti **y)
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
void cvec_round_cvec(int n, cmulti **x, int prec)
{
  int i;
  for(i=0; i<n; i++){ cround(x[i],prec); }
}

/**
 @brief cmulti型のベクトルの値を複成 y=x.
*/
void cvec_clone_cvec(int n, cmulti **y, cmulti **x)
{
  int i;
  for(i=0; i<n; i++){ cclone_c(y[i],x[i]); }
}

/**
 @brief cmulti型のベクトルの値を添字を指定して複成 y=x(I).
*/
void cvec_clone_cvec_index(int n, cmulti **y, cmulti **x, int *I)
{
  int i;
  for(i=0; i<n; i++){ cclone_c(y[i],x[I[i]]); }
}

/**
 @brief cmulti型のベクトルの値を添字を指定して複成 y(I)=x.
*/
void cvec_index_clone_cvec(int n, cmulti **y, int *I, cmulti **x)
{
  int i;
  for(i=0; i<n; i++){ cclone_c(y[I[i]],x[i]); }
}

/**
 @brief cmulti型のベクトルの値をrmulti型のベクトルから複成 y=x.
*/
void cvec_clone_rvec(int n, cmulti **y, rmulti **x)
{
  int i;
  for(i=0; i<n; i++){ cclone_r(y[i],x[i]); }
}

/**
 @brief cmulti型のベクトルの値をrmulti型の2つのベクトルから複成 y=x_r+i*x_i.
*/
void cvec_clone_rvec_rvec(int n, cmulti **y, rmulti **x_r, rmulti **x_i)
{
  int i;
  for(i=0; i<n; i++){ cclone_rr(y[i],x_r[i],x_i[i]); }
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

/**
 @brief cmulti型のベクトルが実数であるかの判定.
*/
int cvec_is_real(int n, cmulti **x)
{
  int i;
  for(i=0; i<n; i++){ if(!cis_real(x[i])){ return 0; } }
  return 1;
}

/** @} */

/////////////////////////////////////////////////////

/** @name cmulti型のベクトルの入出力に関する関数 */
/** @{ */

/**
 @brief cmulti型のベクトルの表示.
*/
void cvec_print(int n, cmulti **x, char *name, char format, int digits)
{
  char **s=NULL;
  if(x==NULL){
    if(name!=NULL){ printf("%s=NULL\n",name); }
    else          { printf("NULL\n"); }
    return;
  }
  s=svec_allocate(n);
  cvec_get_svec(n,s,x,format,digits);
  svec_print(n,s,name);
  s=svec_free(n,s);
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
  char str[100000];
  // file name
  va_start(argp,fmt);
  vsprintf(fname,fmt,argp);
  // open file
  if((fid=fopen(fname,"r"))==0){ ERROR_EXIT("Cant't open file: %s\n",fname); }
  // real
  for(i=0;i<n;i++){
    if(fgets(str,100000,fid)==NULL){ ERROR_EXIT("Error in cvec_read(), fscanf(), fname=%s\n",fname); }
    cset_s(x[i],str);
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
      cx=cvec_allocate_clone_rvec((*n),rx); // allocate
      rx=rvec_free((*n),rx);             // free
    }else if(k==(size_t)l && strncmp(buf,"zvec",l)==0){ /* zvec */
      fclose(fid);                    // close
      zx=zvec_bin_load(n,fname);      // load
      cx=cvec_allocate_prec((*n),53); // allocate
      cvec_set_zvec((*n),cx,zx);         // copy
      zx=zvec_free(zx);               // free
    }else if(k==(size_t)l && strncmp(buf,"dvec",l)==0){ /* dvec */
      fclose(fid);                    // close
      dx=dvec_bin_load(n,fname);      // load
      cx=cvec_allocate_prec((*n),53); // allocate
      cvec_set_dvec((*n),cx,dx);         // copy
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
 @brief cmulti型のベクトルの値をNaNに設定.
*/
void cvec_set_inf(int n, cmulti **x, int rsgn, int isgn)
{
  int i;
  for(i=0; i<n; i++){ cset_inf(x[i],rsgn,isgn); }
}

/**
 @brief cmulti型のベクトルの値を文字列から設定.
*/
void cvec_set_svec(int n, cmulti **x, char **str)
{
  int i;
  for(i=0; i<n; i++){ cset_s(x[i],str[i]); }
}

/**
 @brief cmulti型のベクトルの値のコピー y=x.
 */
void cvec_set_cvec(int n, cmulti **y, cmulti **x)
{
  int i;
  for(i=0; i<n; i++){ cset_c(y[i],x[i]); }
}

/**
 @brief cmulti型のベクトルの値を倍精度複素数から設定.
*/
void cvec_set_zvec(int n, cmulti **y, dcomplex *x)
{
  int i;
  for(i=0; i<n; i++){ cset_z(y[i],x[i]); }
}

/**
 @brief cmulti型のベクトルの値を倍精度実数から設定.
*/
void cvec_set_dvec(int n, cmulti **y, double *x)
{
  int i;
  for(i=0; i<n; i++){ cset_d(y[i],x[i]); }
}

// y=x
void dvec_set_cvec(int n, double *y, cmulti **x)
{
  int i;
  for(i=0; i<n; i++){ y[i]=rget_d(C_R(x[i])); }
}

// y=x
void zvec_set_cvec(int n, dcomplex *y, cmulti **x)
{
  int i;
  for(i=0; i<n; i++){ y[i]=cget_z(x[i]); }
}


/**
 @brief rmulti型のベクトルの値を設定.
 */
void rvec_set_cvec(int n, rmulti **y, cmulti **x)
{
  int i;
  for(i=0; i<n; i++){ rset_r(y[i],C_R(x[i])); }
}


/**
 @brief cmulti型のベクトルの値を設定.
 */
void cvec_set_rvec(int n, cmulti **y, rmulti **x)
{
  int i;
  for(i=0; i<n; i++){ cset_r(y[i],x[i]); }
}

/**
 @brief cmulti型のベクトルの値を倍精度浮実数から設定.
*/
void cvec_set_dvec_dvec(int n, cmulti **y, double *x_r, double *x_i)
{
  int i;
  for(i=0; i<n; i++){ cset_dd(y[i],x_r[i],x_i[i]); }
}

/**
 @brief cmulti型のベクトルの全ての値を倍精度複素数から設定.
*/
void cvec_set_all_z(int n, cmulti **x, dcomplex a)
{
  int i;
  for(i=0; i<n; i++){ cset_z(x[i],a); }
}

/**
 @brief cmulti型のベクトルの全ての値を倍精度実数から設定.
*/
void cvec_set_all_d(int n, cmulti **x, double a)
{
  int i;
  for(i=0; i<n; i++){ cset_d(x[i],a); }
}

/**
 @brief cmulti型のベクトルの全ての値を倍精度実数から設定.
*/
void cvec_set_all_dd(int n, cmulti **x, double a_r, double a_i)
{
  int i;
  for(i=0; i<n; i++){ cset_dd(x[i],a_r,a_i); }
}

/**
 @brief cmulti型のベクトルの値を零に設定.
*/
void cvec_set_zeros(int n, cmulti **x)
{
  return cvec_set_all_d(n,x,0);
}

/**
 @brief cmulti型のベクトルの値を1に設定.
*/
void cvec_set_ones(int n, cmulti **x)
{
  return cvec_set_all_d(n,x,1);
}

/**
 @brief cmulti型のベクトルの値を第k要素が1の基本ベクトルに設定.
*/
void cvec_set_unit(int n, cmulti **x, int k)
{
  cvec_set_zeros(n,x);
  cset_d(x[k],1);
}

/**
 @brief cmulti型のベクトルの値を0,1,2,...に設定.
*/
void cvec_set_grid(int n, cmulti **x)
{
  int i;
  for(i=0; i<n; i++){ cset_d(x[i],i); }
}

/**
 @brief cmulti型のベクトルの値を区間(b,a+b)の疑似乱数値を設定.
*/
void cvec_set_rand(int n, cmulti **x, double a, double b)
{
  dcomplex bz;
  Z_SET(bz,b,b);
  int i;
  for(i=0; i<n; i++){
    cset_rand(x[i]);
    cmul_cd(x[i],x[i],a);
    cadd_cz(x[i],x[i],bz);
  }
}

/** @} */

////////////////////////////////////////////////////////////////////////

/** @name cmulti型のベクトルの型変換に関する関数 */
/** @{ */

/**
 @brief cmulti型のベクトルを倍精度複素数型に変換.
 */
void cvec_get_zvec(int n, dcomplex *y, cmulti **x)
{
  int i;
  for(i=0; i<n; i++){ y[i]=cget_z(x[i]); }
}

/**
 @brief cmulti型のベクトルを倍精度実数型に変換.
 */
void cvec_get_dvec(int n, double *y, cmulti **x)
{
  int i;
  for(i=0; i<n; i++){ y[i]=cget_d(x[i]); }
}

/**
 @brief cmulti型のベクトルを整数型に変換.
 */
void cvec_get_ivec(int n, int *y, cmulti **x)
{
  int i;
  for(i=0; i<n; i++){ y[i]=cget_d(x[i]); }
}

/**
 @brief rmulti型のベクトルを文字列型に変換 y=char(x)
 */
void cvec_get_svec(int n, char **y, cmulti **x, char format, int digits)
{
  char f[1024],buf[1<<13];
  int i;
  sprintf(f,"%%-.%dR%c%%+.%dR%ci",digits,format,digits,format);
  for(i=0; i<n; i++){
    mpfr_sprintf(buf,f,C_R(x[i]),C_I(x[i]));
    y[i]=char_renew(y[i],buf,NULL);
  }
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
  for(i=0; i<n; i++){ cswap(y[i],x[I[i]]); }
  cvec_swap(n,x,y);
  y=cvec_free(n,y);
}

/**
 @brief cmulti型のベクトルの要素を小さい順に入れ替える.
 @details 要素x[i]の値の小さい順に並べ替え，要素番号を配列Iに保存する.
          ただし，I==NULLの場合は保存は行われない.
 @param[in]      n ベクトルのサイズ.
 @param[in,out]  x   [in]初期化済みのrmulti型のベクトル.[out]ソートされたrmulti型のベクトル.
 @param[in,out]  I   [in]入れ替えられた要素の要素番号の保存先.[out]保存された要素番号.
*/
void cvec_sort(int n, cmulti **x, int *I)
{
  if(I!=NULL){ ivec_set_grid(n,I); }
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
    if(lt_cc(x[i],x[left])){
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
  Y=cvec_allocate_clone_cvec(n,X);
  cvec_sort(n,Y,I);
  Y=cvec_free(n,Y);
}


/** @} */

///////////////////////////////////////////////////////////////////////

/** @name cmulti型のベクトルの数学関数に関する関数 */
/** @{ */

/**
 @brief cmulti型のベクトルの値のコピー y=x.
*/
void cvec_copy_cvec(int n, cmulti **y, cmulti **x)
{
  int i;
  if(y==x){ return; }
  for(i=0; i<n; i++){ cset_c(y[i],x[i]); }
}

/**
 @brief cmulti型のベクトルの値のrmulti型からのコピー y=x.
*/
void cvec_set_rvec_rvec(int n, cmulti **y, rmulti **x_r, rmulti **x_i)
{
  int i;
  for(i=0; i<n; i++){ cset_rr(y[i],x_r[i],x_i[i]); }
}

/**
 @brief cmulti型のベクトルの値の要素番号Iに従ったコピー y=x(I).
*/
void cvec_copy_cvec_index(int n, cmulti **y, cmulti **x, int *I)
{
  int i;
  if(y==x){ cvec_swap_index(n,x,I); return; }
  for(i=0; i<n; i++){ cset_c(y[i],x[I[i]]); }
}

/**
 @brief cmulti型のベクトルの値を添字を指定してコピー y(I)=x.
*/
void cvec_index_copy_cvec(int n, cmulti **y, int *I, cmulti **x)
{
  int i;
  for(i=0; i<n; i++){ cset_c(y[I[i]],x[i]); }
}

/**
 @brief cmulti型のベクトルの値を添字を指定してコピー y(I)=x.
*/
void cvec_index_copy_rvec(int n, cmulti **y, int *I, rmulti **x)
{
  int i;
  for(i=0; i<n; i++){ cset_r(y[I[i]],x[i]); }
}

/**
 @brief cmulti型の実部のコピー y=real(x)
 */
void rvec_real_cvec(int n, rmulti **y, cmulti **x)
{
  int i;
  for(i=0; i<n; i++){ rset_r(y[i],C_R(x[i])); }
}

/**
 @brief cmulti型の実部の複製 y=real(x)
 */
void rvec_real_cvec_clone(int n, rmulti **y, cmulti **x)
{
  int i;
  for(i=0; i<n; i++){ rclone_r(y[i],C_R(x[i])); }
}

/**
 @brief cmulti型の虚部のコピー y=imag(x)
 */
void rvec_imag_cvec(int n, rmulti **y, cmulti **x)
{
  int i;
  for(i=0; i<n; i++){ rset_r(y[i],C_I(x[i])); }
}

/**
 @brief cmulti型の虚部の複製 y=imag(x)
 */
void rvec_imag_cvec_clone(int n, rmulti **y, cmulti **x)
{
  int i;
  for(i=0; i<n; i++){ rclone_r(y[i],C_I(x[i])); }
}

/**
 @brief cmulti型の複素共役のコピー y=imag(x)
 */
void cvec_conj_cvec(int n, cmulti **y, cmulti **x)
{
  int i;
  for(i=0; i<n; i++){ cconj_c(y[i],x[i]); }
}

/**
 @brief cmulti型のベクトルの指数部の足し算 y=x*2^p
*/
void cvec_mul_2exp(int n, cmulti **y, cmulti **x, int pr, int pi)
{
  int i;
  for(i=0; i<n; i++){ cmul_2exp(y[i],x[i],pr,pi); }
}

/**
 @brief cmulti型のベクトルの指数部の引き算 y=x/2^p
*/
void cvec_div_2exp(int n, cmulti **y, cmulti **x, int pr, int pi)
{
  int i;
  for(i=0; i<n; i++){ cdiv_2exp(y[i],x[i],pr,pi); }
}

/**
 @brief cmulti型のベクトルの符号反転 y=-x
*/
void cvec_neg_cvec(int n, cmulti **y, cmulti **x)
{
  int i;
  for(i=0; i<n; i++){ cneg_c(y[i],x[i]); }
}

/**
 @brief cmulti型のベクトルの要素の実部，虚部の絶対値 y=abs(real(x)+i*imag(x))
*/
void cvec_absc_cvec(int n, cmulti **y, cmulti **x)
{
  int i;
  for(i=0; i<n; i++){ cabsc_c(y[i],x[i]); }
}

/////////////////////////////////////////////////////

/**
 @brief z=x+y
 */
void cvec_add_cvec(int n, cmulti **z, cmulti **x, cmulti **y)
{
  int i;
  for(i=0; i<n; i++){ cadd_cc(z[i],x[i],y[i]); }
}

/**
 @brief z=x+y
 */
void cvec_add_rvec(int n, cmulti **z, cmulti **x, rmulti **y)
{
  int i;
  for(i=0; i<n; i++){ cadd_cr(z[i],x[i],y[i]); }
}

/**
 @brief z=x+y
 */
void rvec_add_cvec(int n, cmulti **z, rmulti **x, cmulti **y)
{
  int i;
  for(i=0; i<n; i++){ cadd_cr(z[i],y[i],x[i]); }
}

/**
 @brief z=x+y
 */
void cvec_add_zvec(int n, cmulti **z, cmulti **x, dcomplex *y)
{
  int i;
  for(i=0; i<n; i++){ cadd_cz(z[i],x[i],y[i]); }
}

/**
 @brief z=x+y
 */
void zvec_add_cvec(int n, cmulti **z, dcomplex *x, cmulti **y)
{
  int i;
  for(i=0; i<n; i++){ cadd_cz(z[i],y[i],x[i]); }
}

/**
 @brief z=x+y
 */
void cvec_add_dvec(int n, cmulti **z, cmulti **x, double *y)
{
  int i;
  for(i=0; i<n; i++){ cadd_cd(z[i],x[i],y[i]); }
}

/**
 @brief z=x+y
 */
void dvec_add_cvec(int n, cmulti **z, double *x, cmulti **y)
{
  int i;
  for(i=0; i<n; i++){ cadd_cd(z[i],y[i],x[i]); }
}

/**
 @brief z=x+y
 */
void rvec_add_zvec(int n, cmulti **z, rmulti **x, dcomplex *y)
{
  int i;
  for(i=0; i<n; i++){ cadd_rz(z[i],x[i],y[i]); }
}

/**
 @brief z=x+y
 */
void zvec_add_rvec(int n, cmulti **z, dcomplex *x, rmulti **y)
{
  int i;
  for(i=0; i<n; i++){ cadd_zr(z[i],x[i],y[i]); }
}

///////////////////////////////////////////////////////

/**
 @brief z=x+y
 */
void cvec_add_cscalar(int n, cmulti **z, cmulti **x, cmulti *y)
{
  int i;
  for(i=0; i<n; i++){ cadd_cc(z[i],x[i],y); }
}

/**
 @brief z=x+y
 */
void cvec_add_rscalar(int n, cmulti **z, cmulti **x, rmulti *y)
{
  int i;
  for(i=0; i<n; i++){ cadd_cr(z[i],x[i],y); }
}

/**
 @brief z=x+y
 */
void rvec_add_cscalar(int n, cmulti **z, rmulti **x, cmulti *y)
{
  int i;
  for(i=0; i<n; i++){ cadd_cr(z[i],y,x[i]); }
}

/**
 @brief z=x+y
 */
void cvec_add_zscalar(int n, cmulti **z, cmulti **x, dcomplex y)
{
  int i;
  for(i=0; i<n; i++){ cadd_cz(z[i],x[i],y); }
}

/**
 @brief z=x+y
 */
void zvec_add_cscalar(int n, cmulti **z, dcomplex *x, cmulti *y)
{
  int i;
  for(i=0; i<n; i++){ cadd_cz(z[i],y,x[i]); }
}

/**
 @brief z=x+y
 */
void cvec_add_dscalar(int n, cmulti **z, cmulti **x, double y)
{
  int i;
  for(i=0; i<n; i++){ cadd_cd(z[i],x[i],y); }
}

/**
 @brief z=x+y
 */
void dvec_add_cscalar(int n, cmulti **z, double *x, cmulti *y)
{
  int i;
  for(i=0; i<n; i++){ cadd_cd(z[i],y,x[i]); }
}

/**
 @brief z=x+y
 */
void rvec_add_zscalar(int n, cmulti **z, rmulti **x, dcomplex y)
{
  int i;
  for(i=0; i<n; i++){ cadd_rz(z[i],x[i],y); }
}

/**
 @brief z=x+y
 */
void zvec_add_rscalar(int n, cmulti **z, dcomplex *x, rmulti *y)
{
  int i;
  for(i=0; i<n; i++){ cadd_zr(z[i],x[i],y); }
}

///////////////////////////////////////////////////////

/**
 @brief z=x+y
 */
void cscalar_add_cvec(int n, cmulti **z, cmulti *x, cmulti **y)
{
    int i;
    for(i=0; i<n; i++){ cadd_cc(z[i],x,y[i]); }
}

/**
 @brief z=x+y
 */
void cscalar_add_rvec(int n, cmulti **z, cmulti *x, rmulti **y)
{
    int i;
    for(i=0; i<n; i++){ cadd_cr(z[i],x,y[i]); }
}

/**
 @brief z=x+y
 */
void rscalar_add_cvec(int n, cmulti **z, rmulti *x, cmulti **y)
{
    int i;
    for(i=0; i<n; i++){ cadd_rc(z[i],x,y[i]); }
}

/**
 @brief z=x+y
 */
void cscalar_add_zvec(int n, cmulti **z, cmulti *x, dcomplex *y)
{
    int i;
    for(i=0; i<n; i++){ cadd_cz(z[i],x,y[i]); }
}

/**
 @brief z=x+y
 */
void zscalar_add_cvec(int n, cmulti **z, dcomplex x, cmulti **y)
{
    int i;
    for(i=0; i<n; i++){ cadd_zc(z[i],x,y[i]); }
}

/**
 @brief z=x+y
 */
void cscalar_add_dvec(int n, cmulti **z, cmulti *x, double *y)
{
    int i;
    for(i=0; i<n; i++){ cadd_cd(z[i],x,y[i]); }
}

/**
 @brief z=x+y
 */
void dscalar_add_cvec(int n, cmulti **z, double x, cmulti **y)
{
    int i;
    for(i=0; i<n; i++){ cadd_dc(z[i],x,y[i]); }
}

/**
 @brief z=x+y
 */
void rscalar_add_zvec(int n, cmulti **z, rmulti *x, dcomplex *y)
{
    int i;
    for(i=0; i<n; i++){ cadd_rz(z[i],x,y[i]); }
}

/**
 @brief z=x+y
 */
void zscalar_add_rvec(int n, cmulti **z, dcomplex x, rmulti **y)
{
    int i;
    for(i=0; i<n; i++){ cadd_zr(z[i],x,y[i]); }
}

///////////////////////////////////////////////////////

/**
 @brief z=x-y
 */
void cvec_sub_cvec(int n, cmulti **z, cmulti **x, cmulti **y)
{
    int i;
    for(i=0; i<n; i++){ csub_cc(z[i],x[i],y[i]); }
}

/**
 @brief z=x-y
 */
void cvec_sub_rvec(int n, cmulti **z, cmulti **x, rmulti **y)
{
    int i;
    for(i=0; i<n; i++){ csub_cr(z[i],x[i],y[i]); }
}

/**
 @brief z=x-y
 */
void rvec_sub_cvec(int n, cmulti **z, rmulti **x, cmulti **y)
{
    int i;
    for(i=0; i<n; i++){ csub_rc(z[i],x[i],y[i]); }
}

/**
 @brief z=x-y
 */
void cvec_sub_zvec(int n, cmulti **z, cmulti **x, dcomplex *y)
{
  int i;
  for(i=0; i<n; i++){ csub_cz(z[i],x[i],y[i]); }
}

/**
 @brief z=x-y
 */
void zvec_sub_cvec(int n, cmulti **z, dcomplex *x, cmulti **y)
{
  int i;
  for(i=0; i<n; i++){ csub_zc(z[i],x[i],y[i]); }
}

/**
 @brief z=x-y
 */
void cvec_sub_dvec(int n, cmulti **z, cmulti **x, double *y)
{
  int i;
  for(i=0; i<n; i++){ csub_cd(z[i],x[i],y[i]); }
}

/**
 @brief z=x-y
 */
void dvec_sub_cvec(int n, cmulti **z, double *x, cmulti **y)
{
  int i;
  for(i=0; i<n; i++){ csub_dc(z[i],x[i],y[i]); }
}

/**
 @brief z=x-y
 */
void rvec_sub_zvec(int n, cmulti **z, rmulti **x, dcomplex *y)
{
    int i;
    for(i=0; i<n; i++){ csub_rz(z[i],x[i],y[i]); }
}

/**
 @brief z=x-y
 */
void zvec_sub_rvec(int n, cmulti **z, dcomplex *x, rmulti **y)
{
  int i;
  for(i=0; i<n; i++){ csub_zr(z[i],x[i],y[i]); }
}

///////////////////////////////////////////////////////

/**
 @brief z=x-y
 */
void cvec_sub_cscalar(int n, cmulti **z, cmulti **x, cmulti *y)
{
  int i;
  for(i=0; i<n; i++){ csub_cc(z[i],x[i],y); }
}

/**
 @brief z=x-y
 */
void cvec_sub_rscalar(int n, cmulti **z, cmulti **x, rmulti *y)
{
    int i;
    for(i=0; i<n; i++){ csub_cr(z[i],x[i],y); }
}

/**
 @brief z=x-y
 */
void rvec_sub_cscalar(int n, cmulti **z, rmulti **x, cmulti *y)
{
    int i;
    for(i=0; i<n; i++){ csub_rc(z[i],x[i],y); }
}

/**
 @brief z=x-y
 */
void cvec_sub_zscalar(int n, cmulti **z, cmulti **x, dcomplex y)
{
    int i;
    for(i=0; i<n; i++){ csub_cz(z[i],x[i],y); }
}

/**
 @brief z=x-y
 */
void zvec_sub_cscalar(int n, cmulti **z, dcomplex *x, cmulti *y)
{
    int i;
    for(i=0; i<n; i++){ csub_zc(z[i],x[i],y); }
}

/**
 @brief z=x-y
 */
void cvec_sub_dscalar(int n, cmulti **z, cmulti **x, double y)
{
    int i;
    for(i=0; i<n; i++){ csub_cd(z[i],x[i],y); }
}

/**
 @brief z=x-y
 */
void dvec_sub_cscalar(int n, cmulti **z, double *x, cmulti *y)
{
    int i;
    for(i=0; i<n; i++){ csub_dc(z[i],x[i],y); }
}

/**
 @brief z=x-y
 */
void rvec_sub_zscalar(int n, cmulti **z, rmulti **x, dcomplex y)
{
    int i;
    for(i=0; i<n; i++){ csub_rz(z[i],x[i],y); }
}

/**
 @brief z=x-y
 */
void zvec_sub_rscalar(int n, cmulti **z, dcomplex *x, rmulti *y)
{
    int i;
    for(i=0; i<n; i++){ csub_zr(z[i],x[i],y); }
}



///////////////////////////////////////////////////////

/**
 @brief z=x-y
 */
void cscalar_sub_cvec(int n, cmulti **z, cmulti *x, cmulti **y)
{
  int i;
  for(i=0; i<n; i++){ csub_cc(z[i],x,y[i]); }
}

/**
 @brief z=x-y
 */
void cscalar_sub_rvec(int n, cmulti **z, cmulti *x, rmulti **y)
{
  int i;
  for(i=0; i<n; i++){ csub_cr(z[i],x,y[i]); }
}

/**
 @brief z=x-y
 */
void rscalar_sub_cvec(int n, cmulti **z, rmulti *x, cmulti **y)
{
    int i;
    for(i=0; i<n; i++){ csub_rc(z[i],x,y[i]); }
}

/**
 @brief z=x-y
 */
void cscalar_sub_zvec(int n, cmulti **z, cmulti *x, dcomplex *y)
{
  int i;
  for(i=0; i<n; i++){ csub_cz(z[i],x,y[i]); }
}

/**
 @brief z=x-y
 */
void zscalar_sub_cvec(int n, cmulti **z, dcomplex x, cmulti **y)
{
    int i;
    for(i=0; i<n; i++){ csub_zc(z[i],x,y[i]); }
}

/**
 @brief z=x-y
 */
void cscalar_sub_dvec(int n, cmulti **z, cmulti *x, double *y)
{
  int i;
  for(i=0; i<n; i++){ csub_cd(z[i],x,y[i]); }
}

/**
 @brief z=x-y
 */
void dscalar_sub_cvec(int n, cmulti **z, double x, cmulti **y)
{
  int i;
  for(i=0; i<n; i++){ csub_dc(z[i],x,y[i]); }
}

/**
 @brief z=x-y
 */
void rscalar_sub_zvec(int n, cmulti **z, rmulti *x, dcomplex *y)
{
    int i;
    for(i=0; i<n; i++){ csub_rz(z[i],x,y[i]); }
}

/**
 @brief z=x-y
 */
void zscalar_sub_rvec(int n, cmulti **z, dcomplex x, rmulti **y)
{
  int i;
  for(i=0; i<n; i++){ csub_zr(z[i],x,y[i]); }
}

///////////////////////////////////////////////////////

/**
 @brief z=x.*y
 */
void cvec_mul_cvec(int n, cmulti **z, cmulti **x, cmulti **y)
{
    int i;
    for(i=0; i<n; i++){ cmul_cc(z[i],x[i],y[i]); }
}

/**
 @brief z=x.*y
 */
void cvec_mul_rvec(int n, cmulti **z, cmulti **x, rmulti **y)
{
    int i;
    for(i=0; i<n; i++){ cmul_cr(z[i],x[i],y[i]); }
}

/**
 @brief z=x.*y
 */
void rvec_mul_cvec(int n, cmulti **z, rmulti **x, cmulti **y)
{
    int i;
    for(i=0; i<n; i++){ cmul_rc(z[i],x[i],y[i]); }
}

/**
 @brief z=x.*y
 */
void cvec_mul_zvec(int n, cmulti **z, cmulti **x, dcomplex *y)
{
    int i;
    for(i=0; i<n; i++){ cmul_cz(z[i],x[i],y[i]); }
}

/**
 @brief z=x.*y
 */
void zvec_mul_cvec(int n, cmulti **z, dcomplex *x, cmulti **y)
{
    int i;
    for(i=0; i<n; i++){ cmul_zc(z[i],x[i],y[i]); }
}

/**
 @brief z=x.*y
 */
void cvec_mul_dvec(int n, cmulti **z, cmulti **x, double *y)
{
    int i;
    for(i=0; i<n; i++){ cmul_cd(z[i],x[i],y[i]); }
}

/**
 @brief z=x.*y
 */
void dvec_mul_cvec(int n, cmulti **z, double *x, cmulti **y)
{
    int i;
    for(i=0; i<n; i++){ cmul_dc(z[i],x[i],y[i]); }
}

/**
 @brief z=x.*y
 */
void rvec_mul_zvec(int n, cmulti **z, rmulti **x, dcomplex *y)
{
    int i;
    for(i=0; i<n; i++){ cmul_rz(z[i],x[i],y[i]); }
}

/**
 @brief z=x.*y
 */
void zvec_mul_rvec(int n, cmulti **z, dcomplex *x, rmulti **y)
{
    int i;
    for(i=0; i<n; i++){ cmul_zr(z[i],x[i],y[i]); }
}

///////////////////////////////////////////////////////

/**
 @brief z=x.*y
 */
void cvec_mul_cscalar(int n, cmulti **z, cmulti **x, cmulti *y)
{
    int i;
    for(i=0; i<n; i++){ cmul_cc(z[i],x[i],y); }
}

/**
 @brief z=x.*y
 */
void cvec_mul_rscalar(int n, cmulti **z, cmulti **x, rmulti *y)
{
    int i;
    for(i=0; i<n; i++){ cmul_cr(z[i],x[i],y); }
}

/**
 @brief z=x.*y
 */
void rvec_mul_cscalar(int n, cmulti **z, rmulti **x, cmulti *y)
{
    int i;
    for(i=0; i<n; i++){ cmul_rc(z[i],x[i],y); }
}

/**
 @brief z=x.*y
 */
void cvec_mul_zscalar(int n, cmulti **z, cmulti **x, dcomplex y)
{
    int i;
    for(i=0; i<n; i++){ cmul_cz(z[i],x[i],y); }
}

/**
 @brief z=x.*y
 */
void zvec_mul_cscalar(int n, cmulti **z, dcomplex *x, cmulti *y)
{
    int i;
    for(i=0; i<n; i++){ cmul_zc(z[i],x[i],y); }
}

/**
 @brief z=x.*y
 */
void cvec_mul_dscalar(int n, cmulti **z, cmulti **x, double y)
{
    int i;
    for(i=0; i<n; i++){ cmul_cd(z[i],x[i],y); }
}

/**
 @brief z=x.*y
 */
void dvec_mul_cscalar(int n, cmulti **z, double *x, cmulti *y)
{
    int i;
    for(i=0; i<n; i++){ cmul_dc(z[i],x[i],y); }
}

/**
 @brief z=x.*y
 */
void rvec_mul_zscalar(int n, cmulti **z, rmulti **x, dcomplex y)
{
    int i;
    for(i=0; i<n; i++){ cmul_rz(z[i],x[i],y); }
}

/**
 @brief z=x.*y
 */
void zvec_mul_rscalar(int n, cmulti **z, dcomplex *x, rmulti *y)
{
    int i;
    for(i=0; i<n; i++){ cmul_zr(z[i],x[i],y); }
}



///////////////////////////////////////////////////////

/**
 @brief z=x.*y
 */
void cscalar_mul_cvec(int n, cmulti **z, cmulti *x, cmulti **y)
{
    int i;
    for(i=0; i<n; i++){ cmul_cc(z[i],x,y[i]); }
}

/**
 @brief z=x.*y
 */
void cscalar_mul_rvec(int n, cmulti **z, cmulti *x, rmulti **y)
{
    int i;
    for(i=0; i<n; i++){ cmul_cr(z[i],x,y[i]); }
}

/**
 @brief z=x.*y
 */
void rscalar_mul_cvec(int n, cmulti **z, rmulti *x, cmulti **y)
{
    int i;
    for(i=0; i<n; i++){ cmul_rc(z[i],x,y[i]); }
}

/**
 @brief z=x.*y
 */
void cscalar_mul_zvec(int n, cmulti **z, cmulti *x, dcomplex *y)
{
    int i;
    for(i=0; i<n; i++){ cmul_cz(z[i],x,y[i]); }
}

/**
 @brief z=x.*y
 */
void zscalar_mul_cvec(int n, cmulti **z, dcomplex x, cmulti **y)
{
    int i;
    for(i=0; i<n; i++){ cmul_zc(z[i],x,y[i]); }
}

/**
 @brief z=x.*y
 */
void cscalar_mul_dvec(int n, cmulti **z, cmulti *x, double *y)
{
    int i;
    for(i=0; i<n; i++){ cmul_cd(z[i],x,y[i]); }
}

/**
 @brief z=x.*y
 */
void dscalar_mul_cvec(int n, cmulti **z, double x, cmulti **y)
{
    int i;
    for(i=0; i<n; i++){ cmul_dc(z[i],x,y[i]); }
}

/**
 @brief z=x.*y
 */
void rscalar_mul_zvec(int n, cmulti **z, rmulti *x, dcomplex *y)
{
    int i;
    for(i=0; i<n; i++){ cmul_rz(z[i],x,y[i]); }
}

/**
 @brief z=x.*y
 */
void zscalar_mul_rvec(int n, cmulti **z, dcomplex x, rmulti **y)
{
    int i;
    for(i=0; i<n; i++){ cmul_zr(z[i],x,y[i]); }
}

///////////////////////////////////////////////////////

/**
 @brief z=x./y
 */
void cvec_div_cvec(int n, cmulti **z, cmulti **x, cmulti **y)
{
    int i;
    for(i=0; i<n; i++){ cdiv_cc(z[i],x[i],y[i]); }
}

/**
 @brief z=x./y
 */
void cvec_div_rvec(int n, cmulti **z, cmulti **x, rmulti **y)
{
    int i;
    for(i=0; i<n; i++){ cdiv_cr(z[i],x[i],y[i]); }
}

/**
 @brief z=x./y
 */
void rvec_div_cvec(int n, cmulti **z, rmulti **x, cmulti **y)
{
    int i;
    for(i=0; i<n; i++){ cdiv_rc(z[i],x[i],y[i]); }
}

/**
 @brief z=x./y
 */
void cvec_div_zvec(int n, cmulti **z, cmulti **x, dcomplex *y)
{
    int i;
    for(i=0; i<n; i++){ cdiv_cz(z[i],x[i],y[i]); }
}

/**
 @brief z=x./y
 */
void zvec_div_cvec(int n, cmulti **z, dcomplex *x, cmulti **y)
{
    int i;
    for(i=0; i<n; i++){ cdiv_zc(z[i],x[i],y[i]); }
}

/**
 @brief z=x./y
 */
void cvec_div_dvec(int n, cmulti **z, cmulti **x, double *y)
{
    int i;
    for(i=0; i<n; i++){ cdiv_cd(z[i],x[i],y[i]); }
}

/**
 @brief z=x./y
 */
void dvec_div_cvec(int n, cmulti **z, double *x, cmulti **y)
{
    int i;
    for(i=0; i<n; i++){ cdiv_dc(z[i],x[i],y[i]); }
}

/**
 @brief z=x./y
 */
void rvec_div_zvec(int n, cmulti **z, rmulti **x, dcomplex *y)
{
    int i;
    for(i=0; i<n; i++){ cdiv_rz(z[i],x[i],y[i]); }
}

/**
 @brief z=x./y
 */
void zvec_div_rvec(int n, cmulti **z, dcomplex *x, rmulti **y)
{
    int i;
    for(i=0; i<n; i++){ cdiv_zr(z[i],x[i],y[i]); }
}

///////////////////////////////////////////////////////

/**
 @brief z=x./y
 */
void cvec_div_cscalar(int n, cmulti **z, cmulti **x, cmulti *y)
{
    int i;
    for(i=0; i<n; i++){ cdiv_cc(z[i],x[i],y); }
}

/**
 @brief z=x./y
 */
void cvec_div_rscalar(int n, cmulti **z, cmulti **x, rmulti *y)
{
    int i;
    for(i=0; i<n; i++){ cdiv_cr(z[i],x[i],y); }
}

/**
 @brief z=x./y
 */
void rvec_div_cscalar(int n, cmulti **z, rmulti **x, cmulti *y)
{
    int i;
    for(i=0; i<n; i++){ cdiv_rc(z[i],x[i],y); }
}

/**
 @brief z=x./y
 */
void cvec_div_zscalar(int n, cmulti **z, cmulti **x, dcomplex y)
{
    int i;
    for(i=0; i<n; i++){ cdiv_cz(z[i],x[i],y); }
}

/**
 @brief z=x./y
 */
void zvec_div_cscalar(int n, cmulti **z, dcomplex *x, cmulti *y)
{
    int i;
    for(i=0; i<n; i++){ cdiv_zc(z[i],x[i],y); }
}

/**
 @brief z=x./y
 */
void cvec_div_dscalar(int n, cmulti **z, cmulti **x, double y)
{
    int i;
    for(i=0; i<n; i++){ cdiv_cd(z[i],x[i],y); }
}

/**
 @brief z=x./y
 */
void dvec_div_cscalar(int n, cmulti **z, double *x, cmulti *y)
{
    int i;
    for(i=0; i<n; i++){ cdiv_dc(z[i],x[i],y); }
}

/**
 @brief z=x./y
 */
void rvec_div_zscalar(int n, cmulti **z, rmulti **x, dcomplex y)
{
    int i;
    for(i=0; i<n; i++){ cdiv_rz(z[i],x[i],y); }
}

/**
 @brief z=x./y
 */
void zvec_div_rscalar(int n, cmulti **z, dcomplex *x, rmulti *y)
{
    int i;
    for(i=0; i<n; i++){ cdiv_zr(z[i],x[i],y); }
}



///////////////////////////////////////////////////////

/**
 @brief z=x./y
 */
void cscalar_div_cvec(int n, cmulti **z, cmulti *x, cmulti **y)
{
    int i;
    for(i=0; i<n; i++){ cdiv_cc(z[i],x,y[i]); }
}

/**
 @brief z=x./y
 */
void cscalar_div_rvec(int n, cmulti **z, cmulti *x, rmulti **y)
{
    int i;
    for(i=0; i<n; i++){ cdiv_cr(z[i],x,y[i]); }
}

/**
 @brief z=x./y
 */
void rscalar_div_cvec(int n, cmulti **z, rmulti *x, cmulti **y)
{
    int i;
    for(i=0; i<n; i++){ cdiv_rc(z[i],x,y[i]); }
}

/**
 @brief z=x./y
 */
void cscalar_div_zvec(int n, cmulti **z, cmulti *x, dcomplex *y)
{
    int i;
    for(i=0; i<n; i++){ cdiv_cz(z[i],x,y[i]); }
}

/**
 @brief z=x./y
 */
void zscalar_div_cvec(int n, cmulti **z, dcomplex x, cmulti **y)
{
    int i;
    for(i=0; i<n; i++){ cdiv_zc(z[i],x,y[i]); }
}

/**
 @brief z=x./y
 */
void cscalar_div_dvec(int n, cmulti **z, cmulti *x, double *y)
{
    int i;
    for(i=0; i<n; i++){ cdiv_cd(z[i],x,y[i]); }
}

/**
 @brief z=x./y
 */
void dscalar_div_cvec(int n, cmulti **z, double x, cmulti **y)
{
    int i;
    for(i=0; i<n; i++){ cdiv_dc(z[i],x,y[i]); }
}

/**
 @brief z=x./y
 */
void rscalar_div_zvec(int n, cmulti **z, rmulti *x, dcomplex *y)
{
    int i;
    for(i=0; i<n; i++){ cdiv_rz(z[i],x,y[i]); }
}

/**
 @brief z=x./y
 */
void zscalar_div_rvec(int n, cmulti **z, dcomplex x, rmulti **y)
{
    int i;
    for(i=0; i<n; i++){ cdiv_zr(z[i],x,y[i]); }
}

///////////////////////////////////////////////////////



//////////////////////////////////////////////////////////////////////////////////

/**
 @brief cmulti型のベクトルの要素ごとの掛け算の加算 z+=x.*y
*/
void cvec_add_mul_cvec_cvec(int n, cmulti **z, cmulti **x, cmulti **y)
{
  int i;
  for(i=0; i<n; i++){ cadd_mul_cc(z[i],x[i],y[i]); }
}

/**
 @brief cmulti型のベクトルxとスカラーyの掛け算の加算 z+=x*y
*/
void cvec_add_mul_cvec_cscalar(int n, cmulti **z, cmulti **x, cmulti *y)
{
  int i;
  for(i=0; i<n; i++){ cadd_mul_cc(z[i],x[i],y); }
}

/**
 @brief cmulti型のベクトルxとスカラーyの掛け算の加算 z+=x*y
*/
void cvec_add_mul_cvec_rscalar(int n, cmulti **z, cmulti **x, rmulti *y)
{
  int i;
  for(i=0; i<n; i++){ cadd_mul_cr(z[i],x[i],y); }
}

/**
 @brief cmulti型のベクトルxとスカラーyの掛け算の加算 z+=x*y
*/
void cvec_add_mul_cvec_zscalar(int n, cmulti **z, cmulti **x, dcomplex y)
{
  int i;
  for(i=0; i<n; i++){ cadd_mul_cz(z[i],x[i],y); }
}

/**
 @brief cmulti型のベクトルxとスカラーyの掛け算の加算 z+=x*y
*/
void cvec_add_mul_cvec_dscalar(int n, cmulti **z, cmulti **x, double y)
{
  int i;
  for(i=0; i<n; i++){ cadd_mul_cd(z[i],x[i],y); }
}

/**
 @brief cmulti型のベクトルの要素ごとの掛け算の減算 z-=x.*y
*/
void cvec_sub_mul_cvec_cvec(int n, cmulti **z, cmulti **x, cmulti **y)
{
  int i;
  for(i=0; i<n; i++){ csub_mul_cc(z[i],x[i],y[i]); }
}

/**
 @brief cmulti型のベクトルxとスカラーyの掛け算の減算 z-=x*y
*/
void cvec_sub_mul_cvec_cscalar(int n, cmulti **z, cmulti **x, cmulti *y)
{
  int i;
  for(i=0; i<n; i++){ csub_mul_cc(z[i],x[i],y); }
}

/**
 @brief cmulti型のベクトルxとスカラーyの掛け算の減算 z-=x*y
*/
void cvec_sub_mul_r(int n, cmulti **z, cmulti **x, rmulti *y)
{
  int i;
  for(i=0; i<n; i++){ csub_mul_cr(z[i],x[i],y); }
}

/**
 @brief cmulti型のベクトルxとスカラーyの掛け算の減算 z-=x*y
*/
void cvec_sub_mul_cvec_zscalar(int n, cmulti **z, cmulti **x, dcomplex y)
{
  int i;
  for(i=0; i<n; i++){ csub_mul_cz(z[i],x[i],y); }
}

/**
 @brief cmulti型のベクトルxとスカラーyの掛け算の減算 z-=x*y
*/
void cvec_sub_mul_cvec_dscalar(int n, cmulti **z, cmulti **x, double y)
{
  int i;
  for(i=0; i<n; i++){ csub_mul_cd(z[i],x[i],y); }
}

/**
 @brief cmulti型のベクトルの要素ごとの掛け算 z=conj(x).*y
*/
void cvec_dot(int n, cmulti **z, cmulti **x, cmulti **y)
{
  int i;
  for(i=0; i<n; i++){ cdot_cc(z[i],x[i],y[i]); }
}

/**
 @brief cmulti型のベクトルyとスカラーxの掛け算 z=conj(x)*y
*/
void cvec_dot_c1(int n, cmulti **z, cmulti *x, cmulti **y)
{
  int i;
  for(i=0; i<n; i++){ cdot_cc(z[i],x,y[i]); }
}

/**
 @brief cmulti型のベクトルxとスカラーyの掛け算 z=conj(x)*y
*/
void cvec_dot_c2(int n, cmulti **z, cmulti **x, cmulti *y)
{
  int i;
  for(i=0; i<n; i++){ cdot_cc(z[i],x[i],y); }
}

/**
 @brief cmulti型のベクトルyとスカラーxの掛け算 z=conj(x)*y
*/
void cvec_dot_z1(int n, cmulti **z, dcomplex x, cmulti **y)
{
  int i;
  for(i=0; i<n; i++){ cdot_zc(z[i],x,y[i]); }
}

/**
 @brief cmulti型のベクトルxとスカラーyの掛け算 z=conj(x)*y
*/
void cvec_dot_z2(int n, cmulti **z, cmulti **x, dcomplex y)
{
  int i;
  for(i=0; i<n; i++){ cdot_cz(z[i],x[i],y); }
}

/**
 @brief cmulti型のベクトルの要素ごとの掛け算の加算 z+=conj(x).*y
*/
void cvec_add_dot_cvec_cvec(int n, cmulti **z, cmulti **x, cmulti **y)
{
  int i;
  for(i=0; i<n; i++){ cadd_dot_cc(z[i],x[i],y[i]); }
}

/**
 @brief cmulti型のベクトルyとスカラーxの掛け算の加算 z+=conj(x)*y
*/
void cvec_add_dot_cscalar_cvec(int n, cmulti **z, cmulti *x, cmulti **y)
{
  int i;
  for(i=0; i<n; i++){ cadd_dot_cc(z[i],x,y[i]); }
}

/**
 @brief cmulti型のベクトルxとスカラーyの掛け算の加算 z+=conj(x)*y
*/
void cvec_add_dot_cvec_cscalar(int n, cmulti **z, cmulti **x, cmulti *y)
{
  int i;
  for(i=0; i<n; i++){ cadd_dot_cc(z[i],x[i],y); }
}

/**
 @brief cmulti型のベクトルyとスカラーxの掛け算の加算 z+=conj(x)*y
*/
void cvec_add_dot_zscalar_cvec(int n, cmulti **z, dcomplex x, cmulti **y)
{
  int i;
  for(i=0; i<n; i++){ cadd_dot_zc(z[i],x,y[i]); }
}

/**
 @brief cmulti型のベクトルxとスカラーyの掛け算の加算 z+=conj(x)*y
*/
void cvec_add_dot_cvec_zscalar(int n, cmulti **z, cmulti **x, dcomplex y)
{
  int i;
  for(i=0; i<n; i++){ cadd_dot_cz(z[i],x[i],y); }
}

/**
 @brief cmulti型のベクトルの要素ごとの掛け算の減算 z-=conj(x).*y
*/
void cvec_sub_dot_cvec_cvec(int n, cmulti **z, cmulti **x, cmulti **y)
{
  int i;
  for(i=0; i<n; i++){ csub_dot_cc(z[i],x[i],y[i]); }
}

/**
 @brief cmulti型のベクトルyとスカラーxの掛け算の減算 z-=conj(x)*y
*/
void cvec_sub_dot_cscalar_cvec(int n, cmulti **z, cmulti *x, cmulti **y)
{
  int i;
  for(i=0; i<n; i++){ csub_dot_cc(z[i],x,y[i]); }
}

/**
 @brief cmulti型のベクトルxとスカラーyの掛け算の減算 z-=conj(x)*y
*/
void cvec_sub_dot_cvec_cscalar(int n, cmulti **z, cmulti **x, cmulti *y)
{
  int i;
  for(i=0; i<n; i++){ csub_dot_cc(z[i],x[i],y); }
}

/**
 @brief cmulti型のベクトルyとスカラーxの掛け算の減算 z-=conj(x)*y
*/
void cvec_sub_dot_zscalar_cvec(int n, cmulti **z, dcomplex x, cmulti **y)
{
  int i;
  for(i=0; i<n; i++){ csub_dot_zc(z[i],x,y[i]); }
}

/**
 @brief cmulti型のベクトルxとスカラーyの掛け算の減算 z-=conj(x)*y
*/
void cvec_sub_dot_cvec_zscalar(int n, cmulti **z, cmulti **x, dcomplex y)
{
  int i;
  for(i=0; i<n; i++){ csub_dot_cz(z[i],x[i],y); }
}

/**
 @brief cmulti型のベクトルxに関する行列Aによる1次変換 y=A*x
 @param[in]     m   ベクトルyのサイズ.行列Aの行のサイズ.
 @param[in]     n   ベクトルxのサイズ.行列Aの列のサイズ.
 @param[in,out] y   [in]初期化済みのサイズmのベクトル.[out]y=A*xの計算結果.
 @param[in]     A   サイズ(m,n)の行列.
 @param[in]     LDA 行列Aの第1次元.
 @param[in]     x   サイズnのベクトル.
*/
void cvec_lintr(int m, int n, cmulti **y, cmulti **A, int LDA, cmulti **x)
{
  int i,j;
  cmulti **z=NULL;
  z=cvec_allocate_prec(m,cvec_get_prec_max(m,y));
  for(i=0; i<m; i++){
    cset_zero(z[i]); // z=0
    for(j=0; j<n; j++){
      cadd_mul_cc(z[i],MAT(A,i,j,LDA),x[j]); // z+=A*x
    }
  }
  cvec_copy_cvec(m,y,z); // y=z
  z=cvec_free(m,z);
}

/**
 @brief cmulti型のベクトルxに関する行列Aによる1次変換の加算 y+=A*x
 @param[in]     m   ベクトルyのサイズ.行列Aの行のサイズ.
 @param[in]     n   ベクトルxのサイズ.行列Aの列のサイズ.
 @param[in,out] y   [in]初期化済みのサイズmのベクトル.[out]y+=A*xの計算結果.
 @param[in]     A   サイズ(m,n)の行列.
 @param[in]     LDA 行列Aの第1次元.
 @param[in]     x   サイズnのベクトル.
*/
void cvec_add_lintr(int m, int n, cmulti **y, cmulti **A, int LDA, cmulti **x)
{
  int i,j;
  cmulti **z=NULL;
  z=cvec_allocate_prec(m,cvec_get_prec_max(m,y));
  for(i=0; i<m; i++){
    cset_c(z[i],y[i]); // z=y
    for(j=0; j<n; j++){
      cadd_mul_cc(z[i],MAT(A,i,j,LDA),x[j]); // z+=A*x
    }
  }
  cvec_copy_cvec(m,y,z); // y=z
  z=cvec_free(m,z);
}

/**
 @brief cmulti型のベクトルxに関する行列Aによる1次変換の減算 y-=A*x
 @param[in]     m   ベクトルyのサイズ.行列Aの行のサイズ.
 @param[in]     n   ベクトルxのサイズ.行列Aの列のサイズ.
 @param[in,out] y   [in]初期化済みのサイズmのベクトル.[out]y-=A*xの計算結果.
 @param[in]     A   サイズ(m,n)の行列.
 @param[in]     LDA 行列Aの第1次元.
 @param[in]     x   サイズnのベクトル.
*/
void cvec_sub_lintr(int m, int n, cmulti **y, cmulti **A, int LDA, cmulti **x)
{
  int i,j;
  cmulti **z=NULL;
  z=cvec_allocate_prec(m,cvec_get_prec_max(m,y));
  for(i=0; i<m; i++){
    cset_c(z[i],y[i]); // z=y
    for(j=0; j<n; j++){
      csub_mul_cc(z[i],MAT(A,i,j,LDA),x[j]); // z-=A*x
    }
  }
  cvec_copy_cvec(m,y,z); // y=z
  z=cvec_free(m,z);
}

/**
 @brief cmulti型のベクトルxに関する転置行列Aによる1次変換 y=A^T*x
 @param[in]     m   ベクトルxのサイズ.行列Aの行のサイズ.
 @param[in]     n   ベクトルyのサイズ.行列Aの列のサイズ.
 @param[in,out] y   [in]初期化済みのサイズnのベクトル.[out]y=A^T*xの計算結果.
 @param[in]     A   サイズ(m,n)の行列.
 @param[in]     LDA 行列Aの第1次元.
 @param[in]     x   サイズmのベクトル.
*/
void cvec_lintr_t(int m, int n, cmulti **y, cmulti **A, int LDA, cmulti **x)
{
  int i,j;
  cmulti **z=NULL;
  z=cvec_allocate_prec(n,cvec_get_prec_max(n,y));
  for(j=0; j<n; j++){
    cset_zero(z[j]); // z=0
    for(i=0; i<m; i++){
      cadd_mul_cc(z[j],MAT(A,i,j,LDA),x[i]); // z=A^T*x
    }
  }
  cvec_copy_cvec(n,y,z); // y=z
  z=cvec_free(n,z);
}

/**
 @brief cmulti型のベクトルxに関する転置行列Aによる1次変換の加算 y+=A^T*x
 @param[in]     m   ベクトルxのサイズ.行列Aの行のサイズ.
 @param[in]     n   ベクトルyのサイズ.行列Aの列のサイズ.
 @param[in,out] y   [in]初期化済みのサイズnのベクトル.[out]y+=A^T*xの計算結果.
 @param[in]     A   サイズ(m,n)の行列.
 @param[in]     LDA 行列Aの第1次元.
 @param[in]     x   サイズmのベクトル.
*/
void cvec_add_lintr_t(int m, int n, cmulti **y, cmulti **A, int LDA, cmulti **x)
{
  int i,j;
  cmulti **z=NULL;
  z=cvec_allocate_prec(n,cvec_get_prec_max(n,y));
  for(j=0; j<n; j++){
    cset_c(z[j],y[j]); // z=y
    for(i=0; i<m; i++){
      cadd_mul_cc(z[j],MAT(A,i,j,LDA),x[i]); // z+=A^T*x
    }
  }
  cvec_copy_cvec(n,y,z); // y=z
  z=cvec_free(n,z);
}

/**
 @brief cmulti型のベクトルxに関する転置行列Aによる1次変換の減算 y-=A^T*x
 @param[in]     m   ベクトルxのサイズ.行列Aの行のサイズ.
 @param[in]     n   ベクトルyのサイズ.行列Aの列のサイズ.
 @param[in,out] y   [in]初期化済みのサイズnのベクトル.[out]y-=A^T*xの計算結果.
 @param[in]     A   サイズ(m,n)の行列.
 @param[in]     LDA 行列Aの第1次元.
 @param[in]     x   サイズmのベクトル.
*/
void cvec_sub_lintr_t(int m, int n, cmulti **y, cmulti **A, int LDA, cmulti **x)
{
  int i,j;
  cmulti **z=NULL;
  z=cvec_allocate_prec(n,cvec_get_prec_max(n,y));
  for(j=0; j<n; j++){
    cset_c(z[j],y[j]); // z=y
    for(i=0; i<m; i++){
      csub_mul_cc(z[j],MAT(A,i,j,LDA),x[i]); // z-=A^T*x
    }
  }
  cvec_copy_cvec(n,y,z); // y=z
  z=cvec_free(n,z);
}

/**
 @brief cmulti型のベクトルxに関する共役転置行列Aによる1次変換 y=A'*x
 @param[in]     m   ベクトルxのサイズ.行列Aの行のサイズ.
 @param[in]     n   ベクトルyのサイズ.行列Aの列のサイズ.
 @param[in,out] y   [in]初期化済みのサイズnのベクトル.[out]y=A'*xの計算結果.
 @param[in]     A   サイズ(m,n)の行列.
 @param[in]     LDA 行列Aの第1次元.
 @param[in]     x   サイズmのベクトル.
*/
void cvec_lintr_ct(int m, int n, cmulti **y, cmulti **A, int LDA, cmulti **x)
{
  int i,j;
  cmulti **z=NULL;
  z=cvec_allocate_prec(n,cvec_get_prec_max(n,y));
  for(j=0; j<n; j++){
    cset_zero(z[j]); // z=0
    for(i=0; i<m; i++){
      cadd_dot_cc(z[j],MAT(A,i,j,LDA),x[i]); // z=A'*x
    }
  }
  cvec_copy_cvec(n,y,z); // y=z
  z=cvec_free(n,z);
}

/**
 @brief cmulti型のベクトルxに関する共役転置行列Aによる1次変換の加算 y+=A'*x
 @param[in]     m   ベクトルxのサイズ.行列Aの行のサイズ.
 @param[in]     n   ベクトルyのサイズ.行列Aの列のサイズ.
 @param[in,out] y   [in]初期化済みのサイズnのベクトル.[out]y+=A'*xの計算結果.
 @param[in]     A   サイズ(m,n)の行列.
 @param[in]     LDA 行列Aの第1次元.
 @param[in]     x   サイズmのベクトル.
*/
void cvec_add_lintr_ct(int m, int n, cmulti **y, cmulti **A, int LDA, cmulti **x)
{
  int i,j;
  cmulti **z=NULL;
  z=cvec_allocate_prec(n,cvec_get_prec_max(n,y));
  for(j=0; j<n; j++){
    cset_c(z[j],y[j]); // z=y
    for(i=0; i<m; i++){
      cadd_dot_cc(z[j],MAT(A,i,j,LDA),x[i]); // z+=A'*x
    }
  }
  cvec_copy_cvec(n,y,z); // y=z
  z=cvec_free(n,z);
}

/**
 @brief cmulti型のベクトルxに関する共役転置行列Aによる1次変換の減算 y-=A'*x
 @param[in]     m   ベクトルxのサイズ.行列Aの行のサイズ.
 @param[in]     n   ベクトルyのサイズ.行列Aの列のサイズ.
 @param[in,out] y   [in]初期化済みのサイズnのベクトル.[out]y-=A'*xの計算結果.
 @param[in]     A   サイズ(m,n)の行列.
 @param[in]     LDA 行列Aの第1次元.
 @param[in]     x   サイズmのベクトル.
*/
void cvec_sub_lintr_ct(int m, int n, cmulti **y, cmulti **A, int LDA, cmulti **x)
{
  int i,j;
  cmulti **z=NULL;
  z=cvec_allocate_prec(n,cvec_get_prec_max(n,y));
  for(j=0; j<n; j++){
    cset_c(z[j],y[j]); // z=y
    for(i=0; i<m; i++){
      csub_dot_cc(z[j],MAT(A,i,j,LDA),x[i]); // z-=A'*x
    }
  }
  cvec_copy_cvec(n,y,z); // y=z
  z=cvec_free(n,z);
}

/**
 @brief cmulti型のベクトルの要素の総和 value=sum(x)
 */
void csum_cvec(cmulti *value, int n, cmulti **x)
{
  int i;
  cset_zero(value);
  for(i=0; i<n; i++){ cadd_cc(value,value,x[i]); }
}

/**
 @brief cmulti型のベクトルの要素の絶対値の平方の総和 value=sum(abs(x).^2)
*/
void rsum_pow2_abs_cvec(rmulti *value, int n, cmulti **x)
{
  int i;
  rset_zero(value);
  for(i=0; i<n; i++){ radd_abs2_c(value,x[i]); }
}

/**
 @brief cmulti型のベクトルの内積 value=sum(conj(x).*y)
 */
void cvec_sum_dot(cmulti *value, int n, cmulti **x, cmulti **y)
{
  int i;
  cset_zero(value);
  for(i=0; i<n; i++){ cadd_dot_cc(value,x[i],y[i]); }
}

/**
 @brief cmulti型のベクトルの要素の最大値 value=max(x)
*/
void cmax_cvec(cmulti *value, int n, cmulti **x)
{
  int i;
  cset_c(value,x[0]);     // value=x[0]
  for(i=1; i<n; i++){
    if(gt_cc(x[i],value)){ // x[i]<value
      cset_c(value,x[i]); // value=x[i]
    }
  }
}

/**
 @brief rmulti型のベクトルの要素の最小値 value=min(x)
*/
void cmin_cvec(cmulti *value, int n, cmulti **x)
{
  int i;
  cset_c(value,x[0]);     // value=x[0]
  for(i=1; i<n; i++){
    if(lt_cc(x[i],value)){ // x[i]<value
      cset_c(value,x[i]); // value=x[i]
    }
  }
}

/**
 @brief cmulti型のべき乗 z=x.^y
*/
void cvec_pow(int n, cmulti **z, cmulti **x, cmulti **y) 
{
  int i;
  for(i=0; i<n; i++){ cpow_cc(z[i],x[i],y[i]); }
}

/**
 @brief cmulti型のべき乗 z=x^y
*/
void cvec_pow_si(int n, cmulti **z, cmulti **x, int y)
{
  int i;
  for(i=0; i<n; i++){ cpow_c(z[i],x[i],y); }
}

/**
 @brief cmulti型のべき乗 z=x^y
*/
void cvec_pow_r(int n, cmulti **z, cmulti **x, rmulti *y)
{
  int i;
  for(i=0; i<n; i++){ cpow_cr(z[i],x[i],y); }
}

/**
 @brief cmulti型のベクトルの要素ごとの絶対値の対数 y=log2(abs(x))
*/
void rvec_log2_abs_cvec(int n, rmulti **y, cmulti **x)
{
  int i;
  for(i=0; i<n; i++){
    rabs_c(y[i],x[i]);
    rlog2_r(y[i],y[i]);
  }
}

/**
 @brief cmulti型のベクトルの規格化 y=x/sqrt(x'*x)
*/
void cvec_normalize_cvec(int n, cmulti **y, cmulti **x)
{
  rmulti *a=NULL;
  a=rallocate_prec(cvec_get_prec_max(n,y));
  rnorm2_cvec(a,n,x);   // a=sqrt(x'*x)
  rinv_r(a,a);           // a=1/sqrt(x'*x)
  cvec_mul_rscalar(n,y,x,a); // y=x/sqrt(x'*x)
  a=rfree(a);
}

/**
 @brief cmulti型のベクトルの規格化 y=x/sqrt(x'*x)
*/
void cvec_normalize_sgn_cvec(int n, cmulti **y, cmulti **x)
{
  int prec,k=-1;
  rmulti *a=NULL;
  cmulti *b=NULL;
  // allocate
  prec=cvec_get_prec_max(n,y);
  a=rallocate_prec(prec);
  b=callocate_prec(prec);
  // compute
  rmax_abs_cvec_index(a,n,x,&k); // a=abs(x[k])
  cnormalize_c(b,x[k]);     // b=sgn(x[k])
  cconj_c(b,b);             // b=conj(sgn(x[k]))
  rnorm2_cvec(a,n,x);      // a=sqrt(x'*x)
  cdiv_cr(b,b,a);         // b=conj(sgn(x[k]))/sqrt(x'*x)
  cvec_mul_cscalar(n,y,x,b);    // y=x*(conj(sgn(x[k]))/sqrt(x'*x))
  // done
  a=rfree(a);
  b=cfree(b);
}

/**
 @brief cmulti型のベクトルの直交化 y-=(x'*y)*x where x'*x=1
*/
void cvec_orthogonalize(int n, cmulti **y, cmulti **x)
{
  cmulti *a=NULL;
  a=callocate_prec(cvec_get_prec_max(n,y));
  cvec_sum_dot(a,n,x,y);   // a=x'*y
  cvec_sub_mul_cvec_cscalar(n,y,x,a); // y-=x*a
  a=cfree(a);
}

/**
 @brief cmulti型のベクトルの2乗ノルム value=sqrt(sum(abs(x).^2))
*/
void rnorm2_cvec(rmulti *value, int n, cmulti **x)
{
  rsum_pow2_abs_cvec(value,n,x); // value=sum(abs(x).^2)
  rsqrt_r(value,value);       // value=sqrt(value)
}

/**
 @brief cmulti型のベクトルの要素の平均 value=sum(x)/n
*/
void caverage_cvec(cmulti *value, int n, cmulti **x)
{
  csum_cvec(value,n,x);
  cdiv_cd(value,value,n);
}

/**
 @brief cmulti型のベクトルの要素の差の絶対値の平方の最大値 value=max(abs(x-y).^2)
*/
void cvec_max_abs2_sub(rmulti *value, int n, cmulti **x, cmulti **y)
{
  int i;
  cmulti *a=NULL;
  rmulti *b=NULL;
  a=callocate_prec(rget_prec(value));
  b=rallocate_prec(rget_prec(value));
  csub_cc(a,x[0],y[0]);  // value=x[0]-y[0]
  rabs2_c(value,a);     // value=(abs(x[0]-y[0]))^2
  for(i=1; i<n; i++){
    csub_cc(a,x[i],y[i]); // a=x[i]-y[i]
    rabs2_c(b,a);        // b=(abs(x[i]-y[i]))^2
    if(gt_rr(b,value)){
      rset_r(value,b);     // value=b
    }
  }
  a=cfree(a);
  b=rfree(b);
}

/**
 @brief cmulti型のベクトルの要素の絶対値 y=abs(x)
*/
void rvec_abs_cvec(int n, rmulti **y, cmulti **x)
{
  int i;
  for(i=0; i<n; i++){ rabs_c(y[i],x[i]); }
}

/**
 @brief cmulti型のベクトルの要素の差の絶対値 z=abs(x-y)
*/
void cvec_abs_sub(int n, rmulti **z, cmulti **x, cmulti **y)
{
  int i;
  for(i=0; i<n; i++){
    rabs_sub_cc(z[i],x[i],y[i]);
  }
}

/**
 @brief cmulti型のベクトルの要素とスカラーの差の絶対値 z=abs(x-y)
*/
void cvec_abs_sub_c(int n, rmulti **z, cmulti **x, cmulti *y)
{
  int i;
  for(i=0; i<n; i++){
    rabs_sub_cc(z[i],x[i],y);
  }
}

/**
 @brief cmulti型のベクトルの要素とスカラーの差の絶対値 z=abs(x-y)
*/
void cvec_abs_sub_r(int n, rmulti **z, cmulti **x, rmulti *y)
{
  int i;
  for(i=0; i<n; i++){
    rabs_sub_cr(z[i],x[i],y);
  }
}

/**
 @brief rmulti型のベクトルの要素の絶対値の総和 value=sum(abs(x))
 */
void rsum_abs_cvec(rmulti *value, int n, cmulti **x)
{
  int i;
  rset_zero(value); // value=0
  for(i=0; i<n; i++){ radd_abs_c(value,x[i]); } // value+=abs(x[i])
}

/**
 @brief cmulti型のベクトルの要素の差の絶対値の総和 value=sum(abs(x-y))
*/
void cvec_sum_abs_sub(rmulti *value, int n, cmulti **x, cmulti **y)
{
  int i;
  cmulti *a=NULL;
  a=callocate_prec(rget_prec(value));
  rset_zero(value); // value=0
  for(i=0; i<n; i++){
    csub_cc(a,x[i],y[i]); // a=x[i]-y[i]
    radd_abs_c(value,a); // value+=abs(x[i]-y[i])
  }
  a=cfree(a);
}

/**
 @brief cmulti型のベクトルの要素の絶対値の最大値 value=max(abs(x))
*/
void rmax_abs_cvec(rmulti *value, int n, cmulti **x)
{
  int i;
  rmulti *a=NULL;
  a=rallocate_prec(rget_prec(value));
  rabs2_c(value,x[0]);  // value=abs(x[0])^2
  for(i=1; i<n; i++){
    rabs2_c(a,x[i]);    // a=abs(x[i])^2
    if(gt_rr(a,value)){    // a>value
      rset_r(value,a); // value=a
    }
  }
  rsqrt_r(value,value); // value=sqrt(value)
  a=rfree(a);            // free
}

/**
 @brief cmulti型のベクトルの要素の実部，虚部の絶対値の最大値 value=max(abs(x))
*/
void rmax_max_absc_cvec(rmulti *value, int n, cmulti **x)
{
  int i;
  rmulti *a=NULL;
  a=rallocate_prec(rget_prec(value));
  rabs_r(value,C_R(x[0]));
  rabs_r(a,C_I(x[0]));
  if(gt_rr(a,value)){ rset_r(value,a); }
  for(i=1; i<n; i++){
    rabs_r(a,C_R(x[i]));
    if(gt_rr(a,value)){ rset_r(value,a); }
    rabs_r(a,C_I(x[i]));
    if(gt_rr(a,value)){ rset_r(value,a); }
  }
  a=rfree(a);
}

/**
 @brief cmulti型のベクトルの要素の絶対値の最大値 value=max(abs(x)) とその添え字 I
*/
void rmax_abs_cvec_index(rmulti *value, int n, cmulti **x, int *I)
{
  int i;
  rmulti *a=NULL;
  a=rallocate_prec(rget_prec(value));
  rabs2_c(value,x[0]);   // value=abs(x[0])^2
  if(I!=NULL){ (*I)=0; }  // I=i
  for(i=1; i<n; i++){
    rabs2_c(a,x[i]);     // a=abs(x[i])^2
    if(gt_rr(a,value)){     // a>value
      rset_r(value,a);  // value=a
      if(I!=NULL){ (*I)=i; } // I=i
    }
  }
  rsqrt_r(value,value); // value=sqrt(value)
  a=rfree(a);
}


/**
 @brief cmulti型のベクトルの要素の差の絶対値の最大値 value=max(abs(x-y))
*/
void cvec_max_abs_sub(rmulti *value, int n, cmulti **x, cmulti **y)
{
  int i;
  cmulti *a=NULL;
  rmulti *b=NULL;
  a=callocate_prec(rget_prec(value));
  b=rallocate_prec(rget_prec(value));
  csub_cc(a,x[0],y[0]);  // a=x[0]-y[0]
  rabs2_c(value,a);     // value=abs(x[0]-y[0])^2
  for(i=1; i<n; i++){
    csub_cc(a,x[i],y[i]); // a=x[i]-y[i]
    rabs2_c(b,a);        // b=abs(x[i]-y[i])^2
    if(gt_rr(b,value)){
      rset_r(value,b);  // value=b
    }
  }
  rsqrt_r(value,value); // value=sqrt(value)
  a=cfree(a);
  b=rfree(b);
}

/**
 @brief cmulti型のベクトルの要素の絶対値の最小値 value=min(abs(x))
*/
void rmin_abs_cvec(rmulti *value, int n, cmulti **x)
{
  int i;
  rmulti *a=NULL;
  a=rallocate_prec(rget_prec(value));
  rabs2_c(value,x[0]);  // value=abs(x[0])^2
  for(i=1; i<n; i++){
    rabs2_c(a,x[i]);    // a=abs(x[i])^2
    if(lt_rr(a,value)){    // a<value
      rset_r(value,a); // value=a
    }
  }
  rsqrt_r(value,value); // value=sqrt(value)
  a=rfree(a);
}

/**
 @brief cmulti型のベクトルの要素の絶対値の最小値 value=min(abs(x)) とその添え字 I
*/
void rmin_abs_cvec_index(rmulti *value, int n, cmulti **x, int *I)
{
  int i;
  rmulti *a=NULL;
  a=rallocate_prec(rget_prec(value));
  rabs2_c(value,x[0]);   // value=abs(x[0])^2
  if(I!=NULL){ (*I)=0; }  // I=i
  for(i=1; i<n; i++){
    rabs2_c(a,x[i]);     // a=abs(x[i])^2
    if(lt_rr(a,value)){     // a<value
      rset_r(value,a);  // value=a
      if(I!=NULL){ (*I)=i; } // I=i
    }
  }
  rsqrt_r(value,value); // value=sqrt(value)
  a=rfree(a);
}


/**
 @brief cmulti型のベクトルの要素の偏角 y=arg(x)
*/
void rvec_arg_cvec(int n, rmulti **y, cmulti **x)
{
  int i;
  for(i=0; i<n; i++){ rarg_c(y[i],x[i]); }
}



/**
 @brief cmulti型のベクトルの方向余弦 value=x'*y/sqrt(x'*x)/sqrt(y'*y)
*/
void cvec_dcos(cmulti *value, int n, cmulti **x, cmulti **y)
{
  rmulti *a=NULL;
  a=rallocate_prec(cget_prec(value));
  cvec_sum_dot(value,n,x,y); // value=x'*y
  rnorm2_cvec(a,n,x);         // a=sqrt(x'*x)
  cdiv_cr(value,value,a);    // value=x'*y/sqrt(x'*x)
  rnorm2_cvec(a,n,y);         // a=sqrt(y'*y)
  cdiv_cr(value,value,a);    // value=x'*y/sqrt(x'*x)/sqrt(y'*y)
  a=rfree(a);
}

/**
 @brief cmulti型のベクトルの方向余弦の絶対値 value=abs(x'*y)/sqrt(x'*x)/sqrt(y'*y)
*/
void cvec_abs_dcos(rmulti *value, int n, cmulti **x, cmulti **y)
{
  cmulti *a=NULL;
  a=callocate_prec(rget_prec(value));
  cvec_dcos(a,n,x,y);
  rabs_c(value,a);
  a=cfree(a);
}

/**
 @brief cmulti型のベクトルの鋭角の角度 theta=acos(abs(x'*y)/sqrt(x'*x)/sqrt(y'*y))
*/
void cvec_angle(rmulti *theta, int n, cmulti **x, cmulti **y)
{
  rmulti *a=NULL;
  a=rallocate_prec(rget_prec(theta));
  cvec_abs_dcos(a,n,x,y); // a=abs(x'*y)/sqrt(x'*x)/sqrt(y'*y)
  if(gt_rd(a,1)){ rset_d(a,1); }
  racos_r(theta,a); // theta=acos(dcos)
  a=rfree(a);
}

/**
 @brief cmulti型のベクトルの角度[deg] theta=acos(abs(x'*y)/sqrt(x'*x)/sqrt(y'*y))
*/
void cvec_angle_deg(rmulti *theta, int n, cmulti **x, cmulti **y)
{
  cvec_angle(theta,n,x,y);
  rmul_rd(theta,theta,(180.0/M_PI));
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
    value=cmp_cc(x[i],y[i]);
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
  for(i=0; i<n; i++){ if(!eq_cc(x[i],y[i])){ return 0; } }
  return 1;
}

/**
 @brief cmulti型のベクトルの値の比較 X>Y
*/
int cvec_gt(int n, cmulti **x, cmulti **y)
{
  int i;
  for(i=0; i<n; i++){ if(!gt_cc(x[i],y[i])){ return 0; } }
  return 1;
}

/**
 @brief cmulti型のベクトルの値の比較 X>=Y
*/
int cvec_ge(int n, cmulti **x, cmulti **y)
{
  int i;
  for(i=0; i<n; i++){ if(!ge_cc(x[i],y[i])){ return 0; } }
  return 1;
}

/**
 @brief cmulti型のベクトルの値の比較 X<Y
*/
int cvec_lt(int n, cmulti **x, cmulti **y)
{
  int i;
  for(i=0; i<n; i++){ if(!lt_cc(x[i],y[i])){ return 0; } }
  return 1;
}

/**
 @brief cmulti型のベクトルの値の比較 X<=Y
*/
int cvec_le(int n, cmulti **x, cmulti **y)
{
  int i;
  for(i=0; i<n; i++){ if(!le_cc(x[i],y[i])){ return 0; } }
  return 1;
}

/**
 @brief cmulti型のベクトルの値の比較 X==Y
*/
int cvec_eqc(int n, cmulti **x, cmulti **y)
{
  int i;
  for(i=0; i<n; i++){ if(!ceq_cc(x[i],y[i])){ return 0; } }
  return 1;
}

/**
 @brief cmulti型のベクトルの値の比較 X>Y
*/
int cvec_gtc(int n, cmulti **x, cmulti **y)
{
  int i;
  for(i=0; i<n; i++){ if(!cgt_cc(x[i],y[i])){ return 0; } }
  return 1;
}

/**
 @brief cmulti型のベクトルの値の比較 X>=Y
*/
int cvec_gec(int n, cmulti **x, cmulti **y)
{
  int i;
  for(i=0; i<n; i++){ if(!cge_cc(x[i],y[i])){ return 0; } }
  return 1;
}

/**
 @brief cmulti型のベクトルの値の比較 X<Y
*/
int cvec_ltc(int n, cmulti **x, cmulti **y)
{
  int i;
  for(i=0; i<n; i++){ if(!clt_cc(x[i],y[i])){ return 0; } }
  return 1;
}

/**
 @brief cmulti型のベクトルの値の比較 X<=Y
*/
int cvec_lec(int n, cmulti **x, cmulti **y)
{
  int i;
  for(i=0; i<n; i++){ if(!cle_cc(x[i],y[i])){ return 0; } }
  return 1;
}

/** @} */



//EOF
