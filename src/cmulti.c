#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<stdarg.h>
#include"is_macros.h"
#include"is_rmulti.h"
#include"is_cmulti.h"
#include"is_strings.h"

/**
 @file  cmulti.c
 @brief 多倍長精度複素数cmulti型に関する関数の定義.
 @details 多倍長精度複素数cmulti型のベクトルに関しては@link cvec.c@endlink を参照.
          多倍長精度複素数cmulti型の行列に関しては@link cmat.c@endlink を参照.
 */

/////////////////////////////////////////////////////////

/*
 マクロ
 */

#define RAr(X,P) ((X)=rallocate_prec(rget_prec(P)))
#define RAc(X,P) ((X)=rallocate_prec(cget_prec(P)))
#define RF(X) ((X)=rfree(X))
#define CAr(X,P) ((X)=callocate_prec(rget_prec(P)))
#define CAc(X,P) ((X)=callocate_prec(cget_prec(P)))
#define CF(X) ((X)=cfree(X))

/////////////////////////////////////////////////////////

/** @name cmulti型の初期化に関する関数 */
/** @{ */

/**
 @brief cmulti型の新規生成.
*/
cmulti *callocate()
{
  cmulti *x=NULL;
  x=(cmulti*)malloc(sizeof(cmulti));  
  C_R(x)=rallocate();
  C_I(x)=rallocate();
  return x;
}

/**
 @brief cmulti型の精度を指定しての新規生成.
*/
cmulti *callocate_prec(int prec)
{
  cmulti *x=NULL;
  x=(cmulti*)malloc(sizeof(cmulti));  
  C_R(x)=rallocate_prec(prec);
  C_I(x)=rallocate_prec(prec);
  return x;
}

/**
 @brief cmulti型の複製の生成.
*/
cmulti *callocate_clone(cmulti *y)
{
  cmulti *x=NULL;
  x=(cmulti*)malloc(sizeof(cmulti));  
  C_R(x)=rallocate_clone(C_R(y));
  C_I(x)=rallocate_clone(C_I(y));
  return x;
}

/**
 @brief cmulti型のrmulti型からの複製の生成.
*/
cmulti *callocate_clone_r(rmulti *y)
{
  cmulti *x=NULL;
  x=(cmulti*)malloc(sizeof(cmulti));  
  C_R(x)=rallocate_clone(y);
  C_I(x)=rallocate();
  rset_zero(C_I(x));
  return x;
}

/**
 @brief cmulti型の終了処理.
*/
cmulti *cfree(cmulti *x)
{
  if(x==NULL) return NULL; 
  C_R(x)=rfree(C_R(x));
  C_I(x)=rfree(C_I(x));
  free(x);
  x=NULL;
  return x;
}

/**
 @brief 初期化済みのcmulti型の浮動小数点数の精度(ビット数)を変更し再初期化.
*/
int cround(cmulti *x, int prec)
{
  int e=0;
  e+=rround(C_R(x),prec);
  e+=rround(C_I(x),prec);
  return e;
}

/**
 @brief cmulti型の値を複成.
*/
int cclone(cmulti *y, cmulti *x)
{
  int e=0;
  e+=rclone(C_R(y),C_R(x));
  e+=rclone(C_I(y),C_I(x));
  return e;
}

/**
 @brief cmulti型の値をrmulti型から複成.
*/
int cclone_r(cmulti *y, rmulti *x)
{
  int e=0;
  e+=rclone(C_R(y),x);
  e+=rset_zero(C_I(y));
  return e;
}

/**
 @brief cmulti型の値をrmulti型から複成.
*/
int cclone_rr(cmulti *y, rmulti *x_r, rmulti *x_i)
{
  int e=0;
  e+=rclone(C_R(y),x_r);
  e+=rclone(C_I(y),x_i);
  return e;
}

/**
 @brief cmulti型の値の交換.
*/
void cswap(cmulti *x, cmulti *y)
{
  rswap(C_R(x),C_R(y));
  rswap(C_I(x),C_I(y));
}


/** @} */

/////////////////////////////////////////////////////////////

/** @name cmulti型のメンバ変数に関する関数 */
/** @{ */

/**
 @brief cmulti型の浮動小数点数の精度(ビット数)を取得.
*/
int cget_prec(cmulti *x)
{
  int pr,pi;
  pr=rget_prec(C_R(x));
  pi=rget_prec(C_I(x));
  return pr>=pi?pr:pi;
}

/**
 @brief cmulti型の浮動小数点数の指数部の取得.
*/
int cget_exp(cmulti *x)
{
  int pr,pi;
  pr=rget_exp(C_R(x));
  pi=rget_exp(C_I(x));
  return pr>=pi?pr:pi;
}

/**
 @brief cmulti型がNaNであるかの判定.
*/
int cis_nan(cmulti *x){
  return (ris_nan(C_R(x)) || ris_nan(C_I(x)));
}

/**
 @brief cmulti型がInfであるかの判定.
*/
int cis_inf(cmulti *x){
  return (ris_inf(C_R(x)) || ris_inf(C_I(x)));
}

/**
 @brief cmulti型が数であるかの判定.
*/
int cis_number(cmulti *x){
  return (ris_number(C_R(x)) && ris_number(C_I(x)));
}

/**
 @brief cmulti型が零であるかの判定.
*/
int cis_zero(cmulti *x){
  return (cis_number(x) && ris_zero(C_R(x)) && ris_zero(C_I(x)));
}

/**
 @brief cmulti型が実数であるかの判定.
*/
int cis_real(cmulti *x){
  return cis_number(x) && ris_zero(C_I(x));
}

/**
 @brief cmulti型が純虚数であるかの判定.
*/
int cis_pure_imaginary(cmulti *x){
  return cis_number(x) && ris_zero(C_R(x));
}

/**
 @brief cmulti型が虚数であるかの判定.
*/
int cis_imaginary(cmulti *x){
  return cis_number(x) && !ris_zero(C_I(x));
}

/** @} */

/////////////////////////////////////////////////////////////

/** @name cmulti型の入出力に関する関数 */
/** @{ */

/**
 @brief cmulti型のファイルへの保存
*/
void cbin_save(cmulti *x, FILE *fid)
{
  rbin_save(C_R(x),fid);
  rbin_save(C_I(x),fid);
}

/**
 @brief cmulti型のファイルからの読み込み
*/
cmulti *cbin_load(FILE *fid)
{
  cmulti *x=NULL;
  x=(cmulti*)malloc(sizeof(cmulti));
  C_R(x)=rbin_load(fid);
  C_I(x)=rbin_load(fid);
  return x;
}

/** @} */
/////////////////////////////////////////////////////////////

/** @name cmulti型の値の設定に関する関数 */
/** @{ */

/**
 @brief cmulti型の浮動小数点数を文字列から設定.
*/
void cset_ss(cmulti *x, const char *str_real, const char *str_imag)
{
  rset_s(C_R(x),str_real);
  rset_s(C_I(x),str_imag);
}

/**
 @brief cmulti型の浮動小数点数を文字列から設定.
*/
void cset_s(cmulti *x, const char *value)
{
  strings *list=NULL;
  int ret1=0,ret2=0;
  NULL_EXC2(x,value);
  list=strings_split_number(value);
  if(list==NULL){
    rset_nan(C_R(x));
    rset_nan(C_I(x));
  }else{
    cset_zero(x);
    if(strings_size(list)>=1 && strings_at(list,0)!=NULL){ ret1=mpfr_set_str(C_R(x),strings_at(list,0),10,get_round_mode()); }
    if(strings_size(list)>=2 && strings_at(list,1)!=NULL){ ret2=mpfr_set_str(C_I(x),strings_at(list,1),10,get_round_mode()); }
    if(ret1 || ret2){ cset_nan(x); }
  }
  list=strings_del(list);
}

/**
 @brief cmulti型の浮動小数点数を文字列から設定.
*/
void cset_script(cmulti *x, const char *str)
{
  strings *list=NULL;
  list=strings_split(str," \t\n",NULL,NULL,NULL);
  if(list!=NULL && list->n==1){
    cset_s(x,list->str[0]);
  }else if(list!=NULL && list->n==2){
    cset_ss(x,list->str[0],list->str[1]);
  }else{
    printf("Error in void cset_script(cmulti *x, const char *str)\n");
    printf("str=%s\n",str);
    printf("list="); strings_print(list); printf("\n");
    exit(0);
  }
  list=strings_del(list);
}

/**
 @brief cmulti型の浮動小数点数を倍精度浮動小数点数から設定.
*/
int cset_z(cmulti *x, dcomplex value)
{
  int e=0;
  e+=rset_d(C_R(x),Z_R(value));
  e+=rset_d(C_I(x),Z_I(value));
  return e;
}

/**
 @brief cmulti型の浮動小数点数を倍精度浮動小数点数から設定.
*/
int cset_dd(cmulti *x, double real, double imag)
{
  int e=0;
  e+=rset_d(C_R(x),real);
  e+=rset_d(C_I(x),imag);
  return e;
}

/**
 @brief cmulti型の浮動小数点数を倍精度浮動小数点数から設定.
*/
int cset_d(cmulti *x, double real)
{
  int e=0;
  e+=rset_d(C_R(x),real);
  e+=rset_d(C_I(x),0);
  return e;
}

/**
 @brief cmulti型の浮動小数点数を符号なし整数から設定.
*/
int cset_ui(cmulti *x, ulong real)
{
  int e=0;
  e+=rset_ui(C_R(x),real);
  e+=rset_ui(C_I(x),0);
  return e;
}

/**
 @brief cmulti型の浮動小数点数を符号あり整数から設定.
*/
int cset_si(cmulti *x, long real)
{
  int e=0;
  e+=rset_si(C_R(x),real);
  e+=rset_si(C_I(x),0);
  return e;
}

/**
 @brief cmulti型の値をInfに設定.
*/
void cset_inf(cmulti *x, int sgn_r, int sgn_i)
{
  rset_inf(C_R(x),sgn_r);
  rset_inf(C_I(x),sgn_i);
}

/**
 @brief cmulti型の値をNaNに設定.
*/
void cset_nan(cmulti *x)
{
  rset_nan(C_R(x));
  rset_nan(C_I(x));
}

/**
 @brief cmulti型の値を零に設定.
*/
int cset_zero(cmulti *x)
{
  int e=0;
  e+=rset_zero(C_R(x));
  e+=rset_zero(C_I(x));
  return e;
}

/**
 @brief cmulti型の値を1に設定.
*/
int cset_one(cmulti *x)
{
  int e=0;
  e+=rset_one(C_R(x));
  e+=rset_zero(C_I(x));
  return e;
}

/**
 @brief cmulti型の値を-1に設定.
*/
int cset_one_neg(cmulti *x)
{
  int e=0;
  e+=rset_one_neg(C_R(x));
  e+=rset_zero(C_I(x));
  return e;
}

/**
 @brief cmulti型の値を区間(0,1)の疑似乱数値を設定.
*/
void cset_rand(cmulti *x)
{
  rset_rand(C_R(x));
  rset_rand(C_I(x));
}


/** @} */

/////////////////////////////////////////////////////////

/** @name cmulti型の型変換に関する関数 */
/** @{ */

/**
 @brief cmulti型を倍精度複素数型に変換.
 */
dcomplex cget_z(cmulti *x)
{
  dcomplex z;
  Z_R(z)=rget_d(C_R(x));
  Z_I(z)=rget_d(C_I(x));
  return z;
}

/**
 @brief cmulti型を倍精度型に変換.
 */
double cget_d(cmulti *x)
{
  return rget_d(C_R(x));
}

/**
 @brief cmulti型を符号なし整数型に変換.
 */
ulong cget_ui(cmulti *x)
{
  return rget_ui(C_R(x));
}

/**
 @brief cmulti型を符号あり整数型に変換.
 */
long cget_si(cmulti *x)
{
  return rget_si(C_R(x));
}


/** @} */

/////////////////////////////////////////////////////////

/** @name cmulti型の自動精度調整モードが機能する関数 */
/** @{ */

/**
 @brief cmulti型の値のコピー y=x.
*/
int ccopy(cmulti *y, cmulti *x)
{
  int e=0;
  if(y==x){ return 0; }
  e+=rcopy(C_R(y),C_R(x));
  e+=rcopy(C_I(y),C_I(x));
  return e;
}

/**
 @brief cmulti型の値のコピー y=x.
*/
int ccopy_rr(cmulti *y, rmulti *x_r, rmulti *x_i)
{
  int e=0;
  e+=rcopy(C_R(y),x_r);
  e+=rcopy(C_I(y),x_i);
  return e;
}

/**
 @brief cmulti型の値のコピー y=x.
*/
int ccopy_r(cmulti *y, rmulti *x)
{
  int e=0;
  e+=rcopy(C_R(y),x);
  e+=rset_zero(C_I(y));
  return e;
}

/**
 @brief cmulti型の複素共役 y=conj(x).
*/
int cconj(cmulti *y, cmulti *x)
{
  int e=0;
  e+=rcopy(C_R(y),C_R(x));
  e+=rneg(C_I(y),C_I(x));
  return e;
}

/**
 @brief cmulti型の複素共役をとり複製 y=conj(x).
*/
int cconj_clone(cmulti *y, cmulti *x)
{
  int e=0;
  e+=rclone(C_R(y),C_R(x));
  e+=rclone(C_I(y),C_I(x));
  e+=rneg(C_I(y),C_I(y));
  return e;
}

/**
 @brief cmulti型の指数部の足し算 y=x*2^n
*/
int cmul_2exp(cmulti *y, cmulti *x, int nr, int ni)
{
  int e=0;
  e+=rmul_2exp(C_R(y),C_R(x),nr);
  e+=rmul_2exp(C_I(y),C_I(x),ni);
  return e;
}

/**
 @brief cmulti型の指数部の引き算 y=x/2^n
*/
int cdiv_2exp(cmulti *y, cmulti *x, int nr, int ni)
{
  int e=0;
  e+=rdiv_2exp(C_R(y),C_R(x),nr);
  e+=rdiv_2exp(C_I(y),C_I(x),ni);
  return e;
}

/**
 @brief cmulti型の符号反転 y=-x
*/
int cneg(cmulti *y, cmulti *x)
{
  int e=0;
  e+=rneg(C_R(y),C_R(x));
  e+=rneg(C_I(y),C_I(x));
  return e;
}

/**
 @brief cmulti型の絶対値の平方 y=abs(x)^2
*/
int cabs2(rmulti *y, cmulti *x)
{
  int e=0;
  e+=rmul(y,C_R(x),C_R(x));     //  y=x.r*x.r
  e+=radd_mul(y,C_I(x),C_I(x)); // y+=x.i*x.i
  return e;
}

//追加

/**
 @brief cmulti型の絶対値の平方 y=abs(x)^2
*/
int cabs2_ws(rmulti *y, cmulti *x, int *rwss, rmulti **rws)
{
  int e=0;
  e+=rmul(y,C_R(x),C_R(x));     //  y=x.r*x.r
  e+=radd_mul_ws(y,C_I(x),C_I(x),rwss,rws); // y+=x.i*x.i
  return e;
}
//ここまで

/**
 @brief cmulti型の実部と虚部の絶対値 y=abs(x.r)+i*abs(x.i)
*/
int cabsc(cmulti *y, cmulti *x)
{
  int e=0;
  e+=rabs(C_R(y),C_R(x));
  e+=rabs(C_I(y),C_I(x));
  return e;
}

/**
 @brief cmulti型の実部と虚部の絶対値 y=max(abs(x.r),abs(x.i))
*/
int cmax_absc(rmulti *y, cmulti *x)
{
  int e=0;
  rmulti *a=NULL;
  RAc(a,x);
  e+=rabs(y,C_R(x));
  e+=rabs(a,C_I(x));
  if(rgt(a,y)){ e+=rcopy(y,a); }
  RF(a);
  return e;
}

/**
 @brief cmulti型の足し算 z=x+y
*/
int cadd(cmulti *z, cmulti *x, cmulti *y)
{
  int e=0;
  e+=radd(C_R(z),C_R(x),C_R(y)); // z.r=x.r+y.r
  e+=radd(C_I(z),C_I(x),C_I(y)); // z.i=x.i+y.i
  return e;
}

/**
 @brief cmulti型の足し算 z=x+y
*/
int cadd_z(cmulti *z, cmulti *x, dcomplex y)
{
  int e=0;
  e+=radd_d(C_R(z),C_R(x),Z_R(y)); // z.r=x.r+y.r
  e+=radd_d(C_I(z),C_I(x),Z_I(y)); // z.i=x.i+y.i
  return e;
}

/**
 @brief cmulti型の足し算 z=x+y
*/
int cadd_r(cmulti *z, cmulti *x, rmulti *y)
{
  int e=0;
  e+=radd(C_R(z),C_R(x),y); // z.r=x.r+y
  e+=rcopy(C_I(z),C_I(x));  // z.i=x.i
  return e;
}

/**
 @brief cmulti型の足し算 z=x+y
*/
int cadd_d(cmulti *z, cmulti *x, double y)
{
  int e=0;
  e+=radd_d(C_R(z),C_R(x),y); // z.r=x.r+y
  e+=rcopy(C_I(z),C_I(x));    // z.i=x.i
  return e;
}

/**
 @brief cmulti型の足し算 z=x+y
*/
int cadd_ui(cmulti *z, cmulti *x, ulong y)
{
  int e=0;
  e+=radd_ui(C_R(z),C_R(x),y); // z.r=x.r+y
  e+=rcopy(C_I(z),C_I(x));     // z.i=x.i
  return e;
}

/**
 @brief cmulti型の足し算 z=x+y
*/
int cadd_si(cmulti *z, cmulti *x, long y)
{
  int e=0;
  e+=radd_si(C_R(z),C_R(x),y); // z.r=x.r+y
  e+=rcopy(C_I(z),C_I(x));     // z.i=x.i
  return e;
}

/**
 @brief cmulti型の引き算 z=x-y
*/
int csub(cmulti *z, cmulti *x, cmulti *y)
{
  int e=0;
  e+=rsub(C_R(z),C_R(x),C_R(y)); // z.r=x.r-y.r
  e+=rsub(C_I(z),C_I(x),C_I(y)); // z.i=x.i-y.i
  return e;
}

/**
 @brief cmulti型の引き算 z=x-y
*/
int csub_z1(cmulti *z, dcomplex x, cmulti *y)
{
  int e=0;
  e+=rsub_d1(C_R(z),Z_R(x),C_R(y)); // z.r=x.r-y.r
  e+=rsub_d1(C_I(z),Z_I(x),C_I(y)); // z.i=x.i-y.i
  return e;
}

/**
 @brief cmulti型の引き算 z=x-y
*/
int csub_z2(cmulti *z, cmulti *x, dcomplex y)
{
  int e=0;
  e+=rsub_d2(C_R(z),C_R(x),Z_R(y)); // z.r=x.r-y.r
  e+=rsub_d2(C_I(z),C_I(x),Z_I(y)); // z.i=x.i-y.i
  return e;
}

/**
 @brief cmulti型の引き算 z=x-y
*/
int csub_r1(cmulti *z, rmulti *x, cmulti *y)
{
  int e=0;
  e+=rsub(C_R(z),x,C_R(y)); // z.r=x-y.r
  e+=rneg(C_I(z),C_I(y));   // z.i= -y.i
  return e;
}

/**
 @brief cmulti型の引き算 z=x-y
*/
int csub_r2(cmulti *z, cmulti *x, rmulti *y)
{
  int e=0;
  e+=rsub(C_R(z),C_R(x),y); // z.r=x.r-y
  e+=rcopy(C_I(z),C_I(x));  // z.i=x.i
  return e;
}

/**
 @brief cmulti型の引き算 z=x-y
*/
int csub_d1(cmulti *z, double x, cmulti *y)
{
  int e=0;
  e+=rsub_d1(C_R(z),x,C_R(y)); // z.r=x-y.r
  e+=rneg(C_I(z),C_I(y));      // z.i= -y.i
  return e;
}

/**
 @brief cmulti型の引き算 z=x-y
*/
int csub_d2(cmulti *z, cmulti *x, double y)
{
  int e=0;
  e+=rsub_d2(C_R(z),C_R(x),y); // z.r=x.r-y
  e+=rcopy(C_I(z),C_I(x));     // z.i=x.i
  return e;
}

/**
 @brief cmulti型の引き算 z=x-y
*/
int csub_ui1(cmulti *z, ulong x, cmulti *y)
{
  int e=0;
  e+=rsub_ui1(C_R(z),x,C_R(y)); // z.r=x-y.r
  e+=rneg(C_I(z),C_I(y));      // z.i= -y.i
  return e;
}

/**
 @brief cmulti型の引き算 z=x-y
*/
int csub_ui2(cmulti *z, cmulti *x, ulong y)
{
  int e=0;
  e+=rsub_ui2(C_R(z),C_R(x),y); // z.r=x.r-y
  e+=rcopy(C_I(z),C_I(x));     // z.i=x.i
  return e;
}

/**
 @brief cmulti型の引き算 z=x-y
*/
int csub_si1(cmulti *z, long x, cmulti *y)
{
  int e=0;
  e+=rsub_si1(C_R(z),x,C_R(y)); // z.r=x-y.r
  e+=rneg(C_I(z),C_I(y));       // z.i= -y.i
  return e;
}

/**
 @brief cmulti型の引き算 z=x-y
*/
int csub_si2(cmulti *z, cmulti *x, long y)
{
  int e=0;
  e+=rsub_si2(C_R(z),C_R(x),y); // z.r=x.r-y
  e+=rcopy(C_I(z),C_I(x));     // z.i=x.i
  return e;
}


/**
 @brief cmulti型の掛け算 z=x*y
*/
int cmul(cmulti *z, cmulti *x, cmulti *y)
{
  int e=0;
  cmulti *a=NULL;
  CAc(a,z);
  e+=    rmul(C_R(a),C_R(x),C_R(y)); // z.r =x.r*y.r
  e+=rsub_mul(C_R(a),C_I(x),C_I(y)); // z.r-=x.i*y.i
  e+=    rmul(C_I(a),C_I(x),C_R(y)); // z.i =x.i*y.r
  e+=radd_mul(C_I(a),C_R(x),C_I(y)); // z.i+=x.r*y.i
  e+=ccopy(z,a);                     // z=a
  CF(a);
  return e;
}

//追加

/**
 @brief cmulti型の掛け算 z=x*y
*/
int cmul_ws(cmulti *z, cmulti *x, cmulti *y, int *rwss, rmulti **rws, int *cwss, cmulti **cws)
{
  if(*rwss<1){ ERROR_EXIT("Error `rwss=%d<1' in cmul_ws().\n",*rwss); }
  if(*cwss<1){ ERROR_EXIT("Error `cwss=%d<1' in cmul_ws().\n",*cwss); }
  int e=0;
  *cwss=*cwss-1;
  e+=    rmul(C_R(cws[0]),C_R(x),C_R(y)); // z.r =x.r*y.r
  e+=rsub_mul_ws(C_R(cws[0]),C_I(x),C_I(y),rwss,rws); // z.r-=x.i*y.i
  e+=    rmul(C_I(cws[0]),C_I(x),C_R(y)); // z.i =x.i*y.r
  e+=radd_mul_ws(C_I(cws[0]),C_R(x),C_I(y),rwss,rws); // z.i+=x.r*y.i
  e+=ccopy(z,cws[0]);                     // z=cws[0]
  *cwss=*cwss+1;
  return e;
}

//ここまで

/**
 @brief cmulti型の掛け算 z=x*y
*/
int cmul_z(cmulti *z, cmulti *x, dcomplex y)
{
  int e=0;
  cmulti *a=NULL;
  CAc(a,z);
  e+=cset_z(a,y); // a=y
  e+=cmul(z,x,a); // z=x*y
  CF(a);
  return e;
}

/**
 @brief cmulti型の掛け算 z=x*y
*/
int cmul_r(cmulti *z, cmulti *x, rmulti *y)
{
  int e=0;
  e+=rmul(C_R(z),C_R(x),y); // z.r=x.r*y
  e+=rmul(C_I(z),C_I(x),y); // z.i=x.i*y
  return e;
}

/**
 @brief cmulti型の掛け算 z=x*y
*/
int cmul_d(cmulti *z, cmulti *x, double y){
  int e=0;
  e+=rmul_d(C_R(z),C_R(x),y); // z.r=x.r*y
  e+=rmul_d(C_I(z),C_I(x),y); // z.i=x.i*y
  return e;
}

/**
 @brief cmulti型の掛け算 z=x*y
*/
int cmul_ui(cmulti *z, cmulti *x, ulong y){
  int e=0;
  e+=rmul_ui(C_R(z),C_R(x),y); // z.r=x.r*y
  e+=rmul_ui(C_I(z),C_I(x),y); // z.i=x.i*y
  return e;
}

/**
 @brief cmulti型の掛け算 z=x*y
*/
int cmul_si(cmulti *z, cmulti *x, long y){
  int e=0;
  e+=rmul_si(C_R(z),C_R(x),y); // z.r=x.r*y
  e+=rmul_si(C_I(z),C_I(x),y); // z.i=x.i*y
  return e;
}

///////////////////////////////////////////////////////////

/**
 @brief cmulti型の複素共役との掛け算 z=conj(x)*y
*/
int cdot(cmulti *z, cmulti *x, cmulti *y)
{
  int e=0;
  cmulti *a=NULL;
  CAc(a,z);
  e+=    rmul(C_R(a),C_R(x),C_R(y)); // z.r =x.r*y.r
  e+=radd_mul(C_R(a),C_I(x),C_I(y)); // z.r+=x.i*y.i
  e+=    rmul(C_I(a),C_R(x),C_I(y)); // z.i =x.r*y.i
  e+=rsub_mul(C_I(a),C_I(x),C_R(y)); // z.i-=x.i*y.r
  e+=ccopy(z,a);                     // z=a
  CF(a);
  return e;
}

/**
 @brief cmulti型の複素共役との掛け算 z=conj(x)*y
*/
int cdot_z1(cmulti *z, dcomplex x, cmulti *y)
{
  int e=0;
  cmulti *a=NULL;
  CAc(a,z);
  e+=cset_z(a,x); // a=x
  e+=cdot(z,a,y); // z=conj(x)*y  
  CF(a);
  return e;
}

/**
 @brief cmulti型の複素共役との掛け算 z=conj(x)*y
*/
int cdot_z2(cmulti *z, cmulti *x, dcomplex y)
{
  int e=0;
  cmulti *a=NULL;
  CAc(a,z);
  e+=cset_z(a,y); // a=y
  e+=cdot(z,x,a); // z=conj(x)*y 
  CF(a);
  return e;
}


/**
 @brief cmulti型の掛け算の加算 z+=x*y
*/
int cadd_mul(cmulti *z, cmulti *x, cmulti *y)
{
  int e=0;
  cmulti *a=NULL;
  CAc(a,z);
  e+=cmul(a,x,y);   // a=x*y
  e+=cadd(z,z,a);   // z=z+a
  CF(a);
  return e;
}

/**
 @brief cmulti型の掛け算の加算 z+=x*y
*/
int cadd_mul_r(cmulti *z, cmulti *x, rmulti *y)
{
  int e=0;
  cmulti *a=NULL;
  CAc(a,z);
  e+=cmul_r(a,x,y); // a=x*y
  e+=cadd(z,z,a);   // z=z+a
  CF(a);
  return e;
}

/**
 @brief cmulti型の掛け算の加算 z+=x*y
*/
int cadd_mul_z(cmulti *z, cmulti *x, dcomplex y)
{
  int e=0;
  cmulti *a=NULL;
  CAc(a,z);
  e+=cmul_z(a,x,y); // a=x*y
  e+=cadd(z,z,a);   // z=z+a
  CF(a);
  return e;
}

/**
 @brief cmulti型の掛け算の加算 z+=x*y
*/
int cadd_mul_d(cmulti *z, cmulti *x, double y)
{
  int e=0;
  cmulti *a=NULL;
  CAc(a,z);
  e+=cmul_d(a,x,y); // a=x*y
  e+=cadd(z,z,a);   // z=z+a
  CF(a);
  return e;
}

/**
 @brief cmulti型の掛け算の減算 z-=x*y
*/
int csub_mul(cmulti *z, cmulti *x, cmulti *y)
{
  int e=0;
  cmulti *a=NULL;
  CAc(a,z);
  e+=cmul(a,x,y);   // a=x*y
  e+=csub(z,z,a);   // z=z-a
  CF(a);
  return e;
}

//追加

/**
 @brief cmulti型の掛け算の減算 z-=x*y
*/
int csub_mul_ws(cmulti *z, cmulti *x, cmulti *y, int *rwss, rmulti **rws, int *cwss, cmulti **cws)
{
  if(*rwss<1){ ERROR_EXIT("Error `rwss=%d<1' in csub_mul_ws().\n",*rwss); }
  if(*cwss<2){ ERROR_EXIT("Error `cwss=%d<2' in csub_mul_ws().\n",*cwss); }
  int e=0;
  *cwss=*cwss-1;
  e+=cmul_ws(cws[0],x,y,rwss,rws,cwss,&cws[1]);   // cws[0]=x*y
  e+=csub(z,z,cws[0]);   // z=z-cws[0]
  *cwss=*cwss+1;
  return e;
}

//ここまで
/**
 @brief cmulti型の掛け算の減算 z-=x*y
*/
int csub_mul_r(cmulti *z, cmulti *x, rmulti *y)
{
  int e=0;
  cmulti *a=NULL;
  CAc(a,z);
  e+=cmul_r(a,x,y); // a=x*y
  e+=csub(z,z,a);   // z=z-a
  CF(a);
  return e;
}

/**
 @brief cmulti型の掛け算の減算 z-=x*y
*/
int csub_mul_z(cmulti *z, cmulti *x, dcomplex y)
{
  int e=0;
  cmulti *a=NULL;
  CAc(a,z);
  e+=cmul_z(a,x,y); // a=x*y
  e+=csub(z,z,a);   // z=z-a
  CF(a);
  return e;
}

/**
 @brief cmulti型の掛け算の減算 z-=x*y
*/
int csub_mul_d(cmulti *z, cmulti *x, double y)
{
  int e=0;
  cmulti *a=NULL;
  CAc(a,z);
  e+=cmul_d(a,x,y); // a=x*y
  e+=csub(z,z,a);   // z=z-a
  CF(a);
  return e;
}

/**
 @brief cmulti型の複素共役との掛け算の加算 z+=conj(x)*y
*/
int cadd_dot(cmulti *z, cmulti *x, cmulti *y)
{
  int e=0;
  cmulti *a=NULL;
  CAc(a,z);
  e+=cdot(a,x,y);   // a=conj(x)*y
  e+=cadd(z,z,a);   // z=z+a
  CF(a);
  return e;
}

/**
 @brief cmulti型の複素共役との掛け算の加算 z+=conj(x)*y
*/
int cadd_dot_z1(cmulti *z, dcomplex x, cmulti *y)
{
  int e=0;
  cmulti *a=NULL;
  CAc(a,z);
  e+=cdot_z1(a,x,y); // a=conj(x)*y
  e+=cadd(z,z,a);    // z=z+a
  CF(a);
  return e;
}

/**
 @brief cmulti型の複素共役との掛け算の加算 z+=conj(x)*y
*/
int cadd_dot_z2(cmulti *z, cmulti *x, dcomplex y)
{
  int e=0;
  cmulti *a=NULL;
  CAc(a,z);
  e+=cdot_z2(a,x,y); // a=conj(x)*y
  e+=cadd(z,z,a);    // z=z+a
  CF(a);
  return e;
}


/**
 @brief cmulti型の複素共役との掛け算の減算 z-=conj(x)*y
*/
int csub_dot(cmulti *z, cmulti *x, cmulti *y)
{
  int e=0;
  cmulti *a=NULL;
  CAc(a,z);
  e+=cdot(a,x,y);   // a=conj(x)*y
  e+=csub(z,z,a);   // z=z-a
  CF(a);
  return e;
}

/**
 @brief cmulti型の複素共役との掛け算の減算 z-=conj(x)*y
*/
int csub_dot_z1(cmulti *z, dcomplex x, cmulti *y)
{
  int e=0;
  cmulti *a=NULL;
  CAc(a,z);
  e+=cdot_z1(a,x,y);   // a=conj(x)*y
  e+=csub(z,z,a);   // z=z-a
  CF(a);
  return e;
}

/**
 @brief cmulti型の複素共役との掛け算の減算 z-=conj(x)*y
*/
int csub_dot_z2(cmulti *z, cmulti *x, dcomplex y)
{
  int e=0;
  cmulti *a=NULL;
  CAc(a,z);
  e+=cdot_z2(a,x,y);   // a=conj(x)*y
  e+=csub(z,z,a);   // z=z-a
  CF(a);
  return e;
}


/**
 @brief cmulti型の絶対値の平方の加算 y+=abs(x)^2
*/
int cadd_abs2(rmulti *y, cmulti *x)
{
  int e=0;
  rmulti *a=NULL;
  RAr(a,y);
  e+=cabs2(a,x);  // a=abs(x)^2
  e+=radd(y,y,a); // y=y+a
  RF(a);
  return e;
}

/**
 @brief rmulti型のべき乗 y=x^n
*/
int cpow_ui(cmulti *y, cmulti *x, ulong n)
{
  cmulti *a=NULL;
  int e=0;
  ulong i;
  if(cis_real(x)){
    e+=rpow_ui(C_R(y),C_R(x),n);
    e+=rset_si(C_I(y),0);
  }else if(n==0){
    e+=cset_one(y);
  }else{
    CAc(a,y);
    e+=ccopy(a,x);                        // a=x
    for(i=1; i<n; i++){ e+=cmul(a,a,x); } // a*=x
    e+=ccopy(y,a);                        // y=a
    CF(a);
  }
  return e;
}

/** @} */

/////////////////////////////////////////////////////////

/** @name cmulti型の数学関数に関する関数 */
/** @{ */

/**
 @brief cmulti型の逆数 z=1/x
*/
int cinv(cmulti *z, cmulti *x)
{
  int e=0;
  rmulti *den=NULL;
  cmulti *a=NULL;
  RAc(den,z); CAc(a,z);
  e+=cabs2(den,x);     // den=|x|^2
  e+=cconj(a,x);       // z=conj(x)
  e+=cdiv_r2(a,a,den); // z/=den
  e+=ccopy(z,a);       // z=a
  RF(den); CF(a);
  return e;
}

//追加

/**
 @brief cmulti型の逆数 z=1/x
*/
int cinv_ws(cmulti *z, cmulti *x, int *rwss, rmulti **rws, int *cwss, cmulti **cws)
{
  if(*rwss<2){ ERROR_EXIT("Error `rwss=%d<1' in cinv_ws().\n",*rwss); }
  if(*cwss<1){ ERROR_EXIT("Error `cwss=%d<1' in cinv_ws().\n",*cwss); }
  int e=0;
  *rwss=*rwss-1; *cwss=*cwss-1;
  e+=cabs2_ws(rws[0],x,rwss,&rws[1]);     // rws[0]=|x|^2
  e+=cconj(cws[0],x);       // z=conj(x)
  e+=cdiv_r2(cws[0],cws[0],rws[0]); // z/=rws[0]
  e+=ccopy(z,cws[0]);       // z=cws[0]
  *rwss=*rwss+1; *cwss=*cwss+1;
  return e;
}

//ここまで
/**
 @brief cmulti型の割り算 z=x/y
*/
int cdiv(cmulti *z, cmulti *x, cmulti *y)
{
  int e=0;
  rmulti *den=NULL;
  cmulti *a=NULL;
  RAc(den,z); CAc(a,z);
  e+=cabs2(den,y);     // den=|y|^2
  e+=cdot(a,y,x);      // z=conj(y)*x
  e+=cdiv_r2(a,a,den); // z/=den
  e+=ccopy(z,a);       // z=a
  RF(den); CF(a);
  return e;
}

/**
 @brief cmulti型の割り算 z=x/y
*/
int cdiv_z1(cmulti *z, dcomplex x, cmulti *y)
{
  int e=0;
  cmulti *a=NULL;
  CAc(a,z);
  e+=cset_z(a,x);   // a=x
  e+=cdiv(z,a,y);   // z=x/y
  CF(a);
  return e;
}

/**
 @brief cmulti型の割り算 z=x/y
*/
int cdiv_z2(cmulti *z, cmulti *x, dcomplex y)
{
  int e=0;
  cmulti *a=NULL;
  CAc(a,z);
  e+=cset_z(a,y);   // a=y
  e+=cdiv(z,x,a);   // z=x/y
  CF(a);
  return e;
}

/**
 @brief cmulti型の割り算 z=x/y
*/
int cdiv_r1(cmulti *z, rmulti *x, cmulti *y)
{
  int e=0;
  rmulti *den=NULL;
  cmulti *a=NULL;
  RAc(den,z); CAc(a,z);
  e+=cabs2(den,y);       // den=|y|^2
  e+=cconj(a,y);         // z=conj(y)
  e+=cmul_r(a,a,x);      // z*=x
  e+=cdiv_r2(a,a,den);   // z/=den
  e+=ccopy(z,a);         // z=a
  RF(den); CF(a);
  return e;
}

/**
 @brief cmulti型の割り算 z=x/y
*/
int cdiv_r2(cmulti *z, cmulti *x, rmulti *y)
{
  int e=0;
  e+=rdiv(C_R(z),C_R(x),y); // z.r=x.r/y
  e+=rdiv(C_I(z),C_I(x),y); // z.i=x.i/y
  return e;
}

/**
 @brief cmulti型の割り算 z=x/y
*/
int cdiv_d1(cmulti *z, double x, cmulti *y)
{
  int e=0;
  rmulti *a=NULL;
  RAc(a,z);
  e+=rset_d(a,x);    // a=x
  e+=cdiv_r1(z,a,y); // z=x/y
  RF(a);
  return e;
}

/**
 @brief cmulti型の割り算 z=x/y
*/
int cdiv_d2(cmulti *z, cmulti *x, double y)
{
  int e=0;
  e+=rdiv_d2(C_R(z),C_R(x),y); // z.r=x.r/y
  e+=rdiv_d2(C_I(z),C_I(x),y); // z.i=x.i/y
  return e;
}

/**
 @brief cmulti型の割り算 z=x/y
*/
int cdiv_ui1(cmulti *z, ulong x, cmulti *y)
{
  int e=0;
  rmulti *a=NULL;
  RAc(a,z);
  e+=rset_ui(a,x);   // a=x
  e+=cdiv_r1(z,a,y); // z=x/y
  RF(a);
  return e;
}

/**
 @brief cmulti型の割り算 z=x/y
*/
int cdiv_ui2(cmulti *z, cmulti *x, ulong y)
{
  int e=0;
  e+=rdiv_ui2(C_R(z),C_R(x),y); // z.r=x.r/y
  e+=rdiv_ui2(C_I(z),C_I(x),y); // z.i=x.i/y
  return e;
}

/**
 @brief cmulti型の割り算 z=x/y
*/
int cdiv_si1(cmulti *z, long x, cmulti *y)
{
  int e=0;
  rmulti *a=NULL;
  RAc(a,z);
  e+=rset_si(a,x);   // a=x
  e+=cdiv_r1(z,a,y); // z=x/y
  RF(a);
  return e;
}

/**
 @brief cmulti型の割り算 z=x/y
*/
int cdiv_si2(cmulti *z, cmulti *x, long y)
{
  int e=0;
  e+=rdiv_si2(C_R(z),C_R(x),y); // z.r=x.r/y
  e+=rdiv_si2(C_I(z),C_I(x),y); // z.i=x.i/y
  return e;
}

/**
 @brief cmulti型の絶対値 y=abs(x)
*/
int cabsv(rmulti *y, cmulti *x)
{
  int e=0;
  e+=cabs2(y,x); // y=abs(x)^2
  e+=rsqrt(y,y); // y=sqrt(y)
  return e;
}

//追加

/**
 @brief cmulti型の絶対値 y=abs(x)
*/
int cabsv_ws(rmulti *y, cmulti *x, int *rwss, rmulti **rws)
{
  int e=0;
  e+=cabs2_ws(y,x,rwss,rws); // y=abs(x)^2
  e+=rsqrt(y,y); // y=sqrt(y)
  return e;
}
//ここまで

/**
 @brief cmulti型の差の絶対値 z=abs(x-y)
*/
int cabs_sub(rmulti *z, cmulti *x, cmulti *y)
{
  int e=0;
  cmulti *a;
  CAr(a,z);
  e+=csub(a,x,y); // a=x-y
  e+=cabsv(z,a);  // z=abs(x-y)
  CF(a);
  return e;
}

/**
 @brief cmulti型の差の絶対値 z=abs(x-y)
*/
int cabs_sub_r(rmulti *z, cmulti *x, rmulti *y)
{
  int e=0;
  cmulti *a;
  CAr(a,z);
  e+=csub_r2(a,x,y); // a=x-y
  e+=cabsv(z,a);     // z=abs(x-y)
  CF(a);
  return e;
}

/**
 @brief cmulti型の絶対値の割り算 z=x/abs(y)
*/
int cdiv_abs(rmulti *z, rmulti *x, cmulti *y)
{
  int e=0;
  rmulti *a;
  RAr(a,z);
  e+=cabsv(a,y);  // z=abs(y)  
  e+=rdiv(z,x,a); // z=x/abs(y)
  RF(a);
  return e;
}

/**
 @brief cmulti型の絶対値の加算 y+=abs(x)
*/
int cadd_abs(rmulti *y, cmulti *x)
{
  int e=0;
  rmulti *a=NULL;
  RAr(a,y);
  e+=cabsv(a,x);  // a=abs(x)
  e+=radd(y,y,a); // y=y+a
  RF(a);
  return e;
}


/**
 @brief cmulti型の偏角 theta=arg(z)
*/
int cargument(rmulti *theta, cmulti *z)
{
  int e=0;
  e+=ratan2(theta,C_I(z),C_R(z));
  return e;
}


/**
 @brief cmulti型の規格化 y=x/abs(x)
*/
int cnormalize(cmulti *y, cmulti *x)
{
  int e=0;
  rmulti *a=NULL;
  RAc(a,y);
  e+=cabsv(a,x);     // a=abs(x)
  e+=cdiv_r2(y,x,a); // y=x/a
  RF(a);
  return e;
}

/**
 @brief cmulti型の極座標 z=r*exp(i*theta)
*/
int cget_polar(rmulti *r, rmulti *theta, cmulti *z)
{
  int e=0;
  e+=cabsv(r,z);
  e+=ratan2(theta,C_I(z),C_R(z));
  return e;
}

/**
 @brief cmulti型の極座標から実部，虚部へ変換 z=r*exp(i*theta)
*/
int cset_polar(cmulti *z, rmulti *r, rmulti *theta)
{
  int e=0;
  rmulti *s=NULL,*c=NULL;
  RAc(s,z); RAc(c,z);
  e+=rcos(c,theta); e+=rmul(C_R(z),r,c); // z.r=r*cos(theta)
  e+=rsin(s,theta); e+=rmul(C_I(z),r,s); // z.i=r*sin(theta)
  RF(s); RF(c);
  return e;  
}

/**
 @brief cmulti型のべき乗 y=x^n
*/
int cpow_si(cmulti *y, cmulti *x, long n)
{
  cmulti *a=NULL;
  int e=0;
  long i;
  if(n==0){
    e+=cset_one(y);
  }else if(cis_real(x)){
    e+=rpow_si(C_R(y),C_R(x),n);
    e+=rset_zero(C_I(y));
  }else if(n>0){
    CAc(a,y);
    e+=ccopy(a,x);                        // a=x
    for(i=1; i<n; i++){ e+=cmul(a,a,x); } // a*=x
    e+=ccopy(y,a);                        // y=a
    CF(a);
  }else if(n<0){
    CAc(a,y);
    e+=ccopy(a,x);                           // a=x
    for(i=1; i<(-n); i++){ e+=cmul(a,a,x); } // a*=x
    e+=cinv(a,a);                            // a=1/a
    e+=ccopy(y,a);                           // y=a
    CF(a);
  }else{ ERROR_AT; }
  return e;
}

/**
 @brief cmulti型のべき乗 z=x^y
*/
int cpow_r1(cmulti *z, rmulti *x, cmulti *y)
{
  int e=0;
  rmulti *a=NULL;
  cmulti *b=NULL;
  RAc(a,z); CAc(b,z);
  e+=rlog(a,x);     // a=log(x)
  e+=cmul_r(b,y,a); // b=y*log(x)
  e+=cexp_c(z,b);   // z=exp(y*log(x))
  RF(a); CF(b);
  return e;
}

/**
 @brief cmulti型のべき乗 z=x^y
*/
int cpow_r2(cmulti *z, cmulti *x, rmulti *y)
{
  int e=0;
  rmulti *r=NULL,*theta=NULL;
  RAc(r,z); RAc(theta,z);
  e+=cget_polar(r,theta,x); // x=r*exp(i*theta)
  e+=rpow(r,r,y);           // r=r^y
  e+=rmul(theta,theta,y);   // theta=theta*y
  e+=cset_polar(z,r,theta); // z=r^y*cos(theta*y)+i*r^y*sin(theta*y)
  RF(r); RF(theta);
  return e;
}

/**
 @brief cmulti型のべき乗 z=x^y
*/
int cpow_d2(cmulti *z, cmulti *x, double y)
{
  int e=0;
  rmulti *r=NULL,*theta=NULL;
  RAc(r,z); RAc(theta,z);
  e+=cget_polar(r,theta,x); // x=r*exp(i*theta)
  e+=rpow_d2(r,r,y);        // r=r^y
  e+=rmul_d(theta,theta,y); // theta=theta*y
  e+=cset_polar(z,r,theta); // z=r^y*cos(theta*y)+i*r^y*sin(theta*y)
  RF(r); RF(theta);
  return e;
}

/**
 @brief cmulti型のべき乗 z=x^y
*/
int cpow_c(cmulti *z, cmulti *x, cmulti *y)
{
  int e=0;
  cmulti *a=NULL;
  CAc(a,z);
  e+=clog_c(a,x); // a=log(x)
  e+=cmul(a,y,a); // a=y*log(x)
  e+=cexp_c(z,a); // z=exp(y*log(x))
  CF(a);
  return e;
}

/**
 @brief cmulti型の計算 y=sqrt(x)
*/
int csqrt_r(cmulti *y, rmulti *x)
{
  int e=0;
  if(ris_zero(x) || ris_positive(x)){
    e+=rsqrt(C_R(y),x);
    e+=rset_zero(C_I(y));
  }else if(ris_negative(x)){
    e+=rset_zero(C_R(y));
    e+=rneg(C_I(y),x);
    e+=rsqrt(C_I(y),C_I(y));
  }else{ ERROR_AT; }
  return e;
}

/**
 @brief cmulti型の計算 y=sqrt(x)
*/
int csqrt_c(cmulti *y, cmulti *x)
{
  return cpow_d2(y,x,0.5);
}

/**
 @brief cmulti型の計算 y=exp(x)
*/
int cexp_c(cmulti *y, cmulti *x)
{
  int e=0;
  rmulti *r=NULL,*theta=NULL;
  RAc(r,y); RAc(theta,y);
  e+=rexp(r,C_R(x));
  e+=rcopy(theta,C_I(x));
  e+=cset_polar(y,r,theta); // y=exp(x.r)*(cos(x.i)+i*sin(x.i))
  RF(r); RF(theta);
  return e;
}

/**
 @brief cmulti型の計算 y=log(x)
*/
int clog_r(cmulti *y, rmulti *x)
{
  int e=0;
  cmulti *a=NULL;
  CAc(a,y);
  e+=ccopy_r(a,x);
  e+=clog_c(y,a);
  CF(a);
  return e;
}

/**
 @brief cmulti型の計算 y=log(x)
*/
int clog_c(cmulti *y, cmulti *x)
{
  int e=0;
  rmulti *r=NULL,*theta=NULL;
  RAc(r,y); RAc(theta,y);
  e+=cget_polar(r,theta,x); // x=r*exp(i*theta)
  e+=rlog(C_R(y),r);        // y=log(r)+i*theta
  e+=rcopy(C_I(y),theta);   // y=log(r)+i*theta
  RF(r); RF(theta);
  return e;
}

/**
 @brief cmulti型の計算 y=sin(x)
*/
int csin_c(cmulti *y, cmulti *x)
{
  int e=0;
  rmulti *c=NULL,*ch=NULL,*s=NULL,*sh=NULL;
  RAc(c,y); RAc(ch,y); RAc(s,y); RAc(sh,y);
  // z.r=sin(x.r)*cosh(x.i)
  // z.i=cos(x.r)*sinh(x.i)
  e+=rsin(s,C_R(x));
  e+=rcos(c,C_R(x));
  e+=rsinh(sh,C_I(x));
  e+=rcosh(ch,C_I(x));
  e+=rmul(C_R(y),s,ch);
  e+=rmul(C_I(y),c,sh);
  RF(c); RF(ch); RF(s); RF(sh);
  return e;
}

/**
 @brief cmulti型の計算 y=cos(x)
*/
int ccos_c(cmulti *y, cmulti *x)
{
  int e=0;
  rmulti *c=NULL,*ch=NULL,*s=NULL,*sh=NULL;
  RAc(c,y); RAc(ch,y); RAc(s,y); RAc(sh,y);
  // z.r=cos(x.r)*cosh(x.i)
  // z.i=-sin(x.r)*sinh(x.i)
  e+=rsin(s,C_R(x));
  e+=rcos(c,C_R(x));
  e+=rsinh(sh,C_I(x));
  e+=rcosh(ch,C_I(x));
  e+=rmul(C_R(y),c,ch);
  e+=rmul(C_I(y),s,sh);
  e+=rneg(C_I(y),C_I(y));
  RF(c); RF(ch); RF(s); RF(sh);
  return e;
}

/**
 @brief cmulti型の計算 y=tan(x)
*/
int ctan_c(cmulti *y, cmulti *x)
{
  int e=0;
  cmulti *c=NULL,*s=NULL;
  CAc(c,y); CAc(s,y);
  e+=csin_c(s,x);
  e+=ccos_c(c,x);
  e+=cdiv(y,s,c); // y=sin(x)/cos(x)
  CF(c); CF(s);
  return e;
}

/**
 @brief cmulti型の計算 y=arcsin(x)
*/
int casin_r(cmulti *y, rmulti *x)
{
  int e=0;
  cmulti *a=NULL;
  CAc(a,y);
  e+=ccopy_r(a,x);
  e+=casin_c(y,a);
  CF(a);
  return e;
}

/**
 @brief cmulti型の計算 y=arcsin(x)
*/
int casin_c(cmulti *y, cmulti *x)
{
  cmulti *a=NULL,*b=NULL;
  int e=0;
  if(cis_real(x) && rle_si2(C_R(x),1) && rge_si2(C_R(x),-1)){
    e+=rasin(C_R(y),C_R(x));
    e+=rset_zero(C_I(y));
  }else{
    CAc(a,y); CAc(b,y);
    e+=cmul(a,x,x);     // a=x^2
    e+=csub_si1(a,1,a); // a=1-x^2
    e+=csqrt_c(a,a);    // a=sqrt(1-x^2)
    e+=cset_dd(b,0,1);  // b=I
    e+=cmul(b,b,x);     // b=I*x
    e+=cadd(a,a,b);     // a=I*x+sqrt(1-x^2)
    e+=clog_c(a,a);     // a=log(I*x+sqrt(1-x^2))
    e+=cset_dd(b,0,-1); // b=-I
    e+=cmul(a,a,b);     // a=-I*log(I*x+sqrt(1-x^2))
    e+=ccopy(y,a);      // y=-I*log(I*x+sqrt(1-x^2))
    CF(a); CF(b);
  }
  return e;
}

/**
 @brief cmulti型の計算 y=arccos(x)
*/
int cacos_r(cmulti *y, rmulti *x)
{
  int e=0;
  cmulti *a=NULL;
  CAc(a,y);
  e+=ccopy_r(a,x);
  e+=cacos_c(y,a);
  CF(a);
  return e;
}

/**
 @brief cmulti型の計算 y=arccos(x)
*/
int cacos_c(cmulti *y, cmulti *x)
{
  cmulti *a=NULL,*b=NULL;
  int e=0;
  if(cis_real(x) && rle_si2(C_R(x),1) && rge_si2(C_R(x),-1)){
    e+=racos(C_R(y),C_R(x));
    e+=rset_zero(C_I(y));
  }else{
    CAc(a,y); CAc(b,y);
    e+=cmul(a,x,x);     // a=x^2
    e+=csub_si2(a,a,1); // a=x^2-1
    e+=csqrt_c(a,a);    // a=sqrt(x^2-1)
    e+=cadd(a,a,x);     // a=x+sqrt(x^2-1)
    e+=clog_c(a,a);     // a=log(x+sqrt(x^2-1))
    e+=cset_dd(b,0,-1); // b=-I
    e+=cmul(a,a,b);     // a=-I*log(I*x+sqrt(1-x^2))
    e+=ccopy(y,a);      // y=-I*log(I*x+sqrt(1-x^2))
    CF(a); CF(b);
  }
  return e;
}

/**
 @brief cmulti型の計算 y=arctan(x)
*/
int catan_r(cmulti *y, rmulti *x)
{
  int e=0;
  cmulti *a=NULL;
  CAc(a,y);
  e+=ccopy_r(a,x);
  e+=catan_c(y,a);
  CF(a);
  return e;
}

/**
 @brief cmulti型の計算 y=arctan(x)
*/
int catan_c(cmulti *y, cmulti *x)
{
  cmulti *a=NULL,*b=NULL;
  int e=0;
  if(cis_real(x)){
    e+=ratan(C_R(y),C_R(x));
    e+=rset_si(C_I(y),0);
  }else{
    CAc(a,y); CAc(b,y);
    e+=cadd_si(a,x,1);    // a=1+x
    e+=csub_si1(b,1,x);   // b=1-x
    e+=cdiv(a,a,b);       // a=(1+x)/(1-x)
    e+=clog_c(a,a);       // a=log((1+x)/(1-x))
    e+=cset_dd(b,0,-0.5); // b=-I/2
    e+=cmul(a,a,b);       // a=(-I/2)*log((1+x)/(1-x))
    e+=ccopy(y,a);        // y=(-I/2)*log((1+x)/(1-x))
    CF(a); CF(b);
  }
  return e;
}

/**
 @brief cmulti型の計算 y=sinh(x)
*/
int csinh_c(cmulti *y, cmulti *x)
{
  int e=0;
  rmulti *c=NULL,*ch=NULL,*s=NULL,*sh=NULL;
  RAc(c,y); RAc(ch,y); RAc(s,y); RAc(sh,y);
  // z.r=sinh(x.r)*cos(x.i)
  // z.i=cosh(x.r)*sin(x.i)
  e+=rsinh(sh,C_R(x));
  e+=rcosh(ch,C_R(x));
  e+=rsin(s,C_I(x));
  e+=rcos(c,C_I(x));
  e+=rmul(C_R(y),sh,c);
  e+=rmul(C_I(y),ch,s);
  RF(c); RF(ch); RF(s); RF(sh);
  return e;
}

/**
 @brief cmulti型の計算 y=cosh(x)
*/
int ccosh_c(cmulti *y, cmulti *x)
{
  int e=0;
  rmulti *c=NULL,*ch=NULL,*s=NULL,*sh=NULL;
  RAc(c,y); RAc(ch,y); RAc(s,y); RAc(sh,y);
  // z.r=cosh(x.r)*cos(x.i)
  // z.i=sinh(x.r)*sin(x.i)
  e+=rsinh(sh,C_R(x));
  e+=rcosh(ch,C_R(x));
  e+=rsin(s,C_I(x));
  e+=rcos(c,C_I(x));
  e+=rmul(C_R(y),ch,c);
  e+=rmul(C_I(y),sh,s);
  RF(c); RF(ch); RF(s); RF(sh);
  return e;
}

/**
 @brief cmulti型の計算 y=tanh(x)
*/
int ctanh_c(cmulti *y, cmulti *x)
{
  int e=0;
  cmulti *ch=NULL,*sh=NULL;
  CAc(ch,y); CAc(sh,y);
  e+=csinh_c(sh,x);
  e+=ccosh_c(ch,x);
  e+=cdiv(y,sh,ch); // y=sinh(x)/cosh(x)
  CF(ch); CF(sh);
  return e;
}

/**
 @brief cmulti型の計算 y=asinh(x)
*/
int casinh_r(cmulti *y, rmulti *x)
{
  int e=0;
  cmulti *a=NULL;
  CAc(a,y);
  e+=ccopy_r(a,x);
  e+=casinh_c(y,a);
  CF(a);
  return e;
}

/**
 @brief cmulti型の計算 y=asinh(x)
*/
int casinh_c(cmulti *y, cmulti *x)
{
  int e=0;
  cmulti *a=NULL;  
  if(cis_real(x)){
    e+=rasinh(C_R(y),C_R(x));
    e+=rset_zero(C_I(y));
  }else{
    CAc(a,y);
    e+=cmul(a,x,x);
    e+=cadd_si(a,a,1);
    e+=csqrt_c(a,a);
    e+=cadd(a,a,x);
    e+=clog_c(a,a);
    e+=ccopy(y,a);
    CF(a);
  }
  return e;
}

/**
 @brief cmulti型の計算 y=acosh(x)
*/
int cacosh_r(cmulti *y, rmulti *x)
{
  int e=0;
  cmulti *a=NULL;
  CAc(a,y);
  e+=ccopy_r(a,x);
  e+=cacosh_c(y,a);
  CF(a);
  return e;
}

/**
 @brief cmulti型の計算 y=acosh(x)
*/
int cacosh_c(cmulti *y, cmulti *x)
{
  int e=0;
  cmulti *a=NULL;  
  if(cis_real(x) && rge_si2(C_R(x),1)){
    e+=racosh(C_R(y),C_R(x));
    e+=rset_zero(C_I(y));
  }else{
    CAc(a,y);
    e+=cmul(a,x,x);
    e+=csub_si2(a,a,1);
    e+=csqrt_c(a,a);
    e+=cadd(a,a,x);
    e+=clog_c(a,a);
    e+=ccopy(y,a);
    CF(a);
  }
  return e;
}

/**
 @brief cmulti型の計算 y=arctanh(x)
*/
int catanh_r(cmulti *y, rmulti *x)
{
  int e=0;
  cmulti *a=NULL;
  CAc(a,y);
  e+=ccopy_r(a,x);
  e+=catanh_c(y,a);
  CF(a);
  return e;
}

/**
 @brief cmulti型の計算 y=arctanh(x)
*/
int catanh_c(cmulti *y, cmulti *x)
{
  cmulti *a=NULL,*b=NULL;
  int e=0;
  if(cis_real(x) && rle_si2(C_R(x),1) && rge_si2(C_R(x),-1)){
    e+=ratanh(C_R(y),C_R(x));
    e+=rset_zero(C_I(y));
  }else{
    CAc(a,y); CAc(b,y);
    e+=cadd_si(a,x,1);  // a=1+x
    e+=csub_si1(b,1,x); // b=1-x
    e+=cdiv(a,a,b);     // a=(1+x)/(1-x)
    e+=clog_c(a,a);     // a=log((1+x)/(1-x))
    e+=cdiv_si2(a,a,2); // a=(1/2)*log((1+x)/(1-x))
    e+=ccopy(y,a);      // y=(1/2)*log((1+x)/(1-x))
    CF(a); CF(b);
  }
  return e;
}

/** @} */

/////////////////////////////////////////////////////////

/** @name cmulti型の値の比較に関する関数 */
/** @{ */

static int __cmulti_cmp_order=0;

/**
 @brief cmulti型の値の比較 x<=>y を実部，虚部で判定
*/
void ccmp_set_real_imag()
{
  __cmulti_cmp_order=0;
}

/**
 @brief cmulti型の値の比較 x<=>y を絶対値，偏角で判定
*/
void ccmp_set_abs_arg()
{
  __cmulti_cmp_order=1;
}

/**
 @brief cmulti型の値の比較 x<=>y の方法の取得
*/
int ccmp_get_type()
{
  return __cmulti_cmp_order;
}



/**
 @brief cmulti型の値の比較 x<=>y
*/
int ccmp(cmulti *x, cmulti *y)
{
  rmulti *xr=NULL,*yr=NULL,*xt=NULL,*yt=NULL;
  int value;
  if(ccmp_get_type()){
    RAc(xr,x); RAc(yr,y); RAc(xt,x); RAc(yt,y);
    cget_polar(xr,xt,x);
    cget_polar(yr,yt,y);
    value=rcmp(xr,yr);
    if(value==0){ value=rcmp(xt,yt); }
    RF(xr); RF(yr); RF(xt); RF(yt);
    return value;
  }else{
    value=rcmp(C_R(x),C_R(y));
    if(value!=0) return value;
    return rcmp(C_I(x),C_I(y));
  }
}

/**
 @brief cmulti型の値の比較 x<=>y
*/
int ccmp_z(cmulti *x, dcomplex y)
{
  int value;
  value=rcmp_d(C_R(x),Z_R(y));
  if(value!=0) return value;
  return rcmp_d(C_I(x),Z_I(y));
}

/**
 @brief cmulti型の値の比較 x<=>y
*/
int ccmp_r1(rmulti *x, cmulti *y)
{
  return -ccmp_r2(y,x);
}

/**
 @brief cmulti型の値の比較 x<=>y
*/
int ccmp_r2(cmulti *x, rmulti *y)
{
  int value;
  value=rcmp(C_R(x),y);
  if(value!=0) return value;
  return rcmp_d(C_I(x),0);
}

/**
 @brief cmulti型の値の比較 x<=>y
*/
int ccmp_d(cmulti *x, double y)
{
  int value;
  value=rcmp_d(C_R(x),y);
  if(value!=0) return value;
  return rcmp_d(C_I(x),0);
}

/**
 @brief cmulti型の値の比較 x<=>y
*/
int ccmp_ui(cmulti *x, ulong y)
{
  int value;
  value=rcmp_ui(C_R(x),y);
  if(value!=0) return value;
  return rcmp_ui(C_I(x),0);
}

/**
 @brief cmulti型の値の比較 x<=>y
*/
int ccmp_si(cmulti *x, long y)
{
  int value;
  value=rcmp_si(C_R(x),y);
  if(value!=0) return value;
  return rcmp_si(C_I(x),0);
}

/** @brief cmulti型の値の比較 x==y */
int ceq(cmulti *x, cmulti *y)    { return cis_number(x) && cis_number(y) && ccmp   (x,y)==0; }
/** @brief cmulti型の値の比較 x==y */
int ceq_z(cmulti *x, dcomplex y) { return cis_number(x) && ccmp_z (x,y)==0; }
/** @brief cmulti型の値の比較 x==y */
int ceq_r(cmulti *x, rmulti *y)  { return cis_number(x) && ccmp_r2(x,y)==0; }
/** @brief cmulti型の値の比較 x==y */
int ceq_d(cmulti *x, double y)   { return cis_number(x) && ccmp_d (x,y)==0; }
/** @brief cmulti型の値の比較 x==y */
int ceq_ui(cmulti *x, ulong y)   { return cis_number(x) && ccmp_ui(x,y)==0; }
/** @brief cmulti型の値の比較 x==y */
int ceq_si(cmulti *x, long y)    { return cis_number(x) && ccmp_si(x,y)==0; }
/** @brief cmulti型の値の比較 x!=y */
int cne(cmulti *x, cmulti *y)    { return !(cis_number(x) && cis_number(y) && ccmp   (x,y)==0); }
/** @brief cmulti型の値の比較 x!=y */
int cne_z(cmulti *x, dcomplex y) { return !(cis_number(x) && ccmp_z (x,y)==0); }
/** @brief cmulti型の値の比較 x!=y */
int cne_r(cmulti *x, rmulti *y)  { return !(cis_number(x) && ccmp_r2(x,y)==0); }
/** @brief cmulti型の値の比較 x!=y */
int cne_d(cmulti *x, double y)   { return !(cis_number(x) && ccmp_d (x,y)==0); }
/** @brief cmulti型の値の比較 x!=y */
int cne_ui(cmulti *x, ulong y)   { return !(cis_number(x) && ccmp_ui(x,y)==0); }
/** @brief cmulti型の値の比較 x!=y */
int cne_si(cmulti *x, long y)    { return !(cis_number(x) && ccmp_si(x,y)==0); }
/** @brief cmulti型の値の比較 x>y */
int cgt(cmulti *x, cmulti *y)     { return cis_number(x) && cis_number(y) && ccmp   (x,y)>0; }
/** @brief cmulti型の値の比較 x>y */
int cgt_r1(rmulti *x, cmulti *y)  { return ris_number(x) && cis_number(y) && ccmp_r1(x,y)>0; }
/** @brief cmulti型の値の比較 x>y */
int cgt_r2(cmulti *x, rmulti *y)  { return cis_number(x) && ris_number(y) && ccmp_r2(x,y)>0; }
/** @brief cmulti型の値の比較 x>y */
int cgt_z1(dcomplex x, cmulti *y) { return cis_number(y) && ccmp_z (y,x)<0; }
/** @brief cmulti型の値の比較 x>y */
int cgt_z2(cmulti *x, dcomplex y) { return cis_number(x) && ccmp_z (x,y)>0; }
/** @brief cmulti型の値の比較 x>y */
int cgt_d1(double x, cmulti *y)   { return cis_number(y) && ccmp_d (y,x)<0; }
/** @brief cmulti型の値の比較 x>y */
int cgt_d2(cmulti *x, double y)   { return cis_number(x) && ccmp_d (x,y)>0; }
/** @brief cmulti型の値の比較 x>y */
int cgt_ui1(ulong x, cmulti *y)   { return cis_number(y) && ccmp_ui(y,x)<0; }
/** @brief cmulti型の値の比較 x>y */
int cgt_ui2(cmulti *x, ulong y)   { return cis_number(x) && ccmp_ui(x,y)>0; }
/** @brief cmulti型の値の比較 x>y */
int cgt_si1(long x, cmulti *y)    { return cis_number(y) && ccmp_si(y,x)<0; }
/** @brief cmulti型の値の比較 x>y */
int cgt_si2(cmulti *x, long y)    { return cis_number(x) && ccmp_si(x,y)>0; }
/** @brief cmulti型の値の比較 x>=y */
int cge(cmulti *x, cmulti *y)    { return cis_number(x) && cis_number(y) && ccmp   (x,y)>=0; }
/** @brief cmulti型の値の比較 x>=y */
int cge_r1(rmulti *x, cmulti *y) { return ris_number(x) && cis_number(y) && ccmp_r1(x,y)>=0; }
/** @brief cmulti型の値の比較 x>=y */
int cge_r2(cmulti *x, rmulti *y) { return cis_number(x) && ris_number(y) && ccmp_r2(x,y)>=0; }
/** @brief cmulti型の値の比較 x>=y */
int cge_z1(dcomplex x, cmulti *y){ return cis_number(y) && ccmp_z (y,x)<=0; }
/** @brief cmulti型の値の比較 x>=y */
int cge_z2(cmulti *x, dcomplex y){ return cis_number(x) && ccmp_z (x,y)>=0; }
/** @brief cmulti型の値の比較 x>=y */
int cge_d1(double x, cmulti *y)  { return cis_number(y) && ccmp_d (y,x)<=0; }
/** @brief cmulti型の値の比較 x>=y */
int cge_d2(cmulti *x, double y)  { return cis_number(x) && ccmp_d (x,y)>=0; }
/** @brief cmulti型の値の比較 x>=y */
int cge_ui1(ulong x, cmulti *y)  { return cis_number(y) && ccmp_ui(y,x)<=0; }
/** @brief cmulti型の値の比較 x>=y */
int cge_ui2(cmulti *x, ulong y)  { return cis_number(x) && ccmp_ui(x,y)>=0; }
/** @brief cmulti型の値の比較 x>=y */
int cge_si1(long x, cmulti *y)   { return cis_number(y) && ccmp_si(y,x)<=0; }
/** @brief cmulti型の値の比較 x>=y */
int cge_si2(cmulti *x, long y)   { return cis_number(x) && ccmp_si(x,y)>=0; }
/** @brief cmulti型の値の比較 x<y */
int clt(cmulti *x, cmulti *y)    { return cis_number(x) && cis_number(y) && ccmp   (x,y)<0; }
/** @brief cmulti型の値の比較 x<y */
int clt_r1(rmulti *x, cmulti *y) { return ris_number(x) && cis_number(y) && ccmp_r1(x,y)<0; }
/** @brief cmulti型の値の比較 x<y */
int clt_r2(cmulti *x, rmulti *y) { return cis_number(x) && ris_number(y) && ccmp_r2(x,y)<0; }
/** @brief cmulti型の値の比較 x<y */
int clt_z1(dcomplex x, cmulti *y){ return cis_number(y) && ccmp_z (y,x)>0; }
/** @brief cmulti型の値の比較 x<y */
int clt_z2(cmulti *x, dcomplex y){ return cis_number(x) && ccmp_z (x,y)<0; }
/** @brief cmulti型の値の比較 x<y */
int clt_d1(double x, cmulti *y)  { return cis_number(y) && ccmp_d (y,x)>0; }
/** @brief cmulti型の値の比較 x<y */
int clt_d2(cmulti *x, double y)  { return cis_number(x) && ccmp_d (x,y)<0; }
/** @brief cmulti型の値の比較 x<y */
int clt_ui1(ulong x, cmulti *y)  { return cis_number(y) && ccmp_ui(y,x)>0; }
/** @brief cmulti型の値の比較 x<y */
int clt_ui2(cmulti *x, ulong y)  { return cis_number(x) && ccmp_ui(x,y)<0; }
/** @brief cmulti型の値の比較 x<y */
int clt_si1(long x, cmulti *y)   { return cis_number(y) && ccmp_si(y,x)>0; }
/** @brief cmulti型の値の比較 x<y */
int clt_si2(cmulti *x, long y)   { return cis_number(x) && ccmp_si(x,y)<0; }
/** @brief cmulti型の値の比較 x<=y */
int cle(cmulti *x, cmulti *y)    { return cis_number(x) && cis_number(y) && ccmp   (x,y)<=0; }
/** @brief cmulti型の値の比較 x<=y */
int cle_r1(rmulti *x, cmulti *y) { return ris_number(x) && cis_number(y) && ccmp_r1(x,y)<=0; }
/** @brief cmulti型の値の比較 x<=y */
int cle_r2(cmulti *x, rmulti *y) { return cis_number(x) && ris_number(y) && ccmp_r2(x,y)<=0; }
/** @brief cmulti型の値の比較 x<=y */
int cle_z1(dcomplex x, cmulti *y){ return cis_number(y) && ccmp_z (y,x)>=0; }
/** @brief cmulti型の値の比較 x<=y */
int cle_z2(cmulti *x, dcomplex y){ return cis_number(x) && ccmp_z (x,y)<=0; }
/** @brief cmulti型の値の比較 x<=y */
int cle_d1(double x, cmulti *y)  { return cis_number(y) && ccmp_d (y,x)>=0; }
/** @brief cmulti型の値の比較 x<=y */
int cle_d2(cmulti *x, double y)  { return cis_number(x) && ccmp_d (x,y)<=0; }
/** @brief cmulti型の値の比較 x<=y */
int cle_ui1(ulong x, cmulti *y)  { return cis_number(y) && ccmp_ui(y,x)>=0; }
/** @brief cmulti型の値の比較 x<=y */
int cle_ui2(cmulti *x, ulong y)  { return cis_number(x) && ccmp_ui(x,y)<=0; }
/** @brief cmulti型の値の比較 x<=y */
int cle_si1(long x, cmulti *y)   { return cis_number(y) && ccmp_si(y,x)>=0; }
/** @brief cmulti型の値の比較 x<=y */
int cle_si2(cmulti *x, long y)   { return cis_number(x) && ccmp_si(x,y)<=0; }

/** @brief cmulti型の値の比較 abs(x)<=>abs(y) */
int cabs_cmp(cmulti *x, cmulti *y)
{
  int value;  
  rmulti *ax=NULL,*ay=NULL;
  RAc(ax,x); RAc(ay,y);
  cabs2(ax,x); // ax=abs(x)
  cabs2(ay,y); // ay=abs(y)
  value=rcmp(ax,ay);
  RF(ax); RF(ay);  
  return value;
}

/** @brief cmulti型の値の比較 abs(x)==abs(y) */
int cabs_eq(cmulti *x, cmulti *y) { return cis_number(x) && cis_number(y) && cabs_cmp(x,y)==0; }
/** @brief cmulti型の値の比較 abs(x)>abs(y) */
int cabs_gt(cmulti *x, cmulti *y) { return cis_number(x) && cis_number(y) && cabs_cmp(x,y)>0;  }
/** @brief cmulti型の値の比較 abs(x)>=abs(y) */
int cabs_ge(cmulti *x, cmulti *y) { return cis_number(x) && cis_number(y) && cabs_cmp(x,y)>=0; }
/** @brief cmulti型の値の比較 abs(x)<abs(y) */
int cabs_lt(cmulti *x, cmulti *y) { return cis_number(x) && cis_number(y) && cabs_cmp(x,y)<0;  }
/** @brief cmulti型の値の比較 abs(x)<=abs(y) */
int cabs_le(cmulti *x, cmulti *y) { return cis_number(x) && cis_number(y) && cabs_cmp(x,y)<=0; }
/** @brief cmulti型の値の比較 x.r==y.r && x.i==y.i */
int ceqc(cmulti *x, cmulti *y){ return cis_number(x) && cis_number(y) && req(C_R(x),C_R(y)) && req(C_I(x),C_I(y)); }
/** @brief cmulti型の値の比較 x.r>y.r && x.i>y.i */
int cgtc(cmulti *x, cmulti *y){ return cis_number(x) && cis_number(y) && rgt(C_R(x),C_R(y)) && rgt(C_I(x),C_I(y)); }
/** @brief cmulti型の値の比較 x.r>=y.r && x.i>=y.i */
int cgec(cmulti *x, cmulti *y){ return cis_number(x) && cis_number(y) && rge(C_R(x),C_R(y)) && rge(C_I(x),C_I(y)); }
/** @brief cmulti型の値の比較 x.r<y.r && x.i<y.i */
int cltc(cmulti *x, cmulti *y){ return cis_number(x) && cis_number(y) && rlt(C_R(x),C_R(y)) && rlt(C_I(x),C_I(y)); }
/** @brief cmulti型の値の比較 x.r<=y.r && x.i<=y.i */
int clec(cmulti *x, cmulti *y){ return cis_number(x) && cis_number(y) && rle(C_R(x),C_R(y)) && rle(C_I(x),C_I(y)); }

/** @} */

//EOF
