#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<stdarg.h>
#include"is_macros.h"
#include"is_rmulti.h"
#include"is_cmulti.h"
#include"is_irmulti.h"
#include"is_icmulti.h"
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
#define RF(X) ((X)=rmfree(X))
#define CAr(X,P) ((X)=callocate_prec(rget_prec(P)))
#define CAc(X,P) ((X)=callocate_prec(cget_prec(P)))
#define CF(X) ((X)=cmfree(X))

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
cmulti *cmfree(cmulti *x)
{
  if(x==NULL) return NULL; 
  C_R(x)=rmfree(C_R(x));
  C_I(x)=rmfree(C_I(x));
  free(x);
  x=NULL;
  return x;
}

/**
 @brief 初期化済みのcmulti型の浮動小数点数の精度(ビット数)を変更し再初期化.
*/
void cround(cmulti *x, int prec)
{
  rround(C_R(x),prec);
  rround(C_I(x),prec);
}

/**
 @brief cmulti型の値を複成.
*/
void cclone_c(cmulti *y, cmulti *x)
{
  rclone_r(C_R(y),C_R(x));
  rclone_r(C_I(y),C_I(x));
}

/**
 @brief cmulti型の値をrmulti型から複成.
*/
void cclone_r(cmulti *y, rmulti *x)
{
  rclone_r(C_R(y),x);
  rset_zero(C_I(y));
}

/**
 @brief cmulti型の値をrmulti型から複成.
*/
void cclone_rr(cmulti *y, rmulti *x_r, rmulti *x_i)
{
  rclone_r(C_R(y),x_r);
  rclone_r(C_I(y),x_i);
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

/** @} */

/////////////////////////////////////////////////////////////

/** @name cmulti型のメンバ変数に関する関数 */
/** @{ */

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
 @brief cmulti型の値のコピー y=x.
 */
void cset_c(cmulti *y, cmulti *x)
{
  if(y==x){ return; }
  rset_r(C_R(y),C_R(x));
  rset_r(C_I(y),C_I(x));
}

/**
 @brief cmulti型の浮動小数点数を倍精度浮動小数点数から設定.
 */
void cset_dd(cmulti *x, double real, double imag)
{
  rset_d(C_R(x),real);
  rset_d(C_I(x),imag);
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
void cset_zero(cmulti *x)
{
  rset_zero(C_R(x));
  rset_zero(C_I(x));
}

/**
 @brief cmulti型の値を1に設定.
 */
void cset_one(cmulti *x)
{
  rset_one(C_R(x));
  rset_zero(C_I(x));
}

/**
 @brief cmulti型の値を-1に設定.
 */
void cset_one_neg(cmulti *x)
{
  rset_one_neg(C_R(x));
  rset_zero(C_I(x));
}

/**
 @brief cmulti型の値を区間(0,1)の疑似乱数値を設定.
 */
void cset_rand(cmulti *x)
{
  rset_rand(C_R(x));
  rset_rand(C_I(x));
}

/**
 @brief cmulti型の極座標 z=r*exp(i*theta)
 */
void cget_polar(rmulti *r, rmulti *theta, cmulti *z)
{
  rabs_c(r,z);
  ratan2_rr(theta,C_I(z),C_R(z));
}

/**
 @brief cmulti型の極座標から実部，虚部へ変換 z=r*exp(i*theta)
 */
void cset_polar(cmulti *z, rmulti *r, rmulti *theta)
{
  rmulti *s=NULL,*c=NULL;
  RAc(s,z); RAc(c,z);
  rcos_r(c,theta); rmul_rr(C_R(z),r,c); // z.r=r*cos(theta)
  rsin_r(s,theta); rmul_rr(C_I(z),r,s); // z.i=r*sin(theta)
  RF(s); RF(c);
}

/** @} */

/////////////////////////////////////////////////////////

/** @name cmulti型の型変換に関する関数 */
/** @{ */

/**
 @brief cmulti型の浮動小数点数を文字列から設定.
 */
void cset_s(cmulti *x, char *s)
{
  int prec;
  cmulti *a0=NULL,*a1=NULL;
  prec=cget_prec(x);
  a0=callocate_prec(prec);
  a1=callocate_prec(prec);
  icset_s(a0,a1,s);
  irmid(C_R(x),C_R(a0),C_R(a1));
  irmid(C_I(x),C_I(a0),C_I(a1));
}

/**
 @brief cmulti型の浮動小数点数を符号あり整数から設定.
 */
void cset_i(cmulti *x, int real)
{
  rset_i(C_R(x),real);
  rset_i(C_I(x),0);
}

/**
 @brief cmulti型の浮動小数点数を倍精度浮動小数点数から設定.
 */
void cset_d(cmulti *x, double real)
{
  rset_d(C_R(x),real);
  rset_d(C_I(x),0);
}

/**
 @brief y=x.
 */
void cset_rr(cmulti *y, rmulti *x_r, rmulti *x_i)
{
  rset_r(C_R(y),x_r);
  rset_r(C_I(y),x_i);
}

/**
 @brief y=x
 */
void cset_r(cmulti *y, rmulti *x)
{
  rset_r  (C_R(y),x);
  rset_d(C_I(y),0);
}

/**
 @brief cmulti型の浮動小数点数を倍精度浮動小数点数から設定.
 */
void cset_z(cmulti *x, dcomplex value)
{
  rset_d(C_R(x),Z_R(value));
  rset_d(C_I(x),Z_I(value));
}


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
 @brief cmulti型を符号あり整数型に変換.
 */
int cget_i(cmulti *x)
{
  return rget_i(C_R(x));
}


/** @} */

/////////////////////////////////////////////////////////

/** @name cmulti型の関数 y=f(x) */
/** @{ */

/**
 @brief y=-x
 */
void cneg_c(cmulti *y, cmulti *x)
{
  rneg_r(C_R(y),C_R(x));
  rneg_r(C_I(y),C_I(x));
}

/**
 @brief y=-x
 */
void cneg_r(cmulti *y, rmulti *x)
{
  rneg_r(C_R(y),x);
  rset_d(C_I(y),0);
}

/**
 @brief y=-x
 */
void cneg_z(cmulti *y, dcomplex x)
{
  rneg_d(C_R(y),Z_R(x));
  rneg_d(C_I(y),Z_I(x));
}

/**
 @brief y=-x
 */
void cneg_d(cmulti *y, double x)
{
  rneg_d(C_R(y),x);
  rneg_d(C_I(y),0);
}

/**
 @brief y=1/x
 */
void cinv_c_ws(cmulti *z, cmulti *x, int n_rws, rmulti **rws, int n_cws, cmulti **cws)
{
  if(n_rws<2){ ERROR_EXIT("Error `n_rws=%d<1' in cinv_ws().\n",n_rws); }
  if(n_cws<1){ ERROR_EXIT("Error `n_cws=%d<1' in cinv_ws().\n",n_cws); }
  rabs2_c_ws(rws[0],x,n_rws-1,rws+1); // rws[0]=|x|^2
  cconj_c(cws[0],x);                  // z=conj(x)
  cdiv_cr(cws[0],cws[0],rws[0]);     // z/=rws[0]
  cset_c(z,cws[0]);                // z=cws[0]
}

/**
 @brief y=1/x
 */
void cinv_c(cmulti *z, cmulti *x)
{
  rmulti *den=NULL;
  cmulti *a=NULL;
  RAc(den,z); CAc(a,z);
  rabs2_c(den,x);     // den=|x|^2
  cconj_c(a,x);       // z=conj(x)
  cdiv_cr(a,a,den);  // z/=den
  cset_c(z,a);     // z=a
  RF(den); CF(a);
}

/**
 @brief y=1/x
 */
void cinv_r(cmulti *z, rmulti *x)
{
  rinv_r(C_R(z),x);
  rset_d(C_I(z),0);
}

/**
 @brief y=1/x
 */
void cinv_z(cmulti *z, dcomplex x)
{
  cmulti *a=NULL;
  CAc(a,z);
  cset_z(a,x);
  cinv_c(z,a);
  CF(a);
}

/**
 @brief y=1/x
 */
void cinv_d(cmulti *z, double x)
{
  rset_d(C_R(z),1.0/x);
  rset_d(C_I(z),0);
}

/**
 @brief cmulti型の複素共役 y=conj(x).
*/
void cconj_c(cmulti *y, cmulti *x)
{
  rset_r(C_R(y),C_R(x));
  rneg_r(C_I(y),C_I(x));
}

/**
 @brief cmulti型の複素共役をとり複製 y=conj(x).
*/
void cconj_clone(cmulti *y, cmulti *x)
{
  rclone_r(C_R(y),C_R(x));
  rclone_r(C_I(y),C_I(x));
  rneg_r(C_I(y),C_I(y));
}

/**
 @brief cmulti型の絶対値 y=abs(x)
*/
void rabs_c(rmulti *y, cmulti *x)
{
  rabs2_c(y,x); // y=abs(x)^2
  rsqrt_r(y,y); // y=sqrt(y)
}

/**
 @brief cmulti型の絶対値 y=abs(x)
*/
void rabs_c_ws(rmulti *y, cmulti *x, int n_rws, rmulti **rws)
{
  rabs2_c_ws(y,x,n_rws,rws); // y=abs(x)^2
  rsqrt_r(y,y);              // y=sqrt(y)
}

/**
 @brief cmulti型の絶対値の加算 y+=abs(x)
*/
void radd_abs_c(rmulti *y, cmulti *x)
{
  rmulti *a=NULL;
  RAr(a,y);
  rabs_c(a,x);  // a=abs(x)
  radd_rr(y,y,a); // y=y+a
  RF(a);
}

/**
 @brief cmulti型の絶対値の平方 y=abs(x)^2
*/
void rabs2_c(rmulti *y, cmulti *x)
{
  rmul_rr(y,C_R(x),C_R(x));     // y=x.r*x.r
  radd_mul_rr(y,C_I(x),C_I(x)); // y+=x.i*x.i
}

/**
 @brief cmulti型の絶対値の平方 y=abs(x)^2
*/
void rabs2_c_ws(rmulti *y, cmulti *x, int n_rws, rmulti **rws)
{
  rmul_rr(y,C_R(x),C_R(x));                  // y=x.r*x.r
  radd_mul_rr_ws(y,C_I(x),C_I(x),n_rws,rws); // y+=x.i*x.i
}

/**
 @brief cmulti型の絶対値の平方の加算 y+=abs(x)^2
*/
void radd_abs2_c(rmulti *y, cmulti *x)
{
  rmulti *a=NULL;
  RAr(a,y);
  rabs2_c(a,x);  // a=abs(x)^2
  radd_rr(y,y,a); // y=y+a
  RF(a);
}

/**
 @brief cmulti型の実部と虚部の絶対値 y=abs(x.r)+i*abs(x.i)
*/
void cabsc_c(cmulti *y, cmulti *x)
{
  rabs_r(C_R(y),C_R(x));
  rabs_r(C_I(y),C_I(x));
}

/**
 @brief cmulti型の実部と虚部の絶対値 y=max(abs(x.r),abs(x.i))
*/
void rmax_absc_c(rmulti *y, cmulti *x)
{
  rmulti *a=NULL;
  RAc(a,x);
  rabs_r(y,C_R(x));
  rabs_r(a,C_I(x));
  if(gt_rr(a,y)){ rset_r(y,a); }
  RF(a);
}

/**
 @brief cmulti型の規格化 y=x/abs(x)
*/
void cnormalize_c(cmulti *y, cmulti *x)
{
  rmulti *a=NULL;
  RAc(a,y);
  rabs_c(a,x);     // a=abs(x)
  cdiv_cr(y,x,a); // y=x/a
  RF(a);
}

/**
 @brief cmulti型の偏角 theta=arg(z)
*/
void rarg_c(rmulti *theta, cmulti *z)
{
  ratan2_rr(theta,C_I(z),C_R(z));
}

/**
 @brief cmulti型の計算 y=sqrt(x)
*/
void csqrt_r(cmulti *y, rmulti *x)
{
  if(ris_zero(x) || ris_positive(x)){
    rsqrt_r(C_R(y),x);
    rset_zero(C_I(y));
  }else if(ris_negative(x)){
    rset_zero(C_R(y));
    rneg_r(C_I(y),x);
    rsqrt_r(C_I(y),C_I(y));
  }else{ ERROR_AT; }
}

/**
 @brief cmulti型の計算 y=sqrt(x)
*/
void csqrt_c(cmulti *y, cmulti *x)
{
  cpow_cd(y,x,0.5);
}

/**
 @brief cmulti型の計算 y=exp(x)
*/
void cexp_c(cmulti *y, cmulti *x)
{
  rmulti *r=NULL,*theta=NULL;
  RAc(r,y); RAc(theta,y);
  rexp_r(r,C_R(x));
  rset_r(theta,C_I(x));
  // y=exp(x.r)*(cos(x.i)+i*sin(x.i))  
  cset_polar(y,r,theta);
  RF(r); RF(theta);
}

/**
 @brief cmulti型の計算 y=log(x)
*/
void clog_r(cmulti *y, rmulti *x)
{
  cmulti *a=NULL;
  CAc(a,y);
  cset_r(a,x);
  clog_c(y,a);
  CF(a);
}

/**
 @brief cmulti型の計算 y=log(x)
*/
void clog_c(cmulti *y, cmulti *x)
{
  rmulti *r=NULL,*theta=NULL;
  RAc(r,y); RAc(theta,y);
  cget_polar(r,theta,x); // x=r*exp(i*theta)
  rlog_r(C_R(y),r);        // y=log(r)+i*theta
  rset_r(C_I(y),theta);   // y=log(r)+i*theta
  RF(r); RF(theta);
}

/**
 @brief cmulti型の計算 y=sin(x)
*/
void csin_c(cmulti *y, cmulti *x)
{
  rmulti *c=NULL,*ch=NULL,*s=NULL,*sh=NULL;
  RAc(c,y); RAc(ch,y); RAc(s,y); RAc(sh,y);
  // z.r=sin(x.r)*cosh(x.i)
  // z.i=cos(x.r)*sinh(x.i)
  rsin_r(s,C_R(x));
  rcos_r(c,C_R(x));
  rsinh_r(sh,C_I(x));
  rcosh_r(ch,C_I(x));
  rmul_rr(C_R(y),s,ch);
  rmul_rr(C_I(y),c,sh);
  RF(c); RF(ch); RF(s); RF(sh);
}

/**
 @brief cmulti型の計算 y=cos(x)
*/
void ccos_c(cmulti *y, cmulti *x)
{
  rmulti *c=NULL,*ch=NULL,*s=NULL,*sh=NULL;
  RAc(c,y); RAc(ch,y); RAc(s,y); RAc(sh,y);
  // z.r=cos(x.r)*cosh(x.i)
  // z.i=-sin(x.r)*sinh(x.i)
  rsin_r(s,C_R(x));
  rcos_r(c,C_R(x));
  rsinh_r(sh,C_I(x));
  rcosh_r(ch,C_I(x));
  rmul_rr(C_R(y),c,ch);
  rmul_rr(C_I(y),s,sh);
  rneg_r(C_I(y),C_I(y));
  RF(c); RF(ch); RF(s); RF(sh);
}

/**
 @brief cmulti型の計算 y=tan(x)
*/
void ctan_c(cmulti *y, cmulti *x)
{
  cmulti *c=NULL,*s=NULL;
  CAc(c,y); CAc(s,y);
  csin_c(s,x);
  ccos_c(c,x);
  cdiv_cc(y,s,c); // y=sin(x)/cos(x)
  CF(c); CF(s);
}

/**
 @brief cmulti型の計算 y=arcsin(x)
*/
void casin_r(cmulti *y, rmulti *x)
{
  cmulti *a=NULL;
  CAc(a,y);
  cset_r(a,x);
  casin_c(y,a);
  CF(a);
}

/**
 @brief cmulti型の計算 y=arcsin(x)
*/
void casin_c(cmulti *y, cmulti *x)
{
  cmulti *a=NULL,*b=NULL;
  if(cis_real(x) && le_rd(C_R(x),1) && ge_rd(C_R(x),-1)){
    rasin_r(C_R(y),C_R(x));
    rset_zero(C_I(y));
  }else{
    CAc(a,y); CAc(b,y);
    cmul_cc(a,x,x);   // a=x^2
    csub_dc(a,1,a);   // a=1-x^2
    csqrt_c(a,a);    // a=sqrt(1-x^2)
    cset_dd(b,0,1);   // b=I
    cmul_cc(b,b,x);   // b=I*x
    cadd_cc(a,a,b);   // a=I*x+sqrt(1-x^2)
    clog_c(a,a);     // a=log(I*x+sqrt(1-x^2))
    cset_dd(b,0,-1);  // b=-I
    cmul_cc(a,a,b);   // a=-I*log(I*x+sqrt(1-x^2))
    cset_c(y,a);    // y=-I*log(I*x+sqrt(1-x^2))
    CF(a); CF(b);
  }
}

/**
 @brief cmulti型の計算 y=arccos(x)
*/
void cacos_r(cmulti *y, rmulti *x)
{
  cmulti *a=NULL;
  CAc(a,y);
  cset_r(a,x);
  cacos_c(y,a);
  CF(a);
}

/**
 @brief cmulti型の計算 y=arccos(x)
*/
void cacos_c(cmulti *y, cmulti *x)
{
  cmulti *a=NULL,*b=NULL;
  if(cis_real(x) && le_rd(C_R(x),1) && ge_rd(C_R(x),-1)){
    racos_r(C_R(y),C_R(x));
    rset_zero(C_I(y));
  }else{
    CAc(a,y); CAc(b,y);
    cmul_cc(a,x,x);     // a=x^2
    csub_cd(a,a,1);     // a=x^2-1
    csqrt_c(a,a);      // a=sqrt(x^2-1)
    cadd_cc(a,a,x);     // a=x+sqrt(x^2-1)
    clog_c(a,a);       // a=log(x+sqrt(x^2-1))
    cset_dd(b,0,-1);    // b=-I
    cmul_cc(a,a,b);     // a=-I*log(I*x+sqrt(1-x^2))
    cset_c(y,a);      // y=-I*log(I*x+sqrt(1-x^2))
    CF(a); CF(b);
  }
}

/**
 @brief cmulti型の計算 y=arctan(x)
*/
void catan_r(cmulti *y, rmulti *x)
{
  cmulti *a=NULL;
  CAc(a,y);
  cset_r(a,x);
  catan_c(y,a);
  CF(a);
}

/**
 @brief cmulti型の計算 y=arctan(x)
*/
void catan_c(cmulti *y, cmulti *x)
{
  cmulti *a=NULL,*b=NULL;
  if(cis_real(x)){
    ratan_r(C_R(y),C_R(x));
    rset_i(C_I(y),0);
  }else{
    CAc(a,y); CAc(b,y);
    cadd_cd(a,x,1);    // a=1+x
    csub_dc(b,1,x);    // b=1-x
    cdiv_cc(a,a,b);    // a=(1+x)/(1-x)
    clog_c(a,a);      // a=log((1+x)/(1-x))
    cset_dd(b,0,-0.5); // b=-I/2
    cmul_cc(a,a,b);    // a=(-I/2)*log((1+x)/(1-x))
    cset_c(y,a);     // y=(-I/2)*log((1+x)/(1-x))
    CF(a); CF(b);
  }
}

/**
 @brief cmulti型の計算 y=sinh(x)
*/
void csinh_c(cmulti *y, cmulti *x)
{
  rmulti *c=NULL,*ch=NULL,*s=NULL,*sh=NULL;
  RAc(c,y); RAc(ch,y); RAc(s,y); RAc(sh,y);
  // z.r=sinh(x.r)*cos(x.i)
  // z.i=cosh(x.r)*sin(x.i)
  rsinh_r(sh,C_R(x));
  rcosh_r(ch,C_R(x));
  rsin_r(s,C_I(x));
  rcos_r(c,C_I(x));
  rmul_rr(C_R(y),sh,c);
  rmul_rr(C_I(y),ch,s);
  RF(c); RF(ch); RF(s); RF(sh);
}

/**
 @brief cmulti型の計算 y=cosh(x)
*/
void ccosh_c(cmulti *y, cmulti *x)
{
  rmulti *c=NULL,*ch=NULL,*s=NULL,*sh=NULL;
  RAc(c,y); RAc(ch,y); RAc(s,y); RAc(sh,y);
  // z.r=cosh(x.r)*cos(x.i)
  // z.i=sinh(x.r)*sin(x.i)
  rsinh_r(sh,C_R(x));
  rcosh_r(ch,C_R(x));
  rsin_r(s,C_I(x));
  rcos_r(c,C_I(x));
  rmul_rr(C_R(y),ch,c);
  rmul_rr(C_I(y),sh,s);
  RF(c); RF(ch); RF(s); RF(sh);
}

/**
 @brief cmulti型の計算 y=tanh(x)
*/
void ctanh_c(cmulti *y, cmulti *x)
{
  cmulti *ch=NULL,*sh=NULL;
  CAc(ch,y); CAc(sh,y);
  csinh_c(sh,x);
  ccosh_c(ch,x);
  cdiv_cc(y,sh,ch); // y=sinh(x)/cosh(x)
  CF(ch); CF(sh);
}

/**
 @brief cmulti型の計算 y=asinh(x)
*/
void casinh_r(cmulti *y, rmulti *x)
{
  cmulti *a=NULL;
  CAc(a,y);
  cset_r(a,x);
  casinh_c(y,a);
  CF(a);
}

/**
 @brief cmulti型の計算 y=asinh(x)
*/
void casinh_c(cmulti *y, cmulti *x)
{
  cmulti *a=NULL;  
  if(cis_real(x)){
    rasinh_r(C_R(y),C_R(x));
    rset_zero(C_I(y));
  }else{
    CAc(a,y);
    cmul_cc(a,x,x);
    cadd_cd(a,a,1);
    csqrt_c(a,a);
    cadd_cc(a,a,x);
    clog_c(a,a);
    cset_c(y,a);
    CF(a);
  }
}

/**
 @brief cmulti型の計算 y=acosh(x)
*/
void cacosh_r(cmulti *y, rmulti *x)
{
  cmulti *a=NULL;
  CAc(a,y);
  cset_r(a,x);
  cacosh_c(y,a);
  CF(a);
}

/**
 @brief cmulti型の計算 y=acosh(x)
*/
void cacosh_c(cmulti *y, cmulti *x)
{
  cmulti *a=NULL;  
  if(cis_real(x) && ge_rd(C_R(x),1)){
    racosh_r(C_R(y),C_R(x));
    rset_zero(C_I(y));
  }else{
    CAc(a,y);
    cmul_cc(a,x,x);
    csub_cd(a,a,1);
    csqrt_c(a,a);
    cadd_cc(a,a,x);
    clog_c(a,a);
    cset_c(y,a);
    CF(a);
  }
}

/**
 @brief cmulti型の計算 y=arctanh(x)
*/
void catanh_r(cmulti *y, rmulti *x)
{
  cmulti *a=NULL;
  CAc(a,y);
  cset_r(a,x);
  catanh_c(y,a);
  CF(a);
}

/**
 @brief cmulti型の計算 y=arctanh(x)
*/
void catanh_c(cmulti *y, cmulti *x)
{
  cmulti *a=NULL,*b=NULL;
  if(cis_real(x) && le_rd(C_R(x),1) && ge_rd(C_R(x),-1)){
    ratanh_r(C_R(y),C_R(x));
    rset_zero(C_I(y));
  }else{
    CAc(a,y); CAc(b,y);
    cadd_cd(a,x,1);     // a=1+x
    csub_dc(b,1,x);     // b=1-x
    cdiv_cc(a,a,b);     // a=(1+x)/(1-x)
    clog_c(a,a);       // a=log((1+x)/(1-x))
    cdiv_cd(a,a,2);     // a=(1/2)*log((1+x)/(1-x))
    cset_c(y,a);      // y=(1/2)*log((1+x)/(1-x))
    CF(a); CF(b);
  }
}

/** @} */

/////////////////////////////////////////////////////////////////

/** @name cmulti型の関数 z=f(x,y) */
/** @{ */

/**
 @brief z=x+y
 */
void cadd_cc(cmulti *z, cmulti *x, cmulti *y)
{
  NULL_EXC3(x,y,z);
  radd_rr(C_R(z),C_R(x),C_R(y)); // z.r=x.r+y.r
  radd_rr(C_I(z),C_I(x),C_I(y)); // z.i=x.i+y.i
}

/**
 @brief z=x+y
*/
void cadd_cz(cmulti *z, cmulti *x, dcomplex y)
{
  NULL_EXC2(x,z);
  radd_rd(C_R(z),C_R(x),Z_R(y)); // z.r=x.r+y.r
  radd_rd(C_I(z),C_I(x),Z_I(y)); // z.i=x.i+y.i
}

/**
 @brief z=x+y
 */
void cadd_zc(cmulti *z, dcomplex x, cmulti *y)
{
  NULL_EXC2(y,z);
  radd_dr(C_R(z),Z_R(x),C_R(y)); // z.r=x.r+y.r
  radd_dr(C_I(z),Z_I(x),C_I(y)); // z.i=x.i+y.i
}

/**
 @brief z=x+y
 */
void cadd_cr(cmulti *z, cmulti *x, rmulti *y)
{
  NULL_EXC3(x,y,z);
  radd_rr(C_R(z),C_R(x),y); // z.r=x.r+y
  rset_r(C_I(z),C_I(x)  ); // z.i=x.i
}

/**
 @brief z=x+y
 */
void cadd_rc(cmulti *z, rmulti *x, cmulti *y)
{
  NULL_EXC3(x,y,z);
  radd_rr(C_R(z),x,C_R(y)); // z.r=x+y.r
  rset_r(C_I(z),  C_I(y)); // z.i=  y.i
}

/**
 @brief z=x+y
 */
void cadd_rz(cmulti *z, rmulti *x, dcomplex y)
{
  NULL_EXC2(x,z);
  radd_rd(C_R(z),x,Z_R(y)); // z.r=x+y.r
  rset_d(C_I(z),  Z_I(y)); // z.i=  y.i
}

/**
 @brief z=x+y
 */
void cadd_zr(cmulti *z, dcomplex x, rmulti *y)
{
  NULL_EXC2(y,z);
  radd_dr(C_R(z),Z_R(x),y); // z.r=x.r+y
  rset_d(C_I(z),Z_I(x)  ); // z.i=x.i
}

/**
 @brief z=x+y
 */
void cadd_cd(cmulti *z, cmulti *x, double y)
{
  NULL_EXC2(x,z);
  radd_rd(C_R(z),C_R(x),y); // z.r=x.r+y
  rset_r(C_I(z),C_I(x)  ); // z.i=x.i
}

/**
 @brief z=x+y
 */
void cadd_dc(cmulti *z, double x, cmulti *y)
{
  NULL_EXC2(y,z);
  radd_dr(C_R(z),x,C_R(y)); // z.r=x+y.r
  rset_r(C_I(z),  C_I(y)); // z.i=  y.i
}

/**
 @brief z=x-y
 */
void csub_cc(cmulti *z, cmulti *x, cmulti *y)
{
  NULL_EXC3(x,y,z);
  rsub_rr(C_R(z),C_R(x),C_R(y)); // z.r=x.r-y.r
  rsub_rr(C_I(z),C_I(x),C_I(y)); // z.i=x.i-y.i
}

/**
 @brief z=x-y
 */
void csub_zc(cmulti *z, dcomplex x, cmulti *y)
{
  NULL_EXC2(y,z);
  rsub_dr(C_R(z),Z_R(x),C_R(y)); // z.r=x.r-y.r
  rsub_dr(C_I(z),Z_I(x),C_I(y)); // z.i=x.i-y.i
}

/**
 @brief z=x-y
 */
void csub_cz(cmulti *z, cmulti *x, dcomplex y)
{
  NULL_EXC2(x,z);
  rsub_rd(C_R(z),C_R(x),Z_R(y)); // z.r=x.r-y.r
  rsub_rd(C_I(z),C_I(x),Z_I(y)); // z.i=x.i-y.i
}

/**
 @brief z=x-y
 */
void csub_rc(cmulti *z, rmulti *x, cmulti *y)
{
  NULL_EXC3(x,y,z);
  rsub_rr(C_R(z),x,C_R(y)); // z.r=x-y.r
  rneg_r(C_I(z),  C_I(y)); // z.i= -y.i
}

/**
 @brief z=x-y
 */
void csub_cr(cmulti *z, cmulti *x, rmulti *y)
{
  NULL_EXC3(x,y,z);
  rsub_rr(C_R(z),C_R(x),y); // z.r=x.r-y
  rset_r(C_I(z),C_I(x)  ); // z.i=x.i
}

/**
 @brief z=x-y
 */
void csub_dc(cmulti *z, double x, cmulti *y)
{
  NULL_EXC2(y,z);
  rsub_dr(C_R(z),x,C_R(y)); // z.r=x-y.r
  rneg_r(C_I(z),  C_I(y)); // z.i= -y.i
}

/**
 @brief z=x-y
 */
void csub_cd(cmulti *z, cmulti *x, double y)
{
  NULL_EXC2(x,z);
  rsub_rd(C_R(z),C_R(x),y); // z.r=x.r-y
  rset_r(C_I(z),C_I(x)  ); // z.i=x.i
}

/**
 @brief z=x-y
 */
void csub_zr(cmulti *z, dcomplex x, rmulti *y)
{
  NULL_EXC2(y,z);
  rsub_dr(C_R(z),Z_R(x),y); // z.r=x.r-y
  rset_d(C_I(z),Z_I(x)  ); // z.i=x.i
}

/**
 @brief z=x-y
 */
void csub_rz(cmulti *z, rmulti *x, dcomplex y)
{
  NULL_EXC2(x,z);
  rsub_rd(C_R(z),x,Z_R(y)); // z.r=x.r-y.r
  rneg_d(C_I(z),  Z_I(y)); // z.i=   -y.i
}

/**
 @brief z=x*y
 */
void cmul_cc_ws(cmulti *z, cmulti *x, cmulti *y, int n_rws, rmulti **rws, int n_cws, cmulti **cws)
{
  NULL_EXC3(x,y,z);
  if(n_rws<1){ ERROR_EXIT("Error `n_rws=%d<1' in cmul_ws().\n",n_rws); }
  if(n_cws<1){ ERROR_EXIT("Error `n_cws=%d<1' in cmul_ws().\n",n_cws); }
  rmul_rr     (C_R(cws[0]),C_R(x),C_R(y));           // z.r =x.r*y.r
  rsub_mul_rr_ws(C_R(cws[0]),C_I(x),C_I(y),n_rws,rws); // z.r-=x.i*y.i
  rmul_rr     (C_I(cws[0]),C_I(x),C_R(y));           // z.i =x.i*y.r
  radd_mul_rr_ws(C_I(cws[0]),C_R(x),C_I(y),n_rws,rws); // z.i+=x.r*y.i
  cset_c(z,cws[0]);                                // z=cws[0]
}

/**
 @brief z=x*y
 */
void cmul_cc(cmulti *z, cmulti *x, cmulti *y)
{
  cmulti *a=NULL;
  NULL_EXC3(x,y,z);
  CAc(a,z);
  rmul_rr  (C_R(a),C_R(x),C_R(y)); // z.r =x.r*y.r
  rsub_mul_rr(C_R(a),C_I(x),C_I(y)); // z.r-=x.i*y.i
  rmul_rr  (C_I(a),C_I(x),C_R(y)); // z.i =x.i*y.r
  radd_mul_rr(C_I(a),C_R(x),C_I(y)); // z.i+=x.r*y.i
  cset_c(z,a);                     // z=a
  CF(a);
}

/**
 @brief z=x*y
 */
void cmul_cz(cmulti *z, cmulti *x, dcomplex y)
{
  cmulti *a=NULL;
  NULL_EXC2(x,z);
  CAc(a,z);
  cset_z(a,y);   // a=y
  cmul_cc(z,x,a); // z=x*y
  CF(a);
}

/**
 @brief z=x*y
 */
void cmul_zc(cmulti *z, dcomplex x, cmulti *y)
{
  cmulti *a=NULL;
  NULL_EXC2(y,z);
  CAc(a,z);
  cset_z(a,x);   // a=x
  cmul_cc(z,a,y); // z=x*y
  CF(a);
}

/**
 @brief z=x*y
 */
void cmul_cr(cmulti *z, cmulti *x, rmulti *y)
{
  NULL_EXC3(x,y,z);
  rmul_rr(C_R(z),C_R(x),y); // z.r=x.r*y
  rmul_rr(C_I(z),C_I(x),y); // z.i=x.i*y
}

/**
 @brief z=x*y
 */
void cmul_rc(cmulti *z, rmulti *x, cmulti *y)
{
  NULL_EXC3(x,y,z);
  rmul_rr(C_R(z),x,C_R(y)); // z.r=x*y.r
  rmul_rr(C_I(z),x,C_I(y)); // z.i=x*y.i
}

/**
 @brief z=x*y
 */
void cmul_cd(cmulti *z, cmulti *x, double y)
{
  NULL_EXC2(x,z);
  rmul_rd(C_R(z),C_R(x),y); // z.r=x.r*y
  rmul_rd(C_I(z),C_I(x),y); // z.i=x.i*y
}

/**
 @brief z=x*y
 */
void cmul_dc(cmulti *z, double x, cmulti *y)
{
  NULL_EXC2(y,z);
  rmul_dr(C_R(z),x,C_R(y)); // z.r=x*y.r
  rmul_dr(C_I(z),x,C_I(y)); // z.i=x*y.i
}

/**
 @brief z=x*y
 */
void cmul_rz(cmulti *z, rmulti *x, dcomplex y)
{
  NULL_EXC2(x,z);
  rmul_rd(C_R(z),x,Z_R(y)); // z.r=x*y.r
  rmul_rd(C_I(z),x,Z_I(y)); // z.i=x*y.i
}

/**
 @brief z=x*y
 */
void cmul_zr(cmulti *z, dcomplex x, rmulti *y)
{
  NULL_EXC2(y,z);
  rmul_dr(C_R(z),Z_R(x),y); // z.r=x.r*y
  rmul_dr(C_I(z),Z_I(x),y); // z.i=x.i*y
}

/**
 @brief z=conj(x)*y
 */
void cdot_cc(cmulti *z, cmulti *x, cmulti *y)
{
  cmulti *a=NULL;
  NULL_EXC3(x,y,z);
  CAc(a,z);
  rmul_rr  (C_R(a),C_R(x),C_R(y)); // z.r =x.r*y.r
  radd_mul_rr(C_R(a),C_I(x),C_I(y)); // z.r+=x.i*y.i
  rmul_rr  (C_I(a),C_R(x),C_I(y)); // z.i =x.r*y.i
  rsub_mul_rr(C_I(a),C_I(x),C_R(y)); // z.i-=x.i*y.r
  cset_c(z,a);                     // z=a
  CF(a);
}

/**
 @brief z=conj(x)*y
 */
void cdot_cz(cmulti *z, cmulti *x, dcomplex y)
{
  cmulti *a=NULL;
  NULL_EXC2(x,z);
  CAc(a,z);
  cset_z(a,y); // a=y
  cdot_cc(z,x,a); // z=conj(x)*y
  CF(a);
}

/**
 @brief z=conj(x)*y
 */
void cdot_zc(cmulti *z, dcomplex x, cmulti *y)
{
  cmulti *a=NULL;
  NULL_EXC2(y,z);
  CAc(a,z);
  cset_z(a,x); // a=x
  cdot_cc(z,a,y); // z=conj(x)*y  
  CF(a);
}

/**
 @brief z=conj(x)*y
 */
void cdot_cr(cmulti *z, cmulti *x, rmulti *y)
{
  NULL_EXC3(x,y,z);
  rmul_rr(C_R(z),C_R(x),y);                        // z.r= x.r*y
  rmul_rr(C_I(z),C_I(x),y); rneg_r(C_I(z),C_I(z)); // z.i=-x.i*y
}

/**
 @brief z=conj(x)*y
 */
void cdot_rc(cmulti *z, rmulti *x, cmulti *y)
{
  NULL_EXC3(x,y,z);
  rmul_rr(C_R(z),x,C_R(y)); // z.r=x*y.r
  rmul_rr(C_I(z),x,C_I(y)); // z.i=x*y.i
}

/**
 @brief z=conj(x)*y
 */
void cdot_cd(cmulti *z, cmulti *x, double y)
{
  NULL_EXC2(x,z);
  rmul_rd(C_R(z),C_R(x),y);                        // z.r= x.r*y
  rmul_rd(C_I(z),C_I(x),y); rneg_r(C_I(z),C_I(z)); // z.i=-x.i*y
}

/**
 @brief z=conj(x)*y
 */
void cdot_dc(cmulti *z, double x, cmulti *y)
{
  NULL_EXC2(y,z);
  rmul_dr(C_R(z),x,C_R(y)); // z.r=x*y.r
  rmul_dr(C_I(z),x,C_I(y)); // z.i=x*y.i
}

/**
 @brief z=conj(x)*y
 */
void cdot_rz(cmulti *z, rmulti *x, dcomplex y)
{
  NULL_EXC2(x,z);
  rmul_rd(C_R(z),x,Z_R(y)); // z.r=x*y.r
  rmul_rd(C_I(z),x,Z_I(y)); // z.i=x*y.i
}

/**
 @brief z=conj(x)*y
 */
void cdot_zr(cmulti *z, dcomplex x, rmulti *y)
{
  NULL_EXC2(y,z);
  rmul_dr(C_R(z), Z_R(x),y); // z.r= x.r*y
  rmul_dr(C_I(z),-Z_I(x),y); // z.i=-x.i*y
}

/**
 @brief z=x/y
 */
void cdiv_cc(cmulti *z, cmulti *x, cmulti *y)
{
  rmulti *den=NULL;
  cmulti *a=NULL;
  NULL_EXC3(x,y,z);
  RAc(den,z); CAc(a,z);
  rabs2_c(den,y);    // den=|y|^2
  cdot_cc(a,y,x);   // z=conj(y)*x
  cdiv_cr(a,a,den); // z/=den
  cset_c(z,a);    // z=a
  RF(den); CF(a);
}

/**
 @brief z=x/y
 */
void cdiv_zc(cmulti *z, dcomplex x, cmulti *y)
{
  cmulti *a=NULL;
  NULL_EXC2(y,z);
  CAc(a,z);
  cset_z(a,x);   // a=x
  cdiv_cc(z,a,y); // z=x/y
  CF(a);
}

/**
 @brief z=x/y
 */
void cdiv_cz(cmulti *z, cmulti *x, dcomplex y)
{
  cmulti *a=NULL;
  NULL_EXC2(x,z);
  CAc(a,z);
  cset_z(a,y);   // a=y
  cdiv_cc(z,x,a); // z=x/y
  CF(a);
}

/**
 @brief z=x/y
 */
void cdiv_rc(cmulti *z, rmulti *x, cmulti *y)
{
  rmulti *den=NULL;
  cmulti *a=NULL;
  NULL_EXC3(x,y,z);
  RAc(den,z); CAc(a,z);
  rabs2_c(den,y);      // den=|y|^2
  cconj_c(a,y);        // z=conj(y)
  cmul_cr(a,a,x);     // z*=x
  cdiv_cr(a,a,den);   // z/=den
  cset_c(z,a);      // z=a
  RF(den); CF(a);
}

/**
 @brief z=x/y
 */
void cdiv_cr(cmulti *z, cmulti *x, rmulti *y)
{
  NULL_EXC3(x,y,z);
  rdiv_rr(C_R(z),C_R(x),y); // z.r=x.r/y
  rdiv_rr(C_I(z),C_I(x),y); // z.i=x.i/y
}

/**
 @brief z=x/y
 */
void cdiv_dc(cmulti *z, double x, cmulti *y)
{
  rmulti *a=NULL;
  NULL_EXC2(y,z);
  RAc(a,z);
  rset_d(a,x);   // a=x
  cdiv_rc(z,a,y); // z=x/y
  RF(a);
}

/**
 @brief z=x/y
 */
void cdiv_cd(cmulti *z, cmulti *x, double y)
{
  NULL_EXC2(x,z);
  rdiv_rd(C_R(z),C_R(x),y); // z.r=x.r/y
  rdiv_rd(C_I(z),C_I(x),y); // z.i=x.i/y
}


/**
 @brief z=x/y
 */
void cdiv_zr(cmulti *z, dcomplex x, rmulti *y)
{
  NULL_EXC2(z,y);
  rdiv_dr(C_R(z),Z_R(x),y); // z.r=x.r/y
  rdiv_dr(C_I(z),Z_I(x),y); // z.i=x.i/y
}

/**
 @brief z=x/y
 */
void cdiv_rz(cmulti *z, rmulti *x, dcomplex y)
{
  cmulti *a=NULL;
  NULL_EXC2(x,z);
  CAc(a,z);
  cset_z(a,y);   // a=y
  cdiv_rc(z,x,a); // z=x/y
  CF(a);
}

/**
 @brief z=abs(x-y)
 */
void rabs_sub_cc(rmulti *z, cmulti *x, cmulti *y)
{
  cmulti *a;
  NULL_EXC3(x,y,z);
  CAr(a,z);
  csub_cc(a,x,y); // a=x-y
  rabs_c(z,a);    // z=abs(x-y)
  CF(a);
}

/**
 @brief z=abs(x-y)
 */
void rabs_sub_cr(rmulti *z, cmulti *x, rmulti *y)
{
  cmulti *a;
  NULL_EXC3(x,y,z);
  CAr(a,z);
  csub_cr(a,x,y); // a=x-y
  rabs_c(z,a);     // z=abs(x-y)
  CF(a);
}


/**
 @brief z=x^y
 */
void cpow_cc(cmulti *z, cmulti *x, cmulti *y)
{
  cmulti *a=NULL;
  CAc(a,z);
  clog_c(a,x); // a=log(x)
  cmul_cc(a,y,a); // a=y*log(x)
  cexp_c(z,a); // z=exp(y*log(x))
  CF(a);
}

/**
 @brief z=x^y
 */
void cpow_rc(cmulti *z, rmulti *x, cmulti *y)
{
  rmulti *a=NULL;
  cmulti *b=NULL;
  RAc(a,z); CAc(b,z);
  rlog_r(a,x);     // a=log(x)
  cmul_cr(b,y,a); // b=y*log(x)
  cexp_c(z,b);   // z=exp(y*log(x))
  RF(a); CF(b);
}

/**
 @brief z=x^y
 */
void cpow_cr(cmulti *z, cmulti *x, rmulti *y)
{
  rmulti *r=NULL,*theta=NULL;
  RAc(r,z); RAc(theta,z);
  cget_polar(r,theta,x); // x=r*exp(i*theta)
  rpow_rr(r,r,y);           // r=r^y
  rmul_rr(theta,theta,y);   // theta=theta*y
  cset_polar(z,r,theta); // z=r^y*cos(theta*y)+i*r^y*sin(theta*y)
  RF(r); RF(theta);
}

/**
 @brief z=x^y
 */
void cpow_cd(cmulti *z, cmulti *x, double y)
{
  rmulti *r=NULL,*theta=NULL;
  RAc(r,z); RAc(theta,z);
  cget_polar(r,theta,x); // x=r*exp(i*theta)
  rpow_rd(r,r,y);        // r=r^y
  rmul_rd(theta,theta,y); // theta=theta*y
  cset_polar(z,r,theta); // z=r^y*cos(theta*y)+i*r^y*sin(theta*y)
  RF(r); RF(theta);
}

/**
 @brief cmulti型のべき乗 y=x^n
 */
void cpow_c(cmulti *y, cmulti *x, int n)
{
  cmulti *a=NULL;
  long i;
  if(n==0){
    cset_one(y);
  }else if(cis_real(x)){
    rpow_r(C_R(y),C_R(x),n);
    rset_zero(C_I(y));
  }else if(n>0){
    CAc(a,y);
    cset_c(a,x);                        // a=x
    for(i=1; i<n; i++){ cmul_cc(a,a,x); } // a*=x
    cset_c(y,a);                        // y=a
    CF(a);
  }else if(n<0){
    CAc(a,y);
    cset_c(a,x);                           // a=x
    for(i=1; i<(-n); i++){ cmul_cc(a,a,x); } // a*=x
    cinv_c(a,a);                            // a=1/a
    cset_c(y,a);                           // y=a
    CF(a);
  }else{ ERROR_AT; }
}

/**
 @brief y=x*2^n
 */
void cmul_2exp(cmulti *y, cmulti *x, int nr, int ni)
{
  rmul_2exp(C_R(y),C_R(x),nr);
  rmul_2exp(C_I(y),C_I(x),ni);
}

/**
 @brief y=x/2^n
 */
void cdiv_2exp(cmulti *y, cmulti *x, int nr, int ni)
{
  rdiv_2exp(C_R(y),C_R(x),nr);
  rdiv_2exp(C_I(y),C_I(x),ni);
}

/**
 @brief z=10^(floor(log10(abs(x))-y))
 */
void cexp10_floor_log10_abs_sub(cmulti *z, cmulti *x, double y)
{
  rmulti *exp1=NULL,*exp2=NULL;
  RAc(exp1,x); RAc(exp2,x);
  rabs_r(C_R(x),C_R(x));
  rabs_r(C_I(x),C_I(x));
  rlog10_r(C_R(x),C_R(x));
  rlog10_r(C_I(x),C_I(x));
  rsub_rd(C_R(x),C_R(x),y);
  rsub_rd(C_I(x),C_I(x),y);
  rfloor_r(exp1,C_R(x));
  rfloor_r(exp2,C_I(x));
  rexp10_r(C_R(z),exp1);
  rexp10_r(C_I(z),exp2);
  RF(exp1);RF(exp2);
}

/**
 @brief z=2^(floor(log2(abs(x))-y))
 */
void cexp2_floor_log2_abs_sub(cmulti *z, cmulti *x, double y)
{
  rmulti *exp1=NULL,*exp2=NULL;
  RAc(exp1,x); RAc(exp2,x);
  rabs_r(C_R(x),C_R(x));
  rabs_r(C_I(x),C_I(x));
  rlog2_r(C_R(x),C_R(x));
  rlog2_r(C_I(x),C_I(x));
  rsub_rd(C_R(x),C_R(x),y);
  rsub_rd(C_I(x),C_I(x),y);
  rfloor_r(exp1,C_R(x));
  rfloor_r(exp2,C_I(x));
  rexp2_r(C_R(z),exp1);
  rexp2_r(C_I(z),exp2);
  RF(exp1);RF(exp2);
}


/** @} */

/////////////////////////////////////////////////////////////////

/** @name cmulti型の関数 z=z+f(x,y) */
/** @{ */

/**
 @brief z=z+x*y
 */
void cadd_mul_cc(cmulti *z, cmulti *x, cmulti *y)
{
  cmulti *a=NULL;
  NULL_EXC3(x,y,z);
  CAc(a,z);
  cmul_cc(a,x,y);   // a=x*y
  cadd_cc(z,z,a);   // z=z+a
  CF(a);
}

/**
 @brief z=z+x*y
 */
void cadd_mul_cr(cmulti *z, cmulti *x, rmulti *y)
{
  cmulti *a=NULL;
  NULL_EXC3(x,y,z);
  CAc(a,z);
  cmul_cr(a,x,y); // a=x*y
  cadd_cc(z,z,a); // z=z+a
  CF(a);
}

/**
 @brief z=z+x*y
 */
void cadd_mul_cz(cmulti *z, cmulti *x, dcomplex y)
{
  cmulti *a=NULL;
  NULL_EXC2(x,z);
  CAc(a,z);
  cmul_cz(a,x,y); // a=x*y
  cadd_cc(z,z,a); // z=z+a
  CF(a);
}

/**
 @brief z=z+x*y
 */
void cadd_mul_cd(cmulti *z, cmulti *x, double y)
{
  cmulti *a=NULL;
  NULL_EXC2(x,z);
  CAc(a,z);
  cmul_cd(a,x,y); // a=x*y
  cadd_cc(z,z,a); // z=z+a
  CF(a);
}

/**
 @brief z=z-x*y
 */
void csub_mul_cc(cmulti *z, cmulti *x, cmulti *y)
{
  cmulti *a=NULL;
  NULL_EXC3(x,y,z);
  CAc(a,z);
  cmul_cc(a,x,y);   // a=x*y
  csub_cc(z,z,a);   // z=z-a
  CF(a);
}

/**
 @brief z=z-x*y
 */
void csub_mul_cc_ws(cmulti *z, cmulti *x, cmulti *y, int n_rws, rmulti **rws, int n_cws, cmulti **cws)
{
  NULL_EXC3(x,y,z);
  if(n_rws<1){ ERROR_EXIT("Error `n_rws=%d<1' in csub_mul_ws().\n",n_rws); }
  if(n_cws<2){ ERROR_EXIT("Error `n_cws=%d<2' in csub_mul_ws().\n",n_cws); }
  cmul_cc_ws(cws[0],x,y,n_rws,rws,n_cws-1,cws+1); // cws[0]=x*y
  csub_cc(z,z,cws[0]);                            // z=z-cws[0]
}

/**
 @brief z=z-x*y
 */
void csub_mul_cr(cmulti *z, cmulti *x, rmulti *y)
{
  cmulti *a=NULL;
  NULL_EXC3(x,y,z);
  CAc(a,z);
  cmul_cr(a,x,y); // a=x*y
  csub_cc(z,z,a); // z=z-a
  CF(a);
}

/**
 @brief z=z-x*y
 */
void csub_mul_cz(cmulti *z, cmulti *x, dcomplex y)
{
  cmulti *a=NULL;
  NULL_EXC2(x,z);
  CAc(a,z);
  cmul_cz(a,x,y); // a=x*y
  csub_cc(z,z,a); // z=z-a
  CF(a);
}

/**
 @brief z=z-x*y
 */
void csub_mul_cd(cmulti *z, cmulti *x, double y)
{
  cmulti *a=NULL;
  NULL_EXC2(x,z);
  CAc(a,z);
  cmul_cd(a,x,y); // a=x*y
  csub_cc(z,z,a); // z=z-a
  CF(a);
}

/**
 @brief z=z+conj(x)*y
 */
void cadd_dot_cc(cmulti *z, cmulti *x, cmulti *y)
{
  cmulti *a=NULL;
  NULL_EXC3(x,y,z);
  CAc(a,z);
  cdot_cc(a,x,y);   // a=conj(x)*y
  cadd_cc(z,z,a);   // z=z+a
  CF(a);
}

/**
 @brief z=z+conj(x)*y
 */
void cadd_dot_zc(cmulti *z, dcomplex x, cmulti *y)
{
  cmulti *a=NULL;
  NULL_EXC2(y,z);
  CAc(a,z);
  cdot_zc(a,x,y); // a=conj(x)*y
  cadd_cc(z,z,a); // z=z+a
  CF(a);
}

/**
 @brief z=z+conj(x)*y
 */
void cadd_dot_cz(cmulti *z, cmulti *x, dcomplex y)
{
  cmulti *a=NULL;
  NULL_EXC2(x,z);
  CAc(a,z);
  cdot_cz(a,x,y); // a=conj(x)*y
  cadd_cc(z,z,a); // z=z+a
  CF(a);
}


/**
 @brief z=z-conj(x)*y
 */
void csub_dot_cc(cmulti *z, cmulti *x, cmulti *y)
{
  cmulti *a=NULL;
  NULL_EXC3(x,y,z);
  CAc(a,z);
  cdot_cc(a,x,y);   // a=conj(x)*y
  csub_cc(z,z,a);   // z=z-a
  CF(a);
}

/**
 @brief z=z-conj(x)*y
 */
void csub_dot_zc(cmulti *z, dcomplex x, cmulti *y)
{
  cmulti *a=NULL;
  NULL_EXC2(y,z);
  CAc(a,z);
  cdot_zc(a,x,y);   // a=conj(x)*y
  csub_cc(z,z,a);   // z=z-a
  CF(a);
}

/**
 @brief z=z-conj(x)*y
 */
void csub_dot_cz(cmulti *z, cmulti *x, dcomplex y)
{
  cmulti *a=NULL;
  NULL_EXC2(x,z);
  CAc(a,z);
  cdot_cz(a,x,y);   // a=conj(x)*y
  csub_cc(z,z,a);   // z=z-a
  CF(a);
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
int cmp_cc(cmulti *x, cmulti *y)
{
  rmulti *xr=NULL,*yr=NULL,*xt=NULL,*yt=NULL;
  int value;
  if(ccmp_get_type()){
    RAc(xr,x); RAc(yr,y); RAc(xt,x); RAc(yt,y);
    cget_polar(xr,xt,x);
    cget_polar(yr,yt,y);
    value=cmp_rr(xr,yr);
    if(value==0){ value=cmp_rr(xt,yt); }
    RF(xr); RF(yr); RF(xt); RF(yt);
    return value;
  }else{
    value=cmp_rr(C_R(x),C_R(y));
    if(value!=0) return value;
    return cmp_rr(C_I(x),C_I(y));
  }
}

/**
 @brief cmulti型の値の比較 x<=>y
*/
int cmp_cz(cmulti *x, dcomplex y)
{
  int value;
  value=cmp_rd(C_R(x),Z_R(y));
  if(value!=0) return value;
  return cmp_rd(C_I(x),Z_I(y));
}

/**
 @brief cmulti型の値の比較 x<=>y
*/
int cmp_rc(rmulti *x, cmulti *y)
{
  return -cmp_cr(y,x);
}

/**
 @brief cmulti型の値の比較 x<=>y
*/
int cmp_cr(cmulti *x, rmulti *y)
{
  int value;
  value=cmp_rr(C_R(x),y);
  if(value!=0) return value;
  return cmp_rd(C_I(x),0);
}

/**
 @brief cmulti型の値の比較 x<=>y
*/
int cmp_cd(cmulti *x, double y)
{
  int value;
  value=cmp_rd(C_R(x),y);
  if(value!=0) return value;
  return cmp_rd(C_I(x),0);
}

/** @brief cmulti型の値の比較 x==y */
int eq_cc(cmulti *x, cmulti *y)    { return cis_number(x) && cis_number(y) && cmp_cc   (x,y)==0; }
/** @brief cmulti型の値の比較 x==y */
int eq_cz(cmulti *x, dcomplex y) { return cis_number(x) && cmp_cz (x,y)==0; }
/** @brief cmulti型の値の比較 x==y */
int eq_cr(cmulti *x, rmulti *y)  { return cis_number(x) && cmp_cr(x,y)==0; }
/** @brief cmulti型の値の比較 x==y */
int eq_cd(cmulti *x, double y)   { return cis_number(x) && cmp_cd (x,y)==0; }
/** @brief cmulti型の値の比較 x!=y */
int ne_cc(cmulti *x, cmulti *y)    { return !(cis_number(x) && cis_number(y) && cmp_cc   (x,y)==0); }
/** @brief cmulti型の値の比較 x!=y */
int ne_cz(cmulti *x, dcomplex y) { return !(cis_number(x) && cmp_cz (x,y)==0); }
/** @brief cmulti型の値の比較 x!=y */
int ne_cr(cmulti *x, rmulti *y)  { return !(cis_number(x) && cmp_cr(x,y)==0); }
/** @brief cmulti型の値の比較 x!=y */
int ne_cd(cmulti *x, double y)   { return !(cis_number(x) && cmp_cd (x,y)==0); }
/** @brief cmulti型の値の比較 x>y */
int gt_cc(cmulti *x, cmulti *y)     { return cis_number(x) && cis_number(y) && cmp_cc   (x,y)>0; }
/** @brief cmulti型の値の比較 x>y */
int gt_rc(rmulti *x, cmulti *y)  { return ris_number(x) && cis_number(y) && cmp_rc(x,y)>0; }
/** @brief cmulti型の値の比較 x>y */
int gt_cr(cmulti *x, rmulti *y)  { return cis_number(x) && ris_number(y) && cmp_cr(x,y)>0; }
/** @brief cmulti型の値の比較 x>y */
int gt_zc(dcomplex x, cmulti *y) { return cis_number(y) && cmp_cz (y,x)<0; }
/** @brief cmulti型の値の比較 x>y */
int gt_cz(cmulti *x, dcomplex y) { return cis_number(x) && cmp_cz (x,y)>0; }
/** @brief cmulti型の値の比較 x>y */
int gt_dc(double x, cmulti *y)   { return cis_number(y) && cmp_cd (y,x)<0; }
/** @brief cmulti型の値の比較 x>y */
int gt_cd(cmulti *x, double y)   { return cis_number(x) && cmp_cd (x,y)>0; }
/** @brief cmulti型の値の比較 x>=y */
int ge_cc(cmulti *x, cmulti *y)    { return cis_number(x) && cis_number(y) && cmp_cc   (x,y)>=0; }
/** @brief cmulti型の値の比較 x>=y */
int ge_rc(rmulti *x, cmulti *y) { return ris_number(x) && cis_number(y) && cmp_rc(x,y)>=0; }
/** @brief cmulti型の値の比較 x>=y */
int ge_cr(cmulti *x, rmulti *y) { return cis_number(x) && ris_number(y) && cmp_cr(x,y)>=0; }
/** @brief cmulti型の値の比較 x>=y */
int ge_zc(dcomplex x, cmulti *y){ return cis_number(y) && cmp_cz (y,x)<=0; }
/** @brief cmulti型の値の比較 x>=y */
int ge_cz(cmulti *x, dcomplex y){ return cis_number(x) && cmp_cz (x,y)>=0; }
/** @brief cmulti型の値の比較 x>=y */
int ge_dc(double x, cmulti *y)  { return cis_number(y) && cmp_cd (y,x)<=0; }
/** @brief cmulti型の値の比較 x>=y */
int ge_cd(cmulti *x, double y)  { return cis_number(x) && cmp_cd (x,y)>=0; }
/** @brief cmulti型の値の比較 x<y */
int lt_cc(cmulti *x, cmulti *y)    { return cis_number(x) && cis_number(y) && cmp_cc   (x,y)<0; }
/** @brief cmulti型の値の比較 x<y */
int lt_rc(rmulti *x, cmulti *y) { return ris_number(x) && cis_number(y) && cmp_rc(x,y)<0; }
/** @brief cmulti型の値の比較 x<y */
int lt_cr(cmulti *x, rmulti *y) { return cis_number(x) && ris_number(y) && cmp_cr(x,y)<0; }
/** @brief cmulti型の値の比較 x<y */
int lt_zc(dcomplex x, cmulti *y){ return cis_number(y) && cmp_cz (y,x)>0; }
/** @brief cmulti型の値の比較 x<y */
int lt_cz(cmulti *x, dcomplex y){ return cis_number(x) && cmp_cz (x,y)<0; }
/** @brief cmulti型の値の比較 x<y */
int lt_dc(double x, cmulti *y)  { return cis_number(y) && cmp_cd (y,x)>0; }
/** @brief cmulti型の値の比較 x<y */
int lt_cd(cmulti *x, double y)  { return cis_number(x) && cmp_cd (x,y)<0; }
/** @brief cmulti型の値の比較 x<=y */
int le_cc(cmulti *x, cmulti *y)    { return cis_number(x) && cis_number(y) && cmp_cc   (x,y)<=0; }
/** @brief cmulti型の値の比較 x<=y */
int le_rc(rmulti *x, cmulti *y) { return ris_number(x) && cis_number(y) && cmp_rc(x,y)<=0; }
/** @brief cmulti型の値の比較 x<=y */
int le_cr(cmulti *x, rmulti *y) { return cis_number(x) && ris_number(y) && cmp_cr(x,y)<=0; }
/** @brief cmulti型の値の比較 x<=y */
int le_zc(dcomplex x, cmulti *y){ return cis_number(y) && cmp_cz (y,x)>=0; }
/** @brief cmulti型の値の比較 x<=y */
int le_cz(cmulti *x, dcomplex y){ return cis_number(x) && cmp_cz (x,y)<=0; }
/** @brief cmulti型の値の比較 x<=y */
int le_dc(double x, cmulti *y)  { return cis_number(y) && cmp_cd (y,x)>=0; }
/** @brief cmulti型の値の比較 x<=y */
int le_cd(cmulti *x, double y)  { return cis_number(x) && cmp_cd (x,y)<=0; }

/** @brief cmulti型の値の比較 abs(x)<=>abs(y) */
int cmp_abs_cc(cmulti *x, cmulti *y)
{
  int value;  
  rmulti *ax=NULL,*ay=NULL;
  RAc(ax,x); RAc(ay,y);
  rabs2_c(ax,x); // ax=abs(x)
  rabs2_c(ay,y); // ay=abs(y)
  value=cmp_rr(ax,ay);
  RF(ax); RF(ay);  
  return value;
}

/** @brief cmulti型の値の比較 abs(x)==abs(y) */
int eq_abs_cc(cmulti *x, cmulti *y) { return cis_number(x) && cis_number(y) && cmp_abs_cc(x,y)==0; }
/** @brief cmulti型の値の比較 abs(x)>abs(y) */
int gt_abs_cc(cmulti *x, cmulti *y) { return cis_number(x) && cis_number(y) && cmp_abs_cc(x,y)>0;  }
/** @brief cmulti型の値の比較 abs(x)>=abs(y) */
int ge_abs_cc(cmulti *x, cmulti *y) { return cis_number(x) && cis_number(y) && cmp_abs_cc(x,y)>=0; }
/** @brief cmulti型の値の比較 abs(x)<abs(y) */
int lt_abs_cc(cmulti *x, cmulti *y) { return cis_number(x) && cis_number(y) && cmp_abs_cc(x,y)<0;  }
/** @brief cmulti型の値の比較 abs(x)<=abs(y) */
int le_abs_cc(cmulti *x, cmulti *y) { return cis_number(x) && cis_number(y) && cmp_abs_cc(x,y)<=0; }
/** @brief cmulti型の値の比較 x.r==y.r && x.i==y.i */
int ceq_cc(cmulti *x, cmulti *y){ return cis_number(x) && cis_number(y) && eq_rr(C_R(x),C_R(y)) && eq_rr(C_I(x),C_I(y)); }
/** @brief cmulti型の値の比較 x.r>y.r && x.i>y.i */
int cgt_cc(cmulti *x, cmulti *y){ return cis_number(x) && cis_number(y) && gt_rr(C_R(x),C_R(y)) && gt_rr(C_I(x),C_I(y)); }
/** @brief cmulti型の値の比較 x.r>=y.r && x.i>=y.i */
int cge_cc(cmulti *x, cmulti *y){ return cis_number(x) && cis_number(y) && ge_rr(C_R(x),C_R(y)) && ge_rr(C_I(x),C_I(y)); }
/** @brief cmulti型の値の比較 x.r<y.r && x.i<y.i */
int clt_cc(cmulti *x, cmulti *y){ return cis_number(x) && cis_number(y) && lt_rr(C_R(x),C_R(y)) && lt_rr(C_I(x),C_I(y)); }
/** @brief cmulti型の値の比較 x.r<=y.r && x.i<=y.i */
int cle_cc(cmulti *x, cmulti *y){ return cis_number(x) && cis_number(y) && le_rr(C_R(x),C_R(y)) && le_rr(C_I(x),C_I(y)); }

/** @} */

//EOF
