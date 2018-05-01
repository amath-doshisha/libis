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
void cround(cmulti *x, int prec)
{
  rround(C_R(x),prec);
  rround(C_I(x),prec);
}

/**
 @brief cmulti型の値を複成.
*/
void cclone(cmulti *y, cmulti *x)
{
  rclone(C_R(y),C_R(x));
  rclone(C_I(y),C_I(x));
}

/**
 @brief cmulti型の値をrmulti型から複成.
*/
void cclone_r(cmulti *y, rmulti *x)
{
  rclone(C_R(y),x);
  rset_zero(C_I(y));
}

/**
 @brief cmulti型の値をrmulti型から複成.
*/
void cclone_rr(cmulti *y, rmulti *x_r, rmulti *x_i)
{
  rclone(C_R(y),x_r);
  rclone(C_I(y),x_i);
}

/**
 @brief cmulti型の値のコピー y=x.
*/
void ccopy(cmulti *y, cmulti *x)
{
  if(y==x){ return; }
  rcopy(C_R(y),C_R(x));
  rcopy(C_I(y),C_I(x));
}

/**
 @brief cmulti型の値のコピー y=x.
*/
void ccopy_rr(cmulti *y, rmulti *x_r, rmulti *x_i)
{
  rcopy(C_R(y),x_r);
  rcopy(C_I(y),x_i);
}

/**
 @brief cmulti型の値のコピー y=x.
*/
void ccopy_r(cmulti *y, rmulti *x)
{
  rcopy(C_R(y),x);
  rset_zero(C_I(y));
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
void cset_z(cmulti *x, dcomplex value)
{
  rset_d(C_R(x),Z_R(value));
  rset_d(C_I(x),Z_I(value));
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
 @brief cmulti型の浮動小数点数を倍精度浮動小数点数から設定.
*/
void cset_d(cmulti *x, double real)
{
  rset_d(C_R(x),real);
  rset_d(C_I(x),0);
}

/**
 @brief cmulti型の浮動小数点数を符号なし整数から設定.
*/
void cset_ui(cmulti *x, ulong real)
{
  rset_ui(C_R(x),real);
  rset_ui(C_I(x),0);
}

/**
 @brief cmulti型の浮動小数点数を符号あり整数から設定.
*/
void cset_si(cmulti *x, long real)
{
  rset_si(C_R(x),real);
  rset_si(C_I(x),0);
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

/** @name cmulti型の１入力の数学演算子 */
/** @{ */

/**
 @brief cmulti型の符号反転 y=-x
*/
void cneg(cmulti *y, cmulti *x)
{
  rneg(C_R(y),C_R(x));
  rneg(C_I(y),C_I(x));
}


/**
 @brief cmulti型の逆数 z=1/x
*/
void cinv(cmulti *z, cmulti *x)
{
  rmulti *den=NULL;
  cmulti *a=NULL;
  RAc(den,z); CAc(a,z);
  cabs2(den,x);     // den=|x|^2
  cconj(a,x);       // z=conj(x)
  cdiv_r2(a,a,den); // z/=den
  ccopy(z,a);       // z=a
  RF(den); CF(a);
}

//追加

/**
 @brief cmulti型の逆数 z=1/x
*/
void cinv_ws(cmulti *z, cmulti *x, int n_rws, rmulti **rws, int n_cws, cmulti **cws)
{
  if(n_rws<2){ ERROR_EXIT("Error `n_rws=%d<1' in cinv_ws().\n",n_rws); }
  if(n_cws<1){ ERROR_EXIT("Error `n_cws=%d<1' in cinv_ws().\n",n_cws); }
  cabs2_ws(rws[0],x,n_rws-1,rws+1); // rws[0]=|x|^2
  cconj(cws[0],x);                    // z=conj(x)
  cdiv_r2(cws[0],cws[0],rws[0]);      // z/=rws[0]
  ccopy(z,cws[0]);                    // z=cws[0]
}

/**
 @brief cmulti型の複素共役 y=conj(x).
*/
void cconj(cmulti *y, cmulti *x)
{
  rcopy(C_R(y),C_R(x));
  rneg(C_I(y),C_I(x));
}

/**
 @brief cmulti型の複素共役をとり複製 y=conj(x).
*/
void cconj_clone(cmulti *y, cmulti *x)
{
  rclone(C_R(y),C_R(x));
  rclone(C_I(y),C_I(x));
  rneg(C_I(y),C_I(y));
}

/**
 @brief cmulti型の絶対値 y=abs(x)
*/
void cabsv(rmulti *y, cmulti *x)
{
  cabs2(y,x); // y=abs(x)^2
  rsqrt(y,y); // y=sqrt(y)
}

/**
 @brief cmulti型の絶対値 y=abs(x)
*/
void cabsv_ws(rmulti *y, cmulti *x, int n_rws, rmulti **rws)
{
  cabs2_ws(y,x,n_rws,rws); // y=abs(x)^2
  rsqrt(y,y);              // y=sqrt(y)
}

/**
 @brief cmulti型の絶対値の加算 y+=abs(x)
*/
void cadd_abs(rmulti *y, cmulti *x)
{
  rmulti *a=NULL;
  RAr(a,y);
  cabsv(a,x);  // a=abs(x)
  radd(y,y,a); // y=y+a
  RF(a);
}

/**
 @brief cmulti型の絶対値の平方 y=abs(x)^2
*/
void cabs2(rmulti *y, cmulti *x)
{
  rmul(y,C_R(x),C_R(x));     // y=x.r*x.r
  radd_mul(y,C_I(x),C_I(x)); // y+=x.i*x.i
}

/**
 @brief cmulti型の絶対値の平方 y=abs(x)^2
*/
void cabs2_ws(rmulti *y, cmulti *x, int n_rws, rmulti **rws)
{
  rmul(y,C_R(x),C_R(x));                  // y=x.r*x.r
  radd_mul_ws(y,C_I(x),C_I(x),n_rws,rws); // y+=x.i*x.i
}

/**
 @brief cmulti型の絶対値の平方の加算 y+=abs(x)^2
*/
void cadd_abs2(rmulti *y, cmulti *x)
{
  rmulti *a=NULL;
  RAr(a,y);
  cabs2(a,x);  // a=abs(x)^2
  radd(y,y,a); // y=y+a
  RF(a);
}

/**
 @brief cmulti型の実部と虚部の絶対値 y=abs(x.r)+i*abs(x.i)
*/
void cabsc(cmulti *y, cmulti *x)
{
  rabs(C_R(y),C_R(x));
  rabs(C_I(y),C_I(x));
}

/**
 @brief cmulti型の実部と虚部の絶対値 y=max(abs(x.r),abs(x.i))
*/
void cmax_absc(rmulti *y, cmulti *x)
{
  rmulti *a=NULL;
  RAc(a,x);
  rabs(y,C_R(x));
  rabs(a,C_I(x));
  if(rgt(a,y)){ rcopy(y,a); }
  RF(a);
}

/**
 @brief cmulti型の規格化 y=x/abs(x)
*/
void cnormalize(cmulti *y, cmulti *x)
{
  rmulti *a=NULL;
  RAc(a,y);
  cabsv(a,x);     // a=abs(x)
  cdiv_r2(y,x,a); // y=x/a
  RF(a);
}

/**
 @brief cmulti型の偏角 theta=arg(z)
*/
void cargument(rmulti *theta, cmulti *z)
{
  ratan2(theta,C_I(z),C_R(z));
}

/**
 @brief cmulti型の計算 y=sqrt(x)
*/
void csqrt_r(cmulti *y, rmulti *x)
{
  if(ris_zero(x) || ris_positive(x)){
    rsqrt(C_R(y),x);
    rset_zero(C_I(y));
  }else if(ris_negative(x)){
    rset_zero(C_R(y));
    rneg(C_I(y),x);
    rsqrt(C_I(y),C_I(y));
  }else{ ERROR_AT; }
}

/**
 @brief cmulti型の計算 y=sqrt(x)
*/
void csqrt_c(cmulti *y, cmulti *x)
{
  cpow_d2(y,x,0.5);
}

/**
 @brief cmulti型の計算 y=exp(x)
*/
void cexp_c(cmulti *y, cmulti *x)
{
  rmulti *r=NULL,*theta=NULL;
  RAc(r,y); RAc(theta,y);
  rexp(r,C_R(x));
  rcopy(theta,C_I(x));
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
  ccopy_r(a,x);
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
  rlog(C_R(y),r);        // y=log(r)+i*theta
  rcopy(C_I(y),theta);   // y=log(r)+i*theta
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
  rsin(s,C_R(x));
  rcos(c,C_R(x));
  rsinh(sh,C_I(x));
  rcosh(ch,C_I(x));
  rmul(C_R(y),s,ch);
  rmul(C_I(y),c,sh);
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
  rsin(s,C_R(x));
  rcos(c,C_R(x));
  rsinh(sh,C_I(x));
  rcosh(ch,C_I(x));
  rmul(C_R(y),c,ch);
  rmul(C_I(y),s,sh);
  rneg(C_I(y),C_I(y));
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
  cdiv(y,s,c); // y=sin(x)/cos(x)
  CF(c); CF(s);
}

/**
 @brief cmulti型の計算 y=arcsin(x)
*/
void casin_r(cmulti *y, rmulti *x)
{
  cmulti *a=NULL;
  CAc(a,y);
  ccopy_r(a,x);
  casin_c(y,a);
  CF(a);
}

/**
 @brief cmulti型の計算 y=arcsin(x)
*/
void casin_c(cmulti *y, cmulti *x)
{
  cmulti *a=NULL,*b=NULL;
  if(cis_real(x) && rle_si2(C_R(x),1) && rge_si2(C_R(x),-1)){
    rasin(C_R(y),C_R(x));
    rset_zero(C_I(y));
  }else{
    CAc(a,y); CAc(b,y);
    cmul(a,x,x);     // a=x^2
    csub_si1(a,1,a); // a=1-x^2
    csqrt_c(a,a);    // a=sqrt(1-x^2)
    cset_dd(b,0,1);  // b=I
    cmul(b,b,x);     // b=I*x
    cadd(a,a,b);     // a=I*x+sqrt(1-x^2)
    clog_c(a,a);     // a=log(I*x+sqrt(1-x^2))
    cset_dd(b,0,-1); // b=-I
    cmul(a,a,b);     // a=-I*log(I*x+sqrt(1-x^2))
    ccopy(y,a);      // y=-I*log(I*x+sqrt(1-x^2))
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
  ccopy_r(a,x);
  cacos_c(y,a);
  CF(a);
}

/**
 @brief cmulti型の計算 y=arccos(x)
*/
void cacos_c(cmulti *y, cmulti *x)
{
  cmulti *a=NULL,*b=NULL;
  if(cis_real(x) && rle_si2(C_R(x),1) && rge_si2(C_R(x),-1)){
    racos(C_R(y),C_R(x));
    rset_zero(C_I(y));
  }else{
    CAc(a,y); CAc(b,y);
    cmul(a,x,x);     // a=x^2
    csub_si2(a,a,1); // a=x^2-1
    csqrt_c(a,a);    // a=sqrt(x^2-1)
    cadd(a,a,x);     // a=x+sqrt(x^2-1)
    clog_c(a,a);     // a=log(x+sqrt(x^2-1))
    cset_dd(b,0,-1); // b=-I
    cmul(a,a,b);     // a=-I*log(I*x+sqrt(1-x^2))
    ccopy(y,a);      // y=-I*log(I*x+sqrt(1-x^2))
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
  ccopy_r(a,x);
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
    ratan(C_R(y),C_R(x));
    rset_si(C_I(y),0);
  }else{
    CAc(a,y); CAc(b,y);
    cadd_si(a,x,1);    // a=1+x
    csub_si1(b,1,x);   // b=1-x
    cdiv(a,a,b);       // a=(1+x)/(1-x)
    clog_c(a,a);       // a=log((1+x)/(1-x))
    cset_dd(b,0,-0.5); // b=-I/2
    cmul(a,a,b);       // a=(-I/2)*log((1+x)/(1-x))
    ccopy(y,a);        // y=(-I/2)*log((1+x)/(1-x))
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
  rsinh(sh,C_R(x));
  rcosh(ch,C_R(x));
  rsin(s,C_I(x));
  rcos(c,C_I(x));
  rmul(C_R(y),sh,c);
  rmul(C_I(y),ch,s);
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
  rsinh(sh,C_R(x));
  rcosh(ch,C_R(x));
  rsin(s,C_I(x));
  rcos(c,C_I(x));
  rmul(C_R(y),ch,c);
  rmul(C_I(y),sh,s);
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
  cdiv(y,sh,ch); // y=sinh(x)/cosh(x)
  CF(ch); CF(sh);
}

/**
 @brief cmulti型の計算 y=asinh(x)
*/
void casinh_r(cmulti *y, rmulti *x)
{
  cmulti *a=NULL;
  CAc(a,y);
  ccopy_r(a,x);
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
    rasinh(C_R(y),C_R(x));
    rset_zero(C_I(y));
  }else{
    CAc(a,y);
    cmul(a,x,x);
    cadd_si(a,a,1);
    csqrt_c(a,a);
    cadd(a,a,x);
    clog_c(a,a);
    ccopy(y,a);
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
  ccopy_r(a,x);
  cacosh_c(y,a);
  CF(a);
}

/**
 @brief cmulti型の計算 y=acosh(x)
*/
void cacosh_c(cmulti *y, cmulti *x)
{
  cmulti *a=NULL;  
  if(cis_real(x) && rge_si2(C_R(x),1)){
    racosh(C_R(y),C_R(x));
    rset_zero(C_I(y));
  }else{
    CAc(a,y);
    cmul(a,x,x);
    csub_si2(a,a,1);
    csqrt_c(a,a);
    cadd(a,a,x);
    clog_c(a,a);
    ccopy(y,a);
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
  ccopy_r(a,x);
  catanh_c(y,a);
  CF(a);
}

/**
 @brief cmulti型の計算 y=arctanh(x)
*/
void catanh_c(cmulti *y, cmulti *x)
{
  cmulti *a=NULL,*b=NULL;
  if(cis_real(x) && rle_si2(C_R(x),1) && rge_si2(C_R(x),-1)){
    ratanh(C_R(y),C_R(x));
    rset_zero(C_I(y));
  }else{
    CAc(a,y); CAc(b,y);
    cadd_si(a,x,1);  // a=1+x
    csub_si1(b,1,x); // b=1-x
    cdiv(a,a,b);     // a=(1+x)/(1-x)
    clog_c(a,a);     // a=log((1+x)/(1-x))
    cdiv_si2(a,a,2); // a=(1/2)*log((1+x)/(1-x))
    ccopy(y,a);      // y=(1/2)*log((1+x)/(1-x))
    CF(a); CF(b);
  }
}

/** @} */

/////////////////////////////////////////////////////////////////

/** @name cmulti型の２入力の数学演算子 */
/** @{ */


/**
 @brief cmulti型の足し算 z=x+y
*/
void cadd(cmulti *z, cmulti *x, cmulti *y)
{
  radd(C_R(z),C_R(x),C_R(y)); // z.r=x.r+y.r
  radd(C_I(z),C_I(x),C_I(y)); // z.i=x.i+y.i
}

/**
 @brief cmulti型の足し算 z=x+y
*/
void cadd_z(cmulti *z, cmulti *x, dcomplex y)
{
  radd_d(C_R(z),C_R(x),Z_R(y)); // z.r=x.r+y.r
  radd_d(C_I(z),C_I(x),Z_I(y)); // z.i=x.i+y.i
}

/**
 @brief cmulti型の足し算 z=x+y
*/
void cadd_r(cmulti *z, cmulti *x, rmulti *y)
{
  radd(C_R(z),C_R(x),y); // z.r=x.r+y
  rcopy(C_I(z),C_I(x));  // z.i=x.i
}

/**
 @brief cmulti型の足し算 z=x+y
*/
void cadd_d(cmulti *z, cmulti *x, double y)
{
  radd_d(C_R(z),C_R(x),y); // z.r=x.r+y
  rcopy(C_I(z),C_I(x));    // z.i=x.i
}

/**
 @brief cmulti型の足し算 z=x+y
*/
void cadd_ui(cmulti *z, cmulti *x, ulong y)
{
  radd_ui(C_R(z),C_R(x),y); // z.r=x.r+y
  rcopy(C_I(z),C_I(x));     // z.i=x.i
}

/**
 @brief cmulti型の足し算 z=x+y
*/
void cadd_si(cmulti *z, cmulti *x, long y)
{
  radd_si(C_R(z),C_R(x),y); // z.r=x.r+y
  rcopy(C_I(z),C_I(x));     // z.i=x.i
}

/**
 @brief cmulti型の引き算 z=x-y
*/
void csub(cmulti *z, cmulti *x, cmulti *y)
{
  rsub(C_R(z),C_R(x),C_R(y)); // z.r=x.r-y.r
  rsub(C_I(z),C_I(x),C_I(y)); // z.i=x.i-y.i
}

/**
 @brief cmulti型の引き算 z=x-y
*/
void csub_z1(cmulti *z, dcomplex x, cmulti *y)
{
  rsub_d1(C_R(z),Z_R(x),C_R(y)); // z.r=x.r-y.r
  rsub_d1(C_I(z),Z_I(x),C_I(y)); // z.i=x.i-y.i
}

/**
 @brief cmulti型の引き算 z=x-y
*/
void csub_z2(cmulti *z, cmulti *x, dcomplex y)
{
  rsub_d2(C_R(z),C_R(x),Z_R(y)); // z.r=x.r-y.r
  rsub_d2(C_I(z),C_I(x),Z_I(y)); // z.i=x.i-y.i
}

/**
 @brief cmulti型の引き算 z=x-y
*/
void csub_r1(cmulti *z, rmulti *x, cmulti *y)
{
  rsub(C_R(z),x,C_R(y)); // z.r=x-y.r
  rneg(C_I(z),C_I(y));   // z.i= -y.i
}

/**
 @brief cmulti型の引き算 z=x-y
*/
void csub_r2(cmulti *z, cmulti *x, rmulti *y)
{
  rsub(C_R(z),C_R(x),y); // z.r=x.r-y
  rcopy(C_I(z),C_I(x));  // z.i=x.i
}

/**
 @brief cmulti型の引き算 z=x-y
*/
void csub_d1(cmulti *z, double x, cmulti *y)
{
  rsub_d1(C_R(z),x,C_R(y)); // z.r=x-y.r
  rneg(C_I(z),C_I(y));      // z.i= -y.i
}

/**
 @brief cmulti型の引き算 z=x-y
*/
void csub_d2(cmulti *z, cmulti *x, double y)
{
  rsub_d2(C_R(z),C_R(x),y); // z.r=x.r-y
  rcopy(C_I(z),C_I(x));     // z.i=x.i
}

/**
 @brief cmulti型の引き算 z=x-y
*/
void csub_ui1(cmulti *z, ulong x, cmulti *y)
{
  rsub_ui1(C_R(z),x,C_R(y)); // z.r=x-y.r
  rneg(C_I(z),C_I(y));      // z.i= -y.i
}

/**
 @brief cmulti型の引き算 z=x-y
*/
void csub_ui2(cmulti *z, cmulti *x, ulong y)
{
  rsub_ui2(C_R(z),C_R(x),y); // z.r=x.r-y
  rcopy(C_I(z),C_I(x));     // z.i=x.i
}

/**
 @brief cmulti型の引き算 z=x-y
*/
void csub_si1(cmulti *z, long x, cmulti *y)
{
  rsub_si1(C_R(z),x,C_R(y)); // z.r=x-y.r
  rneg(C_I(z),C_I(y));       // z.i= -y.i
}

/**
 @brief cmulti型の引き算 z=x-y
*/
void csub_si2(cmulti *z, cmulti *x, long y)
{
  rsub_si2(C_R(z),C_R(x),y); // z.r=x.r-y
  rcopy(C_I(z),C_I(x));     // z.i=x.i
}

/**
 @brief cmulti型の掛け算 z=x*y
*/
void cmul(cmulti *z, cmulti *x, cmulti *y)
{
  cmulti *a=NULL;
  CAc(a,z);
      rmul(C_R(a),C_R(x),C_R(y)); // z.r =x.r*y.r
  rsub_mul(C_R(a),C_I(x),C_I(y)); // z.r-=x.i*y.i
      rmul(C_I(a),C_I(x),C_R(y)); // z.i =x.i*y.r
  radd_mul(C_I(a),C_R(x),C_I(y)); // z.i+=x.r*y.i
  ccopy(z,a);                     // z=a
  CF(a);
}

/**
 @brief cmulti型の掛け算 z=x*y
*/
void cmul_ws(cmulti *z, cmulti *x, cmulti *y, int n_rws, rmulti **rws, int n_cws, cmulti **cws)
{
  if(n_rws<1){ ERROR_EXIT("Error `n_rws=%d<1' in cmul_ws().\n",n_rws); }
  if(n_cws<1){ ERROR_EXIT("Error `n_cws=%d<1' in cmul_ws().\n",n_cws); }
  rmul(C_R(cws[0]),C_R(x),C_R(y));                  // z.r =x.r*y.r
  rsub_mul_ws(C_R(cws[0]),C_I(x),C_I(y),n_rws,rws); // z.r-=x.i*y.i
  rmul(C_I(cws[0]),C_I(x),C_R(y));                  // z.i =x.i*y.r
  radd_mul_ws(C_I(cws[0]),C_R(x),C_I(y),n_rws,rws); // z.i+=x.r*y.i
  ccopy(z,cws[0]);                                  // z=cws[0]
}

/**
 @brief cmulti型の掛け算 z=x*y
*/
void cmul_z(cmulti *z, cmulti *x, dcomplex y)
{
  cmulti *a=NULL;
  CAc(a,z);
  cset_z(a,y); // a=y
  cmul(z,x,a); // z=x*y
  CF(a);
}

/**
 @brief cmulti型の掛け算 z=x*y
*/
void cmul_r(cmulti *z, cmulti *x, rmulti *y)
{
  rmul(C_R(z),C_R(x),y); // z.r=x.r*y
  rmul(C_I(z),C_I(x),y); // z.i=x.i*y
}

/**
 @brief cmulti型の掛け算 z=x*y
*/
void cmul_d(cmulti *z, cmulti *x, double y){
  rmul_d(C_R(z),C_R(x),y); // z.r=x.r*y
  rmul_d(C_I(z),C_I(x),y); // z.i=x.i*y
}

/**
 @brief cmulti型の掛け算 z=x*y
*/
void cmul_ui(cmulti *z, cmulti *x, ulong y){
  rmul_ui(C_R(z),C_R(x),y); // z.r=x.r*y
  rmul_ui(C_I(z),C_I(x),y); // z.i=x.i*y
}

/**
 @brief cmulti型の掛け算 z=x*y
*/
void cmul_si(cmulti *z, cmulti *x, long y){
  rmul_si(C_R(z),C_R(x),y); // z.r=x.r*y
  rmul_si(C_I(z),C_I(x),y); // z.i=x.i*y
}

/**
 @brief cmulti型の複素共役との掛け算 z=conj(x)*y
*/
void cdot(cmulti *z, cmulti *x, cmulti *y)
{
  cmulti *a=NULL;
  CAc(a,z);
      rmul(C_R(a),C_R(x),C_R(y)); // z.r =x.r*y.r
  radd_mul(C_R(a),C_I(x),C_I(y)); // z.r+=x.i*y.i
      rmul(C_I(a),C_R(x),C_I(y)); // z.i =x.r*y.i
  rsub_mul(C_I(a),C_I(x),C_R(y)); // z.i-=x.i*y.r
  ccopy(z,a);                     // z=a
  CF(a);
}

/**
 @brief cmulti型の複素共役との掛け算 z=conj(x)*y
*/
void cdot_z1(cmulti *z, dcomplex x, cmulti *y)
{
  cmulti *a=NULL;
  CAc(a,z);
  cset_z(a,x); // a=x
  cdot(z,a,y); // z=conj(x)*y  
  CF(a);
}

/**
 @brief cmulti型の複素共役との掛け算 z=conj(x)*y
*/
void cdot_z2(cmulti *z, cmulti *x, dcomplex y)
{
  cmulti *a=NULL;
  CAc(a,z);
  cset_z(a,y); // a=y
  cdot(z,x,a); // z=conj(x)*y 
  CF(a);
}


/**
 @brief cmulti型の掛け算の加算 z+=x*y
*/
void cadd_mul(cmulti *z, cmulti *x, cmulti *y)
{
  cmulti *a=NULL;
  CAc(a,z);
  cmul(a,x,y);   // a=x*y
  cadd(z,z,a);   // z=z+a
  CF(a);
}

/**
 @brief cmulti型の掛け算の加算 z+=x*y
*/
void cadd_mul_r(cmulti *z, cmulti *x, rmulti *y)
{
  cmulti *a=NULL;
  CAc(a,z);
  cmul_r(a,x,y); // a=x*y
  cadd(z,z,a);   // z=z+a
  CF(a);
}

/**
 @brief cmulti型の掛け算の加算 z+=x*y
*/
void cadd_mul_z(cmulti *z, cmulti *x, dcomplex y)
{
  cmulti *a=NULL;
  CAc(a,z);
  cmul_z(a,x,y); // a=x*y
  cadd(z,z,a);   // z=z+a
  CF(a);
}

/**
 @brief cmulti型の掛け算の加算 z+=x*y
*/
void cadd_mul_d(cmulti *z, cmulti *x, double y)
{
  cmulti *a=NULL;
  CAc(a,z);
  cmul_d(a,x,y); // a=x*y
  cadd(z,z,a);   // z=z+a
  CF(a);
}

/**
 @brief cmulti型の掛け算の減算 z-=x*y
*/
void csub_mul(cmulti *z, cmulti *x, cmulti *y)
{
  cmulti *a=NULL;
  CAc(a,z);
  cmul(a,x,y);   // a=x*y
  csub(z,z,a);   // z=z-a
  CF(a);
}

/**
 @brief cmulti型の掛け算の減算 z-=x*y
*/
void csub_mul_ws(cmulti *z, cmulti *x, cmulti *y, int n_rws, rmulti **rws, int n_cws, cmulti **cws)
{
  if(n_rws<1){ ERROR_EXIT("Error `n_rws=%d<1' in csub_mul_ws().\n",n_rws); }
  if(n_cws<2){ ERROR_EXIT("Error `n_cws=%d<2' in csub_mul_ws().\n",n_cws); }
  cmul_ws(cws[0],x,y,n_rws,rws,n_cws-1,cws+1); // cws[0]=x*y
  csub(z,z,cws[0]);                            // z=z-cws[0]
}

/**
 @brief cmulti型の掛け算の減算 z-=x*y
*/
void csub_mul_r(cmulti *z, cmulti *x, rmulti *y)
{
  cmulti *a=NULL;
  CAc(a,z);
  cmul_r(a,x,y); // a=x*y
  csub(z,z,a);   // z=z-a
  CF(a);
}

/**
 @brief cmulti型の掛け算の減算 z-=x*y
*/
void csub_mul_z(cmulti *z, cmulti *x, dcomplex y)
{
  cmulti *a=NULL;
  CAc(a,z);
  cmul_z(a,x,y); // a=x*y
  csub(z,z,a);   // z=z-a
  CF(a);
}

/**
 @brief cmulti型の掛け算の減算 z-=x*y
*/
void csub_mul_d(cmulti *z, cmulti *x, double y)
{
  cmulti *a=NULL;
  CAc(a,z);
  cmul_d(a,x,y); // a=x*y
  csub(z,z,a);   // z=z-a
  CF(a);
}

/**
 @brief cmulti型の複素共役との掛け算の加算 z+=conj(x)*y
*/
void cadd_dot(cmulti *z, cmulti *x, cmulti *y)
{
  cmulti *a=NULL;
  CAc(a,z);
  cdot(a,x,y);   // a=conj(x)*y
  cadd(z,z,a);   // z=z+a
  CF(a);
}

/**
 @brief cmulti型の複素共役との掛け算の加算 z+=conj(x)*y
*/
void cadd_dot_z1(cmulti *z, dcomplex x, cmulti *y)
{
  cmulti *a=NULL;
  CAc(a,z);
  cdot_z1(a,x,y); // a=conj(x)*y
  cadd(z,z,a);    // z=z+a
  CF(a);
}

/**
 @brief cmulti型の複素共役との掛け算の加算 z+=conj(x)*y
*/
void cadd_dot_z2(cmulti *z, cmulti *x, dcomplex y)
{
  cmulti *a=NULL;
  CAc(a,z);
  cdot_z2(a,x,y); // a=conj(x)*y
  cadd(z,z,a);    // z=z+a
  CF(a);
}


/**
 @brief cmulti型の複素共役との掛け算の減算 z-=conj(x)*y
*/
void csub_dot(cmulti *z, cmulti *x, cmulti *y)
{
  cmulti *a=NULL;
  CAc(a,z);
  cdot(a,x,y);   // a=conj(x)*y
  csub(z,z,a);   // z=z-a
  CF(a);
}

/**
 @brief cmulti型の複素共役との掛け算の減算 z-=conj(x)*y
*/
void csub_dot_z1(cmulti *z, dcomplex x, cmulti *y)
{
  cmulti *a=NULL;
  CAc(a,z);
  cdot_z1(a,x,y);   // a=conj(x)*y
  csub(z,z,a);   // z=z-a
  CF(a);
}

/**
 @brief cmulti型の複素共役との掛け算の減算 z-=conj(x)*y
*/
void csub_dot_z2(cmulti *z, cmulti *x, dcomplex y)
{
  cmulti *a=NULL;
  CAc(a,z);
  cdot_z2(a,x,y);   // a=conj(x)*y
  csub(z,z,a);   // z=z-a
  CF(a);
}

/**
 @brief cmulti型の割り算 z=x/y
*/
void cdiv(cmulti *z, cmulti *x, cmulti *y)
{
  rmulti *den=NULL;
  cmulti *a=NULL;
  RAc(den,z); CAc(a,z);
  cabs2(den,y);     // den=|y|^2
  cdot(a,y,x);      // z=conj(y)*x
  cdiv_r2(a,a,den); // z/=den
  ccopy(z,a);       // z=a
  RF(den); CF(a);
}

/**
 @brief cmulti型の割り算 z=x/y
*/
void cdiv_z1(cmulti *z, dcomplex x, cmulti *y)
{
  cmulti *a=NULL;
  CAc(a,z);
  cset_z(a,x);   // a=x
  cdiv(z,a,y);   // z=x/y
  CF(a);
}

/**
 @brief cmulti型の割り算 z=x/y
*/
void cdiv_z2(cmulti *z, cmulti *x, dcomplex y)
{
  cmulti *a=NULL;
  CAc(a,z);
  cset_z(a,y);   // a=y
  cdiv(z,x,a);   // z=x/y
  CF(a);
}

/**
 @brief cmulti型の割り算 z=x/y
*/
void cdiv_r1(cmulti *z, rmulti *x, cmulti *y)
{
  rmulti *den=NULL;
  cmulti *a=NULL;
  RAc(den,z); CAc(a,z);
  cabs2(den,y);       // den=|y|^2
  cconj(a,y);         // z=conj(y)
  cmul_r(a,a,x);      // z*=x
  cdiv_r2(a,a,den);   // z/=den
  ccopy(z,a);         // z=a
  RF(den); CF(a);
}

/**
 @brief cmulti型の割り算 z=x/y
*/
void cdiv_r2(cmulti *z, cmulti *x, rmulti *y)
{
  rdiv(C_R(z),C_R(x),y); // z.r=x.r/y
  rdiv(C_I(z),C_I(x),y); // z.i=x.i/y
}

/**
 @brief cmulti型の割り算 z=x/y
*/
void cdiv_d1(cmulti *z, double x, cmulti *y)
{
  rmulti *a=NULL;
  RAc(a,z);
  rset_d(a,x);    // a=x
  cdiv_r1(z,a,y); // z=x/y
  RF(a);
}

/**
 @brief cmulti型の割り算 z=x/y
*/
void cdiv_d2(cmulti *z, cmulti *x, double y)
{
  rdiv_d2(C_R(z),C_R(x),y); // z.r=x.r/y
  rdiv_d2(C_I(z),C_I(x),y); // z.i=x.i/y
}

/**
 @brief cmulti型の割り算 z=x/y
*/
void cdiv_ui1(cmulti *z, ulong x, cmulti *y)
{
  rmulti *a=NULL;
  RAc(a,z);
  rset_ui(a,x);   // a=x
  cdiv_r1(z,a,y); // z=x/y
  RF(a);
}

/**
 @brief cmulti型の割り算 z=x/y
*/
void cdiv_ui2(cmulti *z, cmulti *x, ulong y)
{
  rdiv_ui2(C_R(z),C_R(x),y); // z.r=x.r/y
  rdiv_ui2(C_I(z),C_I(x),y); // z.i=x.i/y
}

/**
 @brief cmulti型の割り算 z=x/y
*/
void cdiv_si1(cmulti *z, long x, cmulti *y)
{
  rmulti *a=NULL;
  RAc(a,z);
  rset_si(a,x);   // a=x
  cdiv_r1(z,a,y); // z=x/y
  RF(a);
}

/**
 @brief cmulti型の割り算 z=x/y
*/
void cdiv_si2(cmulti *z, cmulti *x, long y)
{
  rdiv_si2(C_R(z),C_R(x),y); // z.r=x.r/y
  rdiv_si2(C_I(z),C_I(x),y); // z.i=x.i/y
}

/**
 @brief cmulti型の差の絶対値 z=abs(x-y)
*/
void cabs_sub(rmulti *z, cmulti *x, cmulti *y)
{
  cmulti *a;
  CAr(a,z);
  csub(a,x,y); // a=x-y
  cabsv(z,a);  // z=abs(x-y)
  CF(a);
}

/**
 @brief cmulti型の差の絶対値 z=abs(x-y)
*/
void cabs_sub_r(rmulti *z, cmulti *x, rmulti *y)
{
  cmulti *a;
  CAr(a,z);
  csub_r2(a,x,y); // a=x-y
  cabsv(z,a);     // z=abs(x-y)
  CF(a);
}

/**
 @brief cmulti型の絶対値の割り算 z=x/abs(y)
*/
void cdiv_abs(rmulti *z, rmulti *x, cmulti *y)
{
  rmulti *a;
  RAr(a,z);
  cabsv(a,y);  // z=abs(y)  
  rdiv(z,x,a); // z=x/abs(y)
  RF(a);
}

/**
 @brief cmulti型のべき乗 z=x^y
*/
void cpow_c(cmulti *z, cmulti *x, cmulti *y)
{
  cmulti *a=NULL;
  CAc(a,z);
  clog_c(a,x); // a=log(x)
  cmul(a,y,a); // a=y*log(x)
  cexp_c(z,a); // z=exp(y*log(x))
  CF(a);
}

/**
 @brief cmulti型のべき乗 z=x^y
*/
void cpow_r1(cmulti *z, rmulti *x, cmulti *y)
{
  rmulti *a=NULL;
  cmulti *b=NULL;
  RAc(a,z); CAc(b,z);
  rlog(a,x);     // a=log(x)
  cmul_r(b,y,a); // b=y*log(x)
  cexp_c(z,b);   // z=exp(y*log(x))
  RF(a); CF(b);
}

/**
 @brief cmulti型のべき乗 z=x^y
*/
void cpow_r2(cmulti *z, cmulti *x, rmulti *y)
{
  rmulti *r=NULL,*theta=NULL;
  RAc(r,z); RAc(theta,z);
  cget_polar(r,theta,x); // x=r*exp(i*theta)
  rpow(r,r,y);           // r=r^y
  rmul(theta,theta,y);   // theta=theta*y
  cset_polar(z,r,theta); // z=r^y*cos(theta*y)+i*r^y*sin(theta*y)
  RF(r); RF(theta);
}

/**
 @brief cmulti型のべき乗 z=x^y
*/
void cpow_d2(cmulti *z, cmulti *x, double y)
{
  rmulti *r=NULL,*theta=NULL;
  RAc(r,z); RAc(theta,z);
  cget_polar(r,theta,x); // x=r*exp(i*theta)
  rpow_d2(r,r,y);        // r=r^y
  rmul_d(theta,theta,y); // theta=theta*y
  cset_polar(z,r,theta); // z=r^y*cos(theta*y)+i*r^y*sin(theta*y)
  RF(r); RF(theta);
}

/**
 @brief cmulti型のべき乗 y=x^n
*/
void cpow_si(cmulti *y, cmulti *x, long n)
{
  cmulti *a=NULL;
  long i;
  if(n==0){
    cset_one(y);
  }else if(cis_real(x)){
    rpow_si(C_R(y),C_R(x),n);
    rset_zero(C_I(y));
  }else if(n>0){
    CAc(a,y);
    ccopy(a,x);                        // a=x
    for(i=1; i<n; i++){ cmul(a,a,x); } // a*=x
    ccopy(y,a);                        // y=a
    CF(a);
  }else if(n<0){
    CAc(a,y);
    ccopy(a,x);                           // a=x
    for(i=1; i<(-n); i++){ cmul(a,a,x); } // a*=x
    cinv(a,a);                            // a=1/a
    ccopy(y,a);                           // y=a
    CF(a);
  }else{ ERROR_AT; }
}

/**
 @brief rmulti型のべき乗 y=x^n
*/
void cpow_ui(cmulti *y, cmulti *x, ulong n)
{
  cmulti *a=NULL;
  ulong i;
  if(cis_real(x)){
    rpow_ui(C_R(y),C_R(x),n);
    rset_si(C_I(y),0);
  }else if(n==0){
    cset_one(y);
  }else{
    CAc(a,y);
    ccopy(a,x);                        // a=x
    for(i=1; i<n; i++){ cmul(a,a,x); } // a*=x
    ccopy(y,a);                        // y=a
    CF(a);
  }
}

/**
 @brief cmulti型の指数部の足し算 y=x*2^n
*/
void cmul_2exp(cmulti *y, cmulti *x, int nr, int ni)
{
  rmul_2exp(C_R(y),C_R(x),nr);
  rmul_2exp(C_I(y),C_I(x),ni);
}

/**
 @brief cmulti型の指数部の引き算 y=x/2^n
*/
void cdiv_2exp(cmulti *y, cmulti *x, int nr, int ni)
{
  rdiv_2exp(C_R(y),C_R(x),nr);
  rdiv_2exp(C_I(y),C_I(x),ni);
}

/**
 @brief cmulti型の指数部で評価 z=10^(floor(log10(abs(x))-y))
*/
void cexp10_floor_log10_abs_sub(cmulti *z, cmulti *x, double y)
{
  rmulti *exp1=NULL,*exp2=NULL;
  RAc(exp1,x); RAc(exp2,x);
  rabs(C_R(x),C_R(x));
  rabs(C_I(x),C_I(x));
  rlog10(C_R(x),C_R(x));
  rlog10(C_I(x),C_I(x));
  rsub_d2(C_R(x),C_R(x),y);
  rsub_d2(C_I(x),C_I(x),y);
  rfloor(exp1,C_R(x));
  rfloor(exp2,C_I(x));
  rexp10(C_R(z),exp1);
  rexp10(C_I(z),exp2);
  RF(exp1);RF(exp2);
}

/**
 @brief cmulti型の指数部で評価 z=2^(floor(log2(abs(x))-y))
*/
void cexp2_floor_log2_abs_sub(cmulti *z, cmulti *x, double y)
{
  rmulti *exp1=NULL,*exp2=NULL;
  RAc(exp1,x); RAc(exp2,x);
  rabs(C_R(x),C_R(x));
  rabs(C_I(x),C_I(x));
  rlog2(C_R(x),C_R(x));
  rlog2(C_I(x),C_I(x));
  rsub_d2(C_R(x),C_R(x),y);
  rsub_d2(C_I(x),C_I(x),y);
  rfloor(exp1,C_R(x));
  rfloor(exp2,C_I(x));
  rexp2(C_R(z),exp1);
  rexp2(C_I(z),exp2);
  RF(exp1);RF(exp2);
}

/**
 @brief cmulti型の極座標 z=r*exp(i*theta)
*/
void cget_polar(rmulti *r, rmulti *theta, cmulti *z)
{
  cabsv(r,z);
  ratan2(theta,C_I(z),C_R(z));
}

/**
 @brief cmulti型の極座標から実部，虚部へ変換 z=r*exp(i*theta)
*/
void cset_polar(cmulti *z, rmulti *r, rmulti *theta)
{
  rmulti *s=NULL,*c=NULL;
  RAc(s,z); RAc(c,z);
  rcos(c,theta); rmul(C_R(z),r,c); // z.r=r*cos(theta)
  rsin(s,theta); rmul(C_I(z),r,s); // z.i=r*sin(theta)
  RF(s); RF(c);
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
