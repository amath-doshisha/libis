#include"is_rmulti.h"
#include"is_rmat.h"
#include"is_rmulti.h"
#include"is_cmulti.h"
#include"is_irmulti.h"
#include"is_icmulti.h"
#include"is_func.h"

/**
 @file  irmulti.c
 @brief 多倍長精度実数型rmultiの機械区間演算に関する関数の定義.
 @details ベクトルに関しては@link irvec.c@endlinkを参照のこと.
          行列に関しては@link irmat.c@endlinkを参照のこと.
 */


#define RA(X,P) ((X)=rallocate_prec(P))
#define CA(X,P) ((X)=callocate_prec(P))
#define RF(X) ((X)=rfree(X))
#define CF(X) ((X)=cfree(X))
#define RAr(X,P) ((X)=rallocate_prec(rget_prec(P)))
#define CAr(X,P) ((X)=callocate_prec(rget_prec(P)))


/** @name irmulti型の値の設定に関する関数 */
/** @{ */

/**
 @brief コピー [y0,y1]=[x0,x1].
 */
void ircopy(rmulti *y0, rmulti *y1, rmulti *x0, rmulti *x1)
{
  mpfr_set(y0,x0,MPFR_RNDD); // lower bound
  mpfr_set(y1,x1,MPFR_RNDU); // upper bound
}

/**
 @brief コピー [y0,y1]=[x0,x1].
 */
void irset(rmulti *y0, rmulti *y1, rmulti *x0, rmulti *x1)
{
  mpfr_set(y0,x0,MPFR_RNDD); // lower bound
  mpfr_set(y1,x1,MPFR_RNDU); // upper bound
}

/**
 @brief 倍精度実数の設定 [y0,y1]=[x,x].
 */
void irset_dr(rmulti *y0, rmulti *y1, double x0, rmulti *x1)
{
  mpfr_set_d(y0,x0,MPFR_RNDD); // lower bound
  mpfr_set  (y1,x1,MPFR_RNDU); // upper bound
}

/**
 @brief 倍精度実数の設定 [y0,y1]=[x,x].
 */
void irset_rd(rmulti *y0, rmulti *y1, rmulti *x0, double x1)
{
  mpfr_set  (y0,x0,MPFR_RNDD); // lower bound
  mpfr_set_d(y1,x1,MPFR_RNDU); // uppper bound
}

/**
 @brief 倍精度実数の設定 [y0,y1]=[x0,x1].
 */
void irset_d(rmulti *y0, rmulti *y1, double x0, double x1)
{
  mpfr_set_d(y0,x0,MPFR_RNDD); // lower bound
  mpfr_set_d(y1,x1,MPFR_RNDU); // upper bound
}


void irset_s(rmulti *x0, rmulti *x1, char *s)
{
  int p0,p1,prec;
  cmulti *a0=NULL,*a1=NULL;
  p0=rget_prec(x0); p1=rget_prec(x1); prec=MAX2(p0,p1);
  a0=callocate_prec(prec);
  a1=callocate_prec(prec);  
  icset_s(a0,a1,s);
  ircopy(x0,x1,C_R(a0),C_R(a1));
}

/** @} */
/** @name irmulti型の型変換に関する関数 */
/** @{ */

/**
 @brief bigint型から[z0,z1]へ型変換.
 */
void irset_bigint(rmulti *z0, rmulti *z1, bigint *x)
{
  mpfr_div(z0,BIGINT_NUM(x),BIGINT_DEN(x),MPFR_RNDD); // lower bound
  mpfr_div(z1,BIGINT_NUM(x),BIGINT_DEN(x),MPFR_RNDU); // upper bound
}

/** @} */
/** @name irmulti型の１入力の数学演算子 */
/** @{ */

/**
 @brief 符号の反転 [y0,y1]=-[x0,x1].
 */
void irneg(rmulti *y0, rmulti *y1, rmulti *x0, rmulti *x1)
{
  int p0,p1,prec;
  rmulti *a0=NULL,*a1=NULL;
  p0=rget_prec(y0); p1=rget_prec(y1); prec=MAX2(p0,p1);
  RA(a0,prec); RA(a1,prec);
  mpfr_neg(a0,x1,MPFR_RNDD); // lower bound
  mpfr_neg(a1,x0,MPFR_RNDU); // upper bound
  ircopy(y0,y1,a0,a1);
  RF(a0); RF(a1);
}

/**
 @brief 符号の正負 [y0,y1]=[-abs(x),abs(x)].
 */
void irpm(rmulti *y0, rmulti *y1, rmulti *x)
{
  mpfr_abs(y1,x, MPFR_RNDU); // uppper bound
  mpfr_neg(y0,y1,MPFR_RNDD); // lower bound
}

/**
 @brief 区間の中心 [m-r,m+r]=[x0,x1]
 */
void irmid(rmulti *mid, rmulti *x0, rmulti *x1)
{
  int prec;
  rmulti *a0=NULL,*a1=NULL;
  prec=rget_prec(mid);
  RA(a0,prec); RA(a1,prec);
  mpfr_sub(a0,x1,x0,MPFR_RNDU);    // a0=(x1-x0) by upper
  mpfr_mul_d(a0,a0,0.5,MPFR_RNDU); // a0=(x1-x0)/2 by upper
  mpfr_add(mid,a0,x0,MPFR_RNDU);   // mid=x0+(x1-x0)/2 by upper
  RF(a0); RF(a1);
}

/**
 @brief 区間の半径 rad=rad([x0,x1])
 */
void irrad(rmulti *rad, rmulti *x0, rmulti *x1)
{
  int prec;
  rmulti *a0=NULL;
  prec=rget_prec(rad);
  RA(a0,prec);
  mpfr_sub(a0,x1,x0,MPFR_RNDU);    // a0=(x1-x0) by upper
  mpfr_mul_d(a0,a0,0.5,MPFR_RNDU); // a0=(x1-x0)/2 by upper
  mpfr_add(a0,a0,x0,MPFR_RNDU);    // a0=x0+(x1-x0)/2 by upper
  mpfr_sub(rad,a0,x0,MPFR_RNDU);   // rad=(x1-x0)/2 by upper
  RF(a0); 
}

/**
 @brief 区間の中心と半径の表示(center-radius form) [m-r,m+r]=[x0,x1]
 */
void irmr(rmulti *mid, rmulti *rad, rmulti *x0, rmulti *x1)
{
  int p0,p1,prec;
  rmulti *a0=NULL,*a1=NULL;
  p0=rget_prec(mid); p1=rget_prec(rad); prec=MAX2(p0,p1);
  RA(a0,prec); RA(a1,prec);
  mpfr_sub(a0,x1,x0,MPFR_RNDU);    // a0=(x1-x0) by upper
  mpfr_mul_d(a0,a0,0.5,MPFR_RNDU); // a0=(x1-x0)/2 by upper
  mpfr_add(a0,a0,x0,MPFR_RNDU);    // a0=x0+(x1-x0)/2 by upper
  mpfr_sub(a1,a0,x0,MPFR_RNDU);    // a1=(x1-x0)/2 by upper
  ircopy(mid,rad,a0,a1);           // [mid,rad]=[a0,a1]
  RF(a0); RF(a1);
}

/**
 @brief 逆数 [z0,z1]=[1,1]/[x0,x1]
*/
void irinv(rmulti *z0, rmulti *z1, rmulti *x0, rmulti *x1)
{
  if(rget_sgn(x0)<0 && rget_sgn(x1)>0){    
    rset_nan(z0); rset_nan(z1); // return nan if divided by 0
  }else{
    mpfr_d_div(z0,1,x1,MPFR_RNDD); // lower bound
    mpfr_d_div(z1,1,x0,MPFR_RNDU); // upper bound
  } 
}

/**
 @brief 絶対値 [y0,y1]=abs([x0,x1])
*/
void irabs(rmulti *y0, rmulti *y1, rmulti *x0, rmulti *x1)
{
  int p0,p1,prec;
  rmulti *a0=NULL,*a1=NULL;
  if(ris_nan(x0) || ris_nan(x1)){ rset_nan(y0); rset_nan(y1); return; }
  // x0>=0 && x1>=0
  if(rget_sgn(x0)>=0 && rget_sgn(x1)>=0){ ircopy(y0,y1,x0,x1); }
  // x0<=0 && x1<=0
  else if(rget_sgn(x0)<=0 && rget_sgn(x1)<=0){ irneg(y0,y1,x0,x1); }
  // x0<0 && x1>0
  else{
    p0=rget_prec(y0); p1=rget_prec(y1); prec=MAX2(p0,p1);
    RA(a0,prec); RA(a1,prec);
    mpfr_neg(a0,x0,MPFR_RNDU);
    rmax2_up(a1,a0,x1);
    mpfr_set_d(y0,0,MPFR_RNDD);
    mpfr_set(y1,a1,MPFR_RNDU);
    RF(a0); RF(a1);
  }
}

/**
 @brief 平方根 [y0,y1]=sqrt([x0,x1])
*/
void irsqrt(rmulti *y0, rmulti *y1, rmulti *x0, rmulti *x1)
{
  mpfr_sqrt(y0,x0,MPFR_RNDD); // lower bound
  mpfr_sqrt(y1,x1,MPFR_RNDU); // upper bound
}


/**
 @brief 指数関数 [y0,y1]=exp([x0,x1])
*/
void irexp(rmulti *y0, rmulti *y1, rmulti *x0, rmulti *x1)
{
  mpfr_exp(y0,x0,MPFR_RNDD); // lower bound
  mpfr_exp(y1,x1,MPFR_RNDU); // upper bound
}

/**
 @brief 対数関数 [y0,y1]=log([x0,x1])
*/
void irlog(rmulti *y0, rmulti *y1, rmulti *x0, rmulti *x1)
{
  mpfr_log(y0,x0,MPFR_RNDD); // lower bound
  mpfr_log(y1,x1,MPFR_RNDU); // upper bound
}

/**
 @brief 対数関数 [y0,y1]=log([x0,x1])
*/
void irlog10(rmulti *y0, rmulti *y1, rmulti *x0, rmulti *x1)
{
  mpfr_log10(y0,x0,MPFR_RNDD); // lower bound
  mpfr_log10(y1,x1,MPFR_RNDU); // upper bound
}

/**
 @brief 三角関数 [y0,y1]=sin([x0,x1])
 */
void irsin(rmulti *y0, rmulti *y1, rmulti *x0, rmulti *x1)
{
  int prec;
  rmulti *p0=NULL,*p1=NULL,*a=NULL,*b0=NULL,*b1=NULL;
  prec=MAX2(rget_prec(y0),rget_prec(y1));
  RA(p0,prec); RA(p1,prec); RA(a,prec); RA(b0,prec); RA(b1,prec);
  mpfr_cos(p0,x0,MPFR_RNDN);
  mpfr_cos(p1,x1,MPFR_RNDN);
  if(rget_sgn(p0)>=0 && rget_sgn(p1)>=0){
    mpfr_sin(y0,x0,MPFR_RNDD);
    mpfr_sin(y1,x1,MPFR_RNDU);
  }else if(rget_sgn(p0)<=0 && rget_sgn(p1)<=0){
    mpfr_sin(y0,x1,MPFR_RNDD);
    mpfr_sin(y1,x0,MPFR_RNDU);
  }else if(rget_sgn(p0)>=0 && rget_sgn(p1)<=0){
    mpfr_sin(b0,x0,MPFR_RNDD);
    mpfr_sin(b1,x1,MPFR_RNDD);
    mpfr_set_d(a,1,MPFR_RNDU);
    if(rle(b0,b1)){ ircopy(y0,y1,b0,a); }else{ ircopy(y0,y1,b1,a); }
  }else if(rget_sgn(p0)<=0 && rget_sgn(p1)>=0){
    mpfr_set_d(a,-1,MPFR_RNDD);
    mpfr_sin(b0,x0,MPFR_RNDU);
    mpfr_sin(b1,x1,MPFR_RNDU);
    if(rge(b0,b1)){ ircopy(y0,y1,a,b0); }else{ ircopy(y0,y1,a,b1); }
  }else{ ERROR_AT; exit(0); }
  RF(p0); RF(p1); RF(a); RF(b0); RF(b1);
}

/**
 @brief 三角関数 [y0,y1]=cos([x0,x1])
 */
void ircos(rmulti *y0, rmulti *y1, rmulti *x0, rmulti *x1)
{
  int prec;
  rmulti *p0=NULL,*p1=NULL,*a=NULL,*b0=NULL,*b1=NULL;
  prec=MAX2(rget_prec(y0),rget_prec(y1));
  RA(p0,prec); RA(p1,prec); RA(a,prec); RA(b0,prec); RA(b1,prec);
  mpfr_sin(p0,x0,MPFR_RNDN); mpfr_neg(p0,p0,MPFR_RNDN);
  mpfr_sin(p1,x1,MPFR_RNDN); mpfr_neg(p1,p1,MPFR_RNDN);
  if(rget_sgn(p0)>=0 && rget_sgn(p1)>=0){
    mpfr_cos(y0,x0,MPFR_RNDD);
    mpfr_cos(y1,x1,MPFR_RNDU);
  }else if(rget_sgn(p0)<=0 && rget_sgn(p1)<=0){
    mpfr_cos(y0,x1,MPFR_RNDD);
    mpfr_cos(y1,x0,MPFR_RNDU);
  }else if(rget_sgn(p0)>=0 && rget_sgn(p1)<=0){
    mpfr_cos(b0,x0,MPFR_RNDD);
    mpfr_cos(b1,x1,MPFR_RNDD);
    mpfr_set_d(a,1,MPFR_RNDU);
    if(rle(b0,b1)){ ircopy(y0,y1,b0,a); }else{ ircopy(y0,y1,b1,a); }
  }else if(rget_sgn(p0)<=0 && rget_sgn(p1)>=0){
    mpfr_set_d(a,-1,MPFR_RNDD);
    mpfr_cos(b0,x0,MPFR_RNDU);
    mpfr_cos(b1,x1,MPFR_RNDU);
    if(rge(b0,b1)){ ircopy(y0,y1,a,b0); }else{ ircopy(y0,y1,a,b1); }
  }else{ ERROR_AT; exit(0); }
  RF(p0); RF(p1); RF(a); RF(b0); RF(b1);
}

/**
 @brief 三角関数 [y0,y1]=tan([x0,x1])
 */
void irtan(rmulti *y0, rmulti *y1, rmulti *x0, rmulti *x1)
{
  mpfr_tan(y0,x0,MPFR_RNDD);
  mpfr_tan(y1,x1,MPFR_RNDU);
}

/**
 @brief 逆三角関数 [y0,y1]=asin([x0,x1])
*/
void irasin(rmulti *y0, rmulti *y1, rmulti *x0, rmulti *x1)
{
  mpfr_asin(y0,x0,MPFR_RNDD);
  mpfr_asin(y1,x1,MPFR_RNDU);
}

/**
 @brief 逆三角関数 [y0,y1]=acos([x0,x1])
*/
void iracos(rmulti *y0, rmulti *y1, rmulti *x0, rmulti *x1)
{
  int prec;
  rmulti *a0=NULL,*a1=NULL;
  prec=MAX2(rget_prec(y0),rget_prec(y1));
  RA(a0,prec); RA(a1,prec);
  mpfr_acos(a0,x1,MPFR_RNDD);
  mpfr_acos(a1,x0,MPFR_RNDU);
  ircopy(y0,y1,a0,a1);
  RF(a0); RF(a1);
}

/**
 @brief 逆三角関数 [y0,y1]=atan([x0,x1])
*/
void iratan(rmulti *y0, rmulti *y1, rmulti *x0, rmulti *x1)
{
  mpfr_atan(y0,x0,MPFR_RNDD);
  mpfr_atan(y1,x1,MPFR_RNDU);
}

/**
 @brief 逆三角関数 [z0,z1]=atan2([y0,y1]/[x0,x1])
*/
void iratan2(rmulti *z0, rmulti *z1, rmulti *y0, rmulti *y1, rmulti *x0, rmulti *x1)
{
  int prec;
  rmulti *a0=NULL,*a1=NULL;
  prec=MAX2(rget_prec(z0),rget_prec(z1));
  RA(a0,prec); RA(a1,prec);
  if(rget_sgn(y0)>=0){
    mpfr_atan2(a0,y0,x1,MPFR_RNDD);
    mpfr_atan2(a1,y1,x0,MPFR_RNDU);
  }else{
    mpfr_atan2(a0,y0,x0,MPFR_RNDD);
    mpfr_atan2(a1,y1,x1,MPFR_RNDU);
  }
  ircopy(z0,z1,a0,a1);
  RF(a0); RF(a1);
}

/**
 @brief 双曲線関数 [y0,y1]=sinh([x0,x1])
*/
void irsinh(rmulti *y0, rmulti *y1, rmulti *x0, rmulti *x1)
{
  mpfr_sinh(y0,x0,MPFR_RNDD);
  mpfr_sinh(y1,x1,MPFR_RNDU);
}

/**
 @brief 双曲線関数 [y0,y1]=cosh([x0,x1])
*/
void ircosh(rmulti *y0, rmulti *y1, rmulti *x0, rmulti *x1)
{
  rmulti *a=NULL,*b0=NULL,*b1=NULL;
  int prec;
  prec=MAX2(rget_prec(y0),rget_prec(y1));
  RA(a,prec); RA(b0,prec); RA(b1,prec);
  if(rget_sgn(x0)>0){
    mpfr_cosh(y0,x0,MPFR_RNDD);
    mpfr_cosh(y1,x1,MPFR_RNDU);
  }else if(rget_sgn(x1)<0){
    mpfr_cosh(b0,x1,MPFR_RNDD);
    mpfr_cosh(b1,x0,MPFR_RNDU);
    ircopy(y0,y1,b0,b1);
  }else{
    mpfr_set_d(a,1,MPFR_RNDD);
    mpfr_cosh(b0,x0,MPFR_RNDU);
    mpfr_cosh(b1,x1,MPFR_RNDU);
    if(rge(b0,b1)){ ircopy(y0,y1,a,b0); }else{ ircopy(y0,y1,a,b1); }
  }
  RF(a); RF(b0); RF(b1);
}

/**
 @brief 双曲線関数 [y0,y1]=tanh([x0,x1])
*/
void irtanh(rmulti *y0, rmulti *y1, rmulti *x0, rmulti *x1)
{
  mpfr_tanh(y0,x0,MPFR_RNDD);
  mpfr_tanh(y1,x1,MPFR_RNDU);
}

/**
 @beif 逆双曲線関数 [y0,y1]=asinh([x0,x1])
*/
void irasinh(rmulti *y0, rmulti *y1, rmulti *x0, rmulti *x1)
{
  mpfr_asinh(y0,x0,MPFR_RNDD);
  mpfr_asinh(y1,x1,MPFR_RNDU);
}

/**
 @brief 逆双曲線関数 [y0,y1]=acosh([x0,x1])
*/
void iracosh(rmulti *y0, rmulti *y1, rmulti *x0, rmulti *x1)
{
  mpfr_acosh(y0,x0,MPFR_RNDD);
  mpfr_acosh(y1,x1,MPFR_RNDU);
}

/**
 @brief 逆双曲線関数 [y0,y1]=atanh([x0,x1])
*/
void iratanh(rmulti *y0, rmulti *y1, rmulti *x0, rmulti *x1)
{
  mpfr_atanh(y0,x0,MPFR_RNDD);
  mpfr_atanh(y1,x1,MPFR_RNDU);
}

/** @} */

/////////////////////////////////////////////////////////////

/** @name irmulti型の２入力の数学演算子 */
/** @{ */

/**
 @brief 足し算 [z0,z1]=[x0,x1]+[y0,y1]
 */
void iradd(rmulti *z0, rmulti *z1, rmulti *x0, rmulti *x1, rmulti *y0, rmulti *y1)
{
  mpfr_add(z0,x0,y0,MPFR_RNDD); // lower bound
  mpfr_add(z1,x1,y1,MPFR_RNDU); // upper bound
}

/**
 @brief 足し算 [z0,z1]=[x0,x1]+[y0,y1]
 */
void idadd_r(rmulti *z0, rmulti *z1, double x0, double x1, rmulti *y0, rmulti *y1)
{
  mpfr_add_d(z0,y0,x0,MPFR_RNDD); // lower bound
  mpfr_add_d(z1,y1,x1,MPFR_RNDU); // upper bound
}
/**
 @brief 足し算 [z0,z1]=[x0,x1]+[y0,y1]
 */
void iradd_d(rmulti *z0, rmulti *z1, rmulti *x0, rmulti *x1, double y0, double y1)
{
  mpfr_add_d(z0,x0,y0,MPFR_RNDD); // lower bound
  mpfr_add_d(z1,x1,y1,MPFR_RNDU); // upper bound
}

///////////////////////////////////

/**
 @brief 引き算 [z0,z1]=[x0,x1]-[y0,y1]
 */
void irsub(rmulti *z0, rmulti *z1, rmulti *x0, rmulti *x1, rmulti *y0, rmulti *y1)
{
  int p0,p1,prec;
  rmulti *a0=NULL,*a1=NULL;
  p0=rget_prec(z0); p1=rget_prec(z1); prec=MAX2(p0,p1);
  RA(a0,prec); RA(a1,prec);
  mpfr_sub(a0,x0,y1,MPFR_RNDD); // lower bound
  mpfr_sub(a1,x1,y0,MPFR_RNDU); // upper bound  
  ircopy(z0,z1,a0,a1);
  RF(a0); RF(a1);
}

/**
 @brief 引き算 [z0,z1]=[x0,x1]-[y0,y1]
 */
void irsub_ws(rmulti *z0, rmulti *z1, rmulti *x0, rmulti *x1, rmulti *y0, rmulti *y1, int n_rws, rmulti **rws)
{
  if(n_rws<2){ ERROR_EXIT("Error `n_rws=%d<2' in irsub_ws().\n",n_rws); }
  mpfr_sub(rws[0],x0,y1,MPFR_RNDD); // lower bound
  mpfr_sub(rws[1],x1,y0,MPFR_RNDU); // upper bound
  ircopy(z0,z1,rws[0],rws[1]);
}

/**
 @brief 引き算 [z0,z1]=[x0,x1]-[y0,y1]
 */
void irsub_d(rmulti *z0, rmulti *z1, rmulti *x0, rmulti *x1, double y0, double y1)
{
  int p0,p1,prec;
  rmulti *a0=NULL,*a1=NULL;
  p0=rget_prec(z0); p1=rget_prec(z1); prec=MAX2(p0,p1);
  RA(a0,prec); RA(a1,prec);
  irset_d(a0,a1,y0,y1);
  irsub(z0,z1,x0,x1,a0,a1);
  RF(a0); RF(a1);
}

/**
 @brief 引き算 [z0,z1]=[x0,x1]-[y,y]
 */
void idsub_r(rmulti *z0, rmulti *z1, double x0, double x1, rmulti *y0, rmulti *y1)
{
  int p0,p1,prec;
  rmulti *a0=NULL,*a1=NULL;
  p0=rget_prec(z0); p1=rget_prec(z1); prec=MAX2(p0,p1);
  RA(a0,prec); RA(a1,prec);
  irset_d(a0,a1,x0,x1);
  irsub(z0,z1,a0,a1,y0,y1);
  RF(a0); RF(a1);
}

///////////////////////////////////


/**
 @brief 区間の拡張 [z0,z1]=[x0,x1]+[-y,y]
 */
void iradd_pm(rmulti *z0, rmulti *z1, rmulti *x0, rmulti *x1, rmulti *y)
{
  mpfr_sub(z0,x0,y,MPFR_RNDD); // lower bound
  mpfr_add(z1,x1,y,MPFR_RNDU); // upper bound
}


/**
 @brief 掛け算 [z0,z1]=[x0,x1]*[y0,y1]
*/
void irmul(rmulti *z0, rmulti *z1, rmulti *x0, rmulti *x1, rmulti *y0, rmulti *y1)
{
  rmulti *Z0=NULL,*Z1=NULL,*a=NULL,*b=NULL;
  int p0,p1,prec;
  p0=rget_prec(z0); p1=rget_prec(z1); prec=MAX2(p0,p1);
  RA(Z0,prec); RA(Z1,prec); RA(a,prec); RA(b,prec);
  if(rget_sgn(x0)>0){
    if     (rget_sgn(y0)>0){ mpfr_mul(Z0,x0,y0,MPFR_RNDD); mpfr_mul(Z1,x1,y1,MPFR_RNDU); }
    else if(rget_sgn(y1)<0){ mpfr_mul(Z0,x1,y0,MPFR_RNDD); mpfr_mul(Z1,x0,y1,MPFR_RNDU); }
    else                   { mpfr_mul(Z0,x1,y0,MPFR_RNDD); mpfr_mul(Z1,x1,y1,MPFR_RNDU); }
  }else if(rget_sgn(x1)<0){
    if     (rget_sgn(y0)>0){ mpfr_mul(Z0,x0,y1,MPFR_RNDD); mpfr_mul(Z1,x1,y0,MPFR_RNDU); }
    else if(rget_sgn(y1)<0){ mpfr_mul(Z0,x1,y1,MPFR_RNDD); mpfr_mul(Z1,x0,y0,MPFR_RNDU); }
    else                   { mpfr_mul(Z0,x0,y1,MPFR_RNDD); mpfr_mul(Z1,x0,y0,MPFR_RNDU); }
  }else{
    if     (rget_sgn(y0)>0){ mpfr_mul(Z0,x0,y1,MPFR_RNDD); mpfr_mul(Z1,x1,y1,MPFR_RNDU); }
    else if(rget_sgn(y1)<0){ mpfr_mul(Z0,x1,y0,MPFR_RNDD); mpfr_mul(Z1,x0,y0,MPFR_RNDU); }
    else{
      mpfr_mul(a,x0,y1,MPFR_RNDD);
      mpfr_mul(b,x1,y0,MPFR_RNDD);
      if(rle(a,b)){ mpfr_set(Z0,a,MPFR_RNDD); }else{ mpfr_set(Z0,b,MPFR_RNDD); }
      mpfr_mul(a,x0,y0,MPFR_RNDU);
      mpfr_mul(b,x1,y1,MPFR_RNDU);
      if(rge(a,b)){ mpfr_set(Z1,a,MPFR_RNDU); }else{ mpfr_set(Z1,b,MPFR_RNDU); }
    }
  }
  ircopy(z0,z1,Z0,Z1);
  RF(Z0); RF(Z1); RF(a); RF(b);
}

/**
 @brief 掛け算 [z0,z1]=[x0,x1]*[y0,y1]
*/
void irmul_ws(rmulti *z0, rmulti *z1, rmulti *x0, rmulti *x1, rmulti *y0, rmulti *y1, int n_rws, rmulti **rws)
{
  if(n_rws<4){ ERROR_EXIT("Error `n_rws=%d<4' in irmul_ws().\n",n_rws); }
  if(rget_sgn(x0)>0){
    if     (rget_sgn(y0)>0){ mpfr_mul(rws[0],x0,y0,MPFR_RNDD); mpfr_mul(rws[1],x1,y1,MPFR_RNDU); }
    else if(rget_sgn(y1)<0){ mpfr_mul(rws[0],x1,y0,MPFR_RNDD); mpfr_mul(rws[1],x0,y1,MPFR_RNDU); }
    else                   { mpfr_mul(rws[0],x1,y0,MPFR_RNDD); mpfr_mul(rws[1],x1,y1,MPFR_RNDU); }
  }else if(rget_sgn(x1)<0){
    if     (rget_sgn(y0)>0){ mpfr_mul(rws[0],x0,y1,MPFR_RNDD); mpfr_mul(rws[1],x1,y0,MPFR_RNDU); }
    else if(rget_sgn(y1)<0){ mpfr_mul(rws[0],x1,y1,MPFR_RNDD); mpfr_mul(rws[1],x0,y0,MPFR_RNDU); }
    else                   { mpfr_mul(rws[0],x0,y1,MPFR_RNDD); mpfr_mul(rws[1],x0,y0,MPFR_RNDU); }
  }else{
    if     (rget_sgn(y0)>0){ mpfr_mul(rws[0],x0,y1,MPFR_RNDD); mpfr_mul(rws[1],x1,y1,MPFR_RNDU); }
    else if(rget_sgn(y1)<0){ mpfr_mul(rws[0],x1,y0,MPFR_RNDD); mpfr_mul(rws[1],x0,y0,MPFR_RNDU); }
    else{
      mpfr_mul(rws[2],x0,y1,MPFR_RNDD);
      mpfr_mul(rws[3],x1,y0,MPFR_RNDD);
      if(rle(rws[2],rws[3])){ mpfr_set(rws[0],rws[2],MPFR_RNDD); }else{ mpfr_set(rws[0],rws[3],MPFR_RNDD); }
      mpfr_mul(rws[2],x0,y0,MPFR_RNDU);
      mpfr_mul(rws[3],x1,y1,MPFR_RNDU);
      if(rge(rws[2],rws[3])){ mpfr_set(rws[1],rws[2],MPFR_RNDU); }else{ mpfr_set(rws[1],rws[3],MPFR_RNDU); }
    }
  }
  ircopy(z0,z1,rws[0],rws[1]);
}

/**
 @brief 掛け算 [z0,z1]=[x0,x1]*[y,y]
*/
void irmul_d(rmulti *z0, rmulti *z1, rmulti *x0, rmulti *x1, double y)
{
  int p0,p1,prec;
  rmulti *a0=NULL,*a1=NULL;
  p0=rget_prec(z0); p1=rget_prec(z1); prec=MAX2(p0,p1);
  RA(a0,prec); RA(a1,prec);
  irset_d(a0,a1,y,y);
  irmul(z0,z1,x0,x1,a0,a1);
  RF(a0); RF(a1);
}

/**
 @brief 割り算 [z0,z1]=[x0,x1]/[y0,y1]
*/
void irdiv(rmulti *z0, rmulti *z1, rmulti *x0, rmulti *x1, rmulti *y0, rmulti *y1)
{
  rmulti *Z0=NULL,*Z1=NULL;
  int p0,p1,prec;
  p0=rget_prec(z0); p1=rget_prec(z1); prec=MAX2(p0,p1);
  RA(Z0,prec); RA(Z1,prec);
  if(rget_sgn(x0)>0){
    if     (rget_sgn(y0)>0){ mpfr_div(Z0,x0,y1,MPFR_RNDD); mpfr_div(Z1,x1,y0,MPFR_RNDU); }
    else if(rget_sgn(y1)<0){ mpfr_div(Z0,x1,y1,MPFR_RNDD); mpfr_div(Z1,x0,y0,MPFR_RNDU); }
    else                   { rset_nan(z0); rset_nan(z1); } // return nan if devided by 0
  }else if(rget_sgn(x1)<0){
    if     (rget_sgn(y0)>0){ mpfr_div(Z0,x0,y0,MPFR_RNDD); mpfr_div(Z1,x1,y1,MPFR_RNDU); }
    else if(rget_sgn(y1)<0){ mpfr_div(Z0,x1,y0,MPFR_RNDD); mpfr_div(Z1,x0,y1,MPFR_RNDU); }
    else                   { rset_nan(z0); rset_nan(z1); } // return nan if devided by 0
  }else{
    if     (rget_sgn(y0)>0){ mpfr_div(Z0,x0,y0,MPFR_RNDD); mpfr_div(Z1,x1,y0,MPFR_RNDU); }
    else if(rget_sgn(y1)<0){ mpfr_div(Z0,x1,y1,MPFR_RNDD); mpfr_div(Z1,x0,y1,MPFR_RNDU); }
    else                   { rset_nan(z0); rset_nan(z1); } // return nan if devided by 0
  }
  ircopy(z0,z1,Z0,Z1);
  RF(Z0); RF(Z1);
}

/**
 @brief 積の加算 [z0,z1]+=[x0,x1]*[y0,y1]
*/
void iradd_mul(rmulti *z0, rmulti *z1, rmulti *x0, rmulti *x1, rmulti *y0, rmulti *y1)
{
  int p0,p1,prec;
  rmulti *a0=NULL,*a1=NULL;
  p0=rget_prec(z0); p1=rget_prec(z1); prec=MAX2(p0,p1);
  RA(a0,prec); RA(a1,prec);
  irmul(a0,a1,x0,x1,y0,y1); // a=x*y
  iradd(z0,z1,z0,z1,a0,a1); // z=z+x*y
  RF(a0); RF(a1);
}

/**
 @brief 積の加算 [z0,z1]+=[x0,x1]*[y0,y1]
*/
void iradd_mul_ws(rmulti *z0, rmulti *z1, rmulti *x0, rmulti *x1, rmulti *y0, rmulti *y1, int n_rws, rmulti **rws)
{
  if(n_rws<6){ ERROR_EXIT("Error `n_rws=%d<6' in iradd_mul_ws().\n",n_rws); }
  irmul_ws(rws[0],rws[1],x0,x1,y0,y1,n_rws-2,rws+2); // a=x*y
  iradd(z0,z1,z0,z1,rws[0],rws[1]); // z=z+x*y
}

/**
 @brief 積の減算 [z0,z1]-=[x0,x1]*[y0,y1]
*/
void irsub_mul(rmulti *z0, rmulti *z1, rmulti *x0, rmulti *x1, rmulti *y0, rmulti *y1)
{
  int p0,p1,prec;
  rmulti *a0=NULL,*a1=NULL;
  p0=rget_prec(z0); p1=rget_prec(z1); prec=MAX2(p0,p1);
  RA(a0,prec); RA(a1,prec);
  irmul(a0,a1,x0,x1,y0,y1); // a=x*y
  irsub(z0,z1,z0,z1,a0,a1); // z=z-x*y
  RF(a0); RF(a1);
}

/**
 @brief 積の減算 [z0,z1]-=[x0,x1]*[y0,y1]
*/
void irsub_mul_ws(rmulti *z0, rmulti *z1, rmulti *x0, rmulti *x1, rmulti *y0, rmulti *y1, int n_rws, rmulti **rws)
{
  if(n_rws<6){ ERROR_EXIT("Error `n_rws=%d<6' in irsub_mul_ws().\n",n_rws); }
  irmul_ws(rws[0],rws[1],x0,x1,y0,y1,n_rws-2,rws+2); // a=x*y
  irsub_ws(z0,z1,z0,z1,rws[0],rws[1],n_rws-2,rws+2); // z=z-x*y
}

/**
 @brief 差の絶対値 [z0,z1]=abs([x0,x1]-[y0,y1])
*/
void irabs_sub(rmulti *z0, rmulti *z1, rmulti *x0, rmulti *x1, rmulti *y0, rmulti *y1)
{
  irsub(z0,z1,x0,x1,y0,y1);
  irabs(z0,z1,z0,z1);
}

/**
 @brief 商の分母が絶対値 [z0,z1]=[x0,x1]/abs([y0,y1])
*/
void irdiv_abs(rmulti *z0, rmulti *z1, rmulti *x0, rmulti *x1, rmulti *y0, rmulti *y1)
{
  irabs(z0,z1,y0,y1);
  irdiv(z0,z1,x0,x1,z0,z1);
}

/**
 @brief 商の分母が絶対値 [z0,z1]=[x0,x1]/abs([y0,y1])
*/
void irdiv_abs_c(rmulti *z0, rmulti *z1, rmulti *x0, rmulti *x1, cmulti *y0, cmulti *y1)
{
  icabs(z0,z1,y0,y1);
  irdiv(z0,z1,x0,x1,z0,z1);
}

/**
 @brief べき乗 [y0,y1]=[x0,x1]^n
 */
void irpow_si(rmulti *y0, rmulti *y1, rmulti *x0, rmulti *x1, long n)
{
  rmulti *a=NULL,*b0=NULL,*b1=NULL;
  int p0,p1,prec;
  p0=rget_prec(y0); p1=rget_prec(y1); prec=MAX2(p0,p1);
  RA(a,prec); RA(b0,prec); RA(b1,prec);
  if(n==0){
    mpfr_pow_si(y0,x0,n,MPFR_RNDD);
    mpfr_pow_si(y1,x1,n,MPFR_RNDU);
  }else if(n>0 && (n%2)>0){
    mpfr_pow_si(y0,x0,n,MPFR_RNDD);
    mpfr_pow_si(y1,x1,n,MPFR_RNDU);
  }else if(n>0){
    if(rget_sgn(x0)>0){
      mpfr_pow_si(y0,x0,n,MPFR_RNDD);
      mpfr_pow_si(y1,x1,n,MPFR_RNDU);
    }else if(rget_sgn(x1)<0){
      mpfr_pow_si(y0,x1,n,MPFR_RNDD);
      mpfr_pow_si(y1,x0,n,MPFR_RNDU);
    }else{
      mpfr_set_d(a,0,MPFR_RNDD);
      mpfr_pow_si(b0,x0,n,MPFR_RNDU);
      mpfr_pow_si(b1,x1,n,MPFR_RNDU);
      if(rge(b0,b1)){ ircopy(y0,y1,a,b0); }else{ ircopy(y0,y1,a,b1); }
    }
  }else{
    rset_d(a,1);
    irpow_si(b0,b1,x0,x1,-n);
    irdiv(y0,y1,a,a,b0,b1);
  }
  RF(a); RF(b0); RF(b1);
}

/** @} */

//////////////////////////////////////////////////////

/** @name irmulti型の判定に関する関数 */
/** @{ */

/**
 @brief 区間の内包の判定 [y0,y1] in [x0,x1], すなわち (x0<=y0 && y1 <=x1).
 */
int irin(rmulti *y0, rmulti *y1, rmulti *x0, rmulti *x1)
{
  return (rle(x0,y0) && rle(y1,x1));
}

/**
 @brief 区間の内包の判定 [y0,y1] in [-abs(x),abs(x)], すなわち (-abs(x)<=y0 && y1 <=abs(x)).
 */
int irin_pm(rmulti *y0, rmulti *y1, rmulti *x)
{
  int value;
  rmulti *z0=NULL,*z1=NULL;
  z0=rallocate_prec(rget_prec(x)); z1=rallocate_prec(rget_prec(x));
  irpm(z0,z1,x);
  value=(rle(z0,y0) && rle(y1,z1));
  z0=rfree(z0); z1=rfree(z1);
  return value;
}

/** @} */

//EOF
