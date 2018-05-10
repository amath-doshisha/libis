#include<isys.h>
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<stdarg.h>
#include"is_cmulti.h"
#include"is_cmat.h"
#include"is_irmulti.h"
#include"is_icmulti.h"
#include"is_func.h"

/**
 @file  icmulti.c
 @brief 多倍長精度実数型cmultiの機械区間演算に関する関数の定義.
 @details ベクトルに関しては@link icvec.c@endlinkを参照のこと.
          行列に関しては@link icmat.c@endlinkを参照のこと.
 */

#define RA(X,P) ((X)=rallocate_prec(P))
#define CA(X,P) ((X)=callocate_prec(P))
#define RF(X) ((X)=rfree(X))
#define CF(X) ((X)=cfree(X))
#define RAr(X,P) ((X)=rallocate_prec(rget_prec(P)))
#define CAr(X,P) ((X)=callocate_prec(rget_prec(P)))


/** @name icmulti型の値の設定に関する関数 */
/** @{ */

/**
 @brief コピー [y0,y1]=[x0,x1].
 */
void icset(cmulti *y0, cmulti *y1, cmulti *x0, cmulti *x1)
{
  ircopy(C_R(y0),C_R(y1),C_R(x0),C_R(x1));
  ircopy(C_I(y0),C_I(y1),C_I(x0),C_I(x1));
}

/**
 @brief 倍精度実数の設定 [y0,y1]=[x,x].
 */
void icset_d(cmulti *y0, cmulti *y1, double x0, double x1)
{
  irset_d(C_R(y0),C_R(y1),x0,x1);
  irset_d(C_I(y0),C_I(y1),0, 0);
}

/**
 @brief 倍精度複素数の設定 [y0,y1]=[x0,x1].
 */
void icset_z(cmulti *y0, cmulti *y1, dcomplex x0, dcomplex x1)
{
  irset_d(C_R(y0),C_R(y1),Z_R(x0),Z_R(x1));
  irset_d(C_I(y0),C_I(y1),Z_I(x0),Z_I(x1));
}

/**
 @brief コピー [y0,y1]=[x0,x1]. rmultiをcmultiにキャスト
 */
void icset_r(cmulti *y0, cmulti *y1, rmulti *x0, rmulti *x1)
{
  irset  (C_R(y0),C_R(y1),x0,x1);
  irset_d(C_I(y0),C_I(y1),0, 0);
}

/**
 @brief 倍精度複素数の設定 [y0,y1]=[x0_r+i*x0_i,xr1_*i*x1_i].
 */
void icset_Z(cmulti *y0, cmulti *y1, double x0_r, double x0_i, double x1_r, double x1_i)
{
  irset_d(C_R(y0),C_R(y1),x0_r,x1_r);
  irset_d(C_I(y0),C_I(y1),x0_i,x1_i);
}

/**
 @brief 倍精度複素数の設定 [y0,y1]=[x0,x1].
 */
void icset_zd(cmulti *y0, cmulti *y1, dcomplex x0, double x1)
{
  irset_d(C_R(y0),C_R(y1),Z_R(x0),x1);
  irset_d(C_I(y0),C_I(y1),Z_I(x0),0);
}

/**
 @brief 倍精度複素数の設定 [y0,y1]=[x0,x1].
 */
void icset_dz(cmulti *y0, cmulti *y1, double x0, dcomplex x1)
{
  irset_d(C_R(y0),C_R(y1),x0,Z_R(x1));
  irset_d(C_I(y0),C_I(y1),0, Z_I(x1));
}

/**
 @brief コピー [y0,y1]=[x0,x1]. rmultiをcmultiにキャスト
 */
void icset_rc(cmulti *y0, cmulti *y1, rmulti *x0, cmulti *x1)
{
  ircopy  (C_R(y0),C_R(y1),x0,C_R(x1));
  irset_dr(C_I(y0),C_I(y1),.0,C_I(x1));
}

/**
 @brief コピー [y0,y1]=[x0,x1]. rmultiをcmultiにキャスト
 */
void icset_cr(cmulti *y0, cmulti *y1, cmulti *x0, rmulti *x1)
{
  ircopy  (C_R(y0),C_R(y1),C_R(x0),x1);
  irset_rd(C_I(y0),C_I(y1),C_I(x0),.0);
}


void icset_s(cmulti *x0, cmulti *x1, char *s)
{
  int ret;
  strings *a=NULL;
  NULL_EXC3(x0,x1,s);
  cset_nan(x0); cset_nan(x1);
  a=strings_split_str_to_irmulti(s);
  if(strings_size(a)!=4){ ERROR_AT; }
  if(strings_at(a,0)!=NULL){ ret=mpfr_set_str(C_R(x0),strings_at(a,0),10,MPFR_RNDD); if(ret){ rset_nan(C_R(x0)); } }
  if(strings_at(a,1)!=NULL){ ret=mpfr_set_str(C_I(x0),strings_at(a,1),10,MPFR_RNDD); if(ret){ rset_nan(C_I(x0)); } }
  if(strings_at(a,2)!=NULL){ ret=mpfr_set_str(C_R(x1),strings_at(a,2),10,MPFR_RNDU); if(ret){ rset_nan(C_R(x1)); } }
  if(strings_at(a,3)!=NULL){ ret=mpfr_set_str(C_I(x1),strings_at(a,3),10,MPFR_RNDU); if(ret){ rset_nan(C_I(x1)); } }
  if(cis_nan(x1)){ mpfr_set(C_R(x1),C_R(x0),MPFR_RNDU); mpfr_set(C_I(x1),C_I(x0),MPFR_RNDU); }
}

/**
 @brief コピー [y0,y1]=[x0,x1].
 */
void iccopy(cmulti *y0, cmulti *y1, cmulti *x0, cmulti *x1)
{
  ircopy(C_R(y0),C_R(y1),C_R(x0),C_R(x1));
  ircopy(C_I(y0),C_I(y1),C_I(x0),C_I(x1));
}

/**
 @brief コピー [y0,y1]=[x0,x1]. rmultiをcmultiにキャスト
 */
void iccopy_r(cmulti *y0, cmulti *y1, rmulti *x0, rmulti *x1)
{
  ircopy (C_R(y0),C_R(y1),x0,x1);
  irset_d(C_I(y0),C_I(y1),0, 0);
}

/**
 @brief コピー [y0,y1]=[x0_r,x1_r]+i*[x0_i,x1_i]
 */
void iccopy_C(cmulti *y0, cmulti *y1, rmulti *x0_r, rmulti *x0_i, rmulti *x1_r, rmulti *x1_i)
{
  ircopy(C_R(y0),C_R(y1),x0_r,x1_r);
  ircopy(C_I(y0),C_I(y1),x0_i,x1_i);
}

/**
 @brief コピー [y0,y1]=[x0,x1]. rmultiをcmultiにキャスト
 */
void iccopy_rc(cmulti *y0, cmulti *y1, rmulti *x0, cmulti *x1)
{
  ircopy  (C_R(y0),C_R(y1),x0,C_R(x1));
  irset_dr(C_I(y0),C_I(y1),.0,C_I(x1));
}

/**
 @brief コピー [y0,y1]=[x0,x1]. rmultiをcmultiにキャスト
 */
void iccopy_cr(cmulti *y0, cmulti *y1, cmulti *x0, rmulti *x1)
{
  ircopy  (C_R(y0),C_R(y1),C_R(x0),x1);
  irset_rd(C_I(y0),C_I(y1),C_I(x0),.0);
}

/** @} */

//////////////////////////////////////////////////////////////

/** @name icmulti型の型変換に関する関数 */
/** @{ */

/**
 @brief bigint型から[z0,z1]へ型変換.
 */
void icset_bigint(cmulti *z0, cmulti *z1, bigint *x)
{
  irset_bigint(C_R(z0),C_R(z1),x);
  irset_d(C_I(z0),C_I(z1),0,0);
}

/** @} */

//////////////////////////////////////////////////////////////

/** @name icmulti型の１入力の数学演算子 */
/** @{ */

/**
 @brief 複素共役 [y0,y1]=conj([x0,x1])
*/
void icconj(cmulti *y0, cmulti *y1, cmulti *x0, cmulti *x1)
{
  ircopy(C_R(y0),C_R(y1),C_R(x0),C_R(x1));
  irneg(C_I(y0),C_I(y1),C_I(x0),C_I(x1));
}

/**
 @brief 符号の反転 [y0,y1]=-[x0,x1].
 */
void icneg(cmulti *y0, cmulti *y1, cmulti *x0, cmulti *x1)
{
  irneg(C_R(y0),C_R(y1),C_R(x0),C_R(x1));
  irneg(C_I(y0),C_I(y1),C_I(x0),C_I(x1));
}

/**
 @brief 符号の正負 [y0,y1]=[-abs(real(x))-i*abs(imag(x))],abs(real(x))+i*abs(imag(x))].
 */
void icpm(cmulti *y0, cmulti *y1, cmulti *x)
{
  irpm(C_R(y0),C_R(y1),C_R(x));
  irpm(C_I(y0),C_I(y1),C_I(x));
}

/**
 @brief 区間の中心 [m-r,m+r]=[x0,x1]
 */
void icmid(cmulti *mid, cmulti *x0, cmulti *x1)
{
  irmid(C_R(mid),C_R(x0),C_R(x1));
  irmid(C_I(mid),C_I(x0),C_I(x1));
}

/**
 @brief 区間の半径 [m-r,m+r]=[x0,x1]
 */
void icrad(cmulti *rad, cmulti *x0, cmulti *x1)
{
  irrad(C_R(rad),C_R(x0),C_R(x1));
  irrad(C_I(rad),C_I(x0),C_I(x1));
}

/**
 @brief 区間の中心と半径 [m-r,m+r]=[x0,x1]
 */
void icmr(cmulti *mid, cmulti *rad, cmulti *x0, cmulti *x1)
{
  irmr(C_R(mid),C_R(rad),C_R(x0),C_R(x1));
  irmr(C_I(mid),C_I(rad),C_I(x0),C_I(x1));
}

/**
 @brief 逆数 [z0,z1]=[1,1]/[x0,x1]
*/
void icinv(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1)
{
  int prec;
  rmulti *den0=NULL,*den1=NULL;
  cmulti *a0=NULL,*a1=NULL;
  prec=MAX2(cget_prec(z0),cget_prec(z1));
  RA(den0,prec); RA(den1,prec); CA(a0,prec); CA(a1,prec);
  icabs2(den0,den1,x0,x1);         // den=|x|^2
  icconj(a0,a1,x0,x1);             // z=conj(x)
  icdiv_r2(a0,a1,a0,a1,den0,den1); // z/=den
  iccopy(z0,z1,a0,a1);
  RF(den0); RF(den1); CF(a0); CF(a1);
}

/**
 @brief 絶対値の平方 [y0,y1]=abs([x0,x1])^2
*/
void icabs2(rmulti *y0, rmulti *y1, cmulti *x0, cmulti *x1)
{
  irmul    (y0,y1,C_R(x0),C_R(x1),C_R(x0),C_R(x1)); // y=x.r*x.r
  iradd_mul(y0,y1,C_I(x0),C_I(x1),C_I(x0),C_I(x1)); // y+=x.i*x.i
}

/**
 @brief 絶対値の平方の加算 [y0,y1]+=abs([x0,x1])^2
*/
void icadd_abs2(rmulti *y0, rmulti *y1, cmulti *x0, cmulti *x1)
{
  int prec;
  rmulti *a0=NULL,*a1=NULL;
  prec=MAX2(rget_prec(y0),rget_prec(y1));
  RA(a0,prec); RA(a1,prec);
  icabs2(a0,a1,x0,x1);      // a=abs(x)^2
  iradd(y0,y1,y0,y1,a0,a1); // y=y+a
  RF(a0); RF(a1);
}

/**
 @brief 絶対値 [y0,y1]=abs([x0,x1])
*/
void icabs(rmulti *y0, rmulti *y1, cmulti *x0, cmulti *x1)
{
  int prec;
  cmulti *a0=NULL,*a1=NULL;
  if(cis_nan(x0) || cis_nan(x1)){rset_nan(y0); rset_nan(y1); return; }
  prec=MAX2(rget_prec(y0),rget_prec(y1));
  CA(a0,prec); CA(a1,prec);
  irabs(C_R(a0),C_R(a1),C_R(x0),C_R(x1)); 
  irabs(C_I(a0),C_I(a1),C_I(x0),C_I(x1));  
  icabs2(y0,y1,a0,a1);
  irsqrt(y0,y1,y0,y1);
  CF(a0); CF(a1);
}

/**
 @brief 絶対値 [y0,y1]=abs([x0,x1])
*/
void icabs_ws(rmulti *y0, rmulti *y1, cmulti *x0, cmulti *x1, int n_cws, cmulti **cws)
{
  if(n_cws<2){ ERROR_EXIT("Error `n_cws=%d<2' in icabs_ws().\n",n_cws); }
  if(cis_nan(x0) || cis_nan(x1)){ rset_nan(y0); rset_nan(y1); return; }
  irabs(C_R(cws[0]),C_R(cws[1]),C_R(x0),C_R(x1)); 
  irabs(C_I(cws[0]),C_I(cws[1]),C_I(x0),C_I(x1));  
  icabs2(y0,y1,cws[0],cws[1]);
  irsqrt(y0,y1,y0,y1);
}

/**
 @brief 平方根 [y0,y1]=sqrt([x0,x1])
*/
void icsqrt(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1)
{
  int prec;
  rmulti *r0=NULL,*r1=NULL,*theta0=NULL,*theta1=NULL;
  prec=MAX2(cget_prec(z0),cget_prec(z1));
  RA(r0,prec); RA(r1,prec); RA(theta0,prec); RA(theta1,prec);
  icget_polar(r0,r1,theta0,theta1,x0,x1);   // x=r*exp(i*theta)
  irsqrt(r0,r1,r0,r1);                      // r=sqrt(r)
  irmul_d(theta0,theta1,theta0,theta1,0.5); // theta=theta/2
  icset_polar(z0,z1,r0,r1,theta0,theta1);   // z=r^y*cos(theta*y)+i*r^y*sin(theta*y)
  RF(r0); RF(r1); RF(theta0); RF(theta1);
}

/**
 @brief 指数関数 [y0,y1]=exp([x0,x1])
*/
void icexp(cmulti *y0, cmulti *y1, cmulti *x0, cmulti *x1)
{
  int prec;
  rmulti *r0=NULL,*r1=NULL,*theta0=NULL,*theta1=NULL;
  prec=MAX2(cget_prec(y0),cget_prec(y1));
  RA(r0,prec); RA(r1,prec); RA(theta0,prec); RA(theta1,prec);
  irexp(r0,r1,C_R(x0),C_R(x1));
  ircopy(theta0,theta1,C_I(x0),C_I(x1));
  icset_polar(y0,y1,r0,r1,theta0,theta1); // y=exp(x.r)*(cos(x.i)+i*sin(x.i))
  RF(r0); RF(r1); RF(theta0); RF(theta1);
}

/**
 @brief 対数関数 [y0,y1]=log([x0,x1])
*/
void iclog(cmulti *y0, cmulti *y1, cmulti *x0, cmulti *x1)
{
  int prec;
  rmulti *r0=NULL,*r1=NULL,*theta0=NULL,*theta1=NULL;
  prec=MAX2(cget_prec(y0),cget_prec(y1));
  RA(r0,prec); RA(r1,prec); RA(theta0,prec); RA(theta1,prec);
  icget_polar(r0,r1,theta0,theta1,x0,x1); // x=r*exp(i*theta)
  irlog(C_R(y0),C_R(y1),r0,r1);           // y=log(r)+i*theta
  ircopy(C_I(y0),C_I(y1),theta0,theta1);  // y=log(r)+i*theta
  RF(r0); RF(r1); RF(theta0); RF(theta1);
}

/**
 @brief 三角関数 [y0,y1]=sin([x0,x1])
 */
void icsin(cmulti *y0, cmulti *y1, cmulti *x0, cmulti *x1)
{
  int prec;
  rmulti *c0=NULL,*c1=NULL,*ch0=NULL,*ch1=NULL,*s0=NULL,*s1=NULL,*sh0=NULL,*sh1=NULL;
  prec=MAX2(cget_prec(y0),cget_prec(y1));
  RA(c0,prec); RA(c1,prec); RA(ch0,prec); RA(ch1,prec); RA(s0,prec); RA(s1,prec); RA(sh0,prec); RA(sh1,prec);
  irsin(s0,s1,C_R(x0),C_R(x1)); ircosh(ch0,ch1,C_I(x0),C_I(x1)); irmul(C_R(y0),C_R(y1),s0,s1,ch0,ch1); // z.r=sin(x.r)*cosh(x.i)
  ircos(c0,c1,C_R(x0),C_R(x1)); irsinh(sh0,sh1,C_I(x0),C_I(x1)); irmul(C_I(y0),C_I(y1),c0,c1,sh0,sh1); // z.i=cos(x.r)*sinh(x.i)
  RF(c0); RF(c1); RF(ch0); RF(ch1); RF(s0); RF(s1); RF(sh0); RF(sh1);
}

/**
 @brief 三角関数 [y0,y1]=cos([x0,x1])
 */
void iccos(cmulti *y0, cmulti *y1, cmulti *x0, cmulti *x1)
{
  int prec;
  rmulti *c0=NULL,*ch0=NULL,*s0=NULL,*sh0=NULL,*c1=NULL,*ch1=NULL,*s1=NULL,*sh1=NULL;
  prec=MAX2(cget_prec(y0),cget_prec(y1));
  RA(c0,prec); RA(ch0,prec); RA(s0,prec); RA(sh0,prec); RA(c1,prec); RA(ch1,prec); RA(s1,prec); RA(sh1,prec);
  ircos(c0,c1,C_R(x0),C_R(x1)); ircosh(ch0,ch1,C_I(x0),C_I(x1)); irmul(C_R(y0),C_R(y1),c0,c1,ch0,ch1);                                            // z.r=cos(x.r)*cosh(x.i)
  irsin(s0,s1,C_R(x0),C_R(x1)); irsinh(sh0,sh1,C_I(x0),C_I(x1)); irmul(C_I(y0),C_I(y1),s0,s1,sh0,sh1); irneg(C_I(y0),C_I(y1),C_I(y0),C_I(y1)); // z.i=-sin(x.r)*sinh(x.i)
  RF(c0); RF(ch0); RF(s0); RF(sh0); RF(c1); RF(ch1); RF(s1); RF(sh1);
}

/**
 @brief 三角関数 [y0,y1]=tan([x0,x1])
 */
void ictan(cmulti *y0, cmulti *y1, cmulti *x0, cmulti *x1)
{
  int prec;
  cmulti *c0=NULL,*s0=NULL,*c1=NULL,*s1=NULL;
  prec=MAX2(cget_prec(y0),cget_prec(y1));
  CA(c0,prec); CA(s0,prec); CA(c1,prec); CA(s1,prec);
  iccos(c0,c1,x0,x1); icsin(s0,s1,x0,x1); icdiv(y0,y1,s0,s1,c0,c1); // y=s/c
  CF(c0); CF(s0); CF(c1); CF(s1);
}

/**
 @brief 逆三角関数 [y0,y1]=asin([x0,x1])
*/
void icasin(cmulti *y0, cmulti *y1, cmulti *x0, cmulti *x1)
{
  int prec;
  cmulti *a0=NULL,*a1=NULL,*b0=NULL,*b1=NULL;
  prec=MAX2(cget_prec(y0),cget_prec(y1));
  CA(a0,prec); CA(a1,prec); CA(b0,prec); CA(b1,prec);
  icmul(a0,a1,x0,x1,x0,x1); // a=x^2 
  icset_d(b0,b1,1,1);       // b=1
  icsub(a0,a1,b0,b1,a0,a1); // a=1-x^2
  icsqrt(a0,a1,a0,a1);      // a=sqrt(1-x^2)
  icset_Z(b0,b1,0,1,0,1);   // b=I
  icmul(b0,b1,b0,b1,x0,x1); // b=I*x
  icadd(a0,a1,a0,a1,b0,b1); // a=I*x+sqrt(1-x^2)
  iclog(a0,a1,a0,a1);       // a=log(I*x+sqrt(1-x^2))
  icset_Z(b0,b1,0,-1,0,01); // b=-I
  icmul(a0,a1,a0,a1,b0,b1); // a=-I*log(I*x+sqrt(1-x^2))
  iccopy(y0,y1,a0,a1);      // y=-I*log(I*x+sqrt(1-x^2))
  CF(a0); CF(a1); CF(b0); CF(b1);
}

/**
 @brief 逆三角関数 [y0,y1]=acos([x0,x1])
*/
void icacos(cmulti *y0, cmulti *y1, cmulti *x0, cmulti *x1)
{
  int prec;
  cmulti *a0=NULL,*a1=NULL,*b0=NULL,*b1=NULL;
  prec=MAX2(cget_prec(y0),cget_prec(y1));
  CA(a0,prec); CA(a1,prec); CA(b0,prec); CA(b1,prec);
  icmul(a0,a1,x0,x1,x0,x1);  // a=x^2 
  icset_d(b0,b1,-1,-1);      // b=-1
  icadd(a0,a1,a0,a1,b0,b1);  // a=x^2-1
  icsqrt(a0,a1,a0,a1);       // a=sqrt(x^2-1)
  icadd(a0,a1,a0,a1,x0,x1);  // a=x+sqrt(x^2-1)
  iclog(a0,a1,a0,a1);        // a=log(x+sqrt(x^2-1))
  icset_Z(b0,b1,0,-1,0,-1);  // b=-I
  icmul(a0,a1,a0,a1,b0,b1);  // a=-I*log(I*x+sqrt(1-x^2))
  iccopy(y0,y1,a0,a1);       // y=-I*log(I*x+sqrt(1-x^2))
  CF(a0); CF(a1); CF(b0); CF(b1);
}

/**
 @brief 逆三角関数 [y0,y1]=atan([x0,x1])
*/
void icatan(cmulti *y0, cmulti *y1, cmulti *x0, cmulti *x1)
{
  int prec;
  cmulti *a0=NULL,*a1=NULL,*b0=NULL,*b1=NULL;
  prec=MAX2(cget_prec(y0),cget_prec(y1));
  CA(a0,prec); CA(a1,prec); CA(b0,prec); CA(b1,prec);
  icset_d(b0,b1,1,1);       // b=1
  icadd(a0,a1,x0,x1,b0,b1); // a=1+x
  icsub(b0,b1,b0,b1,x0,x1); // b=1-x
  icdiv(a0,a1,a0,a1,b0,b1); // a=(1+x)/(1-x)
  iclog(a0,a1,a0,a1);       // a=log((1+x)/(1-x))
  icset_Z(b0,b1,0,-0.5,0,-0.5); // b=-I/2
  icmul(a0,a1,a0,a1,b0,b1); // a=(-I/2)*log((1+x)/(1-x))
  iccopy(y0,y1,a0,a1);      // y=(-I/2)*log((1+x)/(1-x))
  CF(a0); CF(a1); CF(b0); CF(b1);
}

/**
 @brief 双曲線関数 [y0,y1]=sinh([x0,x1])
*/
void icsinh(cmulti *y0, cmulti *y1, cmulti *x0, cmulti *x1)
{
  int prec;
  rmulti *c0=NULL,*ch0=NULL,*s0=NULL,*sh0=NULL,*c1=NULL,*ch1=NULL,*s1=NULL,*sh1=NULL;
  prec=MAX2(cget_prec(y0),cget_prec(y1));
  RA(c0,prec); RA(ch0,prec); RA(s0,prec); RA(sh0,prec); RA(c1,prec); RA(ch1,prec); RA(s1,prec); RA(sh1,prec);
  irsinh(sh0,sh1,C_R(x0),C_R(x1)); ircos(c0,c1,C_I(x0),C_I(x1)); irmul(C_R(y0),C_R(y1),sh0,sh1,c0,c1); // z.r=sinh(x.r)*cos(x.i)
  ircosh(ch0,ch1,C_R(x0),C_R(x1)); irsin(s0,s1,C_I(x0),C_I(x1)); irmul(C_I(y0),C_I(y1),ch0,ch1,s0,s1); // z.i=cosh(x.r)*sin(x.i)
  RF(c0); RF(ch0); RF(s0); RF(sh0); RF(c1); RF(ch1); RF(s1); RF(sh1);
}

/**
 @brief 双曲線関数 [y0,y1]=cosh([x0,x1])
*/
void iccosh(cmulti *y0, cmulti *y1, cmulti *x0, cmulti *x1)
{
  int prec;
  rmulti *c0=NULL,*ch0=NULL,*s0=NULL,*sh0=NULL,*c1=NULL,*ch1=NULL,*s1=NULL,*sh1=NULL;
  prec=MAX2(cget_prec(y0),cget_prec(y1));
  RA(c0,prec); RA(ch0,prec); RA(s0,prec); RA(sh0,prec); RA(c1,prec); RA(ch1,prec); RA(s1,prec); RA(sh1,prec);
  ircosh(ch0,ch1,C_R(x0),C_R(x1)); ircos(c0,c1,C_I(x0),C_I(x1)); irmul(C_R(y0),C_R(y1),ch0,ch1,c0,c1); // z.r=cosh(x.r)*cos(x.i)
  irsinh(sh0,sh1,C_R(x0),C_R(x1)); irsin(s0,s1,C_I(x0),C_I(x1)); irmul(C_I(y0),C_I(y1),sh0,sh1,s0,s1); // z.i=sinh(x.r)*sin(x.i)
  RF(c0); RF(ch0); RF(s0); RF(sh0); RF(c1); RF(ch1); RF(s1); RF(sh1);
}

/**
 @brief 双曲線関数 [y0,y1]=tanh([x0,x1])
*/
void ictanh(cmulti *y0, cmulti *y1, cmulti *x0, cmulti *x1)
{
  int prec;
  cmulti *ch0=NULL,*sh0=NULL,*ch1=NULL,*sh1=NULL;
  prec=MAX2(cget_prec(y0),cget_prec(y1));
  CA(ch0,prec); CA(sh0,prec); CA(ch1,prec); CA(sh1,prec);
  iccosh(ch0,ch1,x0,x1);
  icsinh(sh0,sh1,x0,x1);
  icdiv(y0,y1,sh0,sh1,ch0,ch1); // y=sh/ch
  CF(ch0); CF(sh0); CF(ch1); CF(sh1);
}

/**
 @beif 逆双曲線関数 [y0,y1]=asinh([x0,x1])
*/
void icasinh(cmulti *y0, cmulti *y1, cmulti *x0, cmulti *x1)
{
  int prec;
  cmulti *a0=NULL,*a1=NULL,*b0=NULL,*b1=NULL;
  prec=MAX2(cget_prec(y0),cget_prec(y1));
  CA(a0,prec); CA(a1,prec); CA(b0,prec); CA(b1,prec);
  icset_d(b0,b1,1,1);        // b=1
  icmul(a0,a1,x0,x1,x0,x1);  // a=x^2
  icadd(a0,a1,a0,a1,b0,b1);  // a=x^2+1
  icsqrt(a0,a1,a0,a1);       // a=sqrt(x^2+1)
  icadd(a0,a1,a0,a1,x0,x1);  // a=x+sqrt(x^2+1)
  iclog(a0,a1,a0,a1);        // a=log(x+sqrt(x^2+1))
  iccopy(y0,y1,a0,a1);       // y=log(x+sqrt(x^2+1))
  CF(a0); CF(a1); CF(b0); CF(b1);
}

/**
 @brief 逆双曲線関数 [y0,y1]=acosh([x0,x1])
*/
void icacosh(cmulti *y0, cmulti *y1, cmulti *x0, cmulti *x1)
{
  int prec;
  cmulti *a0=NULL,*a1=NULL,*b0=NULL,*b1=NULL;
  prec=MAX2(cget_prec(y0),cget_prec(y1));
  CA(a0,prec); CA(a1,prec); CA(b0,prec); CA(b1,prec);
  icset_d(b0,b1,-1,-1);      // b=-1
  icmul(a0,a1,x0,x1,x0,x1);  // a=x^2
  icadd(a0,a1,a0,a1,b0,b1);  // a=x^2-1
  icsqrt(a0,a1,a0,a1);       // a=sqrt(x^2-1)
  icadd(a0,a1,a0,a1,x0,x1);  // a=x+sqrt(x^2-1)
  iclog(a0,a1,a0,a1);        // a=log(x+sqrt(x^2-1))
  iccopy(y0,y1,a0,a1);       // y=log(x+sqrt(x^2-1))
  CF(a0); CF(a1); CF(b0); CF(b1);
}

/**
 @brief 逆双曲線関数 [y0,y1]=atanh([x0,x1])
*/
void icatanh(cmulti *y0, cmulti *y1, cmulti *x0, cmulti *x1)
{
  int prec;
  cmulti *a0=NULL,*a1=NULL,*b0=NULL,*b1=NULL;
  prec=MAX2(cget_prec(y0),cget_prec(y1));
  CA(a0,prec); CA(a1,prec); CA(b0,prec); CA(b1,prec);
  icset_d(b0,b1,1,1) ;      // b=1
  icadd(a0,a1,x0,x1,b0,b1); // a=1+x
  icsub(b0,b1,b0,b1,x0,x1); // b=1-x
  icdiv(a0,a1,a0,a1,b0,b1); // a=(1+x)/(1-x)
  iclog(a0,a1,a0,a1);       // a=log((1+x)/(1-x))
  icset_d(b0,b1,0.5,0.5);   // b=0.5
  icmul(a0,a1,a0,a1,b0,b1); // a=(1/2)*log((1+x)/(1-x))
  iccopy(y0,y1,a0,a1);      // y=(1/2)*log((1+x)/(1-x))
  CF(a0); CF(a1); CF(b0); CF(b1);
}

/** @} */

//////////////////////////////////////////////////////////////

/** @name icmulti型の２入力の数学演算子 */
/** @{ */

/**
 @brief 足し算 [z0,z1]=[x0,x1]+[y0,y1]
 */
void icadd(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, cmulti *y0, cmulti *y1)
{
  iradd(C_R(z0),C_R(z1),C_R(x0),C_R(x1),C_R(y0),C_R(y1)); // z.r=x.r+y.r
  iradd(C_I(z0),C_I(z1),C_I(x0),C_I(x1),C_I(y0),C_I(y1)); // z.i=x.i+y.i
}

/**
 @brief 足し算 [z0,z1]=[x0,x1]+[y0,y1]
 */
void icadd_r(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, rmulti *y0, rmulti *y1)
{
  iradd(C_R(z0),C_R(z1),C_R(x0),C_R(x1),y0,y1); // z.r=x.r+y.r
  irset(C_I(z0),C_I(z1),C_I(x0),C_I(x1));       // z.i=x.i
}

/**
 @brief 足し算 [z0,z1]=[x0,x1]+[y0,y1]
 */
void iradd_c(cmulti *z0, cmulti *z1, rmulti *x0, rmulti *x1, cmulti *y0, cmulti *y1)
{
  iradd(C_R(z0),C_R(z1),x0,x1,C_R(y0),C_R(y1)); // z.r=x.r+y.r
  irset(C_I(z0),C_I(z1),      C_I(y0),C_I(y1)); // z.i=    y.i
}

/**
 @brief 足し算 [z0,z1]=[x0,x1]+[y0,y1]
 */
void icadd_z(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, dcomplex y0, dcomplex y1)
{
  iradd_d(C_R(z0),C_R(z1),C_R(x0),C_R(x1),Z_R(y0),Z_R(y1)); // z.r=x.r+y.r
  iradd_d(C_I(z0),C_I(z1),C_I(x0),C_I(x1),Z_I(y0),Z_I(y1)); // z.i=x.i+y.i
}

/**
 @brief 足し算 [z0,z1]=[x0,x1]+[y0,y1]
 */
void izadd_c(cmulti *z0, cmulti *z1, dcomplex x0, dcomplex x1, cmulti *y0, cmulti *y1)
{
  idadd_r(C_R(z0),C_R(z1),Z_R(x0),Z_R(x1),C_R(y0),C_R(y1)); // z.r=x.r+y.r
  idadd_r(C_I(z0),C_I(z1),Z_I(x0),Z_I(x1),C_I(y0),C_I(y1)); // z.i=x.i+y.i
}

/**
 @brief 足し算 [z0,z1]=[x0,x1]+[y0,y1]
 */
void icadd_d(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, double y0, double y1)
{
  iradd_d(C_R(z0),C_R(z1),C_R(x0),C_R(x1),y0,y1); // z.r=x.r+y.r
  irset   (C_I(z0),C_I(z1),C_I(x0),C_I(x1));       // z.i=x.i
}

/**
 @brief 足し算 [z0,z1]=[x0,x1]+[y0,y1]
 */
void idadd_c(cmulti *z0, cmulti *z1, double x0, double x1, cmulti *y0, cmulti *y1)
{
    idadd_r(C_R(z0),C_R(z1),x0,x1,C_R(y0),C_R(y1)); // z.r=x.r+y.r
    irset   (C_I(z0),C_I(z1),      C_I(y0),C_I(y1)); // z.i=    y.i
}

/**
 @brief 足し算 [z0,z1]=[x0,x1]+[y0,y1]
 */
void iradd_z(cmulti *z0, cmulti *z1, rmulti *x0, rmulti *x1, dcomplex y0, dcomplex y1)
{
  iradd_d(C_R(z0),C_R(z1),x0,x1,Z_R(y0),Z_R(y1)); // z.r=x+y.r
  irset_d (C_I(z0),C_I(z1),      Z_I(y0),Z_I(y1)); // z.i=  y.i
}

/**
 @brief 足し算 [z0,z1]=[x0,x1]+[y0,y1]
 */
void izadd_r(cmulti *z0, cmulti *z1, dcomplex x0, dcomplex x1, rmulti *y0, rmulti *y1)
{
  idadd_r(C_R(z0),C_R(z1),Z_R(x0),Z_R(x1),y0,y1); // z.r=x.r+y
  irset_d (C_I(z0),C_I(z1),Z_I(x0),Z_I(x1));       // z.i=x.i
}

////////////////////////////////////////////////////////////

/**
 @brief 引き算 [z0,z1]=[x0,x1]-[y0,y1]
 */
void icsub_ws(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, cmulti *y0, cmulti *y1, int n_rws, rmulti **rws)
{
  if(n_rws<2){ ERROR_EXIT("Error `n_rws=%d<2' in icsub__ws().\n",n_rws); }
  irsub_ws(C_R(z0),C_R(z1),C_R(x0),C_R(x1),C_R(y0),C_R(y1),n_rws,rws); // z.r=x.r-y.r
  irsub_ws(C_I(z0),C_I(z1),C_I(x0),C_I(x1),C_I(y0),C_I(y1),n_rws,rws); // z.i=x.i-y.i
}

////////////////////////////////////////////////////////////

/**
 @brief 引き算 [z0,z1]=[x0,x1]-[y0,y1]
 */
void icsub(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, cmulti *y0, cmulti *y1)
{
  irsub(C_R(z0),C_R(z1),C_R(x0),C_R(x1),C_R(y0),C_R(y1)); // z.r=x.r-y.r
  irsub(C_I(z0),C_I(z1),C_I(x0),C_I(x1),C_I(y0),C_I(y1)); // z.i=x.i-y.i
}


/**
 @brief 引き算 [z0,z1]=[x0,x1]-[y0,y1]
 */
void icsub_r(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, rmulti *y0, rmulti *y1)
{
  irsub (C_R(z0),C_R(z1),C_R(x0),C_R(x1),y0,y1); // z.r=x.r-y.r
  ircopy(C_I(z0),C_I(z1),C_I(x0),C_I(x1));       // z.i=x.i
}

/**
 @brief 引き算 [z0,z1]=[x0,x1]-[y0,y1]
 */
void irsub_c(cmulti *z0, cmulti *z1, rmulti *x0, rmulti *x1, cmulti *y0, cmulti *y1)
{
  irsub(C_R(z0),C_R(z1),x0,x1,C_R(y0),C_R(y1)); // z.r=x.r-y.r
  irneg(C_I(z0),C_I(z1),      C_I(y0),C_I(y1)); // z.i=   -y.i
}

/**
 @brief 引き算 [z0,z1]=[x0,x1]-[y0,y1]
 */
void icsub_z(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, dcomplex y0, dcomplex y1)
{
  irsub_d(C_R(z0),C_R(z1),C_R(x0),C_R(x1),Z_R(y0),Z_R(y1)); // z.r=x.r-y.r
  irsub_d(C_I(z0),C_I(z1),C_I(x0),C_I(x1),Z_I(y0),Z_I(y1)); // z.i=x.i-y.i
}

/**
 @brief 引き算 [z0,z1]=[x0,x1]-[y0,y1]
 */
void izsub_c(cmulti *z0, cmulti *z1, dcomplex x0, dcomplex x1, cmulti *y0, cmulti *y1)
{
  idsub_r(C_R(z0),C_R(z1),Z_R(x0),Z_R(x1),C_R(y0),C_R(y1)); // z.r=x.r-y.r
  idsub_r(C_I(z0),C_I(z1),Z_I(x0),Z_I(x1),C_I(y0),C_I(y1)); // z.i=x.i-y.i
}

/**
 @brief 引き算 [z0,z1]=[x0,x1]-[y0,y1]
 */
void icsub_d(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, double y0, double y1)
{
  irsub_d(C_R(z0),C_R(z1),C_R(x0),C_R(x1),y0,y1); // z.r=x.r-y
  irset  (C_I(z0),C_I(z1),C_I(x0),C_I(x1));       // z.i=x.i
}

/**
 @brief 引き算 [z0,z1]=[x0,x1]-[y0,y1]
 */
void idsub_c(cmulti *z0, cmulti *z1, double x0, double x1, cmulti *y0, cmulti *y1)
{
  idsub_r(C_R(z0),C_R(z1),x0,x1,C_R(y0),C_R(y1)); // z.r=x.r-y.r
  irneg  (C_I(z0),C_I(z1),      C_I(y0),C_I(y1)); // z.i=   -y.i
}

/**
 @brief 引き算 [z0,z1]=[x0,x1]-[y0,y1]
 */
void irsub_z(cmulti *z0, cmulti *z1, rmulti *x0, rmulti *x1, dcomplex y0, dcomplex y1)
{
  irsub_d(C_R(z0),C_R(z1),x0,x1, Z_R(y0), Z_R(y1)); // z.r=x-y.r
  irset_d(C_I(z0),C_I(z1),      -Z_I(y0),-Z_I(y1)); // z.i= -y.i
}

/**
 @brief 引き算 [z0,z1]=[x0,x1]-[y0,y1]
 */
void izsub_r(cmulti *z0, cmulti *z1, dcomplex x0, dcomplex x1, rmulti *y0, rmulti *y1)
{
  idsub_r(C_R(z0),C_R(z1),Z_R(x0),Z_R(x1),y0,y1); // z.r=x.r-y
  irset_d(C_I(z0),C_I(z1),Z_I(x0),Z_I(x1));       // z.i=x.i
}

////////////////////////////////////////////////////////////

/**
 @brief 足し算 [z0,z1]=[x0,x1]+[-y,y]
 */
void icadd_pm(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, cmulti *y)
{
  iradd_pm(C_R(z0),C_R(z1),C_R(x0),C_R(x1),C_R(y)); // z.r=x.r+y.r
  iradd_pm(C_I(z0),C_I(z1),C_I(x0),C_I(x1),C_I(y)); // z.i=x.i+y.i
}

/**
 @brief 掛け算 [z0,z1]=[x0,x1]*[y0,y1]
*/
void icmul(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, cmulti *y0, cmulti *y1)
{
  int prec;
  cmulti *a0=NULL,*a1=NULL;
  prec=MAX2(cget_prec(z0),cget_prec(z1));
  CA(a0,prec); CA(a1,prec);
  irmul    (C_R(a0),C_R(a1),C_R(x0),C_R(x1),C_R(y0),C_R(y1)); // z.r =x.r*y.r
  irsub_mul(C_R(a0),C_R(a1),C_I(x0),C_I(x1),C_I(y0),C_I(y1)); // z.r-=x.i*y.i
  irmul    (C_I(a0),C_I(a1),C_I(x0),C_I(x1),C_R(y0),C_R(y1)); // z.i =x.i*y.r
  iradd_mul(C_I(a0),C_I(a1),C_R(x0),C_R(x1),C_I(y0),C_I(y1)); // z.i+=x.r*y.i
  iccopy(z0,z1,a0,a1);
  CF(a0); CF(a1);
}

/**
 @brief 掛け算 [z0,z1]=[x0,x1]*[y0,y1]
*/
void icmul_ws(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, cmulti *y0, cmulti *y1, int n_rws, rmulti **rws, int n_cws, cmulti **cws)
{
  if(n_rws<8){ ERROR_EXIT("Error `n_rws=%d<8' in icmul_ws().\n",n_rws); }
  if(n_cws<2){ ERROR_EXIT("Error `n_cws=%d<2' in icmul_ws().\n",n_cws); }
  irmul_ws    (C_R(cws[0]),C_R(cws[1]),C_R(x0),C_R(x1),C_R(y0),C_R(y1),n_rws-2,rws+2); // z.r =x.r*y.r
  irsub_mul_ws(C_R(cws[0]),C_R(cws[1]),C_I(x0),C_I(x1),C_I(y0),C_I(y1),n_rws-2,rws+2); // z.r-=x.i*y.i
  irmul_ws    (C_I(cws[0]),C_I(cws[1]),C_I(x0),C_I(x1),C_R(y0),C_R(y1),n_rws-2,rws+2); // z.i =x.i*y.r
  iradd_mul_ws(C_I(cws[0]),C_I(cws[1]),C_R(x0),C_R(x1),C_I(y0),C_I(y1),n_rws-2,rws+2); // z.i+=x.r*y.i
  iccopy(z0,z1,cws[0],cws[1]);
}

/**
 @brief 掛け算 [z0,z1]=[x0,x1]*[y0,y1]
*/
void icmul_r1(cmulti *z0, cmulti *z1, rmulti *x0, rmulti *x1, cmulti *y0, cmulti *y1)
{
  int prec;
  cmulti *a0=NULL,*a1=NULL;
  prec=MAX2(cget_prec(z0),cget_prec(z1));
  CA(a0,prec); CA(a1,prec);
  irmul(C_R(a0),C_R(a1),x0,x1,C_R(y0),C_R(y1)); // z.r =x.r*y.r
  irmul(C_I(a0),C_I(a1),x0,x1,C_I(y0),C_I(y1)); // z.i =x.r*y.i
  iccopy(z0,z1,a0,a1);
  CF(a0); CF(a1);
}

/**
 @brief 掛け算 [z0,z1]=[x0,x1]*[y0,y1]
*/
void icmul_r2(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, rmulti *y0, rmulti *y1)
{
  int prec;
  cmulti *a0=NULL,*a1=NULL;
  prec=MAX2(cget_prec(z0),cget_prec(z1));
  CA(a0,prec); CA(a1,prec);
  irmul(C_R(a0),C_R(a1),C_R(x0),C_R(x1),y0,y1); // z.r =x.r*y.r
  irmul(C_I(a0),C_I(a1),C_I(x0),C_I(x1),y0,y1); // z.i =x.i*y.r
  iccopy(z0,z1,a0,a1);
  CF(a0); CF(a1);
}

/**
 @brief 掛け算 [z0,z1]=[x0,x1]*[y,y]
*/
void icmul_d(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, double y)
{
  irmul_d(C_R(z0),C_R(z1),C_R(x0),C_R(x1),y);
  irmul_d(C_I(z0),C_I(z1),C_I(x0),C_I(x1),y);
}

/**
 @brief 積の加算 [z0,z1]+=[x0,x1]*[y0,y1]
*/
void icadd_mul(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, cmulti *y0, cmulti *y1)
{
  int prec;
  cmulti *a0=NULL,*a1=NULL;
  prec=MAX2(cget_prec(z0),cget_prec(z1));
  CA(a0,prec); CA(a1,prec);
  icmul(a0,a1,x0,x1,y0,y1); // a=x*y
  icadd(z0,z1,z0,z1,a0,a1); // z=z+x*y
  CF(a0); CF(a1);
}

/**
 @brief 積の加算 [z0,z1]+=[x0,x1]*[y0,y1]
*/
void icadd_mul_r1(cmulti *z0, cmulti *z1, rmulti *x0, rmulti *x1, cmulti *y0, cmulti *y1)
{
  int prec;
  cmulti *a0=NULL,*a1=NULL;
  prec=MAX2(cget_prec(z0),cget_prec(z1));
  CA(a0,prec); CA(a1,prec);
  icmul_r1(a0,a1,x0,x1,y0,y1); // a=x*y
  icadd(z0,z1,z0,z1,a0,a1);    // z=z+x*y
  CF(a0); CF(a1);
}

/**
 @brief 積の加算 [z0,z1]+=[x0,x1]*[y0,y1]
*/
void icadd_mul_r2(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, rmulti *y0, rmulti *y1)
{
  int prec;
  cmulti *a0=NULL,*a1=NULL;
  prec=MAX2(cget_prec(z0),cget_prec(z1));
  CA(a0,prec); CA(a1,prec);
  icmul_r2(a0,a1,x0,x1,y0,y1); // a=x*y
  icadd(z0,z1,z0,z1,a0,a1);    // z=z+x*y
  CF(a0); CF(a1);
}

/**
 @brief 積の減算 [z0,z1]-=[x0,x1]*[y0,y1]
*/
void icsub_mul(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, cmulti *y0, cmulti *y1)
{
  int prec;
  cmulti *a0=NULL,*a1=NULL;
  prec=MAX2(cget_prec(z0),cget_prec(z1));
  CA(a0,prec); CA(a1,prec);
  icmul(a0,a1,x0,x1,y0,y1); // a=x*y
  icsub(z0,z1,z0,z1,a0,a1); // z=z-x*y
  CF(a0); CF(a1);
}

/**
 @brief 積の減算 [z0,z1]-=[x0,x1]*[y0,y1]
*/
void icsub_mul_ws(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, cmulti *y0, cmulti *y1, int n_rws, rmulti **rws, int n_cws, cmulti **cws)
{
  if(n_rws<8){ ERROR_EXIT("Error `n_rws=%d<8' in icsub_mul_ws().\n",n_rws); }
  if(n_cws<4){ ERROR_EXIT("Error `n_cws=%d<4' in icsub_mul_ws().\n",n_cws); }
  icmul_ws(cws[0],cws[1],x0,x1,y0,y1,n_rws,rws,n_cws-2,cws+2); // a=x*y
  icsub_ws(z0,z1,z0,z1,cws[0],cws[1],n_rws,rws);               // z=z-x*y
}

/**
 @brief 共役な掛け算 [z0,z1]=conj([x0,x1])*[y0,y1]
*/
void icdot(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, cmulti *y0, cmulti *y1)
{
  int prec;
  cmulti *a0=NULL,*a1=NULL;
  prec=MAX2(cget_prec(z0),cget_prec(z1));
  CA(a0,prec); CA(a1,prec);
  irmul    (C_R(a0),C_R(a1),C_R(x0),C_R(x1),C_R(y0),C_R(y1));  // z.r =x.r*y.r
  iradd_mul(C_R(a0),C_R(a1),C_I(x0),C_I(x1),C_I(y0),C_I(y1));  // z.r+=x.i*y.i
  irmul    (C_I(a0),C_I(a1),C_R(x0),C_R(x1),C_I(y0),C_I(y1));  // z.i =x.r*y.i
  irsub_mul(C_I(a0),C_I(a1),C_I(x0),C_I(x1),C_R(y0),C_R(y1));  // z.i-=x.i*y.r
  iccopy(z0,z1,a0,a1);
  CF(a0); CF(a1);
}

/**
 @brief 共役な積の加算 [z0,z1]+=conj([x0,x1])*[y0,y1]
*/
void icadd_dot(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, cmulti *y0, cmulti *y1)
{
  int prec;
  cmulti *a0=NULL,*a1=NULL;
  prec=MAX2(cget_prec(z0),cget_prec(z1));
  CA(a0,prec); CA(a1,prec);
  icdot(a0,a1,x0,x1,y0,y1); // a=conj(x)*y
  icadd(z0,z1,z0,z1,a0,a1); // z=z+x*y
  CF(a0); CF(a1);
}

/**
 @brief 共役な積の減算 [z0,z1]-=conj([x0,x1])*[y0,y1]
*/
void icsub_dot(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, cmulti *y0, cmulti *y1)
{
  int prec;
  cmulti *a0=NULL,*a1=NULL;
  prec=MAX2(cget_prec(z0),cget_prec(z1));
  CA(a0,prec); CA(a1,prec);
  icdot(a0,a1,x0,x1,y0,y1); // a=conj(x)*y
  icsub(z0,z1,z0,z1,a0,a1); // z=z-x*y
  CF(a0); CF(a1);
}

/**
 @brief 割り算 [z0,z1]=[x0,x1]/[y0,y1]
*/
void icdiv(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, cmulti *y0, cmulti *y1)
{
  int prec;
  rmulti *den0=NULL,*den1=NULL;
  cmulti *a0=NULL,*a1=NULL;
  prec=MAX2(cget_prec(z0),cget_prec(z1));
  RA(den0,prec); RA(den1,prec); CA(a0,prec); CA(a1,prec);
  icabs2(den0,den1,y0,y1);         // den=|y|^2
  icdot(a0,a1,y0,y1,x0,x1);        // z=conj(y)*x
  icdiv_r2(a0,a1,a0,a1,den0,den1); // z/=den
  iccopy(z0,z1,a0,a1);
  RF(den0); RF(den1); CF(a0); CF(a1);
}

/**
 @brief 割り算 [z0,z1]=[x0,x1]/[y0,y1]
*/
void icdiv_r1(cmulti *z0, cmulti *z1, rmulti *x0, rmulti *x1, cmulti *y0, cmulti *y1)
{
  int prec;
  rmulti *den0=NULL,*den1=NULL;
  cmulti *a0=NULL,*a1=NULL;
  prec=MAX2(cget_prec(z0),cget_prec(z1));
  RA(den0,prec); RA(den1,prec); CA(a0,prec); CA(a1,prec);
  icabs2(den0,den1,y0,y1);                       // den=|y|^2
  irmul(C_R(a0),C_R(a1),C_R(y0),C_R(y1),x0,x1);  // z.r =x.r*y.r
  irmul(C_I(a0),C_I(a1),C_I(y0),C_I(y1),x0,x1);  // z.i =x.i*y.r
  icconj(a0,a1,a0,a1);                           // z=-(x.r*y.r+i*x.i*y.r)
  icdiv_r2(a0,a1,a0,a1,den0,den1);               // z=-(x.r*y.r+i*x.i*y.r)/den
  iccopy(z0,z1,a0,a1);
  RF(den0); RF(den1); CF(a0); CF(a1);
}

/**
 @brief 割り算 [z0,z1]=[x0,x1]/[y0,y1]
*/
void icdiv_r2(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, rmulti *y0, rmulti *y1)
{
  irdiv(C_R(z0),C_R(z1),C_R(x0),C_R(x1),y0,y1);
  irdiv(C_I(z0),C_I(z1),C_I(x0),C_I(x1),y0,y1);
}

/**
 @brief 差の絶対値 [z0,z1]=abs([x0,x1]-[y0,y1])
*/
void icabs_sub(rmulti *z0, rmulti *z1, cmulti *x0, cmulti *x1, cmulti *y0, cmulti *y1)
{
  cmulti *a0=NULL,*a1=NULL;
  CAr(a0,z0); CAr(a1,z1);
  icsub(a0,a1,x0,x1,y0,y1);
  icabs(z0,z1,a0,a1);
  CF(a0); CF(a1);
}

/**
 @brief べき乗 [y0,y1]=[x0,x1]^n
 */
void icpow_si(cmulti *y0, cmulti *y1, cmulti *x0, cmulti *x1, long n)
{
  cmulti *a0=NULL,*a1=NULL;
  int prec;
  long i;
  prec=MAX2(cget_prec(y0),cget_prec(y1));
  CA(a0,prec); CA(a1,prec);
  if(n==0){ icset_d(a0,a1,1,1); }
  else if(n>0){
    iccopy(a0,a1,x0,x1);
    for(i=1; i<n; i++){ icmul(a0,a1,a0,a1,x0,x1); }
  }else if(n<0){
    iccopy(a0,a1,x0,x1);
    for(i=1; i<(-n); i++){ icmul(a0,a1,a0,a1,x0,x1); }
    icinv(a0,a1,a0,a1);
  }
  iccopy(y0,y1,a0,a1);
  CF(a0); CF(a1);
}

/**
 @brief 極形式の取得 [z0,z1]=[r0,r1]*exp(i*[theta0,theta1])
 */
void icget_polar(rmulti *r0, rmulti *r1, rmulti *theta0, rmulti *theta1, cmulti *z0, cmulti *z1)
{
  icabs(r0,r1,z0,z1);
  iratan2(theta0,theta1,C_I(z0),C_I(z1),C_R(z0),C_R(z1));
}

/**
 @brief 極形式で値を定め得る [z0,z1]=[r0,r1]*exp(i*[theta0,theta1])
 */
void icset_polar(cmulti *z0, cmulti *z1, rmulti *r0, rmulti *r1, rmulti *theta0, rmulti *theta1)
{
  int prec;
  rmulti *s0=NULL,*s1=NULL,*c0=NULL,*c1=NULL;
  prec=MAX2(cget_prec(z0),cget_prec(z1));
  RA(s0,prec); RA(s1,prec); RA(c0,prec); RA(c1,prec);
  ircos(s0,s1,theta0,theta1); irmul(C_R(z0),C_R(z1),r0,r1,s0,s1); // z.r=r*cos(theta)
  irsin(c0,c1,theta0,theta1); irmul(C_I(z0),C_I(z1),r0,r1,c0,c1); // z.i=r*sin(theta)
  RF(s0); RF(s1); RF(c0); RF(c1);
}

/** @} */

//////////////////////////////////////////////////////

/** @name icmulti型の判定に関する関数 */
/** @{ */

/**
 @brief 区間の内包の判定 [y0,y1] in [x0,x1], すなわち (x0<=y0 && y1 <=x1).
 */
int icin(cmulti *y0, cmulti *y1, cmulti *x0, cmulti *x1)
{
  return (irin(C_R(y0),C_R(y1),C_R(x0),C_R(x1)) && irin(C_I(y0),C_I(y1),C_I(x0),C_I(x1)));
}

/**
 @brief 区間の内包の判定 [y0,y1] in [-abs(x),abs(x)], すなわち (-abs(x)<=y0 && y1 <=abs(x)).
 */
int icin_pm(cmulti *y0, cmulti *y1, cmulti *x)
{
  return (irin_pm(C_R(y0),C_R(y1),C_R(x)) && irin_pm(C_I(y0),C_I(y1),C_I(x)));
}

/** @} */

//EOF
