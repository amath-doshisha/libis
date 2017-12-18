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

/** @name 基本操作 */
/** @{ */


/**
 @brief 倍精度複素数の設定 [y0,y1]=[xr+i*xi,xr*i*xi].
 */
int icset_z(cmulti *y0, cmulti *y1, dcomplex x)
{
  int e=0;
  e+=irset_d(C_R(y0),C_R(y1),Z_R(x));
  e+=irset_d(C_I(y0),C_I(y1),Z_I(x));
  return e;
}

/**
 @brief 倍精度複素数の設定 [y0,y1]=[xr+i*xi,xr*i*xi].
 */
int icset_dd(cmulti *y0, cmulti *y1, double xr, double xi)
{
  int e=0;
  e+=irset_d(C_R(y0),C_R(y1),xr);
  e+=irset_d(C_I(y0),C_I(y1),xi);
  return e;
}

/**
 @brief 倍精度実数の設定 [y0,y1]=[x,x].
 */
int icset_d(cmulti *y0, cmulti *y1, double x)
{
  int e=0;
  e+=irset_d(C_R(y0),C_R(y1),x);
  e+=irset_d(C_I(y0),C_I(y1),0);
  return e;
}

/**
 @brief bigint型から[z0,z1]へ型変換.
 */
int icset_bigint(cmulti *z0, cmulti *z1, bigint *x)
{
  int e=0;
  e+=irset_bigint(C_R(z0),C_R(z1),x);
  e+=irset_d(C_I(z0),C_I(z1),0);
  return e;
}


/**
 @brief コピー [y0,y1]=[x0,x1].
 */
int iccopy(cmulti *y0, cmulti *y1, cmulti *x0, cmulti *x1)
{
  int e=0;
  e+=ircopy(C_R(y0),C_R(y1),C_R(x0),C_R(x1));
  e+=ircopy(C_I(y0),C_I(y1),C_I(x0),C_I(x1));
  return e;
}

/** @} */

//////////////////////////////////////////////////////

/** @name 四則演算 */
/** @{ */

/**
 @brief 複素共役 [y0,y1]=conj([x0,x1])
*/
int icconj(cmulti *y0, cmulti *y1, cmulti *x0, cmulti *x1)
{
  int e=0;
  e+=ircopy(C_R(y0),C_R(y1),C_R(x0),C_R(x1));
  e+=irneg(C_I(y0),C_I(y1),C_I(x0),C_I(x1));
  return e;
}

/**
 @brief 符号の反転 [y0,y1]=-[x0,x1].
 */
int icneg(cmulti *y0, cmulti *y1, cmulti *x0, cmulti *x1)
{
  int e=0;
  e+=irneg(C_R(y0),C_R(y1),C_R(x0),C_R(x1));
  e+=irneg(C_I(y0),C_I(y1),C_I(x0),C_I(x1));
  return e;
}

/**
 @brief 符号の正負 [y0,y1]=[-abs(real(x))-i*abs(imag(x))],abs(real(x))+i*abs(imag(x))].
 */
int icpm(cmulti *y0, cmulti *y1, cmulti *x)
{
  int e=0;
  e+=irpm(C_R(y0),C_R(y1),C_R(x));
  e+=irpm(C_I(y0),C_I(y1),C_I(x));
  return e;
}


/**
 @brief 足し算 [z0,z1]=[x0,x1]+[y0,y1]
 */
int icadd(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, cmulti *y0, cmulti *y1)
{
  int e=0;
  e+=iradd(C_R(z0),C_R(z1),C_R(x0),C_R(x1),C_R(y0),C_R(y1)); // z.r=x.r+y.r
  e+=iradd(C_I(z0),C_I(z1),C_I(x0),C_I(x1),C_I(y0),C_I(y1)); // z.i=x.i+y.i
  return e;
}

/**
 @brief 足し算 [z0,z1]=[x0,x1]+[y0,y1]
 */
int icadd_r(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, rmulti *y0, rmulti *y1)
{
  int e=0;
  e+=iradd (C_R(z0),C_R(z1),C_R(x0),C_R(x1),y0,y1); // z.r=x.r+y.r
  e+=ircopy(C_I(z0),C_I(z1),C_I(x0),C_I(x1));       // z.i=x.i
  return e;
}

/**
 @brief 足し算 [z0,z1]=[x0,x1]+[-y,y]
 */
int icadd_pm(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, cmulti *y)
{
  int e=0;
  e+=iradd_pm(C_R(z0),C_R(z1),C_R(x0),C_R(x1),C_R(y)); // z.r=x.r+y.r
  e+=iradd_pm(C_I(z0),C_I(z1),C_I(x0),C_I(x1),C_I(y)); // z.i=x.i+y.i
  return e;
}

/**
 @brief 区間の中心 [m-r,m+r]=[x0,x1]
 */
int icmid(cmulti *mid, cmulti *x0, cmulti *x1)
{
  int e=0;
  e+=irmid(C_R(mid),C_R(x0),C_R(x1));
  e+=irmid(C_I(mid),C_I(x0),C_I(x1));
  return e;
}

/**
 @brief 区間の半径 [m-r,m+r]=[x0,x1]
 */
int icrad(cmulti *rad, cmulti *x0, cmulti *x1)
{
  int e=0;
  e+=irrad(C_R(rad),C_R(x0),C_R(x1));
  e+=irrad(C_I(rad),C_I(x0),C_I(x1));
  return e;
}

/**
 @brief 区間の中心と半径 [m-r,m+r]=[x0,x1]
 */
int icmr(cmulti *mid, cmulti *rad, cmulti *x0, cmulti *x1)
{
  int e=0;
  e+=irmr(C_R(mid),C_R(rad),C_R(x0),C_R(x1));
  e+=irmr(C_I(mid),C_I(rad),C_I(x0),C_I(x1));
  return e;
}

/**
 @brief 引き算 [z0,z1]=[x0,x1]-[y0,y1]
 */
int icsub(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, cmulti *y0, cmulti *y1)
{
  int e=0;
  e+=irsub(C_R(z0),C_R(z1),C_R(x0),C_R(x1),C_R(y0),C_R(y1)); // z.r=x.r-y.r
  e+=irsub(C_I(z0),C_I(z1),C_I(x0),C_I(x1),C_I(y0),C_I(y1)); // z.i=x.i-y.i
  return e;
}
//編集済み
/**
 @brief 引き算 [z0,z1]=[x0,x1]-[y0,y1]
 */
int icsub_ws(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, cmulti *y0, cmulti *y1, int *rwss, rmulti **rws)
{
  if(*rwss<2){ ERROR_EXIT("Error `rwss=%d<2' in icsub__ws().\n",*rwss); }
  int e=0;
  e+=irsub_ws(C_R(z0),C_R(z1),C_R(x0),C_R(x1),C_R(y0),C_R(y1),rwss,rws); // z.r=x.r-y.r
  e+=irsub_ws(C_I(z0),C_I(z1),C_I(x0),C_I(x1),C_I(y0),C_I(y1),rwss,rws); // z.i=x.i-y.i
  return e;
}

/**
 @brief 引き算 [z0,z1]=[x0,x1]-[y0,y1]
 */
int icsub_r(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, rmulti *y0, rmulti *y1)
{
  int e=0;
  e+=irsub(C_R(z0),C_R(z1),C_R(x0),C_R(x1),y0,y1); // z.r=x.r-y.r
  e+=ircopy(C_I(z0),C_I(z1),C_I(x0),C_I(x1));       // z.i=x.i
  return e;
}

/**
 @brief 引き算 [z0,z1]=[x0,x1]-[y0,y1]
 */
int irsub_c(cmulti *z0, cmulti *z1, rmulti *x0, rmulti *x1, cmulti *y0, cmulti *y1)
{
  int e=0;
  e+=irsub(C_R(z0),C_R(z1),x0,x1,C_R(y0),C_R(y1)); // z.r=x.r-y.r
  e+=irneg(C_I(z0),C_I(z1),C_I(y0),C_I(y1));       //z.i=-(y.i)
  return e;
}
//ここまで
/**
 @brief 引き算 [z0,z1]=[x0,x1]-[y,y]
 */
int icsub_d2(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, double y)
{
  int e=0;
  e+=irsub_d2(C_R(z0),C_R(z1),C_R(x0),C_R(x1),y); // z.r=x.r-y
  e+=ircopy(C_I(z0),C_I(z1),C_I(x0),C_I(x1));     // z.i=x.i
  return e;
}

/**
 @brief 掛け算 [z0,z1]=[x0,x1]*[y0,y1]
*/
int icmul(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, cmulti *y0, cmulti *y1)
{
  int prec,e=0;
  cmulti *a0=NULL,*a1=NULL;
  prec=MAX2(cget_prec(z0),cget_prec(z1));
  CA(a0,prec); CA(a1,prec);
  e+=irmul    (C_R(a0),C_R(a1),C_R(x0),C_R(x1),C_R(y0),C_R(y1)); // z.r =x.r*y.r
  e+=irsub_mul(C_R(a0),C_R(a1),C_I(x0),C_I(x1),C_I(y0),C_I(y1)); // z.r-=x.i*y.i
  e+=irmul    (C_I(a0),C_I(a1),C_I(x0),C_I(x1),C_R(y0),C_R(y1)); // z.i =x.i*y.r
  e+=iradd_mul(C_I(a0),C_I(a1),C_R(x0),C_R(x1),C_I(y0),C_I(y1)); // z.i+=x.r*y.i
  e+=iccopy(z0,z1,a0,a1);
  CF(a0); CF(a1);
  return e;
}
//編集済み

/**
 @brief 掛け算 [z0,z1]=[x0,x1]*[y0,y1]
*/
int icmul_ws(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, cmulti *y0, cmulti *y1, int *rwss, rmulti **rws, int *cwss, cmulti **cws)
{
  if(*rwss<8){ ERROR_EXIT("Error `rwss=%d<8' in icmul_ws().\n",*rwss); }
  if(*cwss<2){ ERROR_EXIT("Error `cwss=%d<2' in icmul_ws().\n",*cwss); }
  int e=0;
  *cwss=*cwss-2;
  e+=irmul_ws    (C_R(cws[0]),C_R(cws[1]),C_R(x0),C_R(x1),C_R(y0),C_R(y1), rwss, &rws[2]); // z.r =x.r*y.r
  e+=irsub_mul_ws(C_R(cws[0]),C_R(cws[1]),C_I(x0),C_I(x1),C_I(y0),C_I(y1), rwss, &rws[2]); // z.r-=x.i*y.i
  e+=irmul_ws    (C_I(cws[0]),C_I(cws[1]),C_I(x0),C_I(x1),C_R(y0),C_R(y1), rwss, &rws[2]); // z.i =x.i*y.r
  e+=iradd_mul_ws(C_I(cws[0]),C_I(cws[1]),C_R(x0),C_R(x1),C_I(y0),C_I(y1), rwss, &rws[2]); // z.i+=x.r*y.i
  e+=iccopy(z0,z1,cws[0],cws[1]);
  *cwss=*cwss+2;
  return e;
}


/**
 @brief 掛け算 [z0,z1]=[x0,x1]*[y0,y1]
*/
int icmul_r1(cmulti *z0, cmulti *z1, rmulti *x0, rmulti *x1, cmulti *y0, cmulti *y1)
{
  int prec,e=0;
  cmulti *a0=NULL,*a1=NULL;
  prec=MAX2(cget_prec(z0),cget_prec(z1));
  CA(a0,prec); CA(a1,prec);
  e+=irmul    (C_R(a0),C_R(a1),x0,x1,C_R(y0),C_R(y1)); // z.r =x.r*y.r
  e+=irmul    (C_I(a0),C_I(a1),x0,x1,C_I(y0),C_I(y1)); // z.i =x.r*y.i
  e+=iccopy(z0,z1,a0,a1);
  CF(a0); CF(a1);
  return e;
}
/**
 @brief 掛け算 [z0,z1]=[x0,x1]*[y0,y1]
*/
int icmul_r2(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, rmulti *y0, rmulti *y1)
{
  int prec,e=0;
  cmulti *a0=NULL,*a1=NULL;
  prec=MAX2(cget_prec(z0),cget_prec(z1));
  CA(a0,prec); CA(a1,prec);
  e+=irmul    (C_R(a0),C_R(a1),C_R(x0),C_R(x1),y0,y1); // z.r =x.r*y.r
  e+=irmul    (C_I(a0),C_I(a1),C_I(x0),C_I(x1),y0,y1); // z.i =x.i*y.r
  e+=iccopy(z0,z1,a0,a1);
  CF(a0); CF(a1);
  return e;
}
//ここまで
/**
 @brief 掛け算 [z0,z1]=[x0,x1]*[y,y]
*/
int icmul_d(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, double y)
{
  int e=0;
  e+=irmul_d(C_R(z0),C_R(z1),C_R(x0),C_R(x1),y);
  e+=irmul_d(C_I(z0),C_I(z1),C_I(x0),C_I(x1),y);
  return e;
}

/**
 @brief 共役な掛け算 [z0,z1]=conj([x0,x1])*[y0,y1]
*/
int icdot(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, cmulti *y0, cmulti *y1)
{
  int prec,e=0;
  cmulti *a0=NULL,*a1=NULL;
  prec=MAX2(cget_prec(z0),cget_prec(z1));
  CA(a0,prec); CA(a1,prec);
  e+=irmul    (C_R(a0),C_R(a1),C_R(x0),C_R(x1),C_R(y0),C_R(y1));  // z.r =x.r*y.r
  e+=iradd_mul(C_R(a0),C_R(a1),C_I(x0),C_I(x1),C_I(y0),C_I(y1));  // z.r+=x.i*y.i
  e+=irmul    (C_I(a0),C_I(a1),C_R(x0),C_R(x1),C_I(y0),C_I(y1));  // z.i =x.r*y.i
  e+=irsub_mul(C_I(a0),C_I(a1),C_I(x0),C_I(x1),C_R(y0),C_R(y1));  // z.i-=x.i*y.r
  e+=iccopy(z0,z1,a0,a1);
  CF(a0); CF(a1);
  return e;
}

/**
 @brief 割り算 [z0,z1]=[x0,x1]/[y0,y1]
*/
int icdiv(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, cmulti *y0, cmulti *y1)
{
  int prec,e=0;
  rmulti *den0=NULL,*den1=NULL;
  cmulti *a0=NULL,*a1=NULL;
  prec=MAX2(cget_prec(z0),cget_prec(z1));
  RA(den0,prec); RA(den1,prec); CA(a0,prec); CA(a1,prec);
  e+=icabs2(den0,den1,y0,y1);         // den=|y|^2
  e+=icdot(a0,a1,y0,y1,x0,x1);        // z=conj(y)*x
  e+=icdiv_r2(a0,a1,a0,a1,den0,den1); // z/=den
  e+=iccopy(z0,z1,a0,a1);
  RF(den0); RF(den1); CF(a0); CF(a1);
  return e;
}

//編集済み
/**
 @brief 割り算 [z0,z1]=[x0,x1]/[y0,y1]
*/
int icdiv_r1(cmulti *z0, cmulti *z1, rmulti *x0, rmulti *x1, cmulti *y0, cmulti *y1)
{
  int prec,e=0;
  rmulti *den0=NULL,*den1=NULL;
  cmulti *a0=NULL,*a1=NULL;
  prec=MAX2(cget_prec(z0),cget_prec(z1));
  RA(den0,prec); RA(den1,prec); CA(a0,prec); CA(a1,prec);
  e+=icabs2(den0,den1,y0,y1);         // den=|y|^2
  e+=irmul(C_R(a0),C_R(a1),C_R(y0),C_R(y1),x0,x1);  // z.r =x.r*y.r
  e+=irmul(C_I(a0),C_I(a1),C_I(y0),C_I(y1),x0,x1);  // z.i =x.i*y.r
  e+=icconj(a0,a1,a0,a1);
  e+=icdiv_r2(a0,a1,a0,a1,den0,den1); // z/=den
  e+=iccopy(z0,z1,a0,a1);
  RF(den0); RF(den1); CF(a0); CF(a1);
  return e;
}
//ここまで
/**
 @brief 割り算 [z0,z1]=[x0,x1]/[y0,y1]
*/
int icdiv_r2(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, rmulti *y0, rmulti *y1)
{
  int e=0;
  e+=irdiv(C_R(z0),C_R(z1),C_R(x0),C_R(x1),y0,y1);
  e+=irdiv(C_I(z0),C_I(z1),C_I(x0),C_I(x1),y0,y1);
  return e;
}

/**
 @brief 逆数 [z0,z1]=[1,1]/[x0,x1]
*/
int icinv(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1)
{
  int prec,e=0;
  rmulti *den0=NULL,*den1=NULL;
  cmulti *a0=NULL,*a1=NULL;
  prec=MAX2(cget_prec(z0),cget_prec(z1));
  RA(den0,prec); RA(den1,prec); CA(a0,prec); CA(a1,prec);
  e+=icabs2(den0,den1,x0,x1);         // den=|x|^2
  e+=icconj(a0,a1,x0,x1);             // z=conj(x)
  e+=icdiv_r2(a0,a1,a0,a1,den0,den1); // z/=den
  e+=iccopy(z0,z1,a0,a1);
  RF(den0); RF(den1); CF(a0); CF(a1);
  return e;
}

/**
 @brief 積の加算 [z0,z1]+=[x0,x1]*[y0,y1]
*/
int icadd_mul(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, cmulti *y0, cmulti *y1)
{
  int prec,e=0;
  cmulti *a0=NULL,*a1=NULL;
  prec=MAX2(cget_prec(z0),cget_prec(z1));
  CA(a0,prec); CA(a1,prec);
  e+=icmul(a0,a1,x0,x1,y0,y1); // a=x*y
  e+=icadd(z0,z1,z0,z1,a0,a1); // z=z+x*y
  CF(a0); CF(a1);
  return e;
}
//編集済み

/**
 @brief 積の加算 [z0,z1]+=[x0,x1]*[y0,y1]
*/
int icadd_mul_r1(cmulti *z0, cmulti *z1, rmulti *x0, rmulti *x1, cmulti *y0, cmulti *y1)
{
  int prec,e=0;
  cmulti *a0=NULL,*a1=NULL;
  prec=MAX2(cget_prec(z0),cget_prec(z1));
  CA(a0,prec); CA(a1,prec);
  e+=icmul_r1(a0,a1,x0,x1,y0,y1); // a=x*y
  e+=icadd(z0,z1,z0,z1,a0,a1); // z=z+x*y
  CF(a0); CF(a1);
  return e;
}

/**
 @brief 積の加算 [z0,z1]+=[x0,x1]*[y0,y1]
*/
int icadd_mul_r2(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, rmulti *y0, rmulti *y1)
{
  int prec,e=0;
  cmulti *a0=NULL,*a1=NULL;
  prec=MAX2(cget_prec(z0),cget_prec(z1));
  CA(a0,prec); CA(a1,prec);
  e+=icmul_r2(a0,a1,x0,x1,y0,y1); // a=x*y
  e+=icadd(z0,z1,z0,z1,a0,a1); // z=z+x*y
  CF(a0); CF(a1);
  return e;
}
//ここまで
/**
 @brief 積の減算 [z0,z1]-=[x0,x1]*[y0,y1]
*/
int icsub_mul(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, cmulti *y0, cmulti *y1)
{
  int prec,e=0;
  cmulti *a0=NULL,*a1=NULL;
  prec=MAX2(cget_prec(z0),cget_prec(z1));
  CA(a0,prec); CA(a1,prec);
  e+=icmul(a0,a1,x0,x1,y0,y1); // a=x*y
  e+=icsub(z0,z1,z0,z1,a0,a1); // z=z-x*y
  CF(a0); CF(a1);
  return e;
}

//追加
/**
 @brief 積の減算 [z0,z1]-=[x0,x1]*[y0,y1]
*/
int icsub_mul_ws(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, cmulti *y0, cmulti *y1, int *rwss, rmulti **rws, int *cwss, cmulti **cws)
{
  int e=0;
  if(*rwss<8){ ERROR_EXIT("Error `rwss=%d<8' in icsub_mul_ws().\n",*rwss); }
  if(*cwss<4){ ERROR_EXIT("Error `cwss=%d<4' in icsub_mul_ws().\n",*cwss); }
  *cwss=*cwss-2;
  e+=icmul_ws(cws[0],cws[1],x0,x1,y0,y1,rwss,rws,cwss,&cws[2]); // a=x*y
  e+=icsub_ws(z0,z1,z0,z1,cws[0],cws[1],rwss,rws); // z=z-x*y
  *cwss=*cwss+2;
  return e;
}

//ここまで

/**
 @brief 共役な積の加算 [z0,z1]+=conj([x0,x1])*[y0,y1]
*/
int icadd_dot(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, cmulti *y0, cmulti *y1)
{
  int prec,e=0;
  cmulti *a0=NULL,*a1=NULL;
  prec=MAX2(cget_prec(z0),cget_prec(z1));
  CA(a0,prec); CA(a1,prec);
  e+=icdot(a0,a1,x0,x1,y0,y1); // a=conj(x)*y
  e+=icadd(z0,z1,z0,z1,a0,a1); // z=z+x*y
  CF(a0); CF(a1);
  return e;
}

/**
 @brief 共役な積の減算 [z0,z1]-=conj([x0,x1])*[y0,y1]
*/
int icsub_dot(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, cmulti *y0, cmulti *y1)
{
  int prec,e=0;
  cmulti *a0=NULL,*a1=NULL;
  prec=MAX2(cget_prec(z0),cget_prec(z1));
  CA(a0,prec); CA(a1,prec);
  e+=icdot(a0,a1,x0,x1,y0,y1); // a=conj(x)*y
  e+=icsub(z0,z1,z0,z1,a0,a1); // z=z-x*y
  CF(a0); CF(a1);
  return e;
}

/**
 @brief 絶対値の平方 [y0,y1]=abs([x0,x1])^2
*/
int icabs2(rmulti *y0, rmulti *y1, cmulti *x0, cmulti *x1)
{
  int e=0;
  e+=irmul    (y0,y1,C_R(x0),C_R(x1),C_R(x0),C_R(x1)); // y=x.r*x.r
  e+=iradd_mul(y0,y1,C_I(x0),C_I(x1),C_I(x0),C_I(x1)); // y+=x.i*x.i
  return e;
}

/**
 @brief 絶対値の平方の加算 [y0,y1]+=abs([x0,x1])^2
*/
int icadd_abs2(rmulti *y0, rmulti *y1, cmulti *x0, cmulti *x1)
{
  int prec,e=0;
  rmulti *a0=NULL,*a1=NULL;
  prec=MAX2(rget_prec(y0),rget_prec(y1));
  RA(a0,prec); RA(a1,prec);
  e+=icabs2(a0,a1,x0,x1);      // a=abs(x)^2
  e+=iradd(y0,y1,y0,y1,a0,a1); // y=y+a
  RF(a0); RF(a1);
  return e;
}

/** @} */

//////////////////////////////////////////////////////

/** @name 数学関数 */
/** @{ */

/**
 @brief 絶対値 [y0,y1]=abs([x0,x1])
*/
int icabs(rmulti *y0, rmulti *y1, cmulti *x0, cmulti *x1)
{
  if(cis_nan(x0) || cis_nan(x1)){rset_nan(y0); rset_nan(y1); return 0;
  }else{
  int prec,e=0;
  cmulti *a0=NULL,*a1=NULL;
  prec=MAX2(rget_prec(y0),rget_prec(y1));
  CA(a0,prec); CA(a1,prec);
  if(rget_sgn(C_R(x0))>=0 && rget_sgn(C_I(x0))>=0 && rget_sgn(C_R(x1))>0 && rget_sgn(C_I(x1))>0){ // 区間が第一象限にある場合
    e+=icabs2(y0,y1,x0,x1);
    e+=irsqrt(y0,y1,y0,y1);
  }else if(rget_sgn(C_R(x0))<0 && rget_sgn(C_I(x0))>=0 && rget_sgn(C_R(x1))<=0 && rget_sgn(C_I(x1))>0){  // 区間が第二象限にある場合
    e+=icabs2(y0,y1,x0,x1);
    e+=irsqrt(y0,y1,y0,y1);
  }else if(rget_sgn(C_R(x0))<0 && rget_sgn(C_I(x0))<0 && rget_sgn(C_R(x1))<=0 && rget_sgn(C_I(x1))<=0){ // 区間が第三象限にある場合
    e+=icabs2(y0,y1,x0,x1);
    e+=irsqrt(y0,y1,y0,y1);
  }else if(rget_sgn(C_R(x0))>=0 && rget_sgn(C_I(x0))<0 && rget_sgn(C_R(x1))>0 && rget_sgn(C_I(x1))<=0){  // 区間が第四象限にある場合
    e+=icabs2(y0,y1,x0,x1);
    e+=irsqrt(y0,y1,y0,y1);
  }else if(rget_sgn(C_R(x0))>=0 && rget_sgn(C_I(x0))<0 && rget_sgn(C_R(x1))>0 && rget_sgn(C_I(x1))>0){  // 区間が正のx軸をまたぐ場合
    if(rabs_gt(C_I(x1),C_I(x0))!=0){
      e+=abs(mpfr_set(C_R(a0),C_R(x0),MPFR_RNDD));
      e+=abs(mpfr_set_d(C_I(a0),0,MPFR_RNDD));
      e+=icabs2(y0,y1,a0,x1);
      e+=irsqrt(y0,y1,y0,y1);
    }else{ 
      e+=abs(mpfr_set(C_R(a1),C_R(x1),MPFR_RNDU));
      e+=abs(mpfr_set_d(C_I(a1),0,MPFR_RNDU));
      e+=icabs2(y0,y1,x0,a1);
      e+=irsqrt(y0,y1,y0,y1);
    }
  }else if(rget_sgn(C_R(x0))<0 && rget_sgn(C_I(x0))>=0 && rget_sgn(C_R(x1))>0 && rget_sgn(C_I(x1))>0){  //  区間が正のy軸をまたぐ場合
    if(rabs_gt(C_R(x1),C_R(x0))!=0){
      e+=abs(mpfr_set_d(C_R(a0),0,MPFR_RNDD));
      e+=abs(mpfr_set(C_I(a0),C_I(x0),MPFR_RNDD));
      e+=icabs2(y0,y1,a0,x1);
      e+=irsqrt(y0,y1,y0,y1);
    }else{
      e+=abs(mpfr_set_d(C_R(a1),0,MPFR_RNDU));
      e+=abs(mpfr_set(C_I(a1),C_I(x1),MPFR_RNDU));
      e+=icabs2(y0,y1,x0,a1);
      e+=irsqrt(y0,y1,y0,y1);
    }
  }else if(rget_sgn(C_R(x0))<0 && rget_sgn(C_I(x0))<0 && rget_sgn(C_R(x1))<=0 && rget_sgn(C_I(x1))>0){  //  区間が負のx軸をまたぐ場合
    if(rabs_gt(C_I(x1),C_I(x0))!=0) { 
      e+=abs(mpfr_set(C_R(a0),C_R(x0),MPFR_RNDD));
      e+=abs(mpfr_set_d(C_I(a0),0,MPFR_RNDD));
      e+=icabs2(y0,y1,a0,x1);
      e+=irsqrt(y0,y1,y0,y1);
    }else{
      e+=abs(mpfr_set(C_R(a1),C_R(x1),MPFR_RNDU));
      e+=abs(mpfr_set_d(C_I(a1),0,MPFR_RNDU));
      e+=icabs2(y0,y1,x0,a1);
      e+=irsqrt(y0,y1,y0,y1);
    }
  }else if(rget_sgn(C_R(x0))<0 && rget_sgn(C_I(x0))<0 && rget_sgn(C_R(x1))>0 && rget_sgn(C_I(x1))<=0){  //  区間が負のy軸をまたぐ場合
    if(rabs_gt(C_R(x1),C_R(x0))!=0){
      e+=abs(mpfr_set_d(C_R(a0),0,MPFR_RNDD));
      e+=abs(mpfr_set(C_I(a0),C_I(x0),MPFR_RNDD));
      e+=icabs2(y0,y1,a0,x1);
      e+=irsqrt(y0,y1,y0,y1);
    }else{
      e+=abs(mpfr_set_d(C_R(a1),0,MPFR_RNDU));
      e+=abs(mpfr_set(C_I(a1),C_I(x1),MPFR_RNDU));
      e+=icabs2(y0,y1,x0,a1);
      e+=irsqrt(y0,y1,y0,y1);
    }
  }else if(rget_sgn(C_R(x0))<0 && rget_sgn(C_I(x0))<0 && rget_sgn(C_R(x1))>0 && rget_sgn(C_I(x1))>0){  // 区間が原点を含む場合
    if(rabs_ge(C_R(x1),C_R(x0))!=0 && rabs_ge(C_I(x1),C_I(x0))!=0){
      e+=abs(mpfr_set_d(C_R(a0),0,MPFR_RNDD));
      e+=abs(mpfr_set_d(C_I(a0),0,MPFR_RNDD));      
      e+=icabs2(y0,y1,a0,x1);
      e+=irsqrt(y0,y1,y0,y1); 
    }else if(rabs_ge(C_R(x0),C_R(x1))!=0 && rabs_ge(C_I(x1),C_I(x0))!=0){
      e+=abs(mpfr_set(C_R(a0),C_R(x0),MPFR_RNDD));
      e+=abs(mpfr_set_d(C_I(a0),0,MPFR_RNDD));      
      e+=abs(mpfr_set_d(C_R(a1),0,MPFR_RNDU));
      e+=abs(mpfr_set(C_I(a1),C_I(x1),MPFR_RNDU));      
      e+=icabs2(y0,y1,a0,a1);
      e+=irsqrt(y0,y1,y0,y1); 
    }else if(rabs_ge(C_R(x0),C_R(x1))!=0 && rabs_ge(C_I(x0),C_I(x1))!=0){
      e+=abs(mpfr_set_d(C_R(a1),0,MPFR_RNDU));
      e+=abs(mpfr_set_d(C_I(a1),0,MPFR_RNDU));      
      e+=icabs2(y0,y1,x0,a1);
      e+=irsqrt(y0,y1,y0,y1); 
    }else if(rabs_ge(C_R(x1),C_R(x0))!=0 && rabs_ge(C_I(x0),C_I(x1))!=0){
      e+=abs(mpfr_set_d(C_R(a0),0,MPFR_RNDD));
      e+=abs(mpfr_set(C_I(a0),C_I(x0),MPFR_RNDD));      
      e+=abs(mpfr_set(C_R(a1),C_R(x1),MPFR_RNDU));
      e+=abs(mpfr_set_d(C_I(a1),0,MPFR_RNDU));      
      e+=icabs2(y0,y1,a0,a1);
      e+=irsqrt(y0,y1,y0,y1); 
    }
  }else{
    CF(a0); CF(a1);
    printf("\x1b[31m"); printf("Error! Special case!!\n");
    mpfr_printf("Input [doun  up] = [%.4Re + %.4Re   %.4Re + %.4Re]\n",C_R(x0),C_I(x0),C_R(x1),C_I(x1));
    printf("\x1b[39m");
    exit(0);
  }
  CF(a0); CF(a1);
  return e;
  }
}

//追加
/**
 @brief 絶対値 [y0,y1]=abs([x0,x1])
*/
int icabs_ws(rmulti *y0, rmulti *y1, cmulti *x0, cmulti *x1, int *cwss, cmulti **cws)
{
  if(*cwss<2){ ERROR_EXIT("Error `cwss=%d<2' in icabs_ws().\n",*cwss); }
  if(cis_nan(x0) || cis_nan(x1)){rset_nan(y0); rset_nan(y1); return 0;
  }else{ int e=0;
  if(rget_sgn(C_R(x0))>=0 && rget_sgn(C_I(x0))>=0 && rget_sgn(C_R(x1))>0 && rget_sgn(C_I(x1))>0){ // 区間が第一象限にある場合
    e+=icabs2(y0,y1,x0,x1);
    e+=irsqrt(y0,y1,y0,y1);
  }else if(rget_sgn(C_R(x0))<0 && rget_sgn(C_I(x0))>=0 && rget_sgn(C_R(x1))<=0 && rget_sgn(C_I(x1))>0){  // 区間が第二象限にある場合
    e+=icabs2(y0,y1,x0,x1);
    e+=irsqrt(y0,y1,y0,y1);
  }else if(rget_sgn(C_R(x0))<0 && rget_sgn(C_I(x0))<0 && rget_sgn(C_R(x1))<=0 && rget_sgn(C_I(x1))<=0){ // 区間が第三象限にある場合
    e+=icabs2(y0,y1,x0,x1);
    e+=irsqrt(y0,y1,y0,y1);
  }else if(rget_sgn(C_R(x0))>=0 && rget_sgn(C_I(x0))<0 && rget_sgn(C_R(x1))>0 && rget_sgn(C_I(x1))<=0){  // 区間が第四象限にある場合
    e+=icabs2(y0,y1,x0,x1);
    e+=irsqrt(y0,y1,y0,y1);
  }else if(rget_sgn(C_R(x0))>=0 && rget_sgn(C_I(x0))<0 && rget_sgn(C_R(x1))>0 && rget_sgn(C_I(x1))>0){  // 区間が正のx軸をまたぐ場合
    if(rabs_gt(C_I(x1),C_I(x0))!=0){
      e+=abs(mpfr_set(C_R(cws[0]),C_R(x0),MPFR_RNDD));
      e+=abs(mpfr_set_d(C_I(cws[0]),0,MPFR_RNDD));
      e+=icabs2(y0,y1,cws[0],x1);
      e+=irsqrt(y0,y1,y0,y1);
    }else{ 
      e+=abs(mpfr_set(C_R(cws[1]),C_R(x1),MPFR_RNDU));
      e+=abs(mpfr_set_d(C_I(cws[1]),0,MPFR_RNDU));
      e+=icabs2(y0,y1,x0,cws[1]);
      e+=irsqrt(y0,y1,y0,y1);
    }
  }else if(rget_sgn(C_R(x0))<0 && rget_sgn(C_I(x0))>=0 && rget_sgn(C_R(x1))>0 && rget_sgn(C_I(x1))>0){  //  区間が正のy軸をまたぐ場合
    if(rabs_gt(C_R(x1),C_R(x0))!=0){
      e+=abs(mpfr_set_d(C_R(cws[0]),0,MPFR_RNDD));
      e+=abs(mpfr_set(C_I(cws[0]),C_I(x0),MPFR_RNDD));
      e+=icabs2(y0,y1,cws[0],x1);
      e+=irsqrt(y0,y1,y0,y1);
    }else{
      e+=abs(mpfr_set_d(C_R(cws[1]),0,MPFR_RNDU));
      e+=abs(mpfr_set(C_I(cws[1]),C_I(x1),MPFR_RNDU));
      e+=icabs2(y0,y1,x0,cws[1]);
      e+=irsqrt(y0,y1,y0,y1);
    }
  }else if(rget_sgn(C_R(x0))<0 && rget_sgn(C_I(x0))<0 && rget_sgn(C_R(x1))<=0 && rget_sgn(C_I(x1))>0){  //  区間が負のx軸をまたぐ場合
    if(rabs_gt(C_I(x1),C_I(x0))!=0) { 
      e+=abs(mpfr_set(C_R(cws[0]),C_R(x0),MPFR_RNDD));
      e+=abs(mpfr_set_d(C_I(cws[0]),0,MPFR_RNDD));
      e+=icabs2(y0,y1,cws[0],x1);
      e+=irsqrt(y0,y1,y0,y1);
    }else{
      e+=abs(mpfr_set(C_R(cws[1]),C_R(x1),MPFR_RNDU));
      e+=abs(mpfr_set_d(C_I(cws[1]),0,MPFR_RNDU));
      e+=icabs2(y0,y1,x0,cws[1]);
      e+=irsqrt(y0,y1,y0,y1);
    }
  }else if(rget_sgn(C_R(x0))<0 && rget_sgn(C_I(x0))<0 && rget_sgn(C_R(x1))>0 && rget_sgn(C_I(x1))<=0){  //  区間が負のy軸をまたぐ場合
    if(rabs_gt(C_R(x1),C_R(x0))!=0){
      e+=abs(mpfr_set_d(C_R(cws[0]),0,MPFR_RNDD));
      e+=abs(mpfr_set(C_I(cws[0]),C_I(x0),MPFR_RNDD));
      e+=icabs2(y0,y1,cws[0],x1);
      e+=irsqrt(y0,y1,y0,y1);
    }else{
      e+=abs(mpfr_set_d(C_R(cws[1]),0,MPFR_RNDU));
      e+=abs(mpfr_set(C_I(cws[1]),C_I(x1),MPFR_RNDU));
      e+=icabs2(y0,y1,x0,cws[1]);
      e+=irsqrt(y0,y1,y0,y1);
    }
  }else if(rget_sgn(C_R(x0))<0 && rget_sgn(C_I(x0))<0 && rget_sgn(C_R(x1))>0 && rget_sgn(C_I(x1))>0){  // 区間が原点を含む場合
    if(rabs_ge(C_R(x1),C_R(x0))!=0 && rabs_ge(C_I(x1),C_I(x0))!=0){
      e+=abs(mpfr_set_d(C_R(cws[0]),0,MPFR_RNDD));
      e+=abs(mpfr_set_d(C_I(cws[0]),0,MPFR_RNDD));      
      e+=icabs2(y0,y1,cws[0],x1);
      e+=irsqrt(y0,y1,y0,y1); 
    }else if(rabs_ge(C_R(x0),C_R(x1))!=0 && rabs_ge(C_I(x1),C_I(x0))!=0){
      e+=abs(mpfr_set(C_R(cws[0]),C_R(x0),MPFR_RNDD));
      e+=abs(mpfr_set_d(C_I(cws[0]),0,MPFR_RNDD));      
      e+=abs(mpfr_set_d(C_R(cws[1]),0,MPFR_RNDU));
      e+=abs(mpfr_set(C_I(cws[1]),C_I(x1),MPFR_RNDU));      
      e+=icabs2(y0,y1,cws[0],cws[1]);
      e+=irsqrt(y0,y1,y0,y1); 
    }else if(rabs_ge(C_R(x0),C_R(x1))!=0 && rabs_ge(C_I(x0),C_I(x1))!=0){
      e+=abs(mpfr_set_d(C_R(cws[1]),0,MPFR_RNDU));
      e+=abs(mpfr_set_d(C_I(cws[1]),0,MPFR_RNDU));      
      e+=icabs2(y0,y1,x0,cws[1]);
      e+=irsqrt(y0,y1,y0,y1); 
    }else if(rabs_ge(C_R(x1),C_R(x0))!=0 && rabs_ge(C_I(x0),C_I(x1))!=0){
      e+=abs(mpfr_set_d(C_R(cws[0]),0,MPFR_RNDD));
      e+=abs(mpfr_set(C_I(cws[0]),C_I(x0),MPFR_RNDD));      
      e+=abs(mpfr_set(C_R(cws[1]),C_R(x1),MPFR_RNDU));
      e+=abs(mpfr_set_d(C_I(cws[1]),0,MPFR_RNDU));      
      e+=icabs2(y0,y1,cws[0],cws[1]);
      e+=irsqrt(y0,y1,y0,y1); 
    }
  }else{
    printf("\x1b[31m"); printf("Error! Special case!!\n");
    mpfr_printf("Input [doun  up] = [%.4Re + %.4Re   %.4Re + %.4Re]\n",C_R(x0),C_I(x0),C_R(x1),C_I(x1));
    printf("\x1b[39m");
    exit(0);
  }
  return e;
  }
}

//ここまで

/**
 @brief 差の絶対値 [z0,z1]=abs([x0,x1]-[y0,y1])
*/
int icabs_sub(rmulti *z0, rmulti *z1, cmulti *x0, cmulti *x1, cmulti *y0, cmulti *y1)
{
  int e=0;
  cmulti *a0=NULL, *a1=NULL;
  CAr(a0,z0); CAr(a1,z1);
  e+=icsub(a0,a1,x0,x1,y0,y1);
  e+=icabs(z0,z1,a0,a1);
  CF(a0); CF(a1);
  return e;
}

/**
 @brief べき乗 [y0,y1]=[x0,x1]^n
 */
int icpow_si(cmulti *y0, cmulti *y1, cmulti *x0, cmulti *x1, long n)
{
  cmulti *a0=NULL,*a1=NULL;
  int prec,e=0;
  long i;
  prec=MAX2(cget_prec(y0),cget_prec(y1));
  CA(a0,prec); CA(a1,prec);
  if(n==0){ e+=icset_d(a0,a1,1); }
  else if(n>0){
    e+=iccopy(a0,a1,x0,x1);
    for(i=1; i<n; i++){ e+=icmul(a0,a1,a0,a1,x0,x1); }
  }else if(n<0){
    e+=iccopy(a0,a1,x0,x1);
    for(i=1; i<(-n); i++){ e+=icmul(a0,a1,a0,a1,x0,x1); }
    e+=icinv(a0,a1,a0,a1);
  }
  e+=iccopy(y0,y1,a0,a1);
  CF(a0); CF(a1);
  return e;
}

/**
 @brief 極形式の取得 [z0,z1]=[r0,r1]*exp(i*[theta0,theta1])
 */
int icget_polar(rmulti *r0, rmulti *r1, rmulti *theta0, rmulti *theta1, cmulti *z0, cmulti *z1)
{
  int e=0;
  e+=icabs(r0,r1,z0,z1);
  e+=iratan2(theta0,theta1,C_I(z0),C_I(z1),C_R(z0),C_R(z1));
  return e;
}

/**
 @brief 極形式で値を定め得る [z0,z1]=[r0,r1]*exp(i*[theta0,theta1])
 */
int icset_polar(cmulti *z0, cmulti *z1, rmulti *r0, rmulti *r1, rmulti *theta0, rmulti *theta1)
{
  int prec,e=0;
  rmulti *s0=NULL,*s1=NULL,*c0=NULL,*c1=NULL;
  prec=MAX2(cget_prec(z0),cget_prec(z1));
  RA(s0,prec); RA(s1,prec); RA(c0,prec); RA(c1,prec);
  e+=ircos(s0,s1,theta0,theta1); e+=irmul(C_R(z0),C_R(z1),r0,r1,s0,s1); // z.r=r*cos(theta)
  e+=irsin(c0,c1,theta0,theta1); e+=irmul(C_I(z0),C_I(z1),r0,r1,c0,c1); // z.i=r*sin(theta)
  RF(s0); RF(s1); RF(c0); RF(c1);
  return e;  
}

/**
 @brief 平方根 [y0,y1]=sqrt([x0,x1])
*/
int icsqrt(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1)
{
  int prec,e=0;
  rmulti *r0=NULL,*r1=NULL,*theta0=NULL,*theta1=NULL;
  prec=MAX2(cget_prec(z0),cget_prec(z1));
  RA(r0,prec); RA(r1,prec); RA(theta0,prec); RA(theta1,prec);
  e+=icget_polar(r0,r1,theta0,theta1,x0,x1);   // x=r*exp(i*theta)
  e+=irsqrt(r0,r1,r0,r1);                      // r=sqrt(r)
  e+=irmul_d(theta0,theta1,theta0,theta1,0.5); // theta=theta/2
  e+=icset_polar(z0,z1,r0,r1,theta0,theta1);   // z=r^y*cos(theta*y)+i*r^y*sin(theta*y)
  RF(r0); RF(r1); RF(theta0); RF(theta1);
  return e;
}

/**
 @brief 指数関数 [y0,y1]=exp([x0,x1])
*/
int icexp(cmulti *y0, cmulti *y1, cmulti *x0, cmulti *x1)
{
  int prec,e=0;
  rmulti *r0=NULL,*r1=NULL,*theta0=NULL,*theta1=NULL;
  prec=MAX2(cget_prec(y0),cget_prec(y1));
  RA(r0,prec); RA(r1,prec); RA(theta0,prec); RA(theta1,prec);
  e+=irexp(r0,r1,C_R(x0),C_R(x1));
  e+=ircopy(theta0,theta1,C_I(x0),C_I(x1));
  e+=icset_polar(y0,y1,r0,r1,theta0,theta1); // y=exp(x.r)*(cos(x.i)+i*sin(x.i))
  RF(r0); RF(r1); RF(theta0); RF(theta1);
  return e;
}

/**
 @brief 対数関数 [y0,y1]=log([x0,x1])
*/
int iclog(cmulti *y0, cmulti *y1, cmulti *x0, cmulti *x1)
{
  int prec,e=0;
  rmulti *r0=NULL,*r1=NULL,*theta0=NULL,*theta1=NULL;
  prec=MAX2(cget_prec(y0),cget_prec(y1));
  RA(r0,prec); RA(r1,prec); RA(theta0,prec); RA(theta1,prec);
  e+=icget_polar(r0,r1,theta0,theta1,x0,x1); // x=r*exp(i*theta)
  e+=irlog(C_R(y0),C_R(y1),r0,r1);           // y=log(r)+i*theta
  e+=ircopy(C_I(y0),C_I(y1),theta0,theta1);  // y=log(r)+i*theta
  RF(r0); RF(r1); RF(theta0); RF(theta1);
  return e;
}

/**
 @brief 三角関数 [y0,y1]=sin([x0,x1])
 */
int icsin(cmulti *y0, cmulti *y1, cmulti *x0, cmulti *x1)
{
  int prec,e=0;
  rmulti *c0=NULL,*c1=NULL,*ch0=NULL,*ch1=NULL,*s0=NULL,*s1=NULL,*sh0=NULL,*sh1=NULL;
  prec=MAX2(cget_prec(y0),cget_prec(y1));
  RA(c0,prec); RA(c1,prec); RA(ch0,prec); RA(ch1,prec); RA(s0,prec); RA(s1,prec); RA(sh0,prec); RA(sh1,prec);
  e+=irsin(s0,s1,C_R(x0),C_R(x1)); e+=ircosh(ch0,ch1,C_I(x0),C_I(x1)); e+=irmul(C_R(y0),C_R(y1),s0,s1,ch0,ch1); // z.r=sin(x.r)*cosh(x.i)
  e+=ircos(c0,c1,C_R(x0),C_R(x1)); e+=irsinh(sh0,sh1,C_I(x0),C_I(x1)); e+=irmul(C_I(y0),C_I(y1),c0,c1,sh0,sh1); // z.i=cos(x.r)*sinh(x.i)
  RF(c0); RF(c1); RF(ch0); RF(ch1); RF(s0); RF(s1); RF(sh0); RF(sh1);
  return e;
}

/**
 @brief 三角関数 [y0,y1]=cos([x0,x1])
 */
int iccos(cmulti *y0, cmulti *y1, cmulti *x0, cmulti *x1)
{
  int prec,e=0;
  rmulti *c0=NULL,*ch0=NULL,*s0=NULL,*sh0=NULL,*c1=NULL,*ch1=NULL,*s1=NULL,*sh1=NULL;
  prec=MAX2(cget_prec(y0),cget_prec(y1));
  RA(c0,prec); RA(ch0,prec); RA(s0,prec); RA(sh0,prec); RA(c1,prec); RA(ch1,prec); RA(s1,prec); RA(sh1,prec);
  e+=ircos(c0,c1,C_R(x0),C_R(x1)); e+=ircosh(ch0,ch1,C_I(x0),C_I(x1)); e+=irmul(C_R(y0),C_R(y1),c0,c1,ch0,ch1);                                            // z.r=cos(x.r)*cosh(x.i)
  e+=irsin(s0,s1,C_R(x0),C_R(x1)); e+=irsinh(sh0,sh1,C_I(x0),C_I(x1)); e+=irmul(C_I(y0),C_I(y1),s0,s1,sh0,sh1); e+=irneg(C_I(y0),C_I(y1),C_I(y0),C_I(y1)); // z.i=-sin(x.r)*sinh(x.i)
  RF(c0); RF(ch0); RF(s0); RF(sh0); RF(c1); RF(ch1); RF(s1); RF(sh1);
  return e;
}

/**
 @brief 三角関数 [y0,y1]=tan([x0,x1])
 */
int ictan(cmulti *y0, cmulti *y1, cmulti *x0, cmulti *x1)
{
  int prec,e=0;
  cmulti *c0=NULL,*s0=NULL,*c1=NULL,*s1=NULL;
  prec=MAX2(cget_prec(y0),cget_prec(y1));
  CA(c0,prec); CA(s0,prec); CA(c1,prec); CA(s1,prec);
  e+=iccos(c0,c1,x0,x1); e+=icsin(s0,s1,x0,x1); e+=icdiv(y0,y1,s0,s1,c0,c1); // y=s/c
  CF(c0); CF(s0); CF(c1); CF(s1);
  return e;
}

/**
 @brief 逆三角関数 [y0,y1]=asin([x0,x1])
*/
int icasin(cmulti *y0, cmulti *y1, cmulti *x0, cmulti *x1)
{
  int prec,e=0;
  cmulti *a0=NULL,*a1=NULL,*b0=NULL,*b1=NULL;
  prec=MAX2(cget_prec(y0),cget_prec(y1));
  CA(a0,prec); CA(a1,prec); CA(b0,prec); CA(b1,prec);
  e+=icmul(a0,a1,x0,x1,x0,x1); // a=x^2 
  e+=icset_dd(b0,b1,1,0);      // b=1
  e+=icsub(a0,a1,b0,b1,a0,a1); // a=1-x^2
  e+=icsqrt(a0,a1,a0,a1);      // a=sqrt(1-x^2)
  e+=icset_dd(b0,b1,0,1);      // b=I
  e+=icmul(b0,b1,b0,b1,x0,x1); // b=I*x
  e+=icadd(a0,a1,a0,a1,b0,b1); // a=I*x+sqrt(1-x^2)
  e+=iclog(a0,a1,a0,a1);       // a=log(I*x+sqrt(1-x^2))
  e+=icset_dd(b0,b1,0,-1);     // b=-I
  e+=icmul(a0,a1,a0,a1,b0,b1); // a=-I*log(I*x+sqrt(1-x^2))
  e+=iccopy(y0,y1,a0,a1);      // y=-I*log(I*x+sqrt(1-x^2))
  CF(a0); CF(a1); CF(b0); CF(b1);
  return e;
}

/**
 @brief 逆三角関数 [y0,y1]=acos([x0,x1])
*/
int icacos(cmulti *y0, cmulti *y1, cmulti *x0, cmulti *x1)
{
  int prec,e=0;
  cmulti *a0=NULL,*a1=NULL,*b0=NULL,*b1=NULL;
  prec=MAX2(cget_prec(y0),cget_prec(y1));
  CA(a0,prec); CA(a1,prec); CA(b0,prec); CA(b1,prec);
  e+=icmul(a0,a1,x0,x1,x0,x1);  // a=x^2 
  e+=icset_dd(b0,b1,-1,0);      // b=-1
  e+=icadd(a0,a1,a0,a1,b0,b1);  // a=x^2-1
  e+=icsqrt(a0,a1,a0,a1);       // a=sqrt(x^2-1)
  e+=icadd(a0,a1,a0,a1,x0,x1);  // a=x+sqrt(x^2-1)
  e+=iclog(a0,a1,a0,a1);        // a=log(x+sqrt(x^2-1))
  e+=icset_dd(b0,b1,0,-1);      // b=-I
  e+=icmul(a0,a1,a0,a1,b0,b1);  // a=-I*log(I*x+sqrt(1-x^2))
  e+=iccopy(y0,y1,a0,a1);       // y=-I*log(I*x+sqrt(1-x^2))
  CF(a0); CF(a1); CF(b0); CF(b1);
  return e;
}

/**
 @brief 逆三角関数 [y0,y1]=atan([x0,x1])
*/
int icatan(cmulti *y0, cmulti *y1, cmulti *x0, cmulti *x1)
{
  int prec,e=0;
  cmulti *a0=NULL,*a1=NULL,*b0=NULL,*b1=NULL;
  prec=MAX2(cget_prec(y0),cget_prec(y1));
  CA(a0,prec); CA(a1,prec); CA(b0,prec); CA(b1,prec);
  e+=icset_dd(b0,b1,1,0);      // b=1
  e+=icadd(a0,a1,x0,x1,b0,b1); // a=1+x
  e+=icsub(b0,b1,b0,b1,x0,x1); // b=1-x
  e+=icdiv(a0,a1,a0,a1,b0,b1); // a=(1+x)/(1-x)
  e+=iclog(a0,a1,a0,a1);       // a=log((1+x)/(1-x))
  e+=icset_dd(b0,b1,0,-0.5);   // b=-I/2
  e+=icmul(a0,a1,a0,a1,b0,b1); // a=(-I/2)*log((1+x)/(1-x))
  e+=iccopy(y0,y1,a0,a1);      // y=(-I/2)*log((1+x)/(1-x))
  CF(a0); CF(a1); CF(b0); CF(b1);
  return e;
}

/**
 @brief 双曲線関数 [y0,y1]=sinh([x0,x1])
*/
int icsinh(cmulti *y0, cmulti *y1, cmulti *x0, cmulti *x1)
{
  int prec,e=0;
  rmulti *c0=NULL,*ch0=NULL,*s0=NULL,*sh0=NULL,*c1=NULL,*ch1=NULL,*s1=NULL,*sh1=NULL;
  prec=MAX2(cget_prec(y0),cget_prec(y1));
  RA(c0,prec); RA(ch0,prec); RA(s0,prec); RA(sh0,prec); RA(c1,prec); RA(ch1,prec); RA(s1,prec); RA(sh1,prec);
  e+=irsinh(sh0,sh1,C_R(x0),C_R(x1)); e+=ircos(c0,c1,C_I(x0),C_I(x1)); e+=irmul(C_R(y0),C_R(y1),sh0,sh1,c0,c1); // z.r=sinh(x.r)*cos(x.i)
  e+=ircosh(ch0,ch1,C_R(x0),C_R(x1)); e+=irsin(s0,s1,C_I(x0),C_I(x1)); e+=irmul(C_I(y0),C_I(y1),ch0,ch1,s0,s1); // z.i=cosh(x.r)*sin(x.i)
  RF(c0); RF(ch0); RF(s0); RF(sh0); RF(c1); RF(ch1); RF(s1); RF(sh1);
  return e;
}

/**
 @brief 双曲線関数 [y0,y1]=cosh([x0,x1])
*/
int iccosh(cmulti *y0, cmulti *y1, cmulti *x0, cmulti *x1)
{
  int prec,e=0;
  rmulti *c0=NULL,*ch0=NULL,*s0=NULL,*sh0=NULL,*c1=NULL,*ch1=NULL,*s1=NULL,*sh1=NULL;
  prec=MAX2(cget_prec(y0),cget_prec(y1));
  RA(c0,prec); RA(ch0,prec); RA(s0,prec); RA(sh0,prec); RA(c1,prec); RA(ch1,prec); RA(s1,prec); RA(sh1,prec);
  e+=ircosh(ch0,ch1,C_R(x0),C_R(x1)); e+=ircos(c0,c1,C_I(x0),C_I(x1)); e+=irmul(C_R(y0),C_R(y1),ch0,ch1,c0,c1); // z.r=cosh(x.r)*cos(x.i)
  e+=irsinh(sh0,sh1,C_R(x0),C_R(x1)); e+=irsin(s0,s1,C_I(x0),C_I(x1)); e+=irmul(C_I(y0),C_I(y1),sh0,sh1,s0,s1); // z.i=sinh(x.r)*sin(x.i)
  RF(c0); RF(ch0); RF(s0); RF(sh0); RF(c1); RF(ch1); RF(s1); RF(sh1);
  return e;
}

/**
 @brief 双曲線関数 [y0,y1]=tanh([x0,x1])
*/
int ictanh(cmulti *y0, cmulti *y1, cmulti *x0, cmulti *x1)
{
  int prec,e=0;
  cmulti *ch0=NULL,*sh0=NULL,*ch1=NULL,*sh1=NULL;
  prec=MAX2(cget_prec(y0),cget_prec(y1));
  CA(ch0,prec); CA(sh0,prec); CA(ch1,prec); CA(sh1,prec);
  e+=iccosh(ch0,ch1,x0,x1);
  e+=icsinh(sh0,sh1,x0,x1);
  e+=icdiv(y0,y1,sh0,sh1,ch0,ch1); // y=sh/ch
  CF(ch0); CF(sh0); CF(ch1); CF(sh1);
  return e;
}

/**
 @beif 逆双曲線関数 [y0,y1]=asinh([x0,x1])
*/
int icasinh(cmulti *y0, cmulti *y1, cmulti *x0, cmulti *x1)
{
  int prec,e=0;
  cmulti *a0=NULL,*a1=NULL,*b0=NULL,*b1=NULL;
  prec=MAX2(cget_prec(y0),cget_prec(y1));
  CA(a0,prec); CA(a1,prec); CA(b0,prec); CA(b1,prec);
  e+=icset_dd(b0,b1,1,0);       // b=1
  e+=icmul(a0,a1,x0,x1,x0,x1);  // a=x^2
  e+=icadd(a0,a1,a0,a1,b0,b1);  // a=x^2+1
  e+=icsqrt(a0,a1,a0,a1);       // a=sqrt(x^2+1)
  e+=icadd(a0,a1,a0,a1,x0,x1);  // a=x+sqrt(x^2+1)
  e+=iclog(a0,a1,a0,a1);        // a=log(x+sqrt(x^2+1))
  e+=iccopy(y0,y1,a0,a1);       // y=log(x+sqrt(x^2+1))
  CF(a0); CF(a1); CF(b0); CF(b1);
  return e;
}

/**
 @brief 逆双曲線関数 [y0,y1]=acosh([x0,x1])
*/
int icacosh(cmulti *y0, cmulti *y1, cmulti *x0, cmulti *x1)
{
  int prec,e=0;
  cmulti *a0=NULL,*a1=NULL,*b0=NULL,*b1=NULL;
  prec=MAX2(cget_prec(y0),cget_prec(y1));
  CA(a0,prec); CA(a1,prec); CA(b0,prec); CA(b1,prec);
  e+=icset_dd(b0,b1,-1,0);      // b=-1
  e+=icmul(a0,a1,x0,x1,x0,x1);  // a=x^2
  e+=icadd(a0,a1,a0,a1,b0,b1);  // a=x^2-1
  e+=icsqrt(a0,a1,a0,a1);       // a=sqrt(x^2-1)
  e+=icadd(a0,a1,a0,a1,x0,x1);  // a=x+sqrt(x^2-1)
  e+=iclog(a0,a1,a0,a1);        // a=log(x+sqrt(x^2-1))
  e+=iccopy(y0,y1,a0,a1);       // y=log(x+sqrt(x^2-1))
  CF(a0); CF(a1); CF(b0); CF(b1);
  return e;
}

/**
 @brief 逆双曲線関数 [y0,y1]=atanh([x0,x1])
*/
int icatanh(cmulti *y0, cmulti *y1, cmulti *x0, cmulti *x1)
{
  int prec,e=0;
  cmulti *a0=NULL,*a1=NULL,*b0=NULL,*b1=NULL;
  prec=MAX2(cget_prec(y0),cget_prec(y1));
  CA(a0,prec); CA(a1,prec); CA(b0,prec); CA(b1,prec);
  e+=icset_dd(b0,b1,1,0);      // b=1
  e+=icadd(a0,a1,x0,x1,b0,b1); // a=1+x
  e+=icsub(b0,b1,b0,b1,x0,x1); // b=1-x
  e+=icdiv(a0,a1,a0,a1,b0,b1); // a=(1+x)/(1-x)
  e+=iclog(a0,a1,a0,a1);       // a=log((1+x)/(1-x))
  e+=icset_dd(b0,b1,0.5,0);    // b=0.5
  e+=icmul(a0,a1,a0,a1,b0,b1); // a=(1/2)*log((1+x)/(1-x))
  e+=iccopy(y0,y1,a0,a1);      // y=(1/2)*log((1+x)/(1-x))
  CF(a0); CF(a1); CF(b0); CF(b1);
  return e;
}

/** @} */

//////////////////////////////////////////////////////

/** @name 判定 */
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
