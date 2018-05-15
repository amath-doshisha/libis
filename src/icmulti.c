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
void icset_c(cmulti *y0, cmulti *y1, cmulti *x0, cmulti *x1)
{
  irset_r(C_R(y0),C_R(y1),C_R(x0),C_R(x1));
  irset_r(C_I(y0),C_I(y1),C_I(x0),C_I(x1));
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
  irset_r  (C_R(y0),C_R(y1),x0,x1);
  irset_d(C_I(y0),C_I(y1),0, 0);
}

/**
 @brief 倍精度複素数の設定 [y0,y1]=[x0_r+i*x0_i,xr1_*i*x1_i].
 */
void icset_dddd(cmulti *y0, cmulti *y1, double x0_r, double x0_i, double x1_r, double x1_i)
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
 @brief コピー [y0,y1]=[x0_r,x1_r]+i*[x0_i,x1_i]
 */
void icset_rrrr(cmulti *y0, cmulti *y1, rmulti *x0_r, rmulti *x0_i, rmulti *x1_r, rmulti *x1_i)
{
  irset_r(C_R(y0),C_R(y1),x0_r,x1_r);
  irset_r(C_I(y0),C_I(y1),x0_i,x1_i);
}

/**
 @brief コピー [y0,y1]=[x0,x1]. rmultiをcmultiにキャスト
 */
void icset_rc(cmulti *y0, cmulti *y1, rmulti *x0, cmulti *x1)
{
  irset_r  (C_R(y0),C_R(y1),x0,C_R(x1));
  irset_dr(C_I(y0),C_I(y1),.0,C_I(x1));
}

/**
 @brief コピー [y0,y1]=[x0,x1]. rmultiをcmultiにキャスト
 */
void icset_cr(cmulti *y0, cmulti *y1, cmulti *x0, rmulti *x1)
{
  irset_r  (C_R(y0),C_R(y1),C_R(x0),x1);
  irset_rd(C_I(y0),C_I(y1),C_I(x0),.0);
}

/** @} */

//////////////////////////////////////////////////////////////

/** @name icmulti型の型変換に関する関数 */
/** @{ */

/** @} */

//////////////////////////////////////////////////////////////

/** @name icmulti型の１入力の数学演算子 */
/** @{ */

/**
 @brief 複素共役 [y0,y1]=conj([x0,x1])
 */
void icconj_c(cmulti *y0, cmulti *y1, cmulti *x0, cmulti *x1)
{
  irset_r(C_R(y0),C_R(y1),C_R(x0),C_R(x1));
  irneg_r(C_I(y0),C_I(y1),C_I(x0),C_I(x1));
}

/**
 @brief 符号の反転 [y0,y1]=-[x0,x1].
 */
void icneg_c(cmulti *y0, cmulti *y1, cmulti *x0, cmulti *x1)
{
  irneg_r(C_R(y0),C_R(y1),C_R(x0),C_R(x1));
  irneg_r(C_I(y0),C_I(y1),C_I(x0),C_I(x1));
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
void icinv_c(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1)
{
  int prec;
  rmulti *den0=NULL,*den1=NULL;
  cmulti *a0=NULL,*a1=NULL;
  prec=MAX2(cget_prec(z0),cget_prec(z1));
  RA(den0,prec); RA(den1,prec); CA(a0,prec); CA(a1,prec);
  irabs2_c(den0,den1,x0,x1);         // den=|x|^2
  icconj_c(a0,a1,x0,x1);             // z=conj(x)
  icdiv_cr(a0,a1,a0,a1,den0,den1); // z/=den
  icset_c(z0,z1,a0,a1);
  RF(den0); RF(den1); CF(a0); CF(a1);
}

/**
 @brief 絶対値の平方 [y0,y1]=abs([x0,x1])^2
 */
void irabs2_c(rmulti *y0, rmulti *y1, cmulti *x0, cmulti *x1)
{
  irmul_rr    (y0,y1,C_R(x0),C_R(x1),C_R(x0),C_R(x1)); // y=x.r*x.r
  iradd_mul_rr(y0,y1,C_I(x0),C_I(x1),C_I(x0),C_I(x1)); // y+=x.i*x.i
}

/**
 @brief 絶対値の平方の加算 [y0,y1]+=abs([x0,x1])^2
 */
void iradd_abs2_c(rmulti *y0, rmulti *y1, cmulti *x0, cmulti *x1)
{
  int prec;
  rmulti *a0=NULL,*a1=NULL;
  prec=MAX2(rget_prec(y0),rget_prec(y1));
  RA(a0,prec); RA(a1,prec);
  irabs2_c(a0,a1,x0,x1);      // a=abs(x)^2
  iradd_rr(y0,y1,y0,y1,a0,a1); // y=y+a
  RF(a0); RF(a1);
}

/**
 @brief 絶対値 [y0,y1]=abs([x0,x1])
 */
void irabs_c(rmulti *y0, rmulti *y1, cmulti *x0, cmulti *x1)
{
  int prec;
  cmulti *a0=NULL,*a1=NULL;
  if(cis_nan(x0) || cis_nan(x1)){rset_nan(y0); rset_nan(y1); return; }
  prec=MAX2(rget_prec(y0),rget_prec(y1));
  CA(a0,prec); CA(a1,prec);
  irabs_r(C_R(a0),C_R(a1),C_R(x0),C_R(x1));
  irabs_r(C_I(a0),C_I(a1),C_I(x0),C_I(x1));
  irabs2_c(y0,y1,a0,a1);
  irsqrt_r(y0,y1,y0,y1);
  CF(a0); CF(a1);
}

/**
 @brief 絶対値 [y0,y1]=abs([x0,x1])
 */
void irabs_c_ws(rmulti *y0, rmulti *y1, cmulti *x0, cmulti *x1, int n_cws, cmulti **cws)
{
  if(n_cws<2){ ERROR_EXIT("Error `n_cws=%d<2' in icabs_ws().\n",n_cws); }
  if(cis_nan(x0) || cis_nan(x1)){ rset_nan(y0); rset_nan(y1); return; }
  irabs_r(C_R(cws[0]),C_R(cws[1]),C_R(x0),C_R(x1));
  irabs_r(C_I(cws[0]),C_I(cws[1]),C_I(x0),C_I(x1));
  irabs2_c(y0,y1,cws[0],cws[1]);
  irsqrt_r(y0,y1,y0,y1);
}

/** @} */

//////////////////////////////////////////////////////////////

/** @name icmulti型の２入力の数学演算子 */
/** @{ */

/**
 @brief [z0,z1]=[x0,x1]+[y0,y1]
 */
void icadd_cc(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, cmulti *y0, cmulti *y1)
{
  iradd_rr(C_R(z0),C_R(z1),C_R(x0),C_R(x1),C_R(y0),C_R(y1)); // z.r=x.r+y.r
  iradd_rr(C_I(z0),C_I(z1),C_I(x0),C_I(x1),C_I(y0),C_I(y1)); // z.i=x.i+y.i
}

/**
 @brief [z0,z1]=[x0,x1]+[y0,y1]
 */
void icadd_cr(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, rmulti *y0, rmulti *y1)
{
  iradd_rr(C_R(z0),C_R(z1),C_R(x0),C_R(x1),y0,y1); // z.r=x.r+y.r
  irset_r(C_I(z0),C_I(z1),C_I(x0),C_I(x1));       // z.i=x.i
}

/**
 @brief [z0,z1]=[x0,x1]+[y0,y1]
 */
void icadd_rc(cmulti *z0, cmulti *z1, rmulti *x0, rmulti *x1, cmulti *y0, cmulti *y1)
{
  iradd_rr(C_R(z0),C_R(z1),x0,x1,C_R(y0),C_R(y1)); // z.r=x.r+y.r
  irset_r(C_I(z0),C_I(z1),      C_I(y0),C_I(y1)); // z.i=    y.i
}

/**
 @brief [z0,z1]=[x0,x1]+[y0,y1]
 */
void icadd_cz(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, dcomplex y0, dcomplex y1)
{
  iradd_rd(C_R(z0),C_R(z1),C_R(x0),C_R(x1),Z_R(y0),Z_R(y1)); // z.r=x.r+y.r
  iradd_rd(C_I(z0),C_I(z1),C_I(x0),C_I(x1),Z_I(y0),Z_I(y1)); // z.i=x.i+y.i
}

/**
 @brief [z0,z1]=[x0,x1]+[y0,y1]
 */
void icadd_zc(cmulti *z0, cmulti *z1, dcomplex x0, dcomplex x1, cmulti *y0, cmulti *y1)
{
  iradd_dr(C_R(z0),C_R(z1),Z_R(x0),Z_R(x1),C_R(y0),C_R(y1)); // z.r=x.r+y.r
  iradd_dr(C_I(z0),C_I(z1),Z_I(x0),Z_I(x1),C_I(y0),C_I(y1)); // z.i=x.i+y.i
}

/**
 @brief [z0,z1]=[x0,x1]+[y0,y1]
 */
void icadd_cd(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, double y0, double y1)
{
  iradd_rd(C_R(z0),C_R(z1),C_R(x0),C_R(x1),y0,y1); // z.r=x.r+y.r
  irset_r(C_I(z0),C_I(z1),C_I(x0),C_I(x1));       // z.i=x.i
}

/**
 @brief [z0,z1]=[x0,x1]+[y0,y1]
 */
void icadd_dc(cmulti *z0, cmulti *z1, double x0, double x1, cmulti *y0, cmulti *y1)
{
  iradd_dr(C_R(z0),C_R(z1),x0,x1,C_R(y0),C_R(y1)); // z.r=x.r+y.r
  irset_r(C_I(z0),C_I(z1),      C_I(y0),C_I(y1)); // z.i=    y.i
}

/**
 @brief [z0,z1]=[x0,x1]+[y0,y1]
 */
void icadd_rz(cmulti *z0, cmulti *z1, rmulti *x0, rmulti *x1, dcomplex y0, dcomplex y1)
{
  iradd_rd(C_R(z0),C_R(z1),x0,x1,Z_R(y0),Z_R(y1)); // z.r=x+y.r
  irset_d(C_I(z0),C_I(z1),      Z_I(y0),Z_I(y1)); // z.i=  y.i
}

/**
 @brief [z0,z1]=[x0,x1]+[y0,y1]
 */
void icadd_zr(cmulti *z0, cmulti *z1, dcomplex x0, dcomplex x1, rmulti *y0, rmulti *y1)
{
  iradd_dr(C_R(z0),C_R(z1),Z_R(x0),Z_R(x1),y0,y1); // z.r=x.r+y
  irset_d(C_I(z0),C_I(z1),Z_I(x0),Z_I(x1));       // z.i=x.i
}

////////////////////////////////////////////////////////////

/**
 @brief [z0,z1]=[x0,x1]-[y0,y1]
 */
void icsub_cc_ws(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, cmulti *y0, cmulti *y1, int n_rws, rmulti **rws)
{
  if(n_rws<2){ ERROR_EXIT("Error `n_rws=%d<2' in icsub__ws().\n",n_rws); }
  irsub_rr_ws(C_R(z0),C_R(z1),C_R(x0),C_R(x1),C_R(y0),C_R(y1),n_rws,rws); // z.r=x.r-y.r
  irsub_rr_ws(C_I(z0),C_I(z1),C_I(x0),C_I(x1),C_I(y0),C_I(y1),n_rws,rws); // z.i=x.i-y.i
}

/**
 @brief [z0,z1]=[x0,x1]-[y0,y1]
 */
void icsub_cc(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, cmulti *y0, cmulti *y1)
{
  irsub_rr(C_R(z0),C_R(z1),C_R(x0),C_R(x1),C_R(y0),C_R(y1)); // z.r=x.r-y.r
  irsub_rr(C_I(z0),C_I(z1),C_I(x0),C_I(x1),C_I(y0),C_I(y1)); // z.i=x.i-y.i
}


/**
 @brief [z0,z1]=[x0,x1]-[y0,y1]
 */
void icsub_cr(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, rmulti *y0, rmulti *y1)
{
  irsub_rr(C_R(z0),C_R(z1),C_R(x0),C_R(x1),y0,y1); // z.r=x.r-y.r
  irset_r(C_I(z0),C_I(z1),C_I(x0),C_I(x1));       // z.i=x.i
}

/**
 @brief [z0,z1]=[x0,x1]-[y0,y1]
 */
void icsub_rc(cmulti *z0, cmulti *z1, rmulti *x0, rmulti *x1, cmulti *y0, cmulti *y1)
{
  irsub_rr(C_R(z0),C_R(z1),x0,x1,C_R(y0),C_R(y1)); // z.r=x.r-y.r
  irneg_r  (C_I(z0),C_I(z1),      C_I(y0),C_I(y1)); // z.i=   -y.i
}

/**
 @brief [z0,z1]=[x0,x1]-[y0,y1]
 */
void icsub_cz(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, dcomplex y0, dcomplex y1)
{
  irsub_rd(C_R(z0),C_R(z1),C_R(x0),C_R(x1),Z_R(y0),Z_R(y1)); // z.r=x.r-y.r
  irsub_rd(C_I(z0),C_I(z1),C_I(x0),C_I(x1),Z_I(y0),Z_I(y1)); // z.i=x.i-y.i
}

/**
 @brief [z0,z1]=[x0,x1]-[y0,y1]
 */
void icsub_zc(cmulti *z0, cmulti *z1, dcomplex x0, dcomplex x1, cmulti *y0, cmulti *y1)
{
  irsub_dr(C_R(z0),C_R(z1),Z_R(x0),Z_R(x1),C_R(y0),C_R(y1)); // z.r=x.r-y.r
  irsub_dr(C_I(z0),C_I(z1),Z_I(x0),Z_I(x1),C_I(y0),C_I(y1)); // z.i=x.i-y.i
}

/**
 @brief [z0,z1]=[x0,x1]-[y0,y1]
 */
void icsub_cd(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, double y0, double y1)
{
  irsub_rd(C_R(z0),C_R(z1),C_R(x0),C_R(x1),y0,y1); // z.r=x.r-y
  irset_r(C_I(z0),C_I(z1),C_I(x0),C_I(x1));       // z.i=x.i
}

/**
 @brief [z0,z1]=[x0,x1]-[y0,y1]
 */
void icsub_dc(cmulti *z0, cmulti *z1, double x0, double x1, cmulti *y0, cmulti *y1)
{
  irsub_dr(C_R(z0),C_R(z1),x0,x1,C_R(y0),C_R(y1)); // z.r=x.r-y.r
  irneg_r  (C_I(z0),C_I(z1),      C_I(y0),C_I(y1)); // z.i=   -y.i
}

/**
 @brief [z0,z1]=[x0,x1]-[y0,y1]
 */
void icsub_rz(cmulti *z0, cmulti *z1, rmulti *x0, rmulti *x1, dcomplex y0, dcomplex y1)
{
  irsub_rd(C_R(z0),C_R(z1),x0,x1, Z_R(y0), Z_R(y1)); // z.r=x-y.r
  irset_d(C_I(z0),C_I(z1),      -Z_I(y0),-Z_I(y1)); // z.i= -y.i
}

/**
 @brief [z0,z1]=[x0,x1]-[y0,y1]
 */
void icsub_zr(cmulti *z0, cmulti *z1, dcomplex x0, dcomplex x1, rmulti *y0, rmulti *y1)
{
  irsub_dr(C_R(z0),C_R(z1),Z_R(x0),Z_R(x1),y0,y1); // z.r=x.r-y
  irset_d(C_I(z0),C_I(z1),Z_I(x0),Z_I(x1));       // z.i=x.i
}

////////////////////////////////////////////////////////////

/**
 @brief [z0,z1]=[x0,x1]*[y0,y1]
 */
void icmul_cc_ws(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, cmulti *y0, cmulti *y1, int n_rws, rmulti **rws, int n_cws, cmulti **cws)
{
  if(n_rws<8){ ERROR_EXIT("Error `n_rws=%d<8' in icmul_ws().\n",n_rws); }
  if(n_cws<2){ ERROR_EXIT("Error `n_cws=%d<2' in icmul_ws().\n",n_cws); }
  irmul_rr_ws  (C_R(cws[0]),C_R(cws[1]),C_R(x0),C_R(x1),C_R(y0),C_R(y1),n_rws-2,rws+2); // z.r =x.r*y.r
  irsub_mul_rr_ws(C_R(cws[0]),C_R(cws[1]),C_I(x0),C_I(x1),C_I(y0),C_I(y1),n_rws-2,rws+2); // z.r-=x.i*y.i
  irmul_rr_ws  (C_I(cws[0]),C_I(cws[1]),C_I(x0),C_I(x1),C_R(y0),C_R(y1),n_rws-2,rws+2); // z.i =x.i*y.r
  iradd_mul_rr_ws(C_I(cws[0]),C_I(cws[1]),C_R(x0),C_R(x1),C_I(y0),C_I(y1),n_rws-2,rws+2); // z.i+=x.r*y.i
  icset_c(z0,z1,cws[0],cws[1]);
}

/**
 @brief [z0,z1]=[x0,x1]*[y0,y1]
 */
void icmul_cc(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, cmulti *y0, cmulti *y1)
{
  int prec;
  cmulti *a0=NULL,*a1=NULL;
  prec=MAX2(cget_prec(z0),cget_prec(z1));
  CA(a0,prec); CA(a1,prec);
  irmul_rr  (C_R(a0),C_R(a1),C_R(x0),C_R(x1),C_R(y0),C_R(y1)); // z.r =x.r*y.r
  irsub_mul_rr(C_R(a0),C_R(a1),C_I(x0),C_I(x1),C_I(y0),C_I(y1)); // z.r-=x.i*y.i
  irmul_rr  (C_I(a0),C_I(a1),C_I(x0),C_I(x1),C_R(y0),C_R(y1)); // z.i =x.i*y.r
  iradd_mul_rr(C_I(a0),C_I(a1),C_R(x0),C_R(x1),C_I(y0),C_I(y1)); // z.i+=x.r*y.i
  icset_c(z0,z1,a0,a1);
  CF(a0); CF(a1);
}

/**
 @brief [z0,z1]=[x0,x1]*[y0,y1]
 */
void icmul_cr(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, rmulti *y0, rmulti *y1)
{
  int prec;
  cmulti *a0=NULL,*a1=NULL;
  prec=MAX2(cget_prec(z0),cget_prec(z1));
  CA(a0,prec); CA(a1,prec);
  irmul_rr(C_R(a0),C_R(a1),C_R(x0),C_R(x1),y0,y1); // z.r =x.r*y.r
  irmul_rr(C_I(a0),C_I(a1),C_I(x0),C_I(x1),y0,y1); // z.i =x.i*y.r
  icset_c(z0,z1,a0,a1);
  CF(a0); CF(a1);
}

/**
 @brief [z0,z1]=[x0,x1]*[y0,y1]
 */
void icmul_rc(cmulti *z0, cmulti *z1, rmulti *x0, rmulti *x1, cmulti *y0, cmulti *y1)
{
  int prec;
  cmulti *a0=NULL,*a1=NULL;
  prec=MAX2(cget_prec(z0),cget_prec(z1));
  CA(a0,prec); CA(a1,prec);
  irmul_rr(C_R(a0),C_R(a1),x0,x1,C_R(y0),C_R(y1)); // z.r =x.r*y.r
  irmul_rr(C_I(a0),C_I(a1),x0,x1,C_I(y0),C_I(y1)); // z.i =x.r*y.i
  icset_c(z0,z1,a0,a1);
  CF(a0); CF(a1);
}

/**
 @brief [z0,z1]=[x0,x1]*[y0,y1]
 */
void icmul_cz(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, dcomplex y0, dcomplex y1)
{
  int prec;
  cmulti *a0=NULL,*a1=NULL;
  prec=MAX2(cget_prec(z0),cget_prec(z1));
  CA(a0,prec); CA(a1,prec);
  icset_z(a0,a1,y0,y1);
  icmul_cc(z0,z1,x0,x1,a0,a1);
  CF(a0); CF(a1);
}

/**
 @brief [z0,z1]=[x0,x1]*[y0,y1]
 */
void icmul_zc(cmulti *z0, cmulti *z1, dcomplex x0, dcomplex x1, cmulti *y0, cmulti *y1)
{
  int prec;
  cmulti *a0=NULL,*a1=NULL;
  prec=MAX2(cget_prec(z0),cget_prec(z1));
  CA(a0,prec); CA(a1,prec);
  icset_z(a0,a1,x0,x1);
  icmul_cc(z0,z1,a0,a1,y0,y1);
  CF(a0); CF(a1);
}

/**
 @brief [z0,z1]=[x0,x1]*[y0,y1]
 */
void icmul_cd(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, double y0, double y1)
{
  int prec;
  rmulti *a0=NULL,*a1=NULL;
  prec=MAX2(cget_prec(z0),cget_prec(z1));
  RA(a0,prec); RA(a1,prec);
  irset_d(a0,a1,y0,y1);
  icmul_cr(z0,z1,x0,x1,a0,a1);
  RF(a0); RF(a1);
}

/**
 @brief [z0,z1]=[x0,x1]*[y0,y1]
 */
void icmul_dc(cmulti *z0, cmulti *z1, double x0, double x1, cmulti *y0, cmulti *y1)
{
  int prec;
  rmulti *a0=NULL,*a1=NULL;
  prec=MAX2(cget_prec(z0),cget_prec(z1));
  RA(a0,prec); RA(a1,prec);
  irset_d(a0,a1,x0,x1);
  icmul_rc(z0,z1,a0,a1,y0,y1);
  RF(a0); RF(a1);
}

/**
 @brief [z0,z1]=[x0,x1]*[y0,y1]
 */
void icmul_rz(cmulti *z0, cmulti *z1, rmulti *x0, rmulti *x1, dcomplex y0, dcomplex y1)
{
  int prec;
  cmulti *a0=NULL,*a1=NULL;
  prec=MAX2(cget_prec(z0),cget_prec(z1));
  CA(a0,prec); CA(a1,prec);
  icset_z(a0,a1,y0,y1);
  icmul_rc(z0,z1,x0,x1,a0,a1);
  CF(a0); CF(a1);
}

/**
 @brief [z0,z1]=[x0,x1]*[y0,y1]
 */
void icmul_zr(cmulti *z0, cmulti *z1, dcomplex x0, dcomplex x1, rmulti *y0, rmulti *y1)
{
  int prec;
  cmulti *a0=NULL,*a1=NULL;
  prec=MAX2(cget_prec(z0),cget_prec(z1));
  CA(a0,prec); CA(a1,prec);
  icset_z(a0,a1,x0,x1);
  icmul_cr(z0,z1,a0,a1,y0,y1);
  CF(a0); CF(a1);
}

////////////////////////////////////////////////////////////

/**
 @brief [z0,z1]=[x0,x1]/[y0,y1]
 */
void icdiv_cc(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, cmulti *y0, cmulti *y1)
{
  int prec;
  rmulti *den0=NULL,*den1=NULL;
  cmulti *a0=NULL,*a1=NULL;
  prec=MAX2(cget_prec(z0),cget_prec(z1));
  RA(den0,prec); RA(den1,prec); CA(a0,prec); CA(a1,prec);
  irabs2_c(den0,den1,y0,y1);        // den=|y|^2
  icdot_cc(a0,a1,y0,y1,x0,x1);       // z=conj(y)*x
  icdiv_cr(a0,a1,a0,a1,den0,den1); // z/=den
  icset_c(z0,z1,a0,a1);
  RF(den0); RF(den1); CF(a0); CF(a1);
}

/**
 @brief [z0,z1]=[x0,x1]/[y0,y1]
 */
void icdiv_cr(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, rmulti *y0, rmulti *y1)
{
  irdiv_rr(C_R(z0),C_R(z1),C_R(x0),C_R(x1),y0,y1);
  irdiv_rr(C_I(z0),C_I(z1),C_I(x0),C_I(x1),y0,y1);
}

/**
 @brief [z0,z1]=[x0,x1]/[y0,y1]
 */
void icdiv_rc(cmulti *z0, cmulti *z1, rmulti *x0, rmulti *x1, cmulti *y0, cmulti *y1)
{
  int prec;
  rmulti *den0=NULL,*den1=NULL;
  cmulti *a0=NULL,*a1=NULL;
  prec=MAX2(cget_prec(z0),cget_prec(z1));
  RA(den0,prec); RA(den1,prec); CA(a0,prec); CA(a1,prec);
  irabs2_c(den0,den1,y0,y1);                         // den=|y|^2
  irmul_rr(C_R(a0),C_R(a1),C_R(y0),C_R(y1),x0,x1);  // z.r =x.r*y.r
  irmul_rr(C_I(a0),C_I(a1),C_I(y0),C_I(y1),x0,x1);  // z.i =x.i*y.r
  icconj_c(a0,a1,a0,a1);                             // z=-(x.r*y.r+i*x.i*y.r)
  icdiv_cr(a0,a1,a0,a1,den0,den1);                  // z=-(x.r*y.r+i*x.i*y.r)/den
  icset_c(z0,z1,a0,a1);
  RF(den0); RF(den1); CF(a0); CF(a1);
}

/**
 @brief [z0,z1]=[x0,x1]/[y0,y1]
 */
void icdiv_cz(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, dcomplex y0, dcomplex y1)
{
  int prec;
  cmulti *a0=NULL,*a1=NULL;
  prec=MAX2(cget_prec(z0),cget_prec(z1));
  CA(a0,prec); CA(a1,prec);
  icset_z(a0,a1,y0,y1);
  icdiv_cc(z0,z1,x0,x1,a0,a1);
  CF(a0); CF(a1);
}

/**
 @brief [z0,z1]=[x0,x1]/[y0,y1]
 */
void icdiv_zc(cmulti *z0, cmulti *z1, dcomplex x0, dcomplex x1, cmulti *y0, cmulti *y1)
{
  int prec;
  cmulti *a0=NULL,*a1=NULL;
  prec=MAX2(cget_prec(z0),cget_prec(z1));
  CA(a0,prec); CA(a1,prec);
  icset_z(a0,a1,x0,x1);
  icdiv_cc(z0,z1,a0,a1,y0,y1);
  CF(a0); CF(a1);
}

/**
 @brief [z0,z1]=[x0,x1]/[y0,y1]
 */
void icdiv_cd(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, double y0, double y1)
{
  int prec;
  rmulti *a0=NULL,*a1=NULL;
  prec=MAX2(cget_prec(z0),cget_prec(z1));
  RA(a0,prec); RA(a1,prec);
  irset_d(a0,a1,y0,y1);
  icdiv_cr(z0,z1,x0,x1,a0,a1);
  RF(a0); RF(a1);
}

/**
 @brief [z0,z1]=[x0,x1]/[y0,y1]
 */
void icdiv_dc(cmulti *z0, cmulti *z1, double x0, double x1, cmulti *y0, cmulti *y1)
{
  int prec;
  rmulti *a0=NULL,*a1=NULL;
  prec=MAX2(cget_prec(z0),cget_prec(z1));
  RA(a0,prec); RA(a1,prec);
  irset_d(a0,a1,x0,x1);
  icdiv_rc(z0,z1,a0,a1,y0,y1);
  RF(a0); RF(a1);
}

/**
 @brief [z0,z1]=[x0,x1]/[y0,y1]
 */
void icdiv_rz(cmulti *z0, cmulti *z1, rmulti *x0, rmulti *x1, dcomplex y0, dcomplex y1)
{
  int prec;
  cmulti *a0=NULL,*a1=NULL;
  prec=MAX2(cget_prec(z0),cget_prec(z1));
  CA(a0,prec); CA(a1,prec);
  icset_z(a0,a1,y0,y1);
  icdiv_rc(z0,z1,x0,x1,a0,a1);
  CF(a0); CF(a1);
}

/**
 @brief [z0,z1]=[x0,x1]/[y0,y1]
 */
void icdiv_zr(cmulti *z0, cmulti *z1, dcomplex x0, dcomplex x1, rmulti *y0, rmulti *y1)
{
  int prec;
  cmulti *a0=NULL,*a1=NULL;
  prec=MAX2(cget_prec(z0),cget_prec(z1));
  CA(a0,prec); CA(a1,prec);
  icset_z(a0,a1,x0,x1);
  icdiv_cr(z0,z1,a0,a1,y0,y1);
  CF(a0); CF(a1);
}

////////////////////////////////////////////////////////////

/**
 @brief 足し算 [z0,z1]=[x0,x1]+[-y,y]
 */
void icadd_pm_cc(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, cmulti *y)
{
  iradd_pm_rr(C_R(z0),C_R(z1),C_R(x0),C_R(x1),C_R(y)); // z.r=x.r+y.r
  iradd_pm_rr(C_I(z0),C_I(z1),C_I(x0),C_I(x1),C_I(y)); // z.i=x.i+y.i
}

/**
 @brief 積の加算 [z0,z1]+=[x0,x1]*[y0,y1]
 */
void icadd_mul_rr(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, cmulti *y0, cmulti *y1)
{
  int prec;
  cmulti *a0=NULL,*a1=NULL;
  prec=MAX2(cget_prec(z0),cget_prec(z1));
  CA(a0,prec); CA(a1,prec);
  icmul_cc(a0,a1,x0,x1,y0,y1); // a=x*y
  icadd_cc(z0,z1,z0,z1,a0,a1); // z=z+x*y
  CF(a0); CF(a1);
}

/**
 @brief 積の加算 [z0,z1]+=[x0,x1]*[y0,y1]
 */
void icadd_mul_rc(cmulti *z0, cmulti *z1, rmulti *x0, rmulti *x1, cmulti *y0, cmulti *y1)
{
  int prec;
  cmulti *a0=NULL,*a1=NULL;
  prec=MAX2(cget_prec(z0),cget_prec(z1));
  CA(a0,prec); CA(a1,prec);
  icmul_rc(a0,a1,x0,x1,y0,y1); // a=x*y
  icadd_cc(z0,z1,z0,z1,a0,a1);    // z=z+x*y
  CF(a0); CF(a1);
}

/**
 @brief 積の加算 [z0,z1]+=[x0,x1]*[y0,y1]
 */
void icadd_mul_cr(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, rmulti *y0, rmulti *y1)
{
  int prec;
  cmulti *a0=NULL,*a1=NULL;
  prec=MAX2(cget_prec(z0),cget_prec(z1));
  CA(a0,prec); CA(a1,prec);
  icmul_cr(a0,a1,x0,x1,y0,y1); // a=x*y
  icadd_cc(z0,z1,z0,z1,a0,a1);    // z=z+x*y
  CF(a0); CF(a1);
}

/**
 @brief 積の減算 [z0,z1]-=[x0,x1]*[y0,y1]
 */
void icsub_mul_cc(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, cmulti *y0, cmulti *y1)
{
  int prec;
  cmulti *a0=NULL,*a1=NULL;
  prec=MAX2(cget_prec(z0),cget_prec(z1));
  CA(a0,prec); CA(a1,prec);
  icmul_cc(a0,a1,x0,x1,y0,y1); // a=x*y
  icsub_cc(z0,z1,z0,z1,a0,a1); // z=z-x*y
  CF(a0); CF(a1);
}

/**
 @brief 積の減算 [z0,z1]-=[x0,x1]*[y0,y1]
 */
void icsub_mul_cc_ws(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, cmulti *y0, cmulti *y1, int n_rws, rmulti **rws, int n_cws, cmulti **cws)
{
  if(n_rws<8){ ERROR_EXIT("Error `n_rws=%d<8' in icsub_mul_ws().\n",n_rws); }
  if(n_cws<4){ ERROR_EXIT("Error `n_cws=%d<4' in icsub_mul_ws().\n",n_cws); }
  icmul_cc_ws(cws[0],cws[1],x0,x1,y0,y1,n_rws,rws,n_cws-2,cws+2); // a=x*y
  icsub_cc_ws(z0,z1,z0,z1,cws[0],cws[1],n_rws,rws);               // z=z-x*y
}

/**
 @brief 共役な掛け算 [z0,z1]=conj([x0,x1])*[y0,y1]
 */
void icdot_cc(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, cmulti *y0, cmulti *y1)
{
  int prec;
  cmulti *a0=NULL,*a1=NULL;
  prec=MAX2(cget_prec(z0),cget_prec(z1));
  CA(a0,prec); CA(a1,prec);
  irmul_rr    (C_R(a0),C_R(a1),C_R(x0),C_R(x1),C_R(y0),C_R(y1));  // z.r =x.r*y.r
  iradd_mul_rr(C_R(a0),C_R(a1),C_I(x0),C_I(x1),C_I(y0),C_I(y1));  // z.r+=x.i*y.i
  irmul_rr    (C_I(a0),C_I(a1),C_R(x0),C_R(x1),C_I(y0),C_I(y1));  // z.i =x.r*y.i
  irsub_mul_rr(C_I(a0),C_I(a1),C_I(x0),C_I(x1),C_R(y0),C_R(y1));  // z.i-=x.i*y.r
  icset_c(z0,z1,a0,a1);
  CF(a0); CF(a1);
}

/**
 @brief 共役な積の加算 [z0,z1]+=conj([x0,x1])*[y0,y1]
 */
void icadd_dot_cc(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, cmulti *y0, cmulti *y1)
{
  int prec;
  cmulti *a0=NULL,*a1=NULL;
  prec=MAX2(cget_prec(z0),cget_prec(z1));
  CA(a0,prec); CA(a1,prec);
  icdot_cc(a0,a1,x0,x1,y0,y1); // a=conj(x)*y
  icadd_cc(z0,z1,z0,z1,a0,a1); // z=z+x*y
  CF(a0); CF(a1);
}

/**
 @brief 共役な積の減算 [z0,z1]-=conj([x0,x1])*[y0,y1]
 */
void icsub_dot_cc(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, cmulti *y0, cmulti *y1)
{
  int prec;
  cmulti *a0=NULL,*a1=NULL;
  prec=MAX2(cget_prec(z0),cget_prec(z1));
  CA(a0,prec); CA(a1,prec);
  icdot_cc(a0,a1,x0,x1,y0,y1); // a=conj(x)*y
  icsub_cc(z0,z1,z0,z1,a0,a1); // z=z-x*y
  CF(a0); CF(a1);
}

/**
 @brief 差の絶対値 [z0,z1]=abs([x0,x1]-[y0,y1])
 */
void irabs_sub_cc(rmulti *z0, rmulti *z1, cmulti *x0, cmulti *x1, cmulti *y0, cmulti *y1)
{
  cmulti *a0=NULL,*a1=NULL;
  CAr(a0,z0); CAr(a1,z1);
  icsub_cc(a0,a1,x0,x1,y0,y1);
  irabs_c(z0,z1,a0,a1);
  CF(a0); CF(a1);
}

/**
 @brief べき乗 [y0,y1]=[x0,x1]^n
 */
void icpow_c(cmulti *y0, cmulti *y1, cmulti *x0, cmulti *x1, int n)
{
  cmulti *a0=NULL,*a1=NULL;
  int prec;
  long i;
  prec=MAX2(cget_prec(y0),cget_prec(y1));
  CA(a0,prec); CA(a1,prec);
  if(n==0){ icset_d(a0,a1,1,1); }
  else if(n>0){
    icset_c(a0,a1,x0,x1);
    for(i=1; i<n; i++){ icmul_cc(a0,a1,a0,a1,x0,x1); }
  }else if(n<0){
    icset_c(a0,a1,x0,x1);
    for(i=1; i<(-n); i++){ icmul_cc(a0,a1,a0,a1,x0,x1); }
    icinv_c(a0,a1,a0,a1);
  }
  icset_c(y0,y1,a0,a1);
  CF(a0); CF(a1);
}


/** @} */



//////////////////////////////////////////////////////

/** @name icmulti型の判定に関する関数 */
/** @{ */

/**
 @brief 区間の内包の判定 [y0,y1] in [x0,x1], すなわち (x0<=y0 && y1 <=x1).
 */
int ic_in_ic(cmulti *y0, cmulti *y1, cmulti *x0, cmulti *x1)
{
  return (ir_in_ir(C_R(y0),C_R(y1),C_R(x0),C_R(x1)) && ir_in_ir(C_I(y0),C_I(y1),C_I(x0),C_I(x1)));
}

/**
 @brief 区間の内包の判定 [y0,y1] in [-abs(x),abs(x)], すなわち (-abs(x)<=y0 && y1 <=abs(x)).
 */
int ic_in_icpm(cmulti *y0, cmulti *y1, cmulti *x)
{
  return (ir_in_irpm(C_R(y0),C_R(y1),C_R(x)) && ir_in_irpm(C_I(y0),C_I(y1),C_I(x)));
}

/** @} */

//EOF
