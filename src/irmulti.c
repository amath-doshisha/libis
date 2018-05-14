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
void irset_r(rmulti *y0, rmulti *y1, rmulti *x0, rmulti *x1)
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
  irset_r(x0,x1,C_R(a0),C_R(a1));
}

/** @} */
/** @name irmulti型の型変換に関する関数 */
/** @{ */

/** @} */
/** @name irmulti型の１入力の数学演算子 */
/** @{ */

/**
 @brief 符号の反転 [y0,y1]=-[x0,x1].
 */
void irneg_r(rmulti *y0, rmulti *y1, rmulti *x0, rmulti *x1)
{
  int p0,p1,prec;
  rmulti *a0=NULL,*a1=NULL;
  p0=rget_prec(y0); p1=rget_prec(y1); prec=MAX2(p0,p1);
  RA(a0,prec); RA(a1,prec);
  mpfr_neg(a0,x1,MPFR_RNDD); // lower bound
  mpfr_neg(a1,x0,MPFR_RNDU); // upper bound
  irset_r(y0,y1,a0,a1);
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
  irset_r(mid,rad,a0,a1);           // [mid,rad]=[a0,a1]
  RF(a0); RF(a1);
}

/**
 @brief 逆数 [z0,z1]=[1,1]/[x0,x1]
*/
void irinv_r(rmulti *z0, rmulti *z1, rmulti *x0, rmulti *x1)
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
void irabs_r(rmulti *y0, rmulti *y1, rmulti *x0, rmulti *x1)
{
  int p0,p1,prec;
  rmulti *a0=NULL,*a1=NULL;
  if(ris_nan(x0) || ris_nan(x1)){ rset_nan(y0); rset_nan(y1); return; }
  // x0>=0 && x1>=0
  if(rget_sgn(x0)>=0 && rget_sgn(x1)>=0){ irset_r(y0,y1,x0,x1); }
  // x0<=0 && x1<=0
  else if(rget_sgn(x0)<=0 && rget_sgn(x1)<=0){ irneg_r(y0,y1,x0,x1); }
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
void irsqrt_r(rmulti *y0, rmulti *y1, rmulti *x0, rmulti *x1)
{
  mpfr_sqrt(y0,x0,MPFR_RNDD); // lower bound
  mpfr_sqrt(y1,x1,MPFR_RNDU); // upper bound
}

/** @} */

/////////////////////////////////////////////////////////////

/** @name irmulti型の２入力の数学演算子 */
/** @{ */

/**
 @brief [z0,z1]=[x0,x1]+[y0,y1]
 */
void iradd_rr(rmulti *z0, rmulti *z1, rmulti *x0, rmulti *x1, rmulti *y0, rmulti *y1)
{
  mpfr_add(z0,x0,y0,MPFR_RNDD); // lower bound
  mpfr_add(z1,x1,y1,MPFR_RNDU); // upper bound
}

/**
 @brief [z0,z1]=[x0,x1]+[y0,y1]
 */
void iradd_dr(rmulti *z0, rmulti *z1, double x0, double x1, rmulti *y0, rmulti *y1)
{
  mpfr_add_d(z0,y0,x0,MPFR_RNDD); // lower bound
  mpfr_add_d(z1,y1,x1,MPFR_RNDU); // upper bound
}
/**
 @brief [z0,z1]=[x0,x1]+[y0,y1]
 */
void iradd_rd(rmulti *z0, rmulti *z1, rmulti *x0, rmulti *x1, double y0, double y1)
{
  mpfr_add_d(z0,x0,y0,MPFR_RNDD); // lower bound
  mpfr_add_d(z1,x1,y1,MPFR_RNDU); // upper bound
}

///////////////////////////////////

/**
 @brief [z0,z1]=[x0,x1]-[y0,y1]
 */
void irsub_rr_ws(rmulti *z0, rmulti *z1, rmulti *x0, rmulti *x1, rmulti *y0, rmulti *y1, int n_rws, rmulti **rws)
{
  if(n_rws<2){ ERROR_EXIT("Error `n_rws=%d<2' in irsub_ws().\n",n_rws); }
  mpfr_sub(rws[0],x0,y1,MPFR_RNDD); // lower bound
  mpfr_sub(rws[1],x1,y0,MPFR_RNDU); // upper bound
  irset_r(z0,z1,rws[0],rws[1]);
}

/**
 @brief [z0,z1]=[x0,x1]-[y0,y1]
 */
void irsub_rr(rmulti *z0, rmulti *z1, rmulti *x0, rmulti *x1, rmulti *y0, rmulti *y1)
{
  int p0,p1,prec;
  rmulti *a0=NULL,*a1=NULL;
  p0=rget_prec(z0); p1=rget_prec(z1); prec=MAX2(p0,p1);
  RA(a0,prec); RA(a1,prec);
  mpfr_sub(a0,x0,y1,MPFR_RNDD); // lower bound
  mpfr_sub(a1,x1,y0,MPFR_RNDU); // upper bound
  irset_r(z0,z1,a0,a1);
  RF(a0); RF(a1);
}

/**
 @brief [z0,z1]=[x0,x1]-[y0,y1]
 */
void irsub_rd(rmulti *z0, rmulti *z1, rmulti *x0, rmulti *x1, double y0, double y1)
{
  int p0,p1,prec;
  rmulti *a0=NULL,*a1=NULL;
  p0=rget_prec(z0); p1=rget_prec(z1); prec=MAX2(p0,p1);
  RA(a0,prec); RA(a1,prec);
  irset_d(a0,a1,y0,y1);
  irsub_rr(z0,z1,x0,x1,a0,a1);
  RF(a0); RF(a1);
}

/**
 @brief [z0,z1]=[x0,x1]-[y,y]
 */
void irsub_dr(rmulti *z0, rmulti *z1, double x0, double x1, rmulti *y0, rmulti *y1)
{
  int p0,p1,prec;
  rmulti *a0=NULL,*a1=NULL;
  p0=rget_prec(z0); p1=rget_prec(z1); prec=MAX2(p0,p1);
  RA(a0,prec); RA(a1,prec);
  irset_d(a0,a1,x0,x1);
  irsub_rr(z0,z1,a0,a1,y0,y1);
  RF(a0); RF(a1);
}

///////////////////////////////////

/**
 @brief [z0,z1]=[x0,x1]*[y0,y1]
 */
void irmul_rr_ws(rmulti *z0, rmulti *z1, rmulti *x0, rmulti *x1, rmulti *y0, rmulti *y1, int n_rws, rmulti **rws)
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
      if(le_rr(rws[2],rws[3])){ mpfr_set(rws[0],rws[2],MPFR_RNDD); }else{ mpfr_set(rws[0],rws[3],MPFR_RNDD); }
      mpfr_mul(rws[2],x0,y0,MPFR_RNDU);
      mpfr_mul(rws[3],x1,y1,MPFR_RNDU);
      if(ge_rr(rws[2],rws[3])){ mpfr_set(rws[1],rws[2],MPFR_RNDU); }else{ mpfr_set(rws[1],rws[3],MPFR_RNDU); }
    }
  }
  irset_r(z0,z1,rws[0],rws[1]);
}

/**
 @brief [z0,z1]=[x0,x1]*[y0,y1]
 */
void irmul_rr(rmulti *z0, rmulti *z1, rmulti *x0, rmulti *x1, rmulti *y0, rmulti *y1)
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
      if(le_rr(a,b)){ mpfr_set(Z0,a,MPFR_RNDD); }else{ mpfr_set(Z0,b,MPFR_RNDD); }
      mpfr_mul(a,x0,y0,MPFR_RNDU);
      mpfr_mul(b,x1,y1,MPFR_RNDU);
      if(ge_rr(a,b)){ mpfr_set(Z1,a,MPFR_RNDU); }else{ mpfr_set(Z1,b,MPFR_RNDU); }
    }
  }
  irset_r(z0,z1,Z0,Z1);
  RF(Z0); RF(Z1); RF(a); RF(b);
}

/**
 @brief [z0,z1]=[x0,x1]*[y0,y1]
 */
void irmul_rd(rmulti *z0, rmulti *z1, rmulti *x0, rmulti *x1, double y0, double y1)
{
  int p0,p1,prec;
  rmulti *a0=NULL,*a1=NULL;
  p0=rget_prec(z0); p1=rget_prec(z1); prec=MAX2(p0,p1);
  RA(a0,prec); RA(a1,prec);
  irset_d(a0,a1,y0,y1);
  irmul_rr(z0,z1,x0,x1,a0,a1);
  RF(a0); RF(a1);
}

/**
 @brief [z0,z1]=[x0,x1]*[y,y]
 */
void irmul_dr(rmulti *z0, rmulti *z1, double x0, double x1, rmulti *y0, rmulti *y1)
{
  int p0,p1,prec;
  rmulti *a0=NULL,*a1=NULL;
  p0=rget_prec(z0); p1=rget_prec(z1); prec=MAX2(p0,p1);
  RA(a0,prec); RA(a1,prec);
  irset_d(a0,a1,x0,x1);
  irmul_rr(z0,z1,a0,a1,y0,y1);
  RF(a0); RF(a1);
}

//////////////////////////////////////////////////////////////

/**
 @brief [z0,z1]=[x0,x1]/[y0,y1]
 */
void irdiv_rr(rmulti *z0, rmulti *z1, rmulti *x0, rmulti *x1, rmulti *y0, rmulti *y1)
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
  irset_r(z0,z1,Z0,Z1);
  RF(Z0); RF(Z1);
}

/**
 @brief [z0,z1]=[x0,x1]/[y0,y1]
 */
void irdiv_rd(rmulti *z0, rmulti *z1, rmulti *x0, rmulti *x1, double y0, double y1)
{
  int p0,p1,prec;
  rmulti *a0=NULL,*a1=NULL;
  p0=rget_prec(z0); p1=rget_prec(z1); prec=MAX2(p0,p1);
  RA(a0,prec); RA(a1,prec);
  irset_d(a0,a1,y0,y1);
  irdiv_rr(z0,z1,x0,x1,a0,a1);
  RF(a0); RF(a1);
}

/**
 @brief [z0,z1]=[x0,x1]/[y,y]
 */
void irdiv_dr(rmulti *z0, rmulti *z1, double x0, double x1, rmulti *y0, rmulti *y1)
{
  int p0,p1,prec;
  rmulti *a0=NULL,*a1=NULL;
  p0=rget_prec(z0); p1=rget_prec(z1); prec=MAX2(p0,p1);
  RA(a0,prec); RA(a1,prec);
  irset_d(a0,a1,x0,x1);
  irdiv_rr(z0,z1,a0,a1,y0,y1);
  RF(a0); RF(a1);
}


///////////////////////////////////


/**
 @brief 区間の拡張 [z0,z1]=[x0,x1]+[-y,y]
 */
void iradd_pm_rr(rmulti *z0, rmulti *z1, rmulti *x0, rmulti *x1, rmulti *y)
{
  mpfr_sub(z0,x0,y,MPFR_RNDD); // lower bound
  mpfr_add(z1,x1,y,MPFR_RNDU); // upper bound
}



/**
 @brief 積の加算 [z0,z1]+=[x0,x1]*[y0,y1]
*/
void iradd_mul_rr(rmulti *z0, rmulti *z1, rmulti *x0, rmulti *x1, rmulti *y0, rmulti *y1)
{
  int p0,p1,prec;
  rmulti *a0=NULL,*a1=NULL;
  p0=rget_prec(z0); p1=rget_prec(z1); prec=MAX2(p0,p1);
  RA(a0,prec); RA(a1,prec);
  irmul_rr(a0,a1,x0,x1,y0,y1); // a=x*y
  iradd_rr(z0,z1,z0,z1,a0,a1); // z=z+x*y
  RF(a0); RF(a1);
}

/**
 @brief 積の加算 [z0,z1]+=[x0,x1]*[y0,y1]
*/
void iradd_mul_rr_ws(rmulti *z0, rmulti *z1, rmulti *x0, rmulti *x1, rmulti *y0, rmulti *y1, int n_rws, rmulti **rws)
{
  if(n_rws<6){ ERROR_EXIT("Error `n_rws=%d<6' in iradd_mul_ws().\n",n_rws); }
  irmul_rr_ws(rws[0],rws[1],x0,x1,y0,y1,n_rws-2,rws+2); // a=x*y
  iradd_rr(z0,z1,z0,z1,rws[0],rws[1]); // z=z+x*y
}

/**
 @brief 積の減算 [z0,z1]-=[x0,x1]*[y0,y1]
*/
void irsub_mul_rr(rmulti *z0, rmulti *z1, rmulti *x0, rmulti *x1, rmulti *y0, rmulti *y1)
{
  int p0,p1,prec;
  rmulti *a0=NULL,*a1=NULL;
  p0=rget_prec(z0); p1=rget_prec(z1); prec=MAX2(p0,p1);
  RA(a0,prec); RA(a1,prec);
  irmul_rr(a0,a1,x0,x1,y0,y1); // a=x*y
  irsub_rr(z0,z1,z0,z1,a0,a1); // z=z-x*y
  RF(a0); RF(a1);
}

/**
 @brief 積の減算 [z0,z1]-=[x0,x1]*[y0,y1]
*/
void irsub_mul_rr_ws(rmulti *z0, rmulti *z1, rmulti *x0, rmulti *x1, rmulti *y0, rmulti *y1, int n_rws, rmulti **rws)
{
  if(n_rws<6){ ERROR_EXIT("Error `n_rws=%d<6' in irsub_mul_ws().\n",n_rws); }
  irmul_rr_ws(rws[0],rws[1],x0,x1,y0,y1,n_rws-2,rws+2); // a=x*y
  irsub_rr_ws(z0,z1,z0,z1,rws[0],rws[1],n_rws-2,rws+2); // z=z-x*y
}

/**
 @brief 差の絶対値 [z0,z1]=abs([x0,x1]-[y0,y1])
*/
void irabs_sub_rr(rmulti *z0, rmulti *z1, rmulti *x0, rmulti *x1, rmulti *y0, rmulti *y1)
{
  irsub_rr(z0,z1,x0,x1,y0,y1);
  irabs_r(z0,z1,z0,z1);
}

/**
 @brief べき乗 [y0,y1]=[x0,x1]^n
 */
void irpow_r(rmulti *y0, rmulti *y1, rmulti *x0, rmulti *x1, int n)
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
      if(ge_rr(b0,b1)){ irset_r(y0,y1,a,b0); }else{ irset_r(y0,y1,a,b1); }
    }
  }else{
    rset_d(a,1);
    irpow_r(b0,b1,x0,x1,-n);
    irdiv_rr(y0,y1,a,a,b0,b1);
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
  return (le_rr(x0,y0) && le_rr(y1,x1));
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
  value=(le_rr(z0,y0) && le_rr(y1,z1));
  z0=rfree(z0); z1=rfree(z1);
  return value;
}

/** @} */

//EOF
