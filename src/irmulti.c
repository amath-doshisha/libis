#include"is_rmulti.h"
#include"is_rmat.h"
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
#define RF(X) ((X)=rfree(X))

/** @name 基本操作 */
/** @{ */

/**
 @brief 倍精度実数の設定 [y0,y1]=[x,x].
 */
int irset_d(rmulti *y0, rmulti *y1, double x)
{
  int e=0;
  e+=abs(mpfr_set_d(y0,x,MPFR_RNDD));
  e+=abs(mpfr_set_d(y1,x,MPFR_RNDU));
  return e;
}

//追加

/**
 @brief 倍精度実数の設定 [y0,y1]=[x0,x1].
 */
int irset_dd(rmulti *y0, rmulti *y1, double x0, double x1)
{
  int e=0;
  e+=abs(mpfr_set_d(y0,x0,MPFR_RNDD));
  e+=abs(mpfr_set_d(y1,x1,MPFR_RNDU));
  return e;
}

//ここまで

/**
 @brief bigint型から[z0,z1]へ型変換.
 */
int irset_bigint(rmulti *z0, rmulti *z1, bigint *x)
{
  int e=0;
  e+=abs(mpfr_div(z0,BIGINT_NUM(x),BIGINT_DEN(x),MPFR_RNDD));
  e+=abs(mpfr_div(z1,BIGINT_NUM(x),BIGINT_DEN(x),MPFR_RNDU));
  return e;
}

/**
 @brief コピー [y0,y1]=[x0,x1].
 */
int ircopy(rmulti *y0, rmulti *y1, rmulti *x0, rmulti *x1)
{
  int e=0;
  e+=abs(mpfr_set(y0,x0,MPFR_RNDD));
  e+=abs(mpfr_set(y1,x1,MPFR_RNDU));
  return e;
}


/** @} */

//////////////////////////////////////////////////////

/** @name 四則演算 */
/** @{ */

/**
 @brief 符号の反転 [y0,y1]=-[x0,x1].
 */
int irneg(rmulti *y0, rmulti *y1, rmulti *x0, rmulti *x1)
{
  int e=0,p0,p1,prec;
  rmulti *a0=NULL,*a1=NULL;
  p0=rget_prec(y0); p1=rget_prec(y1); prec=MAX2(p0,p1);
  RA(a0,prec); RA(a1,prec);
  e+=abs(mpfr_neg(a0,x1,MPFR_RNDD));
  e+=abs(mpfr_neg(a1,x0,MPFR_RNDU));
  e+=ircopy(y0,y1,a0,a1);
  RF(a0); RF(a1);
  return e;
}

/**
 @brief 符号の正負 [y0,y1]=[-abs(x),abs(x)].
 */
int irpm(rmulti *y0, rmulti *y1, rmulti *x)
{
  int e=0;
  e+=abs(mpfr_abs(y1,x, MPFR_RNDU));
  e+=abs(mpfr_neg(y0,y1,MPFR_RNDD));
  return e;
}

/**
 @brief 足し算 [z0,z1]=[x0,x1]+[y0,y1]
 */
int iradd(rmulti *z0, rmulti *z1, rmulti *x0, rmulti *x1, rmulti *y0, rmulti *y1)
{
  int e=0;
  e+=abs(mpfr_add(z0,x0,y0,MPFR_RNDD)); // lower bound
  e+=abs(mpfr_add(z1,x1,y1,MPFR_RNDU)); // upper bound
  return e;
}

/**
 @brief 区間の拡張 [z0,z1]=[x0,x1]+[-y,y]
 */
int iradd_pm(rmulti *z0, rmulti *z1, rmulti *x0, rmulti *x1, rmulti *y)
{
  int e=0;
  e+=abs(mpfr_sub(z0,x0,y,MPFR_RNDD)); // lower bound
  e+=abs(mpfr_add(z1,x1,y,MPFR_RNDU)); // upper bound
  return e;
}

/**
 @brief 区間の中心 [m-r,m+r]=[x0,x1]
 */
int irmid(rmulti *mid, rmulti *x0, rmulti *x1)
{
  int prec,e=0;
  rmulti *a0=NULL,*a1=NULL;
  prec=rget_prec(mid);
  RA(a0,prec); RA(a1,prec);
  e+=abs(mpfr_sub(a0,x1,x0,MPFR_RNDU));    // a0=(x1-x0) by upper
  e+=abs(mpfr_mul_d(a0,a0,0.5,MPFR_RNDU)); // a0=(x1-x0)/2 by upper
  e+=abs(mpfr_add(mid,a0,x0,MPFR_RNDU));   // mid=x0+(x1-x0)/2 by upper
  RF(a0); RF(a1);
  return e;
}

//変更

/**
 @brief 区間の半径 rad=rad([x0,x1])
 */
int irrad(rmulti *rad, rmulti *x0, rmulti *x1)
{
  int prec,e=0;
  rmulti *a0=NULL;
  prec=rget_prec(rad);
  RA(a0,prec);
  e+=abs(mpfr_sub(a0,x1,x0,MPFR_RNDU));    // a0=(x1-x0) by upper
  e+=abs(mpfr_mul_d(a0,a0,0.5,MPFR_RNDU)); // a0=(x1-x0)/2 by upper
  e+=abs(mpfr_add(a0,a0,x0,MPFR_RNDU));    // a0=x0+(x1-x0)/2 by upper
  e+=abs(mpfr_sub(rad,a0,x0,MPFR_RNDU));   // rad=(x1-x0)/2 by upper
  RF(a0); 
  return e;
}

//ここまで



/**
 @brief 区間の中心と半径の表示(center-radius form) [m-r,m+r]=[x0,x1]
 */
int irmr(rmulti *mid, rmulti *rad, rmulti *x0, rmulti *x1)
{
  int p0,p1,prec,e=0;
  rmulti *a0=NULL,*a1=NULL;
  p0=rget_prec(mid); p1=rget_prec(rad); prec=MAX2(p0,p1);
  RA(a0,prec); RA(a1,prec);
  e+=abs(mpfr_sub(a0,x1,x0,MPFR_RNDU));    // a0=(x1-x0) by upper
  e+=abs(mpfr_mul_d(a0,a0,0.5,MPFR_RNDU)); // a0=(x1-x0)/2 by upper
  e+=abs(mpfr_add(a0,a0,x0,MPFR_RNDU));    // a0=x0+(x1-x0)/2 by upper
  e+=abs(mpfr_sub(a1,a0,x0,MPFR_RNDU));    // a1=(x1-x0)/2 by upper
  e+=ircopy(mid,rad,a0,a1);                // [mid,rad]=[a0,a1]
  RF(a0); RF(a1);
  return e;
}

/**
 @brief 引き算 [z0,z1]=[x0,x1]-[y0,y1]
 */
int irsub(rmulti *z0, rmulti *z1, rmulti *x0, rmulti *x1, rmulti *y0, rmulti *y1)
{
  int p0,p1,prec,e=0;
  rmulti *a0=NULL,*a1=NULL;
  p0=rget_prec(z0); p1=rget_prec(z1); prec=MAX2(p0,p1);
  RA(a0,prec); RA(a1,prec);
  e+=abs(mpfr_sub(a0,x0,y1,MPFR_RNDD)); // lower bound
  e+=abs(mpfr_sub(a1,x1,y0,MPFR_RNDU)); // upper bound
  e+=ircopy(z0,z1,a0,a1);
  RF(a0); RF(a1);
  return e;
}
//追加
/**
 @brief 引き算 [z0,z1]=[x0,x1]-[y0,y1]
 */
int irsub_ws(rmulti *z0, rmulti *z1, rmulti *x0, rmulti *x1, rmulti *y0, rmulti *y1, int *rwss, rmulti **rws)
{
  if(*rwss<2){ ERROR_EXIT("Error `rwss=%d<2' in irsub_ws().\n",*rwss); }
  int e=0;
  e+=abs(mpfr_sub(rws[0],x0,y1,MPFR_RNDD)); // lower bound
  e+=abs(mpfr_sub(rws[1],x1,y0,MPFR_RNDU)); // upper bound
  e+=ircopy(z0,z1,rws[0],rws[1]);
  return e;
}
//ここまで
/**
 @brief 引き算 [z0,z1]=[x0,x1]-[y,y]
 */
int irsub_d2(rmulti *z0, rmulti *z1, rmulti *x0, rmulti *x1, double y)
{
  int p0,p1,prec,e=0;
  rmulti *a0=NULL,*a1=NULL;
  p0=rget_prec(z0); p1=rget_prec(z1); prec=MAX2(p0,p1);
  RA(a0,prec); RA(a1,prec);
  e+=irset_d(a0,a1,y);
  e+=irsub(z0,z1,x0,x1,a0,a1);
  RF(a0); RF(a1);
  return e;
}

/**
 @brief 掛け算 [z0,z1]=[x0,x1]*[y0,y1]
*/
int irmul(rmulti *z0, rmulti *z1, rmulti *x0, rmulti *x1, rmulti *y0, rmulti *y1)
{
  rmulti *Z0=NULL,*Z1=NULL,*a=NULL,*b=NULL;
  int p0,p1,prec,e=0;
  p0=rget_prec(z0); p1=rget_prec(z1); prec=MAX2(p0,p1);
  RA(Z0,prec); RA(Z1,prec); RA(a,prec); RA(b,prec);
  if(rget_sgn(x0)>0){
    if     (rget_sgn(y0)>0){ e+=abs(mpfr_mul(Z0,x0,y0,MPFR_RNDD)); e+=abs(mpfr_mul(Z1,x1,y1,MPFR_RNDU)); }
    else if(rget_sgn(y1)<0){ e+=abs(mpfr_mul(Z0,x1,y0,MPFR_RNDD)); e+=abs(mpfr_mul(Z1,x0,y1,MPFR_RNDU)); }
    else                   { e+=abs(mpfr_mul(Z0,x1,y0,MPFR_RNDD)); e+=abs(mpfr_mul(Z1,x1,y1,MPFR_RNDU)); }
  }else if(rget_sgn(x1)<0){
    if     (rget_sgn(y0)>0){ e+=abs(mpfr_mul(Z0,x0,y1,MPFR_RNDD)); e+=abs(mpfr_mul(Z1,x1,y0,MPFR_RNDU)); }
    else if(rget_sgn(y1)<0){ e+=abs(mpfr_mul(Z0,x1,y1,MPFR_RNDD)); e+=abs(mpfr_mul(Z1,x0,y0,MPFR_RNDU)); }
    else                   { e+=abs(mpfr_mul(Z0,x0,y1,MPFR_RNDD)); e+=abs(mpfr_mul(Z1,x0,y0,MPFR_RNDU)); }
  }else{
    if     (rget_sgn(y0)>0){ e+=abs(mpfr_mul(Z0,x0,y1,MPFR_RNDD)); e+=abs(mpfr_mul(Z1,x1,y1,MPFR_RNDU)); }
    else if(rget_sgn(y1)<0){ e+=abs(mpfr_mul(Z0,x1,y0,MPFR_RNDD)); e+=abs(mpfr_mul(Z1,x0,y0,MPFR_RNDU)); }
    else{
      e+=abs(mpfr_mul(a,x0,y1,MPFR_RNDD));
      e+=abs(mpfr_mul(b,x1,y0,MPFR_RNDD));
      if(rle(a,b)){ e+=abs(mpfr_set(Z0,a,MPFR_RNDD)); }else{ e+=abs(mpfr_set(Z0,b,MPFR_RNDD)); }
      e+=abs(mpfr_mul(a,x0,y0,MPFR_RNDU));
      e+=abs(mpfr_mul(b,x1,y1,MPFR_RNDU));
      if(rge(a,b)){ e+=abs(mpfr_set(Z1,a,MPFR_RNDU)); }else{ e+=abs(mpfr_set(Z1,b,MPFR_RNDU)); }
    }
  }
  e+=ircopy(z0,z1,Z0,Z1);
  RF(Z0); RF(Z1); RF(a); RF(b);
  return e;
}


//追加
/**
 @brief 掛け算 [z0,z1]=[x0,x1]*[y0,y1]
*/
int irmul_ws(rmulti *z0, rmulti *z1, rmulti *x0, rmulti *x1, rmulti *y0, rmulti *y1, int *rwss, rmulti **rws)
{
  int e=0;
    if(*rwss<4){ ERROR_EXIT("Error `rwss=%d<4' in irmul_ws().\n",*rwss); }
  if(rget_sgn(x0)>0){
    if     (rget_sgn(y0)>0){ e+=abs(mpfr_mul(rws[0],x0,y0,MPFR_RNDD)); e+=abs(mpfr_mul(rws[1],x1,y1,MPFR_RNDU)); }
    else if(rget_sgn(y1)<0){ e+=abs(mpfr_mul(rws[0],x1,y0,MPFR_RNDD)); e+=abs(mpfr_mul(rws[1],x0,y1,MPFR_RNDU)); }
    else                   { e+=abs(mpfr_mul(rws[0],x1,y0,MPFR_RNDD)); e+=abs(mpfr_mul(rws[1],x1,y1,MPFR_RNDU)); }
  }else if(rget_sgn(x1)<0){
    if     (rget_sgn(y0)>0){ e+=abs(mpfr_mul(rws[0],x0,y1,MPFR_RNDD)); e+=abs(mpfr_mul(rws[1],x1,y0,MPFR_RNDU)); }
    else if(rget_sgn(y1)<0){ e+=abs(mpfr_mul(rws[0],x1,y1,MPFR_RNDD)); e+=abs(mpfr_mul(rws[1],x0,y0,MPFR_RNDU)); }
    else                   { e+=abs(mpfr_mul(rws[0],x0,y1,MPFR_RNDD)); e+=abs(mpfr_mul(rws[1],x0,y0,MPFR_RNDU)); }
  }else{
    if     (rget_sgn(y0)>0){ e+=abs(mpfr_mul(rws[0],x0,y1,MPFR_RNDD)); e+=abs(mpfr_mul(rws[1],x1,y1,MPFR_RNDU)); }
    else if(rget_sgn(y1)<0){ e+=abs(mpfr_mul(rws[0],x1,y0,MPFR_RNDD)); e+=abs(mpfr_mul(rws[1],x0,y0,MPFR_RNDU)); }
    else{
      e+=abs(mpfr_mul(rws[2],x0,y1,MPFR_RNDD));
      e+=abs(mpfr_mul(rws[3],x1,y0,MPFR_RNDD));
      if(rle(rws[2],rws[3])){ e+=abs(mpfr_set(rws[0],rws[2],MPFR_RNDD)); }else{ e+=abs(mpfr_set(rws[0],rws[3],MPFR_RNDD)); }
      e+=abs(mpfr_mul(rws[2],x0,y0,MPFR_RNDU));
      e+=abs(mpfr_mul(rws[3],x1,y1,MPFR_RNDU));
      if(rge(rws[2],rws[3])){ e+=abs(mpfr_set(rws[1],rws[2],MPFR_RNDU)); }else{ e+=abs(mpfr_set(rws[1],rws[3],MPFR_RNDU)); }
    }
  }
  e+=ircopy(z0,z1,rws[0],rws[1]);
  return e;
}
//ここまで
/**
 @brief 掛け算 [z0,z1]=[x0,x1]*[y,y]
*/
int irmul_d(rmulti *z0, rmulti *z1, rmulti *x0, rmulti *x1, double y)
{
  int p0,p1,prec,e=0;
  rmulti *a0=NULL,*a1=NULL;
  p0=rget_prec(z0); p1=rget_prec(z1); prec=MAX2(p0,p1);
  RA(a0,prec); RA(a1,prec);
  e+=irset_d(a0,a1,y);
  e+=irmul(z0,z1,x0,x1,a0,a1);
  RF(a0); RF(a1);
  return e;
}

/**
 @brief 割り算 [z0,z1]=[x0,x1]/[y0,y1]
*/
int irdiv(rmulti *z0, rmulti *z1, rmulti *x0, rmulti *x1, rmulti *y0, rmulti *y1)
{
  rmulti *Z0=NULL,*Z1=NULL;
  int p0,p1,prec,e=0;
  p0=rget_prec(z0); p1=rget_prec(z1); prec=MAX2(p0,p1);
  RA(Z0,prec); RA(Z1,prec);
  if(rget_sgn(x0)>0){
    if     (rget_sgn(y0)>0){ e+=abs(mpfr_div(Z0,x0,y1,MPFR_RNDD)); e+=abs(mpfr_div(Z1,x1,y0,MPFR_RNDU)); }
    else if(rget_sgn(y1)<0){ e+=abs(mpfr_div(Z0,x1,y1,MPFR_RNDD)); e+=abs(mpfr_div(Z1,x0,y0,MPFR_RNDU)); }
    else                   { rset_nan(z0); rset_nan(z1); }//{ERROR_AT; mpfr_printf("Devided by 0, y=[%.3Re %.3Re].\n",y0,y1); }
  }else if(rget_sgn(x1)<0){
    if     (rget_sgn(y0)>0){ e+=abs(mpfr_div(Z0,x0,y0,MPFR_RNDD)); e+=abs(mpfr_div(Z1,x1,y1,MPFR_RNDU)); }
    else if(rget_sgn(y1)<0){ e+=abs(mpfr_div(Z0,x1,y0,MPFR_RNDD)); e+=abs(mpfr_div(Z1,x0,y1,MPFR_RNDU)); }
    else                   { rset_nan(z0); rset_nan(z1); }//{ERROR_AT; mpfr_printf("Devided by 0, y=[%.3Re %.3Re].\n",y0,y1); }
  }else{
    if     (rget_sgn(y0)>0){ e+=abs(mpfr_div(Z0,x0,y0,MPFR_RNDD)); e+=abs(mpfr_div(Z1,x1,y0,MPFR_RNDU)); }
    else if(rget_sgn(y1)<0){ e+=abs(mpfr_div(Z0,x1,y1,MPFR_RNDD)); e+=abs(mpfr_div(Z1,x0,y1,MPFR_RNDU)); }
    else                   { rset_nan(z0); rset_nan(z1); }//{ERROR_AT; mpfr_printf("Devided by 0, y=[%.3Re %.3Re].\n",y0,y1); }
  }
  e+=ircopy(z0,z1,Z0,Z1);
  RF(Z0); RF(Z1);
  return e;
}

//編集済み

/**
 @brief 逆数 [z0,z1]=[1,1]/[x0,x1]
*/
int irinv(rmulti *z0, rmulti *z1, rmulti *x0, rmulti *x1)
{
  int e=0;
  if (rget_sgn(x0)<0 && rget_sgn(x1)>0){ rset_nan(z0); rset_nan(z1); }//{ERROR_AT; mpfr_printf("Divided by 0, y=[%.3Re %.3Re].\n",x0,x1); }
  else {
    e+=abs(mpfr_d_div(z0,1,x1,MPFR_RNDD));
    e+=abs(mpfr_d_div(z1,1,x0,MPFR_RNDU));
  } 
  return e;
}
//ここまで

/**
 @brief 積の加算 [z0,z1]+=[x0,x1]*[y0,y1]
*/
int iradd_mul(rmulti *z0, rmulti *z1, rmulti *x0, rmulti *x1, rmulti *y0, rmulti *y1)
{
  int p0,p1,prec,e=0;
  rmulti *a0=NULL,*a1=NULL;
  p0=rget_prec(z0); p1=rget_prec(z1); prec=MAX2(p0,p1);
  RA(a0,prec); RA(a1,prec);
  e+=irmul(a0,a1,x0,x1,y0,y1); // a=x*y
  e+=iradd(z0,z1,z0,z1,a0,a1); // z=z+x*y
  RF(a0); RF(a1);
  return e;
}

//追加
/**
 @brief 積の加算 [z0,z1]+=[x0,x1]*[y0,y1]
*/
int iradd_mul_ws(rmulti *z0, rmulti *z1, rmulti *x0, rmulti *x1, rmulti *y0, rmulti *y1, int *rwss, rmulti **rws)
{
  if(*rwss<6){ ERROR_EXIT("Error `rwss=%d<6' in iradd_mul_ws().\n",*rwss); }
  int e=0;
  *rwss=*rwss-2;
  e+=irmul_ws(rws[0],rws[1],x0,x1,y0,y1,rwss,&rws[2]); // a=x*y
  e+=iradd(z0,z1,z0,z1,rws[0],rws[1]); // z=z+x*y
  *rwss=*rwss+2;
  return e;
}

//ここまで
/**
 @brief 積の減算 [z0,z1]-=[x0,x1]*[y0,y1]
*/
int irsub_mul(rmulti *z0, rmulti *z1, rmulti *x0, rmulti *x1, rmulti *y0, rmulti *y1)
{
  int p0,p1,prec,e=0;
  rmulti *a0=NULL,*a1=NULL;
  p0=rget_prec(z0); p1=rget_prec(z1); prec=MAX2(p0,p1);
  RA(a0,prec); RA(a1,prec);
  e+=irmul(a0,a1,x0,x1,y0,y1); // a=x*y
  e+=irsub(z0,z1,z0,z1,a0,a1); // z=z-x*y
  RF(a0); RF(a1);
  return e;
}

//追加
/**
 @brief 積の減算 [z0,z1]-=[x0,x1]*[y0,y1]
*/
int irsub_mul_ws(rmulti *z0, rmulti *z1, rmulti *x0, rmulti *x1, rmulti *y0, rmulti *y1, int *rwss, rmulti **rws)
{
  int e=0;
  if(*rwss<6){ ERROR_EXIT("Error `rwss=%d<8' in irsub_mul_ws().\n",*rwss); }
  *rwss=*rwss-2;
  e+=irmul_ws(rws[0],rws[1],x0,x1,y0,y1,rwss,&rws[2]); // a=x*y
  e+=irsub_ws(z0,z1,z0,z1,rws[0],rws[1],rwss,&rws[2]); // z=z-x*y
  *rwss=*rwss+2;
  return e;
}
//ここまで

/** @} */

//////////////////////////////////////////////////////

/** @name 数学関数 */
/** @{ */
//変更
/**
 @brief 絶対値 [y0,y1]=abs([x0,x1])
*/
int irabs(rmulti *y0, rmulti *y1, rmulti *x0, rmulti *x1)
{
  if(ris_nan(x0) || ris_nan(x1)){rset_nan(y0); rset_nan(y1); return 0;
  }else{
  int e=0;
  if(rget_sgn(x0)>=0 && rget_sgn(x1)>=0){         //x0>=0 && x1>=0
    e+=ircopy(y0,y1,x0,x1);   
  }else if(rget_sgn(x0)<=0 && rget_sgn(x1)<=0){   //x0<=0 && x1<=0
    e+=irneg(y0,y1,x0,x1);   
  }else{                                          //x0<0 && x1>0
    int p0,p1,prec;
    rmulti *a0=NULL,*a1=NULL;
    p0=rget_prec(y0); p1=rget_prec(y1); prec=MAX2(p0,p1);
    RA(a0,prec);RA(a1,prec);
    e+=abs(mpfr_neg(a0,x0,MPFR_RNDU));
    e+=rmax2(a1,a0,x1);
    e+=abs(mpfr_set_d(y0,0,MPFR_RNDD));
    e+=abs(mpfr_set(y1,a1, MPFR_RNDU));
    RF(a0);RF(a1);
  }
  return e;
  }
}
//ここまで

/**
 @brief 差の絶対値 [z0,z1]=abs([x0,x1]-[y0,y1])
*/
int irabs_sub(rmulti *z0, rmulti *z1, rmulti *x0, rmulti *x1, rmulti *y0, rmulti *y1)
{
  int e=0;
  e+=irsub(z0,z1,x0,x1,y0,y1);
  e+=irabs(z0,z1,z0,z1);
  return e;
}

/**
 @brief 商の分母が絶対値 [z0,z1]=[x0,x1]/abs([y0,y1])
*/
int irdiv_abs(rmulti *z0, rmulti *z1, rmulti *x0, rmulti *x1, rmulti *y0, rmulti *y1)
{
  int e=0;
  e+=irabs(z0,z1,y0,y1);
  e+=irdiv(z0,z1,x0,x1,z0,z1);
  return e; 
}

/**
 @brief 商の分母が絶対値 [z0,z1]=[x0,x1]/abs([y0,y1])
*/
int irdiv_abs_c(rmulti *z0, rmulti *z1, rmulti *x0, rmulti *x1, cmulti *y0, cmulti *y1)
{
  int e=0;
  e+=icabs(z0,z1,y0,y1);
  e+=irdiv(z0,z1,x0,x1,z0,z1);
  return e; 
}

/**
 @brief べき乗 [y0,y1]=[x0,x1]^n
 */
int irpow_si(rmulti *y0, rmulti *y1, rmulti *x0, rmulti *x1, long n)
{
  rmulti *a=NULL,*b0=NULL,*b1=NULL;
  int p0,p1,prec,e=0;
  p0=rget_prec(y0); p1=rget_prec(y1); prec=MAX2(p0,p1);
  RA(a,prec); RA(b0,prec); RA(b1,prec);
  if(n==0){
    e+=abs(mpfr_pow_si(y0,x0,n,MPFR_RNDD));
    e+=abs(mpfr_pow_si(y1,x1,n,MPFR_RNDU));
  }else if(n>0 && (n%2)>0){
    e+=abs(mpfr_pow_si(y0,x0,n,MPFR_RNDD));
    e+=abs(mpfr_pow_si(y1,x1,n,MPFR_RNDU));
  }else if(n>0){
    if(rget_sgn(x0)>0){
      e+=abs(mpfr_pow_si(y0,x0,n,MPFR_RNDD));
      e+=abs(mpfr_pow_si(y1,x1,n,MPFR_RNDU));
    }else if(rget_sgn(x1)<0){
      e+=abs(mpfr_pow_si(y0,x1,n,MPFR_RNDD));
      e+=abs(mpfr_pow_si(y1,x0,n,MPFR_RNDU));
    }else{
      e+=abs(mpfr_set_d(a,0,MPFR_RNDD));
      e+=abs(mpfr_pow_si(b0,x0,n,MPFR_RNDU));
      e+=abs(mpfr_pow_si(b1,x1,n,MPFR_RNDU));
      if(rge(b0,b1)){ e+=ircopy(y0,y1,a,b0); }else{ e+=ircopy(y0,y1,a,b1); }
    }
  }else{
    e+=rset_d(a,1);
    e+=irpow_si(b0,b1,x0,x1,-n);
    e+=irdiv(y0,y1,a,a,b0,b1);
  }
  RF(a); RF(b0); RF(b1);
  return e;
}

/**
 @brief 平方根 [y0,y1]=sqrt([x0,x1])
*/
int irsqrt(rmulti *y0, rmulti *y1, rmulti *x0, rmulti *x1)
{
  int e=0;
  e+=abs(mpfr_sqrt(y0,x0,MPFR_RNDD));
  e+=abs(mpfr_sqrt(y1,x1,MPFR_RNDU));
  return e;
}


/**
 @brief 指数関数 [y0,y1]=exp([x0,x1])
*/
int irexp(rmulti *y0, rmulti *y1, rmulti *x0, rmulti *x1)
{
  int e=0;
  e+=abs(mpfr_exp(y0,x0,MPFR_RNDD));
  e+=abs(mpfr_exp(y1,x1,MPFR_RNDU));
  return e;
}

/**
 @brief 対数関数 [y0,y1]=log([x0,x1])
*/
int irlog(rmulti *y0, rmulti *y1, rmulti *x0, rmulti *x1)
{
  int e=0;
  e+=abs(mpfr_log(y0,x0,MPFR_RNDD));
  e+=abs(mpfr_log(y1,x1,MPFR_RNDU));
  return e;
}

//追加

/**
 @brief 対数関数 [y0,y1]=log([x0,x1])
*/
int irlog10(rmulti *y0, rmulti *y1, rmulti *x0, rmulti *x1)
{
  int e=0;
  e+=abs(mpfr_log10(y0,x0,MPFR_RNDD));
  e+=abs(mpfr_log10(y1,x1,MPFR_RNDU));
  return e;
}

//ここまで

/**
 @brief 三角関数 [y0,y1]=sin([x0,x1])
 */
int irsin(rmulti *y0, rmulti *y1, rmulti *x0, rmulti *x1)
{
  int prec,e=0;
  rmulti *p0=NULL,*p1=NULL,*a=NULL,*b0=NULL,*b1=NULL;
  prec=MAX2(rget_prec(y0),rget_prec(y1));
  RA(p0,prec); RA(p1,prec); RA(a,prec); RA(b0,prec); RA(b1,prec);
  mpfr_cos(p0,x0,MPFR_RNDN);
  mpfr_cos(p1,x1,MPFR_RNDN);
  if(rget_sgn(p0)>=0 && rget_sgn(p1)>=0){
    e+=abs(mpfr_sin(y0,x0,MPFR_RNDD));
    e+=abs(mpfr_sin(y1,x1,MPFR_RNDU));
  }else if(rget_sgn(p0)<=0 && rget_sgn(p1)<=0){
    e+=abs(mpfr_sin(y0,x1,MPFR_RNDD));
    e+=abs(mpfr_sin(y1,x0,MPFR_RNDU));
  }else if(rget_sgn(p0)>=0 && rget_sgn(p1)<=0){
    e+=abs(mpfr_sin(b0,x0,MPFR_RNDD));
    e+=abs(mpfr_sin(b1,x1,MPFR_RNDD));
    e+=abs(mpfr_set_d(a,1,MPFR_RNDU));
    if(rle(b0,b1)){ e+=ircopy(y0,y1,b0,a); }else{ e+=ircopy(y0,y1,b1,a); }
  }else if(rget_sgn(p0)<=0 && rget_sgn(p1)>=0){
    e+=abs(mpfr_set_d(a,-1,MPFR_RNDD));
    e+=abs(mpfr_sin(b0,x0,MPFR_RNDU));
    e+=abs(mpfr_sin(b1,x1,MPFR_RNDU));
    if(rge(b0,b1)){ e+=ircopy(y0,y1,a,b0); }else{ e+=ircopy(y0,y1,a,b1); }
  }else{ ERROR_AT; exit(0); }
  RF(p0); RF(p1); RF(a); RF(b0); RF(b1);
  return e;
}

/**
 @brief 三角関数 [y0,y1]=cos([x0,x1])
 */
int ircos(rmulti *y0, rmulti *y1, rmulti *x0, rmulti *x1)
{
  int prec,e=0;
  rmulti *p0=NULL,*p1=NULL,*a=NULL,*b0=NULL,*b1=NULL;
  prec=MAX2(rget_prec(y0),rget_prec(y1));
  RA(p0,prec); RA(p1,prec); RA(a,prec); RA(b0,prec); RA(b1,prec);
  mpfr_sin(p0,x0,MPFR_RNDN); mpfr_neg(p0,p0,MPFR_RNDN);
  mpfr_sin(p1,x1,MPFR_RNDN); mpfr_neg(p1,p1,MPFR_RNDN);
  if(rget_sgn(p0)>=0 && rget_sgn(p1)>=0){
    e+=abs(mpfr_cos(y0,x0,MPFR_RNDD));
    e+=abs(mpfr_cos(y1,x1,MPFR_RNDU));
  }else if(rget_sgn(p0)<=0 && rget_sgn(p1)<=0){
    e+=abs(mpfr_cos(y0,x1,MPFR_RNDD));
    e+=abs(mpfr_cos(y1,x0,MPFR_RNDU));
  }else if(rget_sgn(p0)>=0 && rget_sgn(p1)<=0){
    e+=abs(mpfr_cos(b0,x0,MPFR_RNDD));
    e+=abs(mpfr_cos(b1,x1,MPFR_RNDD));
    e+=abs(mpfr_set_d(a,1,MPFR_RNDU));
    if(rle(b0,b1)){ e+=ircopy(y0,y1,b0,a); }else{ e+=ircopy(y0,y1,b1,a); }
  }else if(rget_sgn(p0)<=0 && rget_sgn(p1)>=0){
    e+=abs(mpfr_set_d(a,-1,MPFR_RNDD));
    e+=abs(mpfr_cos(b0,x0,MPFR_RNDU));
    e+=abs(mpfr_cos(b1,x1,MPFR_RNDU));
    if(rge(b0,b1)){ e+=ircopy(y0,y1,a,b0); }else{ e+=ircopy(y0,y1,a,b1); }
  }else{ ERROR_AT; exit(0); }
  RF(p0); RF(p1); RF(a); RF(b0); RF(b1);
  return e;
}

/**
 @brief 三角関数 [y0,y1]=tan([x0,x1])
 */
int irtan(rmulti *y0, rmulti *y1, rmulti *x0, rmulti *x1)
{
  int e=0;
  e+=abs(mpfr_tan(y0,x0,MPFR_RNDD));
  e+=abs(mpfr_tan(y1,x1,MPFR_RNDU));
  return e;
}

/**
 @brief 逆三角関数 [y0,y1]=asin([x0,x1])
*/
int irasin(rmulti *y0, rmulti *y1, rmulti *x0, rmulti *x1)
{
  int e=0;
  e+=abs(mpfr_asin(y0,x0,MPFR_RNDD));
  e+=abs(mpfr_asin(y1,x1,MPFR_RNDU));
  return e;
}

/**
 @brief 逆三角関数 [y0,y1]=acos([x0,x1])
*/
int iracos(rmulti *y0, rmulti *y1, rmulti *x0, rmulti *x1)
{
  int prec,e=0;
  rmulti *a0=NULL,*a1=NULL;
  prec=MAX2(rget_prec(y0),rget_prec(y1));
  RA(a0,prec); RA(a1,prec);
  e+=abs(mpfr_acos(a0,x1,MPFR_RNDD));
  e+=abs(mpfr_acos(a1,x0,MPFR_RNDU));
  e+=ircopy(y0,y1,a0,a1);
  RF(a0); RF(a1);
  return e;
}

/**
 @brief 逆三角関数 [y0,y1]=atan([x0,x1])
*/
int iratan(rmulti *y0, rmulti *y1, rmulti *x0, rmulti *x1)
{
  int e=0;
  e+=abs(mpfr_atan(y0,x0,MPFR_RNDD));
  e+=abs(mpfr_atan(y1,x1,MPFR_RNDU));
  return e;
}

/**
 @brief 逆三角関数 [z0,z1]=atan2([y0,y1]/[x0,x1])
*/
int iratan2(rmulti *z0, rmulti *z1, rmulti *y0, rmulti *y1, rmulti *x0, rmulti *x1)
{
  int prec,e=0;
  rmulti *a0=NULL,*a1=NULL;
  prec=MAX2(rget_prec(z0),rget_prec(z1));
  RA(a0,prec); RA(a1,prec);
  if(rget_sgn(y0)>=0){
    e+=abs(mpfr_atan2(a0,y0,x1,MPFR_RNDD));
    e+=abs(mpfr_atan2(a1,y1,x0,MPFR_RNDU));
  }else{
    e+=abs(mpfr_atan2(a0,y0,x0,MPFR_RNDD));
    e+=abs(mpfr_atan2(a1,y1,x1,MPFR_RNDU));
  }
  e+=ircopy(z0,z1,a0,a1);
  RF(a0); RF(a1);
  return e;
}

/**
 @brief 双曲線関数 [y0,y1]=sinh([x0,x1])
*/
int irsinh(rmulti *y0, rmulti *y1, rmulti *x0, rmulti *x1)
{
  int e=0;
  e+=abs(mpfr_sinh(y0,x0,MPFR_RNDD));
  e+=abs(mpfr_sinh(y1,x1,MPFR_RNDU));
  return e;
}

/**
 @brief 双曲線関数 [y0,y1]=cosh([x0,x1])
*/
int ircosh(rmulti *y0, rmulti *y1, rmulti *x0, rmulti *x1)
{
  rmulti *a=NULL,*b0=NULL,*b1=NULL;
  int prec,e=0;
  prec=MAX2(rget_prec(y0),rget_prec(y1));
  RA(a,prec); RA(b0,prec); RA(b1,prec);
  if(rget_sgn(x0)>0){
    e+=abs(mpfr_cosh(y0,x0,MPFR_RNDD));
    e+=abs(mpfr_cosh(y1,x1,MPFR_RNDU));
  }else if(rget_sgn(x1)<0){
    e+=abs(mpfr_cosh(b0,x1,MPFR_RNDD));
    e+=abs(mpfr_cosh(b1,x0,MPFR_RNDU));
    e+=ircopy(y0,y1,b0,b1);
  }else{
    e+=abs(mpfr_set_d(a,1,MPFR_RNDD));
    e+=abs(mpfr_cosh(b0,x0,MPFR_RNDU));
    e+=abs(mpfr_cosh(b1,x1,MPFR_RNDU));
    if(rge(b0,b1)){ e+=ircopy(y0,y1,a,b0); }else{ e+=ircopy(y0,y1,a,b1); }
  }
  RF(a); RF(b0); RF(b1);
  return e;
}

/**
 @brief 双曲線関数 [y0,y1]=tanh([x0,x1])
*/
int irtanh(rmulti *y0, rmulti *y1, rmulti *x0, rmulti *x1)
{
  int e=0;
  e+=abs(mpfr_tanh(y0,x0,MPFR_RNDD));
  e+=abs(mpfr_tanh(y1,x1,MPFR_RNDU));
  return e;
}

/**
 @beif 逆双曲線関数 [y0,y1]=asinh([x0,x1])
*/
int irasinh(rmulti *y0, rmulti *y1, rmulti *x0, rmulti *x1)
{
  int e=0;
  e+=abs(mpfr_asinh(y0,x0,MPFR_RNDD));
  e+=abs(mpfr_asinh(y1,x1,MPFR_RNDU));
  return e;
}

/**
 @brief 逆双曲線関数 [y0,y1]=acosh([x0,x1])
*/
int iracosh(rmulti *y0, rmulti *y1, rmulti *x0, rmulti *x1)
{
  int e=0;
  e+=abs(mpfr_acosh(y0,x0,MPFR_RNDD));
  e+=abs(mpfr_acosh(y1,x1,MPFR_RNDU));
  return e;
}

/**
 @brief 逆双曲線関数 [y0,y1]=atanh([x0,x1])
*/
int iratanh(rmulti *y0, rmulti *y1, rmulti *x0, rmulti *x1)
{
  int e=0;
  e+=abs(mpfr_atanh(y0,x0,MPFR_RNDD));
  e+=abs(mpfr_atanh(y1,x1,MPFR_RNDU));
  return e;
}

/** @} */

//////////////////////////////////////////////////////

/** @name 判定 */
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
