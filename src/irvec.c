#include"is_svec.h"
#include"is_rmulti.h"
#include"is_cmulti.h"
#include"is_rvec.h"
#include"is_cvec.h"
#include"is_rmat.h"
#include"is_cmat.h"
#include"is_irmulti.h"
#include"is_icmulti.h"
#include"is_irvec.h"
#include"is_icvec.h"
#include"is_func.h"

/**
 @file  irvec.c
 @brief 多倍長精度実数型rmultiの機械区間演算のベクトルに関する関数の定義.
 @details スカラーに関しては@link irmulti.c@endlinkを参照のこと.
          ベクトルに関しては@link irvec.c@endlinkを参照のこと.
          行列に関しては@link irmat.c@endlinkを参照のこと.
 */

#define RA(X,Y)      ((X)=rallocate_prec(rget_prec(Y)))
#define RA2(X,Y0,Y1) ((X)=rallocate_prec(MAX2(rget_prec(Y0),rget_prec(Y1))))
#define RAP(X,Z,N)   ((X)=rallocate_prec(rvec_get_prec_max(N,Z)))
#define RF(X)        ((X)=rmfree(X))
#define RVA(X,Y,N)   { X=rvec_allocate_prec(N,rvec_get_prec_max(N,Y)); }
#define RVF(X,N)     { X=rvec_free(N,X); }


/** @name irmulti型ベクトルの値の設定に関する関数 */
/** @{ */

/**
 @brief irmulti型のベクトルの値の設定
 */
void irvec_set_svec(int n, rmulti **y0, rmulti **y1, char **x)
{
  int i;
  for(i=0; i<n; i++){ irset_s(y0[i],y1[i],x[i]); }
}

/**
 @brief irmulti型のベクトルの値の設定
 */
void irvec_set_ivec(int n, rmulti **y0, rmulti **y1, int *x0, int *x1)
{
  int i;
  for(i=0; i<n; i++){ irset_d(y0[i],y1[i],x0[i],x1[i]); }
}

/**
 @brief irmulti型のベクトルの値の設定
 */
void irvec_set_rvec(int n, rmulti **y0, rmulti **y1, rmulti **x0, rmulti **x1)
{
  int i;
  for(i=0; i<n; i++){ irset_r(y0[i],y1[i],x0[i],x1[i]); }
}

/**
 @brief rmulti型のベクトルの値を設定.
 */
void rvec_set_irvec(int n, rmulti **y, rmulti **x0, rmulti **x1)
{
  int i;
  for(i=0; i<n; i++){ irmid(y[i],x0[i],x1[i]); }
}

/**
 @brief cmulti型のベクトルの値を設定.
 */
void cvec_set_irvec(int n, cmulti **y, rmulti **x0, rmulti **x1)
{
  int i;
  for(i=0; i<n; i++){ irmid(C_R(y[i]),x0[i],x1[i]); rset_d(C_I(y[i]),0); }
}

/**
 @brief irmulti型のベクトルの値の設定
*/
void irvec_set_dvec(int n, rmulti **y0, rmulti **y1, double *x0, double *x1)
{
  int i;
  for(i=0; i<n; i++){ irset_d(y0[i],y1[i],x0[i],x1[i]); }
}

/**
 @brief irmulti型のベクトルの値の設定
*/
void irvec_set_zvec(int n, rmulti **y0, rmulti **y1, dcomplex *x0, dcomplex *x1)
{
  int i;
  for(i=0; i<n; i++){ irset_d(y0[i],y1[i],Z_R(x0[i]),Z_R(x1[i])); }
}

/**
 @brief irmulti型のベクトルの値の設定
*/
void irvec_set_cvec(int n, rmulti **y0, rmulti **y1, cmulti **x0, cmulti **x1)
{
  int i;
  for(i=0; i<n; i++){ irset_r(y0[i],y1[i],C_R(x0[i]),C_R(x1[i])); }
}

// y=x
void dvec_set_irvec(int n, double *y, rmulti **x0, rmulti **x1)
{
  int i;
  for(i=0; i<n; i++){ y[i]=0.5*(rget_d(x0[i])+rget_d(x1[i])); }
}

// y=x
void zvec_set_irvec(int n, dcomplex *y, rmulti **x0, rmulti **x1)
{
  double a;
  int i;
  for(i=0; i<n; i++){
    a=0.5*(rget_d(x0[i])+rget_d(x1[i]));
    Z_SET(y[i],a,0);
  }
}

///////////////////////////////////////////////////////////////////////

/** @* */
/** @name irmulti型ベクトルの型変換に関する関数 */
/** @{ */

/**
 @brief irmulti型のベクトルを整数型に変換.
 */
void irvec_get_ivec(int n, int *y, rmulti **x0, rmulti **x1)
{
  int i;
  for(i=0; i<n; i++){ y[i]=0.5*(rget_d(x0[i])+rget_d(x1[i])); }
}

/**
 @brief irmulti型のベクトルを倍精度実数型に変換.
 */
void irvec_get_dvec(int n, double *y, rmulti **x0, rmulti **x1)
{
  int i;
  for(i=0; i<n; i++){ y[i]=0.5*(rget_d(x0[i])+rget_d(x1[i])); }
}

/**
 @brief irmulti型のベクトルを文字列型に変換 y=char(x)
 */
void irvec_get_svec(int n, char **y, rmulti **x0, rmulti **x1, char format, int digits)
{
  char f[1024],buf[1<<13];
  int i;
  if(format=='e'){ mpfr_sprintf(f,"[%%+.%dR%c, %%+.%dR%c]",digits,format,digits,format); }
  else           { mpfr_sprintf(f,"[%%.%dR%c, %%.%dR%c]",digits,format,digits,format); }
  for(i=0; i<n; i++){
    mpfr_sprintf(buf,f,x0[i],x1[i]);
    y[i]=char_renew(y[i],buf,NULL);
  }
}

/** @} */
/** @name irmulti型ベクトルの入出力に関する関数 */
/** @{ */

/**
 @brief irmulti型のベクトルの表示
 */
void irvec_print(int n, rmulti **x0, rmulti **x1, char *name, char format, int digits)
{
  char **s=NULL;
  if(x0==NULL || x1==NULL){
    if(name!=NULL){ printf("%s=NULL\n",name); }
    else          { printf("NULL\n"); }
    return;
  }
  s=svec_allocate(n);
  irvec_get_svec(n,s,x0,x1,format,digits);
  svec_print(n,s,name);
  s=svec_free(n,s);
}


/** @} */
/** @name irmulti型ベクトルの型変換に関する関数 */
/** @{ */

/** @} */
/** @name irmulti型ベクトルに関する関数子 */
/** @{ */

/**
 @brief irmulti型のベクトルのコピー [y0,y1]=[x0,x1]
 */
void irvec_copy_rvec(int n, rmulti **y0, rmulti **y1, rmulti **x0, rmulti **x1)
{
  int i;
  for(i=0; i<n; i++){ irset_r(y0[i],y1[i],x0[i],x1[i]); }
}


/**
 @brief irmulti型のベクトルの区間の中心 xc=(x1+x0)/2 と半径 xr=x1-x0
*/
void irvec_mid_rad_rvec(int n, rmulti **xc, rmulti **xr, rmulti **x0, rmulti **x1)
{
  mpfr_rnd_t mode;
  mode=get_round_mode();
  set_round_mode(MPFR_RNDU);  // up
  rvec_sub_rvec_rvec(n,xc,x1,x0);    // xc=x1-x0
  rvec_mul_rvec_dscalar(n,xc,xc,0.5); // xc=(x1-x0)/2
  rvec_add_rvec_rvec(n,xc,xc,x0);    // xc=(x1-x0)/2+x0
  rvec_sub_rvec_rvec(n,xr,xc,x0);    // xr=xc-x0
  set_round_mode(mode);       // back
}

/**
 @brief irmulti型のベクトルの符号の反転 [y0,y1]=-[x0,x1]
*/
void irvec_neg_rvec(int n, rmulti **y0, rmulti **y1, rmulti **x0, rmulti **x1)
{
  int i;
  for(i=0; i<n; i++){ irneg_r(y0[i],y1[i],x0[i],x1[i]); }
}

/**
 @brief irmulti型のベクトルの符号の正負 [y0,y1]=[-abs(x),abs(x)]
*/
void irvec_pm_rvec(int n, rmulti **y0, rmulti **y1, rmulti **x)
{
  int i;
  for(i=0; i<n; i++){ irpm(y0[i],y1[i],x[i]); }
}

/**
 @brief irmulti型のベクトルの区間の拡張 [z0,z1]=[x0,x1]+[-y,y]
*/
void irvec_add_pm_rvec(int n, rmulti **z0, rmulti **z1, rmulti **x0, rmulti **x1, rmulti **y)
{
  int i;
  for(i=0; i<n; i++){ iradd_pm_rr(z0[i],z1[i],x0[i],x1[i],y[i]); }
}

/**
 @brief irmulti型のベクトルの区間の中心 y=mid([x0,x1])
 */
void irvec_get_rvec(int n, rmulti **y, rmulti **x0, rmulti **x1)
{
  int i;
  for(i=0; i<n; i++){ irmid(y[i],x0[i],x1[i]); }
}

/**
 @brief irmulti型のベクトルの区間の中心 y=mid([x0,x1])
 */
void irvec_get_cvec(int n, cmulti **y, rmulti **x0, rmulti **x1)
{
  int i;
  for(i=0; i<n; i++){ irmid(C_R(y[i]),x0[i],x1[i]); rset_d(C_I(y[i]),0); }
}

/**
 @brief irmulti型のベクトルの区間の中心 [m-r,m+r]=[x0,x1]
 */
void irvec_mid_rvec(int n, rmulti **mid, rmulti **x0, rmulti **x1)
{
  int i;
  for(i=0; i<n; i++){ irmid(mid[i],x0[i],x1[i]); }
}

/**
 @brief irmulti型のベクトルの区間の半径 [m-r,m+r]=[x0,x1]
 */
void irvec_rad_rvec(int n, rmulti **rad, rmulti **x0, rmulti **x1)
{
  int i;
  for(i=0; i<n; i++){ irrad(rad[i],x0[i],x1[i]); }
}

/**
 @brief irmulti型のベクトルの区間の中心と半径 [m-r,m+r]=[x0,x1]
 */
void irvec_mr_rvec(int n, rmulti **mid, rmulti **rad, rmulti **x0, rmulti **x1)
{
  int i;
  for(i=0; i<n; i++){ irmr(mid[i],rad[i],x0[i],x1[i]); }
}

///////////////////////////////////////////////////////

/**
 @brief [z0,z1]=[x0,z1]+[y0,y1]
 */
void irvec_add_rvec_rvec(int n, rmulti **z0, rmulti **z1, rmulti **x0, rmulti **x1, rmulti **y0, rmulti **y1)
{
  int i;
  for(i=0; i<n; i++){ iradd_rr(z0[i],z1[i],x0[i],x1[i],y0[i],y1[i]); }
}

/**
 @brief [z0,z1]=[x0,z1]+[y0,y1]
 */
void irvec_add_rvec_dvec(int n, rmulti **z0, rmulti **z1, rmulti **x0, rmulti **x1, double *y0, double *y1)
{
  int i;
  for(i=0; i<n; i++){ iradd_rd(z0[i],z1[i],x0[i],x1[i],y0[i],y1[i]); }
}

/**
 @brief [z0,z1]=[x0,z1]+[y0,y1]
 */
void irvec_add_dvec_rvec(int n, rmulti **z0, rmulti **z1, double *x0, double *x1, rmulti **y0, rmulti **y1)
{
    int i;
    for(i=0; i<n; i++){ iradd_dr(z0[i],z1[i],x0[i],x1[i],y0[i],y1[i]); }
}

/**
 @brief [z0,z1]=[x0,z1]+[y0,y1]
 */
void irvec_add_rvec_rscalar(int n, rmulti **z0, rmulti **z1, rmulti **x0, rmulti **x1, rmulti *y0, rmulti *y1)
{
  int i;
  for(i=0; i<n; i++){ iradd_rr(z0[i],z1[i],x0[i],x1[i],y0,y1); }
}

/**
 @brief [z0,z1]=[x0,z1]+[y0,y1]
 */
void irvec_add_rvec_dscalar(int n, rmulti **z0, rmulti **z1, rmulti **x0, rmulti **x1, double y0, double y1)
{
  int i;
  for(i=0; i<n; i++){ iradd_rd(z0[i],z1[i],x0[i],x1[i],y0,y1); }
}

/**
 @brief [z0,z1]=[x0,z1]+[y0,y1]
 */
void irvec_add_dvec_rscalar(int n, rmulti **z0, rmulti **z1, double *x0, double *x1, rmulti *y0, rmulti *y1)
{
  int i;
  for(i=0; i<n; i++){ iradd_dr(z0[i],z1[i],x0[i],x1[i],y0,y1); }
}

/**
 @brief [z0,z1]=[x0,z1]+[y0,y1]
 */
void irvec_add_rscalar_rvec(int n, rmulti **z0, rmulti **z1, rmulti *x0, rmulti *x1, rmulti **y0, rmulti **y1)
{
  int i;
  for(i=0; i<n; i++){ iradd_rr(z0[i],z1[i],x0,x1,y0[i],y1[i]); }
}

/**
 @brief [z0,z1]=[x0,z1]+[y0,y1]
 */
void irvec_add_rscalar_dvec(int n, rmulti **z0, rmulti **z1, rmulti *x0, rmulti *x1, double *y0, double *y1)
{
  int i;
  for(i=0; i<n; i++){ iradd_rd(z0[i],z1[i],x0,x1,y0[i],y1[i]); }
}

/**
 @brief [z0,z1]=[x0,z1]+[y0,y1]
 */
void irvec_add_dscalar_rvec(int n, rmulti **z0, rmulti **z1, double x0, double x1, rmulti **y0, rmulti **y1)
{
  int i;
  for(i=0; i<n; i++){ iradd_dr(z0[i],z1[i],x0,x1,y0[i],y1[i]); }
}

/////////////////////////////////////////////////////////////////////////////

/**
 @brief [z0,z1]=[x0,z1]-[y0,y1]
 */
void irvec_sub_rvec_rvec(int n, rmulti **z0, rmulti **z1, rmulti **x0, rmulti **x1, rmulti **y0, rmulti **y1)
{
  int i;
  for(i=0; i<n; i++){ irsub_rr(z0[i],z1[i],x0[i],x1[i],y0[i],y1[i]); }
}

/**
 @brief [z0,z1]=[x0,z1]-[y0,y1]
 */
void irvec_sub_rvec_dvec(int n, rmulti **z0, rmulti **z1, rmulti **x0, rmulti **x1, double *y0, double *y1)
{
  int i;
  for(i=0; i<n; i++){ irsub_rd(z0[i],z1[i],x0[i],x1[i],y0[i],y1[i]); }
}

/**
 @brief [z0,z1]=[x0,z1]-[y0,y1]
 */
void irvec_sub_dvec_rvec(int n, rmulti **z0, rmulti **z1, double *x0, double *x1, rmulti **y0, rmulti **y1)
{
  int i;
  for(i=0; i<n; i++){ irsub_dr(z0[i],z1[i],x0[i],x1[i],y0[i],y1[i]); }
}

/**
 @brief [z0,z1]=[x0,z1]-[y0,y1]
 */
void irvec_sub_rvec_rscalar(int n, rmulti **z0, rmulti **z1, rmulti **x0, rmulti **x1, rmulti *y0, rmulti *y1)
{
  int i;
  for(i=0; i<n; i++){ irsub_rr(z0[i],z1[i],x0[i],x1[i],y0,y1); }
}

/**
 @brief [z0,z1]=[x0,z1]-[y0,y1]
 */
void irvec_sub_rvec_dscalar(int n, rmulti **z0, rmulti **z1, rmulti **x0, rmulti **x1, double y0, double y1)
{
  int i;
  for(i=0; i<n; i++){ irsub_rd(z0[i],z1[i],x0[i],x1[i],y0,y1); }
}

/**
 @brief [z0,z1]=[x0,z1]-[y0,y1]
 */
void irvec_sub_dvec_rscalar(int n, rmulti **z0, rmulti **z1, double *x0, double *x1, rmulti *y0, rmulti *y1)
{
  int i;
  for(i=0; i<n; i++){ irsub_dr(z0[i],z1[i],x0[i],x1[i],y0,y1); }
}

/**
 @brief [z0,z1]=[x0,z1]-[y0,y1]
 */
void irvec_sub_rscalar_rvec(int n, rmulti **z0, rmulti **z1, rmulti *x0, rmulti *x1, rmulti **y0, rmulti **y1)
{
  int i;
  for(i=0; i<n; i++){ irsub_rr(z0[i],z1[i],x0,x1,y0[i],y1[i]); }
}

/**
 @brief [z0,z1]=[x0,z1]-[y0,y1]
 */
void irvec_sub_rscalar_dvec(int n, rmulti **z0, rmulti **z1, rmulti *x0, rmulti *x1, double *y0, double *y1)
{
  int i;
  for(i=0; i<n; i++){ irsub_rd(z0[i],z1[i],x0,x1,y0[i],y1[i]); }
}

/**
 @brief [z0,z1]=[x0,z1]-[y0,y1]
 */
void irvec_sub_dscalar_rvec(int n, rmulti **z0, rmulti **z1, double x0, double x1, rmulti **y0, rmulti **y1)
{
  int i;
  for(i=0; i<n; i++){ irsub_dr(z0[i],z1[i],x0,x1,y0[i],y1[i]); }
}

/////////////////////////////////////////////////////////////////////////////

/**
 @brief [z0,z1]=[x0,x1]*[y0,y1]
 */
void irvec_mul_rvec_rvec(int n, rmulti **z0, rmulti **z1, rmulti **x0, rmulti **x1, rmulti **y0, rmulti **y1)
{
  int i;
  for(i=0; i<n; i++){ irmul_rr(z0[i],z1[i],x0[i],x1[i],y0[i],y1[i]); }
}

/**
 @brief [z0,z1]=[x0,x1]*[y0,y1]
 */
void irvec_mul_rvec_dvec(int n, rmulti **z0, rmulti **z1, rmulti **x0, rmulti **x1, double *y0, double *y1)
{
  int i;
  for(i=0; i<n; i++){ irmul_rd(z0[i],z1[i],x0[i],x1[i],y0[i],y1[i]); }
}

/**
 @brief [z0,z1]=[x0,x1]*[y0,y1]
 */
void irvec_mul_dvec_rvec(int n, rmulti **z0, rmulti **z1, double *x0, double *x1, rmulti **y0, rmulti **y1)
{
  int i;
  for(i=0; i<n; i++){ irmul_dr(z0[i],z1[i],x0[i],x1[i],y0[i],y1[i]); }
}

/**
 @brief [z0,z1]=[x0,x1]*[y0,y1]
 */
void irvec_mul_rvec_rscalar(int n, rmulti **z0, rmulti **z1, rmulti **x0, rmulti **x1, rmulti *y0, rmulti *y1)
{
  int i;
  for(i=0; i<n; i++){ irmul_rr(z0[i],z1[i],x0[i],x1[i],y0,y1); }
}

/**
 @brief [z0,z1]=[x0,x1]*[y0,y1]
 */
void irvec_mul_rvec_dscalar(int n, rmulti **z0, rmulti **z1, rmulti **x0, rmulti **x1, double y0, double y1)
{
  int i;
  for(i=0; i<n; i++){ irmul_rd(z0[i],z1[i],x0[i],x1[i],y0,y1); }
}

/**
 @brief [z0,z1]=[x0,x1]*[y0,y1]
 */
void irvec_mul_dvec_rscalar(int n, rmulti **z0, rmulti **z1, double *x0, double *x1, rmulti *y0, rmulti *y1)
{
  int i;
  for(i=0; i<n; i++){ irmul_dr(z0[i],z1[i],x0[i],x1[i],y0,y1); }
}

/**
 @brief [z0,z1]=[x0,x1]*[y0,y1]
 */
void irvec_mul_rscalar_rvec(int n, rmulti **z0, rmulti **z1, rmulti *x0, rmulti *x1, rmulti **y0, rmulti **y1)
{
  int i;
  for(i=0; i<n; i++){ irmul_rr(z0[i],z1[i],x0,x1,y0[i],y1[i]); }
}

/**
 @brief [z0,z1]=[x0,x1]*[y0,y1]
 */
void irvec_mul_rscalar_dvec(int n, rmulti **z0, rmulti **z1, rmulti *x0, rmulti *x1, double *y0, double *y1)
{
  int i;
  for(i=0; i<n; i++){ irmul_rd(z0[i],z1[i],x0,x1,y0[i],y1[i]); }
}

/**
 @brief [z0,z1]=[x0,x1]*[y0,y1]
 */
void irvec_mul_dscalar_rvec(int n, rmulti **z0, rmulti **z1, double x0, double x1, rmulti **y0, rmulti **y1)
{
  int i;
  for(i=0; i<n; i++){ irmul_dr(z0[i],z1[i],x0,x1,y0[i],y1[i]); }
}

/////////////////////////////////////////////////////////////////////////////

/**
 @brief [z0,z1]=[x0,x1]/[y0,y1]
 */
void irvec_div_rvec_rvec(int n, rmulti **z0, rmulti **z1, rmulti **x0, rmulti **x1, rmulti **y0, rmulti **y1)
{
  int i;
  for(i=0; i<n; i++){ irdiv_rr(z0[i],z1[i],x0[i],x1[i],y0[i],y1[i]); }
}

/**
 @brief [z0,z1]=[x0,x1]/[y0,y1]
 */
void irvec_div_rvec_dvec(int n, rmulti **z0, rmulti **z1, rmulti **x0, rmulti **x1, double *y0, double *y1)
{
  int i;
  for(i=0; i<n; i++){ irdiv_rd(z0[i],z1[i],x0[i],x1[i],y0[i],y1[i]); }
}

/**
 @brief [z0,z1]=[x0,x1]/[y0,y1]
 */
void irvec_div_dvec_rvec(int n, rmulti **z0, rmulti **z1, double *x0, double *x1, rmulti **y0, rmulti **y1)
{
  int i;
  for(i=0; i<n; i++){ irdiv_dr(z0[i],z1[i],x0[i],x1[i],y0[i],y1[i]); }
}

/**
 @brief [z0,z1]=[x0,x1]/[y0,y1]
 */
void irvec_div_rvec_rscalar(int n, rmulti **z0, rmulti **z1, rmulti **x0, rmulti **x1, rmulti *y0, rmulti *y1)
{
  int i;
  for(i=0; i<n; i++){ irdiv_rr(z0[i],z1[i],x0[i],x1[i],y0,y1); }
}

/**
 @brief [z0,z1]=[x0,x1]/[y0,y1]
 */
void irvec_div_rvec_dscalar(int n, rmulti **z0, rmulti **z1, rmulti **x0, rmulti **x1, double y0, double y1)
{
  int i;
  for(i=0; i<n; i++){ irdiv_rd(z0[i],z1[i],x0[i],x1[i],y0,y1); }
}

/**
 @brief [z0,z1]=[x0,x1]/[y0,y1]
 */
void irvec_div_dvec_rscalar(int n, rmulti **z0, rmulti **z1, double *x0, double *x1, rmulti *y0, rmulti *y1)
{
  int i;
  for(i=0; i<n; i++){ irdiv_dr(z0[i],z1[i],x0[i],x1[i],y0,y1); }
}

/**
 @brief [z0,z1]=[x0,x1]/[y0,y1]
 */
void irvec_div_rscalar_rvec(int n, rmulti **z0, rmulti **z1, rmulti *x0, rmulti *x1, rmulti **y0, rmulti **y1)
{
  int i;
  for(i=0; i<n; i++){ irdiv_rr(z0[i],z1[i],x0,x1,y0[i],y1[i]); }
}

/**
 @brief [z0,z1]=[x0,x1]/[y0,y1]
 */
void irvec_div_rscalar_dvec(int n, rmulti **z0, rmulti **z1, rmulti *x0, rmulti *x1, double *y0, double *y1)
{
  int i;
  for(i=0; i<n; i++){ irdiv_rd(z0[i],z1[i],x0,x1,y0[i],y1[i]); }
}

/**
 @brief [z0,z1]=[x0,x1]/[y0,y1]
 */
void irvec_div_dscalar_rvec(int n, rmulti **z0, rmulti **z1, double x0, double x1, rmulti **y0, rmulti **y1)
{
  int i;
  for(i=0; i<n; i++){ irdiv_dr(z0[i],z1[i],x0,x1,y0[i],y1[i]); }
}

/////////////////////////////////////////////////////////////////////////////


/**
 @brief irmulti型のベクトルの平方の総和 [y0,y1]=sum([x0,x1].^2)
*/
void irsum_pow2_abs_rvec(rmulti *y0, rmulti *y1, int n, rmulti **x0, rmulti **x1)
{
  int i;
  irset_d(y0,y1,0,0);
  for(i=0; i<n; i++){ iradd_mul_rr(y0,y1,x0[i],x1[i],x0[i],x1[i]); }
}

/**
 @brief irmulti型のベクトルの絶対値 [y0,y1]=abs([x0,x1])
*/
void irvec_abs_rvec(int n, rmulti **y0, rmulti **y1, rmulti **x0, rmulti **x1)
{
  int i;
  for(i=0; i<n; i++){ irabs_r(y0[i],y1[i],x0[i],x1[i]); }
}

/**
 @brief irmulti型のベクトルの差の絶対値 [z0,z1]=abs([x0,x1]-[y0,y1])
*/
void irvec_abs_sub_rvec_rvec(int n, rmulti **z0, rmulti **z1, rmulti **x0, rmulti **x1, rmulti **y0, rmulti **y1)
{
  int i;
  for(i=0; i<n; i++){ irabs_sub_rr(z0[i],z1[i],x0[i],x1[i],y0[i],y1[i]); }
}

/**
 @brief irmulti型のベクトルの最大値 [y0,y1]=[max(x0),max(x1)]
*/
void irmax_rvec(rmulti *y0, rmulti *y1, int n, rmulti **x0, rmulti **x1)
{
  int i;
  irset_r(y0,y1,x0[0],x1[0]);        // y0=x0[0], y1=x1[0]
  for(i=1; i<n; i++){
    if(gt_rr(x0[i],y0)){              // x0[i]>y0
      mpfr_set(y0,x0[i],MPFR_RNDD); // y0=x0[i]            
    }
    if(gt_rr(x1[i],y1)){              // x1[i]>y1
      mpfr_set(y1,x1[i],MPFR_RNDU); // y1=x1[i]              
    }
  }
}

/**
 @brief irmulti型のベクトルの最大値 [y0,y1]=[x0,max(x1)] 上限で比較
*/
void irmax_up_rvec(rmulti *y0, rmulti *y1, int n, rmulti **x0, rmulti **x1)
{
  int i;
  irset_r(y0,y1,x0[0],x1[0]);      // y0=x0[0], y1=x1[0]
  for(i=1; i<n; i++){
    if(gt_rr(x1[i],y1)){            // x1[i]>y1
      irset_r(y0,y1,x0[i],x1[i]);  // y0=x0[i], y1=x1[i]
    }
  }
}

/**
 @brief irmulti型のベクトルの最大値 [y0,y1]=[max(x0),x1] 下限で比較
*/
void irmax_down_rvec(rmulti *y0, rmulti *y1, int n, rmulti **x0, rmulti **x1)
{
  int i;
  irset_r(y0,y1,x0[0],x1[0]);      // y0=x0[0], y1=x1[0]
  for(i=1; i<n; i++){
    if(gt_rr(x0[i],y0)){            // x0[i]>y0
      irset_r(y0,y1,x0[i],x1[i]);  // y0=x0[i], y1=x1[i]
    }
  }
}

/**
 @brief irmulti型のベクトルの最小値 [y0,y1]=[min(x0),min(x1)]
*/
void irmin_rvec(rmulti *y0, rmulti *y1, int n, rmulti **x0, rmulti **x1)
{
  int i;
  irset_r(y0,y1,x0[0],x1[0]);        // y0=x0[0], y1=x1[0]
  for(i=1; i<n; i++){
    if(lt_rr(x0[i],y0)){              // x0[i]<y0
      mpfr_set(y0,x0[i],MPFR_RNDD); // y0=x0[i]        
    }
    if(lt_rr(x1[i],y1)){              // x1[i]<y1
      mpfr_set(y1,x1[i],MPFR_RNDD); // y0=x0[i]       
    }
  }
}

/**
 @brief irmulti型のベクトルの最小値 [y0,y1]=[x0,min(x1)]
*/
void irmin_up_rvec(rmulti *y0, rmulti *y1, int n, rmulti **x0, rmulti **x1)
{
  int i;
  irset_r(y0,y1,x0[0],x1[0]);     // y0=x0[0], y1=x1[0]
  for(i=1; i<n; i++){
    if(lt_rr(x1[i],y1)){           // x0[i]<y0
      irset_r(y0,y1,x0[i],x1[i]); // y0=x0[i], y1=x1[i]
    }
  }
}

/**
 @brief irmulti型のベクトルの最小値 [y0,y1]=[min(x0),x1]
*/
void irmin_down_rvec(rmulti *y0, rmulti *y1, int n, rmulti **x0, rmulti **x1)
{
  int i;
  irset_r(y0,y1,x0[0],x1[0]);     // y0=x0[0], y1=x1[0]
  for(i=1; i<n; i++){
    if(lt_rr(x0[i],y0)){           // x0[i]<y0
      irset_r(y0,y1,x0[i],x1[i]); // y0=x0[i], y1=x1[i]
    }
  }
}

/**
 @brief irmulti型のベクトルの要素の絶対値の最大値 value=max(abs(x))
*/
void irmax_up_abs_rvec(rmulti *y0, rmulti *y1, int n, rmulti **x0, rmulti **x1)
{
  int prec;
  rmulti **ax0=NULL,**ax1=NULL;
  prec=MAX2(rget_prec(y0),rget_prec(y1));
  ax0=rvec_allocate_prec(n,prec);
  ax1=rvec_allocate_prec(n,prec);
  irvec_abs_rvec(n,ax0,ax1,x0,x1);
  irmax_up_rvec(y0,y1,n,ax0,ax1);
  RVF(ax0,n); RVF(ax1,n);
}

/**
 @brief irmulti型のベクトルの要素の絶対値の最小値 value=min(abs(x))
*/
void irmin_down_abs_rvec(rmulti *y0, rmulti *y1, int n, rmulti **x0, rmulti **x1)
{
  int prec;
  rmulti **ax0=NULL,**ax1=NULL;
  prec=MAX2(rget_prec(y0),rget_prec(y1));
  ax0=rvec_allocate_prec(n,prec);
  ax1=rvec_allocate_prec(n,prec);
  irvec_abs_rvec(n,ax0,ax1,x0,x1);
  irmin_down_rvec(y0,y1,n,ax0,ax1);
  RVF(ax0,n); RVF(ax1,n);
}

/**
 @brief irmulti型のベクトルの要素の総和 value=sum(x)
 */
void irsum_rvec(rmulti *y0, rmulti *y1, int n, rmulti **x0, rmulti **x1)
{
  int i;
  irset_d(y0,y1,0,0);
  for(i=0; i<n; i++){ iradd_rr(y0,y1,y0,y1,x0[i],x1[i]); }
}

/**
 @brief irmulti型のベクトルの線形変換 [y0,y1]=[A0,A1]*[x0,x1]
*/
void irvec_mul_rmat_rvec(int m, int n, rmulti **y0, rmulti **y1, rmulti **A0, int LDA0, rmulti **A1, int LDA1, rmulti **x0, rmulti **x1)
{
  int i,j;
  rmulti **z0=NULL,**z1=NULL;
  RVA(z0,y0,m); RVA(z1,y1,m);
  for(i=0; i<m; i++){
    irset_d(z0[i],z1[i],0,0);
    for(j=0; j<n; j++){
      iradd_mul_rr(z0[i],z1[i],MAT(A0,i,j,LDA0),MAT(A1,i,j,LDA1),x0[j],x1[j]);
    }
  }
  irvec_copy_rvec(m,y0,y1,z0,z1);
  RVF(z0,m); RVF(z1,m);
}

/**
 @brief irmulti型のベクトルの線形変換の加算 [y0,y1]+=[A0,A1]*[x0,x1]
*/
void irvec_add_lintr(int m, int n, rmulti **y0, rmulti **y1, rmulti **A0, int LDA0, rmulti **A1, int LDA1, rmulti **x0, rmulti **x1)
{
  int i,j;
  rmulti **z0=NULL,**z1=NULL;
  RVA(z0,y0,m); RVA(z1,y1,m);
  for(i=0; i<m; i++){
    irset_r(z0[i],z1[i],y0[i],y1[i]);
    for(j=0; j<n; j++){
      iradd_mul_rr(z0[i],z1[i],MAT(A0,i,j,LDA0),MAT(A1,i,j,LDA1),x0[j],x1[j]);
    }
  }
  irvec_copy_rvec(m,y0,y1,z0,z1);
  RVF(z0,m); RVF(z1,m);
}

/**
 @brief irmulti型のベクトルの線形変換の減算 [y0,y1]-=[A0,A1]*[x0,x1]
*/
void irvec_sub_lintr(int m, int n, rmulti **y0, rmulti **y1, rmulti **A0, int LDA0, rmulti **A1, int LDA1, rmulti **x0, rmulti **x1)
{
  int i,j;
  rmulti **z0=NULL,**z1=NULL;
  RVA(z0,y0,m); RVA(z1,y1,m);
  for(i=0; i<m; i++){
    irset_r(z0[i],z1[i],y0[i],y1[i]);
    for(j=0; j<n; j++){
      irsub_mul_rr(z0[i],z1[i],MAT(A0,i,j,LDA0),MAT(A1,i,j,LDA1),x0[j],x1[j]);
    }
  }
  irvec_copy_rvec(m,y0,y1,z0,z1);
  RVF(z0,m); RVF(z1,m);
}

/**
 @brief irmulti型のベクトルの転置行列の線形変換 [y0,y1]=[A0,A1]'*[x0,x1]
*/
void irvec_lintr_t(int m, int n, rmulti **y0, rmulti **y1, rmulti **A0, int LDA0, rmulti **A1, int LDA1, rmulti **x0, rmulti **x1)
{
  int i,j;
  rmulti **z0=NULL,**z1=NULL;
  RVA(z0,y0,n); RVA(z1,y1,n);
  for(j=0; j<n; j++){
    irset_d(z0[j],z1[j],0,0);
    for(i=0; i<m; i++){
      iradd_mul_rr(z0[j],z1[j],MAT(A0,i,j,LDA0),MAT(A1,i,j,LDA1),x0[i],x1[i]);
    }
  }
  irvec_copy_rvec(n,y0,y1,z0,z1);
  RVF(z0,n); RVF(z1,n);
}

/**
 @brief irmulti型のベクトルの転置行列の線形変換の加算 [y0,y1]+=[A0,A1]'*[x0,x1]
*/
void irvec_add_lintr_t(int m, int n, rmulti **y0, rmulti **y1, rmulti **A0, int LDA0, rmulti **A1, int LDA1, rmulti **x0, rmulti **x1)
{
  int i,j;
  rmulti **z0=NULL,**z1=NULL;
  RVA(z0,y0,n); RVA(z1,y1,n);
  for(j=0; j<n; j++){
    irset_r(z0[j],z1[j],y0[j],y1[j]);
    for(i=0; i<m; i++){
      iradd_mul_rr(z0[j],z1[j],MAT(A0,i,j,LDA0),MAT(A1,i,j,LDA1),x0[i],x1[i]);
    }
  }
  irvec_copy_rvec(n,y0,y1,z0,z1);
  RVF(z0,n); RVF(z1,n);
}

/**
 @brief irmulti型のベクトルの転置行列の線形変換の減算 [y0,y1]-=[A0,A1]'*[x0,x1]
*/
void irvec_sub_lintr_t(int m, int n, rmulti **y0, rmulti **y1, rmulti **A0, int LDA0, rmulti **A1, int LDA1, rmulti **x0, rmulti **x1)
{
  int i,j;
  rmulti **z0=NULL,**z1=NULL;
  RVA(z0,y0,n); RVA(z1,y1,n);
  for(j=0; j<n; j++){
    irset_r(z0[j],z1[j],y0[j],y1[j]);
    for(i=0; i<m; i++){
      irsub_mul_rr(z0[j],z1[j],MAT(A0,i,j,LDA0),MAT(A1,i,j,LDA1),x0[i],x1[i]);
    }
  }
  irvec_copy_rvec(n,y0,y1,z0,z1);
  RVF(z0,n); RVF(z1,n);
}

/** @} */

//EOF
