#include"is_strings.h"
#include"is_svec.h"
#include"is_zvec.h"
#include"is_cmulti.h"
#include"is_rvec.h"
#include"is_cvec.h"
#include"is_cmat.h"
#include"is_irmulti.h"
#include"is_icmulti.h"
#include"is_irvec.h"
#include"is_icvec.h"
#include"is_func.h"

/**
 @file  icvec.c
 @brief 多倍長精度実数型cmultiの機械区間演算のベクトルに関する関数の定義.
 @details スカラーに関しては@link icmulti.c@endlinkを参照のこと.
          ベクトルに関しては@link icvec.c@endlinkを参照のこと.
          行列に関しては@link icmat.c@endlinkを参照のこと.
 */

#define RA(X,P)    ((X)=rallocate_prec(P))
#define RF(X)      ((X)=rfree(X))
#define CA(X,P)    ((X)=callocate_prec(P))
#define CF(X)      ((X)=cfree(X))
#define CVA(X,N,P) { X=cvec_allocate_prec(N,P); }
#define CVF(X,N)   { X=cvec_free(N,X); }
#define RVA(X,Y,N) { X=rvec_allocate_prec(N,rvec_get_prec_max(N,Y)); }
#define RVF(X,N)   { X=rvec_free(N,X); }

/** @name icmulti型ベクトルの値の設定に関する関数 */
/** @{ */


/**
 @brief コピー [y0,y1]=x
 */
void icvec_set_svec(int n, cmulti **y0, cmulti **y1, char **x)
{
  int i;
  for(i=0; i<n; i++){ icset_s(y0[i],y1[i],x[i]); }
}

/**
 @brief コピー [y0,y1]=x
 */
void icvec_set_ivec(int n, cmulti **y0, cmulti **y1, int *x)
{
  int i;
  for(i=0; i<n; i++){ icset_d(y0[i],y1[i],x[i],x[i]); }
}

/**
 @brief コピー [y0,y1]=[x0,x1]
 */
void icvec_set_cvec(int n, cmulti **y0, cmulti **y1, cmulti **x0, cmulti **x1)
{
  int i;
  for(i=0; i<n; i++){ icset_c(y0[i],y1[i],x0[i],x1[i]); }
}

/**
 @brief コピー [y0,y1]=x
 */
void icvec_set_dvec(int n, cmulti **y0, cmulti **y1, double *x0, double *x1)
{
  int i;
  for(i=0; i<n; i++){ icset_d(y0[i],y1[i],x0[i],x1[i]); }
}

/**
 @brief コピー [y0,y1]=x
 */
void icvec_set_zvec(int n, cmulti **y0, cmulti **y1, dcomplex *x0, dcomplex *x1)
{
  int i;
  for(i=0; i<n; i++){ icset_z(y0[i],y1[i],x0[i],x1[i]); }
}

/**
 @brief コピー [y0,y1]=x
 */
void icvec_set_rvec(int n, cmulti **y0, cmulti **y1, rmulti **x0, rmulti **x1)
{
  int i;
  for(i=0; i<n; i++){ icset_r(y0[i],y1[i],x0[i],x1[i]); }
}

/**
 @brief rmulti型のベクトルの値を設定.
 */
void rvec_set_icvec(int n, rmulti **y, cmulti **x0, cmulti **x1)
{
  int i;
  for(i=0; i<n; i++){ irmid(y[i],C_R(x0[i]),C_R(x1[i])); }
}

/**
 @brief cmulti型のベクトルの値を設定.
 */
void cvec_set_icvec(int n, cmulti **y, cmulti **x0, cmulti **x1)
{
  int i;
  for(i=0; i<n; i++){ icmid(y[i],x0[i],x1[i]); }
}


/**
 @brief コピー [y0,y1]=x
 */
void icvec_set_dvec_zvec(int n, cmulti **y0, cmulti **y1, double *x0, dcomplex *x1)
{
  int i;
  for(i=0; i<n; i++){ icset_dz(y0[i],y1[i],x0[i],x1[i]); }
}

/**
 @brief コピー [y0,y1]=x
 */
void icvec_set_zvec_dvec(int n, cmulti **y0, cmulti **y1, dcomplex *x0, double *x1)
{
  int i;
  for(i=0; i<n; i++){ icset_zd(y0[i],y1[i],x0[i],x1[i]); }
}

/**
 @brief コピー [y0,y1]=[x0,x1]
 */
void icvec_set_cvec_rvec(int n, cmulti **y0, cmulti **y1, cmulti **x0, rmulti **x1)
{
  int i;
  for(i=0; i<n; i++){ icset_cr(y0[i],y1[i],x0[i],x1[i]); }
}

/**
 @brief コピー [y0,y1]=[x0,x1]
 */
void icvec_set_rvec_cvec(int n, cmulti **y0, cmulti **y1, rmulti **x0, cmulti **x1)
{
  int i;
  for(i=0; i<n; i++){ icset_rc(y0[i],y1[i],x0[i],x1[i]); }
}

// y=x
void dvec_set_icvec(int n, double *y, cmulti **x0, cmulti **x1)
{
  dcomplex z0,z1;
  int i;
  for(i=0; i<n; i++){
    z0=cget_z(x0[i]);
    z1=cget_z(x1[i]);
    Z_ADD(z0,z1);
    Z_SCALE(z0,0.5);
    y[i]=Z_R(z0);
  }
}

// y=x
void zvec_set_icvec(int n, dcomplex *y, cmulti **x0, cmulti **x1)
{
  dcomplex z0,z1;
  int i;
  for(i=0; i<n; i++){
    z0=cget_z(x0[i]);
    z1=cget_z(x1[i]);
    Z_ADD(z0,z1);
    Z_SCALE(z0,0.5);
    y[i]=z0;
  }
}

//////////////////////////////////////////////////////////////////////

/** @} */
/** @name icmulti型ベクトルの入出力に関する関数 */
/** @{ */

/**
 @brief 表示
 */
void icvec_print(int n, cmulti **x0, cmulti **x1, char *name, char format, int digits)
{
  char **s=NULL;
  if(x0==NULL || x1==NULL){
    if(name!=NULL){ printf("%s=NULL\n",name); }
    else          { printf("NULL\n"); }
    return;
  }
  s=svec_allocate(n);
  icvec_get_svec(n,s,x0,x1,format,digits);
  svec_print(n,s,name);
  s=svec_free(n,s);
}

/** @} */
/** @name icmulti型ベクトルの型変換に関する関数 */
/** @{ */

/**
 @brief cmulti型のベクトルを整数型に変換.
 */
void icvec_get_ivec(int n, int *y, cmulti **x0, cmulti **x1)
{
  dcomplex z0,z1;
  int i;
  for(i=0; i<n; i++){
    z0=cget_z(x0[i]);
    z1=cget_z(x1[i]);
    Z_ADD(z0,z1);
    Z_SCALE(z0,0.5);
    y[i]=Z_R(z0);
  }
}

/**
 @brief cmulti型のベクトルを倍精度複素数型に変換.
 */
void icvec_get_zvec(int n, dcomplex *y, cmulti **x0, cmulti **x1)
{
  dcomplex z0,z1;
  int i;
  for(i=0; i<n; i++){
    z0=cget_z(x0[i]);
    z1=cget_z(x1[i]);
    Z_ADD(z0,z1);
    Z_SCALE(z0,0.5);
    y[i]=z0;
  }
}

/**
 @brief icmulti型のベクトルを文字列型に変換 y=char(x)
 */
void icvec_get_svec(int n, char **y, cmulti **x0, cmulti **x1, char format, int digits)
{
  char f[1024],buf[1<<13];
  int i;
  if(format=='e'){ sprintf(f,"[%%+.%dR%c%%+.%dR%ci, %%+.%dR%c%%+.%dR%ci]",digits,format,digits,format,digits,format,digits,format); }
  else           { sprintf(f,"[%%-.%dR%c%%+.%dR%ci, %%-.%dR%c%%+.%dR%ci]",digits,format,digits,format,digits,format,digits,format); }
  for(i=0; i<n; i++){
    mpfr_sprintf(buf,f,C_R(x0[i]),C_I(x0[i]),C_R(x1[i]),C_I(x1[i]));
    y[i]=char_renew(y[i],buf,NULL);
  }
}

/**
 @brief 区間の中心 y=mid([x0,x1])
 */
void icvec_get_cvec(int n, cmulti **y, cmulti **x0, cmulti **x1)
{
  int i;
  for(i=0; i<n; i++){ icmid(y[i],x0[i],x1[i]); }
}


/** @} */
/** @name icmulti型ベクトルに関する関数子 */
/** @{ */

/**
 @brief コピー [y0,y1]=[x0,x1]
 */
void icvec_copy_cr(int n, cmulti **y0, cmulti **y1, cmulti **x0, rmulti **x1)
{
  int i;
  for(i=0; i<n; i++){ icset_cr(y0[i],y1[i],x0[i],x1[i]); }
}

/**
 @brief コピー [y0,y1]=[x0,x1]
 */
void icvec_copy_rc(int n, cmulti **y0, cmulti **y1, rmulti **x0, cmulti **x1)
{
  int i;
  for(i=0; i<n; i++){ icset_rc(y0[i],y1[i],x0[i],x1[i]); }
}

/**
 @brief コピー [y0,y1]=[x0,x1]
 */
void icvec_copy(int n, cmulti **y0, cmulti **y1, cmulti **x0, cmulti **x1)
{
  int i;
  for(i=0; i<n; i++){ icset_c(y0[i],y1[i],x0[i],x1[i]); }
}

/**
 @brief コピー [y0,y1]=[x0,x1]
 */
void icvec_copy_r(int n, cmulti **y0, cmulti **y1, rmulti **x0, rmulti **x1)
{
  int i;
  for(i=0; i<n; i++){ icset_r(y0[i],y1[i],x0[i],x1[i]); }
}

/**
 @brief コピー [y0,y1]=[x0,x1]
 */
void icvec_set_rvec_rvec(int n, cmulti **y0, cmulti **y1, rmulti **x0_r, rmulti **x0_i, rmulti **x1_r, rmulti **x1_i)
{
  int i;
  for(i=0; i<n; i++){ icset_rrrr(y0[i],y1[i],x0_r[i],x0_i[i],x1_r[i],x1_i[i]); }
}

/**
 @brief 区間の中心 xc=(x1+x0)/2 と半径 xr=x1-x0
*/
void icvec_center_radius(int n, cmulti **xc, cmulti **xr, cmulti **x0, cmulti **x1)
{
  mpfr_rnd_t mode;
  mode=get_round_mode();
  set_round_mode(MPFR_RNDU);  // up
  cvec_sub_cvec(n,xc,x1,x0);    // xc=x1-x0
  cvec_mul_dscalar(n,xc,xc,0.5); // xc=(x1-x0)/2
  cvec_add_cvec(n,xc,xc,x0);    // xc=(x1-x0)/2+x0
  cvec_sub_cvec(n,xr,xc,x0);    // xr=xc-x0
  set_round_mode(mode);       // back
}


/**
 @brief icmulti型の実部のコピー y=real(x)
 */
void icvec_real(int n, rmulti **y0, rmulti **y1, cmulti **x0, cmulti **x1)
{
  int i;
  for(i=0; i<n; i++){ irset_r(y0[i],y1[i],C_R(x0[i]),C_R(x1[i])); }
}

/**
 @brief [y0,y1]=conj([x0,x1])
*/
void icvec_conj(int n, cmulti **y0, cmulti **y1, cmulti **x0, cmulti **x1)
{
  int i;
  for(i=0; i<n; i++){ icconj_c(y0[i],y1[i],x0[i],x1[i]); }
}

/**
 @brief 符号の反転 [y0,y1]=-[x0,x1]
*/
void icvec_neg(int n, cmulti **y0, cmulti **y1, cmulti **x0, cmulti **x1)
{
  int i;
  for(i=0; i<n; i++){ icneg_c(y0[i],y1[i],x0[i],x1[i]); }
}

/**
 @brief 符号の正負 [y0,y1]=[-abs(x),abs(x)]
*/
void icvec_pm(int n, cmulti **y0, cmulti **y1, cmulti **x)
{
  int i;
  for(i=0; i<n; i++){ icpm(y0[i],y1[i],x[i]); }
}

/**
 @brief 区間の拡張 [z0,z1]=[x0,x1]+[-y,y]
*/
void icvec_add_pm(int n, cmulti **z0, cmulti **z1, cmulti **x0, cmulti **x1, cmulti **y)
{
  int i;
  for(i=0; i<n; i++){ icadd_pm_cc(z0[i],z1[i],x0[i],x1[i],y[i]); }
}

/**
 @brief 区間の中心 [m-r,m+r]=[x0,x1]
 */
void icvec_mid(int n, cmulti **mid, cmulti **x0, cmulti **x1)
{
  int i;
  for(i=0; i<n; i++){ icmid(mid[i],x0[i],x1[i]); }
}

/**
 @brief 区間の半径 [m-r,m+r]=[x0,x1]
 */
void icvec_rad(int n, cmulti **rad, cmulti **x0, cmulti **x1)
{
  int i;
  for(i=0; i<n; i++){ icrad(rad[i],x0[i],x1[i]); }
}

/**
 @brief 区間の中心と半径 [m-r,m+r]=[x0,x1]
 */
void icvec_mr(int n, cmulti **mid, cmulti **rad, cmulti **x0, cmulti **x1)
{
  int i;
  for(i=0; i<n; i++){ icmr(mid[i],rad[i],x0[i],x1[i]); }
}

////////////////////////////////////////////////

/**
 @brief [z0,z1]=[x0,x1]+[y0,y1]
 */
void icvec_add_cvec(int n, cmulti **z0, cmulti **z1, cmulti **x0, cmulti **x1, cmulti **y0, cmulti **y1)
{
  int i;
  for(i=0; i<n; i++){ icadd_cc(z0[i],z1[i],x0[i],x1[i],y0[i],y1[i]); }
}

/**
 @brief [z0,z1]=[x0,x1]+[y0,y1]
 */
void icvec_add_rvec(int n, cmulti **z0, cmulti **z1, cmulti **x0, cmulti **x1, rmulti **y0, rmulti **y1)
{
  int i;
  for(i=0; i<n; i++){ icadd_cr(z0[i],z1[i],x0[i],x1[i],y0[i],y1[i]); }
}

/**
 @brief [z0,z1]=[x0,x1]+[y0,y1]
 */
void irvec_add_cvec(int n, cmulti **z0, cmulti **z1, rmulti **x0, rmulti **x1, cmulti **y0, cmulti **y1)
{
  int i;
  for(i=0; i<n; i++){ icadd_rc(z0[i],z1[i],x0[i],x1[i],y0[i],y1[i]); }
}

/**
 @brief [z0,z1]=[x0,x1]+[y0,y1]
 */
void icvec_add_zvec(int n, cmulti **z0, cmulti **z1, cmulti **x0, cmulti **x1, dcomplex *y0, dcomplex *y1)
{
  int i;
  for(i=0; i<n; i++){ icadd_cz(z0[i],z1[i],x0[i],x1[i],y0[i],y1[i]); }
}

/**
 @brief [z0,z1]=[x0,x1]+[y0,y1]
 */
void izvec_add_cvec(int n, cmulti **z0, cmulti **z1, dcomplex *x0, dcomplex *x1, cmulti **y0, cmulti **y1)
{
  int i;
  for(i=0; i<n; i++){ icadd_zc(z0[i],z1[i],x0[i],x1[i],y0[i],y1[i]); }
}

/**
 @brief [z0,z1]=[x0,x1]+[y0,y1]
 */
void icvec_add_dvec(int n, cmulti **z0, cmulti **z1, cmulti **x0, cmulti **x1, double *y0, double *y1)
{
  int i;
  for(i=0; i<n; i++){ icadd_cd(z0[i],z1[i],x0[i],x1[i],y0[i],y1[i]); }
}

/**
 @brief [z0,z1]=[x0,x1]+[y0,y1]
 */
void idvec_add_cvec(int n, cmulti **z0, cmulti **z1, double *x0, double *x1, cmulti **y0, cmulti **y1)
{
  int i;
  for(i=0; i<n; i++){ icadd_dc(z0[i],z1[i],x0[i],x1[i],y0[i],y1[i]); }
}

/**
 @brief [z0,z1]=[x0,x1]+[y0,y1]
 */
void irvec_add_zvec(int n, cmulti **z0, cmulti **z1, rmulti **x0, rmulti **x1, dcomplex *y0, dcomplex *y1)
{
  int i;
  for(i=0; i<n; i++){ icadd_rz(z0[i],z1[i],x0[i],x1[i],y0[i],y1[i]); }
}

/**
 @brief [z0,z1]=[x0,x1]+[y0,y1]
 */
void izvec_add_rvec(int n, cmulti **z0, cmulti **z1, dcomplex *x0, dcomplex *x1, rmulti **y0, rmulti **y1)
{
  int i;
  for(i=0; i<n; i++){ icadd_zr(z0[i],z1[i],x0[i],x1[i],y0[i],y1[i]); }
}

/**
 @brief [z0,z1]=[x0,x1]+[y0,y1]
 */
void icvec_add_cscalar(int n, cmulti **z0, cmulti **z1, cmulti **x0, cmulti **x1, cmulti *y0, cmulti *y1)
{
  int i;
  for(i=0; i<n; i++){ icadd_cc(z0[i],z1[i],x0[i],x1[i],y0,y1); }
}

/**
 @brief [z0,z1]=[x0,x1]+[y0,y1]
 */
void icvec_add_rscalar(int n, cmulti **z0, cmulti **z1, cmulti **x0, cmulti **x1, rmulti *y0, rmulti *y1)
{
  int i;
  for(i=0; i<n; i++){ icadd_cr(z0[i],z1[i],x0[i],x1[i],y0,y1); }
}

/**
 @brief [z0,z1]=[x0,x1]+[y0,y1]
 */
void irvec_add_cscalar(int n, cmulti **z0, cmulti **z1, rmulti **x0, rmulti **x1, cmulti *y0, cmulti *y1)
{
  int i;
  for(i=0; i<n; i++){ icadd_rc(z0[i],z1[i],x0[i],x1[i],y0,y1); }
}

/**
 @brief [z0,z1]=[x0,x1]+[y0,y1]
 */
void icvec_add_zscalar(int n, cmulti **z0, cmulti **z1, cmulti **x0, cmulti **x1, dcomplex y0, dcomplex y1)
{
  int i;
  for(i=0; i<n; i++){ icadd_cz(z0[i],z1[i],x0[i],x1[i],y0,y1); }
}

/**
 @brief [z0,z1]=[x0,x1]+[y0,y1]
 */
void izvec_add_cscalar(int n, cmulti **z0, cmulti **z1, dcomplex *x0, dcomplex *x1, cmulti *y0, cmulti *y1)
{
  int i;
  for(i=0; i<n; i++){ icadd_zc(z0[i],z1[i],x0[i],x1[i],y0,y1); }
}

/**
 @brief [z0,z1]=[x0,x1]+[y0,y1]
 */
void icvec_add_dscalar(int n, cmulti **z0, cmulti **z1, cmulti **x0, cmulti **x1, double y0, double y1)
{
  int i;
  for(i=0; i<n; i++){ icadd_cd(z0[i],z1[i],x0[i],x1[i],y0,y1); }
}

/**
 @brief [z0,z1]=[x0,x1]+[y0,y1]
 */
void idvec_add_cscalar(int n, cmulti **z0, cmulti **z1, double *x0, double *x1, cmulti *y0, cmulti *y1)
{
  int i;
  for(i=0; i<n; i++){ icadd_dc(z0[i],z1[i],x0[i],x1[i],y0,y1); }
}

/**
 @brief [z0,z1]=[x0,x1]+[y0,y1]
 */
void irvec_add_zscalar(int n, cmulti **z0, cmulti **z1, rmulti **x0, rmulti **x1, dcomplex y0, dcomplex y1)
{
  int i;
  for(i=0; i<n; i++){ icadd_rz(z0[i],z1[i],x0[i],x1[i],y0,y1); }
}

/**
 @brief [z0,z1]=[x0,x1]+[y0,y1]
 */
void izvec_add_rscalar(int n, cmulti **z0, cmulti **z1, dcomplex *x0, dcomplex *x1, rmulti *y0, rmulti *y1)
{
  int i;
  for(i=0; i<n; i++){ icadd_zr(z0[i],z1[i],x0[i],x1[i],y0,y1); }
}

/**
 @brief [z0,z1]=[x0,x1]+[y0,y1]
 */
void icscalar_add_cvec(int n, cmulti **z0, cmulti **z1, cmulti *x0, cmulti *x1, cmulti **y0, cmulti **y1)
{
  int i;
  for(i=0; i<n; i++){ icadd_cc(z0[i],z1[i],x0,x1,y0[i],y1[i]); }
}

/**
 @brief [z0,z1]=[x0,x1]+[y0,y1]
 */
void icscalar_add_rvec(int n, cmulti **z0, cmulti **z1, cmulti *x0, cmulti *x1, rmulti **y0, rmulti **y1)
{
  int i;
  for(i=0; i<n; i++){ icadd_cr(z0[i],z1[i],x0,x1,y0[i],y1[i]); }
}

/**
 @brief [z0,z1]=[x0,x1]-[y0,y1]
 */
void irscalar_add_cvec(int n, cmulti **z0, cmulti **z1, rmulti *x0, rmulti *x1, cmulti **y0, cmulti **y1)
{
  int i;
  for(i=0; i<n; i++){ icadd_rc(z0[i],z1[i],x0,x1,y0[i],y1[i]); }
}

/**
 @brief [z0,z1]=[x0,x1]+[y0,y1]
 */
void icscalar_add_zvec(int n, cmulti **z0, cmulti **z1, cmulti *x0, cmulti *x1, dcomplex *y0, dcomplex *y1)
{
  int i;
  for(i=0; i<n; i++){ icadd_cz(z0[i],z1[i],x0,x1,y0[i],y1[i]); }
}

/**
 @brief [z0,z1]=[x0,x1]+[y0,y1]
 */
void izscalar_add_cvec(int n, cmulti **z0, cmulti **z1, dcomplex x0, dcomplex x1, cmulti **y0, cmulti **y1)
{
  int i;
  for(i=0; i<n; i++){ icadd_zc(z0[i],z1[i],x0,x1,y0[i],y1[i]); }
}

/**
 @brief [z0,z1]=[x0,x1]+[y0,y1]
 */
void icscalar_add_dvec(int n, cmulti **z0, cmulti **z1, cmulti *x0, cmulti *x1, double *y0, double *y1)
{
  int i;
  for(i=0; i<n; i++){ icadd_cd(z0[i],z1[i],x0,x1,y0[i],y1[i]); }
}

/**
 @brief [z0,z1]=[x0,x1]+[y0,y1]
 */
void idscalar_add_cvec(int n, cmulti **z0, cmulti **z1, double x0, double x1, cmulti **y0, cmulti **y1)
{
  int i;
  for(i=0; i<n; i++){ icadd_dc(z0[i],z1[i],x0,x1,y0[i],y1[i]); }
}

/**
 @brief [z0,z1]=[x0,x1]+[y0,y1]
 */
void irscalar_add_zvec(int n, cmulti **z0, cmulti **z1, rmulti *x0, rmulti *x1, dcomplex *y0, dcomplex *y1)
{
  int i;
  for(i=0; i<n; i++){ icadd_rz(z0[i],z1[i],x0,x1,y0[i],y1[i]); }
}

/**
 @brief [z0,z1]=[x0,x1]+[y0,y1]
 */
void izscalar_add_rvec(int n, cmulti **z0, cmulti **z1, dcomplex x0, dcomplex x1, rmulti **y0, rmulti **y1)
{
  int i;
  for(i=0; i<n; i++){ icadd_zr(z0[i],z1[i],x0,x1,y0[i],y1[i]); }
}

///////////////////////////////////////////////////

/**
 @brief [z0,z1]=[x0,x1]-[y0,y1]
 */
void icvec_sub_cvec(int n, cmulti **z0, cmulti **z1, cmulti **x0, cmulti **x1, cmulti **y0, cmulti **y1)
{
  int i;
  for(i=0; i<n; i++){ icsub_cc(z0[i],z1[i],x0[i],x1[i],y0[i],y1[i]); }
}

/**
 @brief [z0,z1]=[x0,x1]-[y0,y1]
 */
void icvec_sub_rvec(int n, cmulti **z0, cmulti **z1, cmulti **x0, cmulti **x1, rmulti **y0, rmulti **y1)
{
  int i;
  for(i=0; i<n; i++){ icsub_cr(z0[i],z1[i],x0[i],x1[i],y0[i],y1[i]); }
}

/**
 @brief [z0,z1]=[x0,x1]-[y0,y1]
 */
void irvec_sub_cvec(int n, cmulti **z0, cmulti **z1, rmulti **x0, rmulti **x1, cmulti **y0, cmulti **y1)
{
  int i;
  for(i=0; i<n; i++){ icsub_rc(z0[i],z1[i],x0[i],x1[i],y0[i],y1[i]); }
}

/**
 @brief [z0,z1]=[x0,x1]-[y0,y1]
 */
void icvec_sub_zvec(int n, cmulti **z0, cmulti **z1, cmulti **x0, cmulti **x1, dcomplex *y0, dcomplex *y1)
{
  int i;
  for(i=0; i<n; i++){ icsub_cz(z0[i],z1[i],x0[i],x1[i],y0[i],y1[i]); }
}

/**
 @brief [z0,z1]=[x0,x1]-[y0,y1]
 */
void izvec_sub_cvec(int n, cmulti **z0, cmulti **z1, dcomplex *x0, dcomplex *x1, cmulti **y0, cmulti **y1)
{
  int i;
  for(i=0; i<n; i++){ icsub_zc(z0[i],z1[i],x0[i],x1[i],y0[i],y1[i]); }
}

/**
 @brief [z0,z1]=[x0,x1]-[y0,y1]
 */
void icvec_sub_dvec(int n, cmulti **z0, cmulti **z1, cmulti **x0, cmulti **x1, double *y0, double *y1)
{
  int i;
  for(i=0; i<n; i++){ icsub_cd(z0[i],z1[i],x0[i],x1[i],y0[i],y1[i]); }
}

/**
 @brief [z0,z1]=[x0,x1]-[y0,y1]
 */
void idvec_sub_cvec(int n, cmulti **z0, cmulti **z1, double *x0, double *x1, cmulti **y0, cmulti **y1)
{
  int i;
  for(i=0; i<n; i++){ icsub_dc(z0[i],z1[i],x0[i],x1[i],y0[i],y1[i]); }
}

/**
 @brief [z0,z1]=[x0,x1]-[y0,y1]
 */
void irvec_sub_zvec(int n, cmulti **z0, cmulti **z1, rmulti **x0, rmulti **x1, dcomplex *y0, dcomplex *y1)
{
  int i;
  for(i=0; i<n; i++){ icsub_rz(z0[i],z1[i],x0[i],x1[i],y0[i],y1[i]); }
}

/**
 @brief [z0,z1]=[x0,x1]-[y0,y1]
 */
void izvec_sub_rvec(int n, cmulti **z0, cmulti **z1, dcomplex *x0, dcomplex *x1, rmulti **y0, rmulti **y1)
{
  int i;
  for(i=0; i<n; i++){ icsub_zr(z0[i],z1[i],x0[i],x1[i],y0[i],y1[i]); }
}

/**
 @brief [z0,z1]=[x0,x1]-[y0,y1]
 */
void icvec_sub_cscalar(int n, cmulti **z0, cmulti **z1, cmulti **x0, cmulti **x1, cmulti *y0, cmulti *y1)
{
  int i;
  for(i=0; i<n; i++){ icsub_cc(z0[i],z1[i],x0[i],x1[i],y0,y1); }
}

/**
 @brief [z0,z1]=[x0,x1]-[y0,y1]
 */
void icvec_sub_rscalar(int n, cmulti **z0, cmulti **z1, cmulti **x0, cmulti **x1, rmulti *y0, rmulti *y1)
{
  int i;
  for(i=0; i<n; i++){ icsub_cr(z0[i],z1[i],x0[i],x1[i],y0,y1); }
}

/**
 @brief [z0,z1]=[x0,x1]-[y0,y1]
 */
void irvec_sub_cscalar(int n, cmulti **z0, cmulti **z1, rmulti **x0, rmulti **x1, cmulti *y0, cmulti *y1)
{
  int i;
  for(i=0; i<n; i++){ icsub_rc(z0[i],z1[i],x0[i],x1[i],y0,y1); }
}

/**
 @brief [z0,z1]=[x0,x1]-[y0,y1]
 */
void icvec_sub_zscalar(int n, cmulti **z0, cmulti **z1, cmulti **x0, cmulti **x1, dcomplex y0, dcomplex y1)
{
  int i;
  for(i=0; i<n; i++){ icsub_cz(z0[i],z1[i],x0[i],x1[i],y0,y1); }
}

/**
 @brief [z0,z1]=[x0,x1]-[y0,y1]
 */
void izvec_sub_cscalar(int n, cmulti **z0, cmulti **z1, dcomplex *x0, dcomplex *x1, cmulti *y0, cmulti *y1)
{
  int i;
  for(i=0; i<n; i++){ icsub_zc(z0[i],z1[i],x0[i],x1[i],y0,y1); }
}

/**
 @brief [z0,z1]=[x0,x1]-[y0,y1]
 */
void icvec_sub_dscalar(int n, cmulti **z0, cmulti **z1, cmulti **x0, cmulti **x1, double y0, double y1)
{
  int i;
  for(i=0; i<n; i++){ icsub_cd(z0[i],z1[i],x0[i],x1[i],y0,y1); }
}

/**
 @brief [z0,z1]=[x0,x1]-[y0,y1]
 */
void idvec_sub_cscalar(int n, cmulti **z0, cmulti **z1, double *x0, double *x1, cmulti *y0, cmulti *y1)
{
  int i;
  for(i=0; i<n; i++){ icsub_dc(z0[i],z1[i],x0[i],x1[i],y0,y1); }
}

/**
 @brief [z0,z1]=[x0,x1]-[y0,y1]
 */
void irvec_sub_zscalar(int n, cmulti **z0, cmulti **z1, rmulti **x0, rmulti **x1, dcomplex y0, dcomplex y1)
{
  int i;
  for(i=0; i<n; i++){ icsub_rz(z0[i],z1[i],x0[i],x1[i],y0,y1); }
}

/**
 @brief [z0,z1]=[x0,x1]-[y0,y1]
 */
void izvec_sub_rscalar(int n, cmulti **z0, cmulti **z1, dcomplex *x0, dcomplex *x1, rmulti *y0, rmulti *y1)
{
  int i;
  for(i=0; i<n; i++){ icsub_zr(z0[i],z1[i],x0[i],x1[i],y0,y1); }
}

/**
 @brief [z0,z1]=[x0,x1]-[y0,y1]
 */
void icscalar_sub_cvec(int n, cmulti **z0, cmulti **z1, cmulti *x0, cmulti *x1, cmulti **y0, cmulti **y1)
{
  int i;
  for(i=0; i<n; i++){ icsub_cc(z0[i],z1[i],x0,x1,y0[i],y1[i]); }
}

/**
 @brief [z0,z1]=[x0,x1]-[y0,y1]
 */
void icscalar_sub_rvec(int n, cmulti **z0, cmulti **z1, cmulti *x0, cmulti *x1, rmulti **y0, rmulti **y1)
{
  int i;
  for(i=0; i<n; i++){ icsub_cr(z0[i],z1[i],x0,x1,y0[i],y1[i]); }
}

/**
 @brief [z0,z1]=[x0,x1]-[y0,y1]
 */
void irscalar_sub_cvec(int n, cmulti **z0, cmulti **z1, rmulti *x0, rmulti *x1, cmulti **y0, cmulti **y1)
{
  int i;
  for(i=0; i<n; i++){ icsub_rc(z0[i],z1[i],x0,x1,y0[i],y1[i]); }
}

/**
 @brief [z0,z1]=[x0,x1]-[y0,y1]
 */
void icscalar_sub_zvec(int n, cmulti **z0, cmulti **z1, cmulti *x0, cmulti *x1, dcomplex *y0, dcomplex *y1)
{
  int i;
  for(i=0; i<n; i++){ icsub_cz(z0[i],z1[i],x0,x1,y0[i],y1[i]); }
}

/**
 @brief [z0,z1]=[x0,x1]-[y0,y1]
 */
void izscalar_sub_cvec(int n, cmulti **z0, cmulti **z1, dcomplex x0, dcomplex x1, cmulti **y0, cmulti **y1)
{
  int i;
  for(i=0; i<n; i++){ icsub_zc(z0[i],z1[i],x0,x1,y0[i],y1[i]); }
}

/**
 @brief [z0,z1]=[x0,x1]-[y0,y1]
 */
void icscalar_sub_dvec(int n, cmulti **z0, cmulti **z1, cmulti *x0, cmulti *x1, double *y0, double *y1)
{
  int i;
  for(i=0; i<n; i++){ icsub_cd(z0[i],z1[i],x0,x1,y0[i],y1[i]); }
}

/**
 @brief [z0,z1]=[x0,x1]-[y0,y1]
 */
void idscalar_sub_cvec(int n, cmulti **z0, cmulti **z1, double x0, double x1, cmulti **y0, cmulti **y1)
{
  int i;
  for(i=0; i<n; i++){ icsub_dc(z0[i],z1[i],x0,x1,y0[i],y1[i]); }
}

/**
 @brief [z0,z1]=[x0,x1]-[y0,y1]
 */
void irscalar_sub_zvec(int n, cmulti **z0, cmulti **z1, rmulti *x0, rmulti *x1, dcomplex *y0, dcomplex *y1)
{
  int i;
  for(i=0; i<n; i++){ icsub_rz(z0[i],z1[i],x0,x1,y0[i],y1[i]); }
}

/**
 @brief [z0,z1]=[x0,x1]-[y0,y1]
 */
void izscalar_sub_rvec(int n, cmulti **z0, cmulti **z1, dcomplex x0, dcomplex x1, rmulti **y0, rmulti **y1)
{
  int i;
  for(i=0; i<n; i++){ icsub_zr(z0[i],z1[i],x0,x1,y0[i],y1[i]); }
}

///////////////////////////////////////////////////

/**
 @brief [z0,z1]=[x0,x1]*[y0,y1]
 */
void icvec_mul_cvec(int n, cmulti **z0, cmulti **z1, cmulti **x0, cmulti **x1, cmulti **y0, cmulti **y1)
{
  int i;
  for(i=0; i<n; i++){ icmul_cc(z0[i],z1[i],x0[i],x1[i],y0[i],y1[i]); }
}

/**
 @brief [z0,z1]=[x0,x1]*[y0,y1]
 */
void icvec_mul_rvec(int n, cmulti **z0, cmulti **z1, cmulti **x0, cmulti **x1, rmulti **y0, rmulti **y1)
{
  int i;
  for(i=0; i<n; i++){ icmul_cr(z0[i],z1[i],x0[i],x1[i],y0[i],y1[i]); }
}

/**
 @brief [z0,z1]=[x0,x1]*[y0,y1]
 */
void irvec_mul_cvec(int n, cmulti **z0, cmulti **z1, rmulti **x0, rmulti **x1, cmulti **y0, cmulti **y1)
{
  int i;
  for(i=0; i<n; i++){ icmul_rc(z0[i],z1[i],x0[i],x1[i],y0[i],y1[i]); }
}

/**
 @brief [z0,z1]=[x0,x1]*[y0,y1]
 */
void icvec_mul_zvec(int n, cmulti **z0, cmulti **z1, cmulti **x0, cmulti **x1, dcomplex *y0, dcomplex *y1)
{
  int i;
  for(i=0; i<n; i++){ icmul_cz(z0[i],z1[i],x0[i],x1[i],y0[i],y1[i]); }
}

/**
 @brief [z0,z1]=[x0,x1]*[y0,y1]
 */
void izvec_mul_cvec(int n, cmulti **z0, cmulti **z1, dcomplex *x0, dcomplex *x1, cmulti **y0, cmulti **y1)
{
  int i;
  for(i=0; i<n; i++){ icmul_zc(z0[i],z1[i],x0[i],x1[i],y0[i],y1[i]); }
}

/**
 @brief [z0,z1]=[x0,x1]*[y0,y1]
 */
void icvec_mul_dvec(int n, cmulti **z0, cmulti **z1, cmulti **x0, cmulti **x1, double *y0, double *y1)
{
  int i;
  for(i=0; i<n; i++){ icmul_cd(z0[i],z1[i],x0[i],x1[i],y0[i],y1[i]); }
}

/**
 @brief [z0,z1]=[x0,x1]*[y0,y1]
 */
void idvec_mul_cvec(int n, cmulti **z0, cmulti **z1, double *x0, double *x1, cmulti **y0, cmulti **y1)
{
  int i;
  for(i=0; i<n; i++){ icmul_dc(z0[i],z1[i],x0[i],x1[i],y0[i],y1[i]); }
}

/**
 @brief [z0,z1]=[x0,x1]*[y0,y1]
 */
void irvec_mul_zvec(int n, cmulti **z0, cmulti **z1, rmulti **x0, rmulti **x1, dcomplex *y0, dcomplex *y1)
{
  int i;
  for(i=0; i<n; i++){ icmul_rz(z0[i],z1[i],x0[i],x1[i],y0[i],y1[i]); }
}

/**
 @brief [z0,z1]=[x0,x1]*[y0,y1]
 */
void izvec_mul_rvec(int n, cmulti **z0, cmulti **z1, dcomplex *x0, dcomplex *x1, rmulti **y0, rmulti **y1)
{
  int i;
  for(i=0; i<n; i++){ icmul_zr(z0[i],z1[i],x0[i],x1[i],y0[i],y1[i]); }
}

/**
 @brief [z0,z1]=[x0,x1]*[y0,y1]
 */
void icvec_mul_cscalar(int n, cmulti **z0, cmulti **z1, cmulti **x0, cmulti **x1, cmulti *y0, cmulti *y1)
{
  int i;
  for(i=0; i<n; i++){ icmul_cc(z0[i],z1[i],x0[i],x1[i],y0,y1); }
}

/**
 @brief [z0,z1]=[x0,x1]*[y0,y1]
 */
void icvec_mul_rscalar(int n, cmulti **z0, cmulti **z1, cmulti **x0, cmulti **x1, rmulti *y0, rmulti *y1)
{
  int i;
  for(i=0; i<n; i++){ icmul_cr(z0[i],z1[i],x0[i],x1[i],y0,y1); }
}

/**
 @brief [z0,z1]=[x0,x1]*[y0,y1]
 */
void irvec_mul_cscalar(int n, cmulti **z0, cmulti **z1, rmulti **x0, rmulti **x1, cmulti *y0, cmulti *y1)
{
  int i;
  for(i=0; i<n; i++){ icmul_rc(z0[i],z1[i],x0[i],x1[i],y0,y1); }
}

/**
 @brief [z0,z1]=[x0,x1]*[y0,y1]
 */
void icvec_mul_zscalar(int n, cmulti **z0, cmulti **z1, cmulti **x0, cmulti **x1, dcomplex y0, dcomplex y1)
{
  int i;
  for(i=0; i<n; i++){ icmul_cz(z0[i],z1[i],x0[i],x1[i],y0,y1); }
}

/**
 @brief [z0,z1]=[x0,x1]*[y0,y1]
 */
void izvec_mul_cscalar(int n, cmulti **z0, cmulti **z1, dcomplex *x0, dcomplex *x1, cmulti *y0, cmulti *y1)
{
  int i;
  for(i=0; i<n; i++){ icmul_zc(z0[i],z1[i],x0[i],x1[i],y0,y1); }
}

/**
 @brief [z0,z1]=[x0,x1]*[y0,y1]
 */
void icvec_mul_dscalar(int n, cmulti **z0, cmulti **z1, cmulti **x0, cmulti **x1, double y0, double y1)
{
  int i;
  for(i=0; i<n; i++){ icmul_cd(z0[i],z1[i],x0[i],x1[i],y0,y1); }
}

/**
 @brief [z0,z1]=[x0,x1]*[y0,y1]
 */
void idvec_mul_cscalar(int n, cmulti **z0, cmulti **z1, double *x0, double *x1, cmulti *y0, cmulti *y1)
{
  int i;
  for(i=0; i<n; i++){ icmul_dc(z0[i],z1[i],x0[i],x1[i],y0,y1); }
}

/**
 @brief [z0,z1]=[x0,x1]*[y0,y1]
 */
void irvec_mul_zscalar(int n, cmulti **z0, cmulti **z1, rmulti **x0, rmulti **x1, dcomplex y0, dcomplex y1)
{
  int i;
  for(i=0; i<n; i++){ icmul_rz(z0[i],z1[i],x0[i],x1[i],y0,y1); }
}

/**
 @brief [z0,z1]=[x0,x1]*[y0,y1]
 */
void izvec_mul_rscalar(int n, cmulti **z0, cmulti **z1, dcomplex *x0, dcomplex *x1, rmulti *y0, rmulti *y1)
{
  int i;
  for(i=0; i<n; i++){ icmul_zr(z0[i],z1[i],x0[i],x1[i],y0,y1); }
}

/**
 @brief [z0,z1]=[x0,x1]*[y0,y1]
 */
void icscalar_mul_cvec(int n, cmulti **z0, cmulti **z1, cmulti *x0, cmulti *x1, cmulti **y0, cmulti **y1)
{
  int i;
  for(i=0; i<n; i++){ icmul_cc(z0[i],z1[i],x0,x1,y0[i],y1[i]); }
}

/**
 @brief [z0,z1]=[x0,x1]*[y0,y1]
 */
void icscalar_mul_rvec(int n, cmulti **z0, cmulti **z1, cmulti *x0, cmulti *x1, rmulti **y0, rmulti **y1)
{
  int i;
  for(i=0; i<n; i++){ icmul_cr(z0[i],z1[i],x0,x1,y0[i],y1[i]); }
}

/**
 @brief [z0,z1]=[x0,x1]*[y0,y1]
 */
void irscalar_mul_cvec(int n, cmulti **z0, cmulti **z1, rmulti *x0, rmulti *x1, cmulti **y0, cmulti **y1)
{
  int i;
  for(i=0; i<n; i++){ icmul_rc(z0[i],z1[i],x0,x1,y0[i],y1[i]); }
}

/**
 @brief [z0,z1]=[x0,x1]*[y0,y1]
 */
void icscalar_mul_zvec(int n, cmulti **z0, cmulti **z1, cmulti *x0, cmulti *x1, dcomplex *y0, dcomplex *y1)
{
  int i;
  for(i=0; i<n; i++){ icmul_cz(z0[i],z1[i],x0,x1,y0[i],y1[i]); }
}

/**
 @brief [z0,z1]=[x0,x1]*[y0,y1]
 */
void izscalar_mul_cvec(int n, cmulti **z0, cmulti **z1, dcomplex x0, dcomplex x1, cmulti **y0, cmulti **y1)
{
  int i;
  for(i=0; i<n; i++){ icmul_zc(z0[i],z1[i],x0,x1,y0[i],y1[i]); }
}

/**
 @brief [z0,z1]=[x0,x1]*[y0,y1]
 */
void icscalar_mul_dvec(int n, cmulti **z0, cmulti **z1, cmulti *x0, cmulti *x1, double *y0, double *y1)
{
  int i;
  for(i=0; i<n; i++){ icmul_cd(z0[i],z1[i],x0,x1,y0[i],y1[i]); }
}

/**
 @brief [z0,z1]=[x0,x1]*[y0,y1]
 */
void idscalar_mul_cvec(int n, cmulti **z0, cmulti **z1, double x0, double x1, cmulti **y0, cmulti **y1)
{
  int i;
  for(i=0; i<n; i++){ icmul_dc(z0[i],z1[i],x0,x1,y0[i],y1[i]); }
}

/**
 @brief [z0,z1]=[x0,x1]*[y0,y1]
 */
void irscalar_mul_zvec(int n, cmulti **z0, cmulti **z1, rmulti *x0, rmulti *x1, dcomplex *y0, dcomplex *y1)
{
  int i;
  for(i=0; i<n; i++){ icmul_rz(z0[i],z1[i],x0,x1,y0[i],y1[i]); }
}

/**
 @brief [z0,z1]=[x0,x1]*[y0,y1]
 */
void izscalar_mul_rvec(int n, cmulti **z0, cmulti **z1, dcomplex x0, dcomplex x1, rmulti **y0, rmulti **y1)
{
  int i;
  for(i=0; i<n; i++){ icmul_zr(z0[i],z1[i],x0,x1,y0[i],y1[i]); }
}

///////////////////////////////////////////////////


/**
 @brief [z0,z1]=[x0,x1]/[y0,y1]
 */
void icvec_div_cvec(int n, cmulti **z0, cmulti **z1, cmulti **x0, cmulti **x1, cmulti **y0, cmulti **y1)
{
  int i;
  for(i=0; i<n; i++){ icdiv_cc(z0[i],z1[i],x0[i],x1[i],y0[i],y1[i]); }
}

/**
 @brief [z0,z1]=[x0,x1]/[y0,y1]
 */
void icvec_div_rvec(int n, cmulti **z0, cmulti **z1, cmulti **x0, cmulti **x1, rmulti **y0, rmulti **y1)
{
  int i;
  for(i=0; i<n; i++){ icdiv_cr(z0[i],z1[i],x0[i],x1[i],y0[i],y1[i]); }
}

/**
 @brief [z0,z1]=[x0,x1]/[y0,y1]
 */
void irvec_div_cvec(int n, cmulti **z0, cmulti **z1, rmulti **x0, rmulti **x1, cmulti **y0, cmulti **y1)
{
  int i;
  for(i=0; i<n; i++){ icdiv_rc(z0[i],z1[i],x0[i],x1[i],y0[i],y1[i]); }
}

/**
 @brief [z0,z1]=[x0,x1]/[y0,y1]
 */
void icvec_div_zvec(int n, cmulti **z0, cmulti **z1, cmulti **x0, cmulti **x1, dcomplex *y0, dcomplex *y1)
{
  int i;
  for(i=0; i<n; i++){ icdiv_cz(z0[i],z1[i],x0[i],x1[i],y0[i],y1[i]); }
}

/**
 @brief [z0,z1]=[x0,x1]/[y0,y1]
 */
void izvec_div_cvec(int n, cmulti **z0, cmulti **z1, dcomplex *x0, dcomplex *x1, cmulti **y0, cmulti **y1)
{
  int i;
  for(i=0; i<n; i++){ icdiv_zc(z0[i],z1[i],x0[i],x1[i],y0[i],y1[i]); }
}

/**
 @brief [z0,z1]=[x0,x1]/[y0,y1]
 */
void icvec_div_dvec(int n, cmulti **z0, cmulti **z1, cmulti **x0, cmulti **x1, double *y0, double *y1)
{
  int i;
  for(i=0; i<n; i++){ icdiv_cd(z0[i],z1[i],x0[i],x1[i],y0[i],y1[i]); }
}

/**
 @brief [z0,z1]=[x0,x1]/[y0,y1]
 */
void idvec_div_cvec(int n, cmulti **z0, cmulti **z1, double *x0, double *x1, cmulti **y0, cmulti **y1)
{
  int i;
  for(i=0; i<n; i++){ icdiv_dc(z0[i],z1[i],x0[i],x1[i],y0[i],y1[i]); }
}

/**
 @brief [z0,z1]=[x0,x1]/[y0,y1]
 */
void irvec_div_zvec(int n, cmulti **z0, cmulti **z1, rmulti **x0, rmulti **x1, dcomplex *y0, dcomplex *y1)
{
  int i;
  for(i=0; i<n; i++){ icdiv_rz(z0[i],z1[i],x0[i],x1[i],y0[i],y1[i]); }
}

/**
 @brief [z0,z1]=[x0,x1]/[y0,y1]
 */
void izvec_div_rvec(int n, cmulti **z0, cmulti **z1, dcomplex *x0, dcomplex *x1, rmulti **y0, rmulti **y1)
{
  int i;
  for(i=0; i<n; i++){ icdiv_zr(z0[i],z1[i],x0[i],x1[i],y0[i],y1[i]); }
}

/**
 @brief [z0,z1]=[x0,x1]/[y0,y1]
 */
void icvec_div_cscalar(int n, cmulti **z0, cmulti **z1, cmulti **x0, cmulti **x1, cmulti *y0, cmulti *y1)
{
  int i;
  for(i=0; i<n; i++){ icdiv_cc(z0[i],z1[i],x0[i],x1[i],y0,y1); }
}

/**
 @brief [z0,z1]=[x0,x1]/[y0,y1]
 */
void icvec_div_rscalar(int n, cmulti **z0, cmulti **z1, cmulti **x0, cmulti **x1, rmulti *y0, rmulti *y1)
{
  int i;
  for(i=0; i<n; i++){ icdiv_cr(z0[i],z1[i],x0[i],x1[i],y0,y1); }
}

/**
 @brief [z0,z1]=[x0,x1]/[y0,y1]
 */
void irvec_div_cscalar(int n, cmulti **z0, cmulti **z1, rmulti **x0, rmulti **x1, cmulti *y0, cmulti *y1)
{
  int i;
  for(i=0; i<n; i++){ icdiv_rc(z0[i],z1[i],x0[i],x1[i],y0,y1); }
}

/**
 @brief [z0,z1]=[x0,x1]/[y0,y1]
 */
void icvec_div_zscalar(int n, cmulti **z0, cmulti **z1, cmulti **x0, cmulti **x1, dcomplex y0, dcomplex y1)
{
  int i;
  for(i=0; i<n; i++){ icdiv_cz(z0[i],z1[i],x0[i],x1[i],y0,y1); }
}

/**
 @brief [z0,z1]=[x0,x1]/[y0,y1]
 */
void izvec_div_cscalar(int n, cmulti **z0, cmulti **z1, dcomplex *x0, dcomplex *x1, cmulti *y0, cmulti *y1)
{
  int i;
  for(i=0; i<n; i++){ icdiv_zc(z0[i],z1[i],x0[i],x1[i],y0,y1); }
}

/**
 @brief [z0,z1]=[x0,x1]/[y0,y1]
 */
void icvec_div_dscalar(int n, cmulti **z0, cmulti **z1, cmulti **x0, cmulti **x1, double y0, double y1)
{
  int i;
  for(i=0; i<n; i++){ icdiv_cd(z0[i],z1[i],x0[i],x1[i],y0,y1); }
}

/**
 @brief [z0,z1]=[x0,x1]/[y0,y1]
 */
void idvec_div_cscalar(int n, cmulti **z0, cmulti **z1, double *x0, double *x1, cmulti *y0, cmulti *y1)
{
  int i;
  for(i=0; i<n; i++){ icdiv_dc(z0[i],z1[i],x0[i],x1[i],y0,y1); }
}

/**
 @brief [z0,z1]=[x0,x1]/[y0,y1]
 */
void irvec_div_zscalar(int n, cmulti **z0, cmulti **z1, rmulti **x0, rmulti **x1, dcomplex y0, dcomplex y1)
{
  int i;
  for(i=0; i<n; i++){ icdiv_rz(z0[i],z1[i],x0[i],x1[i],y0,y1); }
}

/**
 @brief [z0,z1]=[x0,x1]/[y0,y1]
 */
void izvec_div_rscalar(int n, cmulti **z0, cmulti **z1, dcomplex *x0, dcomplex *x1, rmulti *y0, rmulti *y1)
{
  int i;
  for(i=0; i<n; i++){ icdiv_zr(z0[i],z1[i],x0[i],x1[i],y0,y1); }
}

/**
 @brief [z0,z1]=[x0,x1]/[y0,y1]
 */
void icscalar_div_cvec(int n, cmulti **z0, cmulti **z1, cmulti *x0, cmulti *x1, cmulti **y0, cmulti **y1)
{
  int i;
  for(i=0; i<n; i++){ icdiv_cc(z0[i],z1[i],x0,x1,y0[i],y1[i]); }
}

/**
 @brief [z0,z1]=[x0,x1]/[y0,y1]
 */
void icscalar_div_rvec(int n, cmulti **z0, cmulti **z1, cmulti *x0, cmulti *x1, rmulti **y0, rmulti **y1)
{
  int i;
  for(i=0; i<n; i++){ icdiv_cr(z0[i],z1[i],x0,x1,y0[i],y1[i]); }
}

/**
 @brief [z0,z1]=[x0,x1]/[y0,y1]
 */
void irscalar_div_cvec(int n, cmulti **z0, cmulti **z1, rmulti *x0, rmulti *x1, cmulti **y0, cmulti **y1)
{
  int i;
  for(i=0; i<n; i++){ icdiv_rc(z0[i],z1[i],x0,x1,y0[i],y1[i]); }
}

/**
 @brief [z0,z1]=[x0,x1]/[y0,y1]
 */
void icscalar_div_zvec(int n, cmulti **z0, cmulti **z1, cmulti *x0, cmulti *x1, dcomplex *y0, dcomplex *y1)
{
  int i;
  for(i=0; i<n; i++){ icdiv_cz(z0[i],z1[i],x0,x1,y0[i],y1[i]); }
}

/**
 @brief [z0,z1]=[x0,x1]/[y0,y1]
 */
void izscalar_div_cvec(int n, cmulti **z0, cmulti **z1, dcomplex x0, dcomplex x1, cmulti **y0, cmulti **y1)
{
  int i;
  for(i=0; i<n; i++){ icdiv_zc(z0[i],z1[i],x0,x1,y0[i],y1[i]); }
}

/**
 @brief [z0,z1]=[x0,x1]/[y0,y1]
 */
void icscalar_div_dvec(int n, cmulti **z0, cmulti **z1, cmulti *x0, cmulti *x1, double *y0, double *y1)
{
  int i;
  for(i=0; i<n; i++){ icdiv_cd(z0[i],z1[i],x0,x1,y0[i],y1[i]); }
}

/**
 @brief [z0,z1]=[x0,x1]/[y0,y1]
 */
void idscalar_div_cvec(int n, cmulti **z0, cmulti **z1, double x0, double x1, cmulti **y0, cmulti **y1)
{
  int i;
  for(i=0; i<n; i++){ icdiv_dc(z0[i],z1[i],x0,x1,y0[i],y1[i]); }
}

/**
 @brief [z0,z1]=[x0,x1]/[y0,y1]
 */
void irscalar_div_zvec(int n, cmulti **z0, cmulti **z1, rmulti *x0, rmulti *x1, dcomplex *y0, dcomplex *y1)
{
  int i;
  for(i=0; i<n; i++){ icdiv_rz(z0[i],z1[i],x0,x1,y0[i],y1[i]); }
}

/**
 @brief [z0,z1]=[x0,x1]/[y0,y1]
 */
void izscalar_div_rvec(int n, cmulti **z0, cmulti **z1, dcomplex x0, dcomplex x1, rmulti **y0, rmulti **y1)
{
  int i;
  for(i=0; i<n; i++){ icdiv_zr(z0[i],z1[i],x0,x1,y0[i],y1[i]); }
}

///////////////////////////////////////////////////

/**
 @brief 絶対値の平方の総和 [y0,y1]=sum(abs([x0,x1]).^2)
*/
void icvec_sum_abs2(rmulti *y0, rmulti *y1, int n, cmulti **x0, cmulti **x1)
{
  int i;
  irset_d(y0,y1,0,0);
  for(i=0; i<n; i++){ iradd_abs2_c(y0,y1,x0[i],x1[i]); }
}

/**
 @brief 絶対値 [y0,y1]=abs([x0,x1])
*/
void icvec_abs(int n, rmulti **y0, rmulti **y1, cmulti **x0, cmulti **x1)
{
  int i;
  for(i=0; i<n; i++){ irabs_c(y0[i],y1[i],x0[i],x1[i]); }
}

/**
 @brief 差の絶対値 [z0,z1]=abs([x0,x1]-[y0,y1])
*/
void icvec_abs_sub(int n, rmulti **z0, rmulti **z1, cmulti **x0, cmulti **x1, cmulti **y0, cmulti **y1)
{
  int i;
  for(i=0; i<n; i++){ irabs_sub_cc(z0[i],z1[i],x0[i],x1[i],y0[i],y1[i]); }
}

/**
 @brief 最大値 [y0,y1]=[max(x0),max(x1)]
*/
void icvec_max(cmulti *y0, cmulti *y1, int n, cmulti **x0, cmulti **x1)
{
  int i;
  icset_c(y0,y1,x0[0],x1[0]);      // y0=x0[0], y1=x1[0]
  for(i=1; i<n; i++){
    if(cgt_c(x0[i],y0)){            // x0[i]>y0
      cset_c(y0,x0[i]);            // y0=x0[i]
    }
    if(cgt_c(x1[i],y1)){            // x1[i]>y1
      cset_c(y1,x1[i]);            // y1=x1[i]
    }
  }
}

/**
 @brief 最大値 [y0,y1]=[x0,max(x1)]
*/
void icvec_umax(cmulti *y0, cmulti *y1, int n, cmulti **x0, cmulti **x1)
{
  int i;
  icset_c(y0,y1,x0[0],x1[0]);      // y0=x0[0], y1=x1[0]
  for(i=1; i<n; i++){
    if(cgt_c(x1[i],y1)){            // x1[i]>y1
      icset_c(y0,y1,x0[i],x1[i]);  // y0=x0[i], y1=x1[i]
    }
  }
}

/**
 @brief 最大値 [y0,y1]=[max(x0),x1]
*/
void icvec_dmax(cmulti *y0, cmulti *y1, int n, cmulti **x0, cmulti **x1)
{
  int i;
  icset_c(y0,y1,x0[0],x1[0]);      // y0=x0[0], y1=x1[0]
  for(i=1; i<n; i++){
    if(cgt_c(x0[i],y0)){            // x0[i]>y0
      icset_c(y0,y1,x0[i],x1[i]);  // y0=x0[i], y1=x1[i]
    }
  }
}

/**
 @brief 最小値 [y0,y1]=[min(x0),min(x1)]
*/
void icvec_min(cmulti *y0, cmulti *y1, int n, cmulti **x0, cmulti **x1)
{
  int i;
  icset_c(y0,y1,x0[0],x1[0]);   // y0=x0[0], y1=x1[0]
  for(i=1; i<n; i++){
    if(clt_c(x0[i],y0)){         // x0[i]<y1
      cset_c(y0,x0[i]);         // y0=x0[i]
    }
    if(clt_c(x1[i],y1)){         // x1[i]<y1
      cset_c(y1,x1[i]);         // y1=x1[i]
    }
  }
}

//ここまで

/**
 @brief 最小値 [y0,y1]=[x0,min(x1)]
*/
void icvec_umin(cmulti *y0, cmulti *y1, int n, cmulti **x0, cmulti **x1)
{
  int i;
  icset_c(y0,y1,x0[0],x1[0]);      // y0=x0[0], y1=x1[0]
  for(i=1; i<n; i++){
    if(clt_c(x1[i],y1)){            // x1[i]<y1
      icset_c(y0,y1,x0[i],x1[i]);  // y0=x0[i], y1=x1[i]
    }
  }
}

/**
 @brief 最小値 [y0,y1]=[min(x0),x1]
*/
void icvec_dmin(cmulti *y0, cmulti *y1, int n, cmulti **x0, cmulti **x1)
{
  int i;
  icset_c(y0,y1,x0[0],x1[0]);      // y0=x0[0], y1=x1[0]
  for(i=1; i<n; i++){
    if(clt_c(x0[i],y0)){            // x0[i]<y0
      icset_c(y0,y1,x0[i],x1[i]);  // y0=x0[i], y1=x1[i]
    }
  }
}

/**
 @brief irmulti型のベクトルの要素の絶対値の最大値 value=max(abs(x))
*/
void icvec_umax_abs(rmulti *y0, rmulti *y1, int n, cmulti **x0, cmulti **x1)
{
  int prec;
  rmulti **ax0=NULL,**ax1=NULL;
  prec=MAX2(rget_prec(y0),rget_prec(y1));
  ax0=rvec_allocate_prec(n,prec);
  ax1=rvec_allocate_prec(n,prec);
  icvec_abs(n,ax0,ax1,x0,x1);
  irvec_umax(y0,y1,n,ax0,ax1);
  RVF(ax0,n); RVF(ax1,n);
}

/**
 @brief irmulti型のベクトルの要素の絶対値の最小値 value=min(abs(x))
*/
void icvec_dmin_abs(rmulti *y0, rmulti *y1, int n, cmulti **x0, cmulti **x1)
{
  int prec;
  rmulti **ax0=NULL,**ax1=NULL;
  prec=MAX2(rget_prec(y0),rget_prec(y1));
  ax0=rvec_allocate_prec(n,prec);
  ax1=rvec_allocate_prec(n,prec);
  icvec_abs(n,ax0,ax1,x0,x1);
  irvec_dmin(y0,y1,n,ax0,ax1);
  RVF(ax0,n); RVF(ax1,n);
}

/**
 @brief icmulti型のベクトルの要素の総和 value=sum(x)
 */
void icvec_sum(cmulti *y0, cmulti *y1, int n, cmulti **x0, cmulti **x1)
{
  int i;
  icset_d(y0,y1,0,0);
  for(i=0; i<n; i++){ icadd_cc(y0,y1,y0,y1,x0[i],x1[i]); }
}

/**
 @brief 線形変換 [y0,y1]=[A0,A1]*[x0,x1]
*/
void icvec_lintr(int m, int n, cmulti **y0, cmulti **y1, cmulti **A0, int LDA0, cmulti **A1, int LDA1, cmulti **x0, cmulti **x1)
{
  int prec,i,j;
  cmulti **z0=NULL,**z1=NULL;
  prec=MAX2(cvec_get_prec_max(m,y0),cvec_get_prec_max(m,y1));
  CVA(z0,m,prec); CVA(z1,m,prec);
  for(i=0; i<m; i++){
    icset_d(z0[i],z1[i],0,0);
    for(j=0; j<n; j++){
      icadd_mul_rr(z0[i],z1[i],MAT(A0,i,j,LDA0),MAT(A1,i,j,LDA1),x0[j],x1[j]);
    }
  }
  icvec_copy(m,y0,y1,z0,z1);
  CVF(z0,m); CVF(z1,m);
}

/**
 @brief 線形変換の加算 [y0,y1]+=[A0,A1]*[x0,x1]
*/
void icvec_add_lintr(int m, int n, cmulti **y0, cmulti **y1, cmulti **A0, int LDA0, cmulti **A1, int LDA1, cmulti **x0, cmulti **x1)
{
  int prec,i,j;
  cmulti **z0=NULL,**z1=NULL;
  prec=MAX2(cvec_get_prec_max(m,y0),cvec_get_prec_max(m,y1));
  CVA(z0,m,prec); CVA(z1,m,prec);
  for(i=0; i<m; i++){
    icset_c(z0[i],z1[i],y0[i],y1[i]);
    for(j=0; j<n; j++){
      icadd_mul_rr(z0[i],z1[i],MAT(A0,i,j,LDA0),MAT(A1,i,j,LDA1),x0[j],x1[j]);
    }
  }
  icvec_copy(m,y0,y1,z0,z1);
  CVF(z0,m); CVF(z1,m);
}

/**
 @brief 線形変換の減算 [y0,y1]-=[A0,A1]*[x0,x1]
*/
void icvec_sub_lintr(int m, int n, cmulti **y0, cmulti **y1, cmulti **A0, int LDA0, cmulti **A1, int LDA1, cmulti **x0, cmulti **x1)
{
  int prec,i,j;
  cmulti **z0=NULL,**z1=NULL;
  prec=MAX2(cvec_get_prec_max(m,y0),cvec_get_prec_max(m,y1));
  CVA(z0,m,prec); CVA(z1,m,prec);
  for(i=0; i<m; i++){
    icset_c(z0[i],z1[i],y0[i],y1[i]);
    for(j=0; j<n; j++){
      icsub_mul_cc(z0[i],z1[i],MAT(A0,i,j,LDA0),MAT(A1,i,j,LDA1),x0[j],x1[j]);
    }
  }
  icvec_copy(m,y0,y1,z0,z1);
  CVF(z0,m); CVF(z1,m);
}

/**
 @brief 転置行列の線形変換 [y0,y1]=[A0,A1]'*[x0,x1]
*/
void icvec_lintr_t(int m, int n, cmulti **y0, cmulti **y1, cmulti **A0, int LDA0, cmulti **A1, int LDA1, cmulti **x0, cmulti **x1)
{
  int prec,i,j;
  cmulti **z0=NULL,**z1=NULL;
  prec=MAX2(cvec_get_prec_max(n,y0),cvec_get_prec_max(n,y1));
  CVA(z0,n,prec); CVA(z1,n,prec);
  for(j=0; j<n; j++){
    icset_d(z0[j],z1[j],0,0);
    for(i=0; i<m; i++){
      icadd_mul_rr(z0[j],z1[j],MAT(A0,i,j,LDA0),MAT(A1,i,j,LDA1),x0[i],x1[i]);
    }
  }
  icvec_copy(n,y0,y1,z0,z1);
  CVF(z0,n); CVF(z1,n);
}

/**
 @brief 転置行列の線形変換の加算 [y0,y1]+=[A0,A1]'*[x0,x1]
*/
void icvec_add_lintr_t(int m, int n, cmulti **y0, cmulti **y1, cmulti **A0, int LDA0, cmulti **A1, int LDA1, cmulti **x0, cmulti **x1)
{
  int prec,i,j;
  cmulti **z0=NULL,**z1=NULL;
  prec=MAX2(cvec_get_prec_max(n,y0),cvec_get_prec_max(n,y1));
  CVA(z0,n,prec); CVA(z1,n,prec);
  for(j=0; j<n; j++){
    icset_c(z0[j],z1[j],y0[j],y1[j]);
    for(i=0; i<m; i++){
      icadd_mul_rr(z0[j],z1[j],MAT(A0,i,j,LDA0),MAT(A1,i,j,LDA1),x0[i],x1[i]);
    }
  }
  icvec_copy(n,y0,y1,z0,z1);
  CVF(z0,n); CVF(z1,n);
}

/**
 @brief 転置行列の線形変換の減算 [y0,y1]-=[A0,A1]'*[x0,x1]
*/
void icvec_sub_lintr_t(int m, int n, cmulti **y0, cmulti **y1, cmulti **A0, int LDA0, cmulti **A1, int LDA1, cmulti **x0, cmulti **x1)
{
  int prec,i,j;
  cmulti **z0=NULL,**z1=NULL;
  prec=MAX2(cvec_get_prec_max(n,y0),cvec_get_prec_max(n,y1));
  CVA(z0,n,prec); CVA(z1,n,prec);
  for(j=0; j<n; j++){
    icset_c(z0[j],z1[j],y0[j],y1[j]);
    for(i=0; i<m; i++){
      icsub_mul_cc(z0[j],z1[j],MAT(A0,i,j,LDA0),MAT(A1,i,j,LDA1),x0[i],x1[i]);
    }
  }
  icvec_copy(n,y0,y1,z0,z1);
  CVF(z0,n); CVF(z1,n);
}

/**
 @brief 共役転置行列の線形変換 [y0,y1]=[A0,A1]'*[x0,x1]
*/
void icvec_lintr_ct(int m, int n, cmulti **y0, cmulti **y1, cmulti **A0, int LDA0, cmulti **A1, int LDA1, cmulti **x0, cmulti **x1)
{
  int prec,i,j;
  cmulti **z0=NULL,**z1=NULL;
  prec=MAX2(cvec_get_prec_max(n,y0),cvec_get_prec_max(n,y1));
  CVA(z0,n,prec); CVA(z1,n,prec);
  for(j=0; j<n; j++){
    icset_d(z0[j],z1[j],0,0);
    for(i=0; i<m; i++){
      icadd_dot_cc(z0[j],z1[j],MAT(A0,i,j,LDA0),MAT(A1,i,j,LDA1),x0[i],x1[i]);
    }
  }
  icvec_copy(n,y0,y1,z0,z1);
  CVF(z0,n); CVF(z1,n);
}

/**
 @brief 共役転置行列の線形変換の加算 [y0,y1]+=[A0,A1]'*[x0,x1]
*/
void icvec_add_lintr_ct(int m, int n, cmulti **y0, cmulti **y1, cmulti **A0, int LDA0, cmulti **A1, int LDA1, cmulti **x0, cmulti **x1)
{
  int prec,i,j;
  cmulti **z0=NULL,**z1=NULL;
  prec=MAX2(cvec_get_prec_max(n,y0),cvec_get_prec_max(n,y1));
  CVA(z0,n,prec); CVA(z1,n,prec);
  for(j=0; j<n; j++){
    icset_c(z0[j],z1[j],y0[j],y1[j]);
    for(i=0; i<m; i++){
      icadd_dot_cc(z0[j],z1[j],MAT(A0,i,j,LDA0),MAT(A1,i,j,LDA1),x0[i],x1[i]);
    }
  }
  icvec_copy(n,y0,y1,z0,z1);
  CVF(z0,n); CVF(z1,n);
}

/**
 @brief 共役転置行列の線形変換の減算 [y0,y1]-=[A0,A1]'*[x0,x1]
*/
void icvec_sub_lintr_ct(int m, int n, cmulti **y0, cmulti **y1, cmulti **A0, int LDA0, cmulti **A1, int LDA1, cmulti **x0, cmulti **x1)
{
  int prec,i,j;
  cmulti **z0=NULL,**z1=NULL;
  prec=MAX2(cvec_get_prec_max(n,y0),cvec_get_prec_max(n,y1));
  CVA(z0,n,prec); CVA(z1,n,prec);
  for(j=0; j<n; j++){
    icset_c(z0[j],z1[j],y0[j],y1[j]);
    for(i=0; i<m; i++){
      icsub_dot_cc(z0[j],z1[j],MAT(A0,i,j,LDA0),MAT(A1,i,j,LDA1),x0[i],x1[i]);
    }
  }
  icvec_copy(n,y0,y1,z0,z1);
  CVF(z0,n); CVF(z1,n);
}

/** @} */

//EOF
