#include"is_cmulti.h"
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


/** @name 基本操作 */
/** @{ */

/**
 @brief コピー [y0,y1]=[x0,x1]
*/
int icvec_copy(int n, cmulti **y0, cmulti **y1, cmulti **x0, cmulti **x1)
{
  int i,e=0;
  for(i=0; i<n; i++){ e+=iccopy(y0[i],y1[i],x0[i],x1[i]); }
  return e;
}

/**
 @brief 表示
 */
void icvec_print(int n, cmulti **x0, cmulti **x1, const char *name, const char *f, int digits)
{
  int i,k;
  char format[128];
  if(STR_EQ(f,"f") || STR_EQ(f,"F")){ k=3; }
  else if((strcmp(f,"e")==0) || (strcmp(f,"E")==0)){ k=7; }
  else{ k=3; }
  sprintf(format,"[%%%d.%dR%s, %%%d.%dR%s] [%%%d.%dR%s, %%%d.%dR%s]\n",digits+k,digits,f,digits+k,digits,f,digits+k,digits,f,digits+k,digits,f);
  if(name!=NULL){ printf("%s\n",name); }
  if(x0==NULL || x1==NULL){ printf("NULL\n"); return; }
  for(i=0; i<n; i++){ mpfr_printf(format,C_R(x0[i]),C_R(x1[i]),C_I(x0[i]),C_I(x1[i])); }
}

/**
 @brief 区間の中心 xc=(x1+x0)/2 と半径 xr=x1-x0
*/
int icvec_center_radius(int n, cmulti **xc, cmulti **xr, cmulti **x0, cmulti **x1)
{
  int e=0;
  mpfr_rnd_t mode;
  mode=get_round_mode();
  set_round_mode(MPFR_RNDU);  // up
  e+=cvec_sub(n,xc,x1,x0);    // xc=x1-x0
  e+=cvec_mul_d(n,xc,xc,0.5); // xc=(x1-x0)/2
  e+=cvec_add(n,xc,xc,x0);    // xc=(x1-x0)/2+x0
  e+=cvec_sub(n,xr,xc,x0);    // xr=xc-x0
  set_round_mode(mode);       // back
  return e; 
}


/** @} */

//////////////////////////////////////////////////////

/** @name 四則演算 */
/** @{ */

/**
 @brief 符号の反転 [y0,y1]=-[x0,x1]
*/
int icvec_neg(int n, cmulti **y0, cmulti **y1, cmulti **x0, cmulti **x1)
{
  int i,e=0;
  for(i=0; i<n; i++){ e+=icneg(y0[i],y1[i],x0[i],x1[i]); }
  return e;
}

/**
 @brief 符号の正負 [y0,y1]=[-abs(x),abs(x)]
*/
int icvec_pm(int n, cmulti **y0, cmulti **y1, cmulti **x)
{
  int i,e=0;
  for(i=0; i<n; i++){ e+=icpm(y0[i],y1[i],x[i]); }
  return e;
}

/**
 @brief 区間の拡張 [z0,z1]=[x0,x1]+[-y,y]
*/
int icvec_add_pm(int n, cmulti **z0, cmulti **z1, cmulti **x0, cmulti **x1, cmulti **y)
{
  int i,e=0;
  for(i=0; i<n; i++){ e+=icadd_pm(z0[i],z1[i],x0[i],x1[i],y[i]); }
  return e;
}

/**
 @brief 区間の中心と半径 [m-r,m+r]=[x0,x1]
 */
int icvec_mr(int n, cmulti **mid, cmulti **rad, cmulti **x0, cmulti **x1)
{
  int i,e=0;
  for(i=0; i<n; i++){ e+=icmr(mid[i],rad[i],x0[i],x1[i]); }
  return e;
}

/**
 @brief 足し算 [z0,z1]=[x0,x1]+[y0,y1]
*/
int icvec_add(int n, cmulti **z0, cmulti **z1, cmulti **x0, cmulti **x1, cmulti **y0, cmulti **y1)
{
  int i,e=0;
  for(i=0; i<n; i++){ e+=icadd(z0[i],z1[i],x0[i],x1[i],y0[i],y1[i]); }
  return e;
}

/**
 @brief 足し算 [z0,z1]=[x0,x1]+[y0,y1]
*/
int icvec_add_c(int n, cmulti **z0, cmulti **z1, cmulti **x0, cmulti **x1, cmulti *y0, cmulti *y1)
{
  int i,e=0;
  for(i=0; i<n; i++){ e+=icadd(z0[i],z1[i],x0[i],x1[i],y0,y1); }
  return e;
}

/**
 @brief 引き算 [z0,z1]=[x0,x1]-[y0,y1]
*/
int icvec_sub(int n, cmulti **z0, cmulti **z1, cmulti **x0, cmulti **x1, cmulti **y0, cmulti **y1)
{
  int i,e=0;
  for(i=0; i<n; i++){ e+=icsub(z0[i],z1[i],x0[i],x1[i],y0[i],y1[i]); }
  return e;
}

/**
 @brief 掛け算 [z0,z1]=[x0,x1]*[y0,y1]
 */
int icvec_mul_c(int n, cmulti **z0, cmulti **z1, cmulti **x0, cmulti **x1, cmulti *y0, cmulti *y1)
{
  int i,e=0;
  for(i=0; i<n; i++){ e+=icmul(z0[i],z1[i],x0[i],x1[i],y0,y1); }
  return e;
}

/**
 @brief 絶対値の平方の総和 [y0,y1]=sum(abs([x0,x1]).^2)
*/
int icvec_sum_abs2(rmulti *y0, rmulti *y1, int n, cmulti **x0, cmulti **x1)
{
  int i,e=0;
  e+=irset_d(y0,y1,0);
  for(i=0; i<n; i++){ e+=icadd_abs2(y0,y1,x0[i],x1[i]); }
  return e;
}

/** @} */

//////////////////////////////////////////////////////

/** @name 数学関数 */
/** @{ */

/**
 @brief 絶対値 [y0,y1]=abs([x0,x1])
*/
int icvec_abs(int n, rmulti **y0, rmulti **y1, cmulti **x0, cmulti **x1)
{
  int i,e=0;
  for(i=0; i<n; i++){
    e+=icabs(y0[i],y1[i],x0[i],x1[i]);
  }
  return e;
}

/**
 @brief 差の絶対値 [z0,z1]=abs([x0,x1]-[y0,y1])
*/
int icvec_abs_sub(int n, rmulti **z0, rmulti **z1, cmulti **x0, cmulti **x1, cmulti **y0, cmulti **y1)
{
  int i,e=0;
  for(i=0; i<n; i++){
    e+=icabs_sub(z0[i],z1[i],x0[i],x1[i],y0[i],y1[i]);
  }
  return e;
}

/**
 @brief 最大値 [y0,y1]=[x0,max(x1)]
*/
int icvec_umax(cmulti *y0, cmulti *y1, int n, cmulti **x0, cmulti **x1)
{
  int i,e=0;
  e+=iccopy(y0,y1,x0[0],x1[0]);      // y0=x0[0], y1=x1[0]
  for(i=1; i<n; i++){
    if(cgt(x1[i],y1)){               // x1[i]>y1
      e+=iccopy(y0,y1,x0[i],x1[i]);  // y0=x0[i], y1=x1[i]
    }
  }
  return e;
}

/**
 @brief 最大値 [y0,y1]=[max(x0),x1]
*/
int icvec_dmax(cmulti *y0, cmulti *y1, int n, cmulti **x0, cmulti **x1)
{
  int i,e=0;
  e+=iccopy(y0,y1,x0[0],x1[0]);      // y0=x0[0], y1=x1[0]
  for(i=1; i<n; i++){
    if(cgt(x0[i],y0)){               // x0[i]>y0
      e+=iccopy(y0,y1,x0[i],x1[i]);  // y0=x0[i], y1=x1[i]
    }
  }
  return e;
}

/**
 @brief 最小値 [y0,y1]=[x0,min(x1)]
*/
int icvec_umin(rmulti *y0, rmulti *y1, int n, rmulti **x0, rmulti **x1)
{
  int i,e=0;
  e+=ircopy(y0,y1,x0[0],x1[0]);      // y0=x0[0], y1=x1[0]
  for(i=1; i<n; i++){
    if(rlt(x1[i],y1)){               // x0[i]<y0
      e+=ircopy(y0,y1,x0[i],x1[i]);  // y0=x0[i], y1=x1[i]
    }
  }
  return e;
}

/**
 @brief 最小値 [y0,y1]=[min(x0),x1]
*/
int icvec_dmin(cmulti *y0, cmulti *y1, int n, cmulti **x0, cmulti **x1)
{
  int i,e=0;
  e+=iccopy(y0,y1,x0[0],x1[0]);      // y0=x0[0], y1=x1[0]
  for(i=1; i<n; i++){
    if(clt(x0[i],y0)){               // x0[i]<y0
      e+=iccopy(y0,y1,x0[i],x1[i]);  // y0=x0[i], y1=x1[i]
    }
  }
  return e;
}

/**
 @brief irmulti型のベクトルの要素の絶対値の最大値 value=max(abs(x))
*/
int icvec_umax_abs(rmulti *y0, rmulti *y1, int n, cmulti **x0, cmulti **x1)
{
  int prec,e=0;
  rmulti **ax0=NULL,**ax1=NULL;
  prec=MAX2(rget_prec(y0),rget_prec(y1));
  ax0=rvec_allocate_prec(n,prec);
  ax1=rvec_allocate_prec(n,prec);
  icvec_abs(n,ax0,ax1,x0,x1);
  irvec_umax(y0,y1,n,ax0,ax1);
  RVF(ax0,n); RVF(ax1,n);
  return e;
}

/**
 @brief irmulti型のベクトルの要素の絶対値の最小値 value=min(abs(x))
*/
int icvec_dmin_abs(rmulti *y0, rmulti *y1, int n, cmulti **x0, cmulti **x1)
{
  int prec,e=0;
  rmulti **ax0=NULL,**ax1=NULL;
  prec=MAX2(rget_prec(y0),rget_prec(y1));
  ax0=rvec_allocate_prec(n,prec);
  ax1=rvec_allocate_prec(n,prec);
  icvec_abs(n,ax0,ax1,x0,x1);
  irvec_dmin(y0,y1,n,ax0,ax1);
  RVF(ax0,n); RVF(ax1,n);
  return e;
}

/** @} */

//////////////////////////////////////////////////////////////////////////

/** @name 線形変換 */
/** @{ */

/**
 @brief 線形変換 [y0,y1]=[A0,A1]*[x0,x1]
*/
int icvec_lintr(int m, int n, cmulti **y0, cmulti **y1, cmulti **A0, int LDA0, cmulti **A1, int LDA1, cmulti **x0, cmulti **x1)
{
  int prec,i,j,e=0;
  cmulti **z0=NULL,**z1=NULL;
  prec=MAX2(cvec_get_prec_max(m,y0),cvec_get_prec_max(m,y1));
  CVA(z0,m,prec); CVA(z1,m,prec);
  for(i=0; i<m; i++){
    e+=icset_d(z0[i],z1[i],0);
    for(j=0; j<n; j++){
      e+=icadd_mul(z0[i],z1[i],MAT(A0,i,j,LDA0),MAT(A1,i,j,LDA1),x0[j],x1[j]);
    }
  }
  e+=icvec_copy(m,y0,y1,z0,z1);
  CVF(z0,m); CVF(z1,m);
  return e;
}

/**
 @brief 線形変換の加算 [y0,y1]+=[A0,A1]*[x0,x1]
*/
int icvec_add_lintr(int m, int n, cmulti **y0, cmulti **y1, cmulti **A0, int LDA0, cmulti **A1, int LDA1, cmulti **x0, cmulti **x1)
{
  int prec,i,j,e=0;
  cmulti **z0=NULL,**z1=NULL;
  prec=MAX2(cvec_get_prec_max(m,y0),cvec_get_prec_max(m,y1));
  CVA(z0,m,prec); CVA(z1,m,prec);
  for(i=0; i<m; i++){
    e+=iccopy(z0[i],z1[i],y0[i],y1[i]);
    for(j=0; j<n; j++){
      e+=icadd_mul(z0[i],z1[i],MAT(A0,i,j,LDA0),MAT(A1,i,j,LDA1),x0[j],x1[j]);
    }
  }
  e+=icvec_copy(m,y0,y1,z0,z1);
  CVF(z0,m); CVF(z1,m);
  return e;
}

/**
 @brief 線形変換の減算 [y0,y1]-=[A0,A1]*[x0,x1]
*/
int icvec_sub_lintr(int m, int n, cmulti **y0, cmulti **y1, cmulti **A0, int LDA0, cmulti **A1, int LDA1, cmulti **x0, cmulti **x1)
{
  int prec,i,j,e=0;
  cmulti **z0=NULL,**z1=NULL;
  prec=MAX2(cvec_get_prec_max(m,y0),cvec_get_prec_max(m,y1));
  CVA(z0,m,prec); CVA(z1,m,prec);
  for(i=0; i<m; i++){
    e+=iccopy(z0[i],z1[i],y0[i],y1[i]);
    for(j=0; j<n; j++){
      e+=icsub_mul(z0[i],z1[i],MAT(A0,i,j,LDA0),MAT(A1,i,j,LDA1),x0[j],x1[j]);
    }
  }
  e+=icvec_copy(m,y0,y1,z0,z1);
  CVF(z0,m); CVF(z1,m);
  return e;
}

/**
 @brief 転置行列の線形変換 [y0,y1]=[A0,A1]'*[x0,x1]
*/
int icvec_lintr_t(int m, int n, cmulti **y0, cmulti **y1, cmulti **A0, int LDA0, cmulti **A1, int LDA1, cmulti **x0, cmulti **x1)
{
  int prec,i,j,e=0;
  cmulti **z0=NULL,**z1=NULL;
  prec=MAX2(cvec_get_prec_max(n,y0),cvec_get_prec_max(n,y1));
  CVA(z0,n,prec); CVA(z1,n,prec);
  for(j=0; j<n; j++){
    e+=icset_d(z0[j],z1[j],0);
    for(i=0; i<m; i++){
      e+=icadd_mul(z0[j],z1[j],MAT(A0,i,j,LDA0),MAT(A1,i,j,LDA1),x0[i],x1[i]);
    }
  }
  e+=icvec_copy(n,y0,y1,z0,z1);
  CVF(z0,n); CVF(z1,n);
  return e;
}

/**
 @brief 転置行列の線形変換の加算 [y0,y1]+=[A0,A1]'*[x0,x1]
*/
int icvec_add_lintr_t(int m, int n, cmulti **y0, cmulti **y1, cmulti **A0, int LDA0, cmulti **A1, int LDA1, cmulti **x0, cmulti **x1)
{
  int prec,i,j,e=0;
  cmulti **z0=NULL,**z1=NULL;
  prec=MAX2(cvec_get_prec_max(n,y0),cvec_get_prec_max(n,y1));
  CVA(z0,n,prec); CVA(z1,n,prec);
  for(j=0; j<n; j++){
    e+=iccopy(z0[j],z1[j],y0[j],y1[j]);
    for(i=0; i<m; i++){
      e+=icadd_mul(z0[j],z1[j],MAT(A0,i,j,LDA0),MAT(A1,i,j,LDA1),x0[i],x1[i]);
    }
  }
  e+=icvec_copy(n,y0,y1,z0,z1);
  CVF(z0,n); CVF(z1,n);
  return e;
}

/**
 @brief 転置行列の線形変換の減算 [y0,y1]-=[A0,A1]'*[x0,x1]
*/
int icvec_sub_lintr_t(int m, int n, cmulti **y0, cmulti **y1, cmulti **A0, int LDA0, cmulti **A1, int LDA1, cmulti **x0, cmulti **x1)
{
  int prec,i,j,e=0;
  cmulti **z0=NULL,**z1=NULL;
  prec=MAX2(cvec_get_prec_max(n,y0),cvec_get_prec_max(n,y1));
  CVA(z0,n,prec); CVA(z1,n,prec);
  for(j=0; j<n; j++){
    e+=iccopy(z0[j],z1[j],y0[j],y1[j]);
    for(i=0; i<m; i++){
      e+=icsub_mul(z0[j],z1[j],MAT(A0,i,j,LDA0),MAT(A1,i,j,LDA1),x0[i],x1[i]);
    }
  }
  e+=icvec_copy(n,y0,y1,z0,z1);
  CVF(z0,n); CVF(z1,n);
  return e;
}

/**
 @brief 共役転置行列の線形変換 [y0,y1]=[A0,A1]'*[x0,x1]
*/
int icvec_lintr_ct(int m, int n, cmulti **y0, cmulti **y1, cmulti **A0, int LDA0, cmulti **A1, int LDA1, cmulti **x0, cmulti **x1)
{
  int prec,i,j,e=0;
  cmulti **z0=NULL,**z1=NULL;
  prec=MAX2(cvec_get_prec_max(n,y0),cvec_get_prec_max(n,y1));
  CVA(z0,n,prec); CVA(z1,n,prec);
  for(j=0; j<n; j++){
    e+=icset_d(z0[j],z1[j],0);
    for(i=0; i<m; i++){
      e+=icadd_dot(z0[j],z1[j],MAT(A0,i,j,LDA0),MAT(A1,i,j,LDA1),x0[i],x1[i]);
    }
  }
  e+=icvec_copy(n,y0,y1,z0,z1);
  CVF(z0,n); CVF(z1,n);
  return e;
}

/**
 @brief 共役転置行列の線形変換の加算 [y0,y1]+=[A0,A1]'*[x0,x1]
*/
int icvec_add_lintr_ct(int m, int n, cmulti **y0, cmulti **y1, cmulti **A0, int LDA0, cmulti **A1, int LDA1, cmulti **x0, cmulti **x1)
{
  int prec,i,j,e=0;
  cmulti **z0=NULL,**z1=NULL;
  prec=MAX2(cvec_get_prec_max(n,y0),cvec_get_prec_max(n,y1));
  CVA(z0,n,prec); CVA(z1,n,prec);
  for(j=0; j<n; j++){
    e+=iccopy(z0[j],z1[j],y0[j],y1[j]);
    for(i=0; i<m; i++){
      e+=icadd_dot(z0[j],z1[j],MAT(A0,i,j,LDA0),MAT(A1,i,j,LDA1),x0[i],x1[i]);
    }
  }
  e+=icvec_copy(n,y0,y1,z0,z1);
  CVF(z0,n); CVF(z1,n);
  return e;
}

/**
 @brief 共役転置行列の線形変換の減算 [y0,y1]-=[A0,A1]'*[x0,x1]
*/
int icvec_sub_lintr_ct(int m, int n, cmulti **y0, cmulti **y1, cmulti **A0, int LDA0, cmulti **A1, int LDA1, cmulti **x0, cmulti **x1)
{
  int prec,i,j,e=0;
  cmulti **z0=NULL,**z1=NULL;
  prec=MAX2(cvec_get_prec_max(n,y0),cvec_get_prec_max(n,y1));
  CVA(z0,n,prec); CVA(z1,n,prec);
  for(j=0; j<n; j++){
    e+=iccopy(z0[j],z1[j],y0[j],y1[j]);
    for(i=0; i<m; i++){
      e+=icsub_dot(z0[j],z1[j],MAT(A0,i,j,LDA0),MAT(A1,i,j,LDA1),x0[i],x1[i]);
    }
  }
  e+=icvec_copy(n,y0,y1,z0,z1);
  CVF(z0,n); CVF(z1,n);
  return e;
}


/** @} */

//////////////////////////////////////////////////////////

/** @name 写像 */
/** @{ */

/**
 @brief 写像 [y0,y1]=f([x0,x1])
*/
// [y0,y1]=f([x0,x1])
int icvec_func(cmulti *y0, cmulti *y1, func_t *f, int n, cmulti **x0, cmulti **x1)
{
  int prec,i,e=0;
  cmulti *a0=NULL,*a1=NULL,*b0=NULL,*b1=NULL,*z0=NULL,*z1=NULL;
  prec=MAX2(cget_prec(y0),cget_prec(y1));
  CA(a0,prec); CA(a1,prec); CA(b0,prec); CA(b1,prec); CA(z0,prec); CA(z1,prec);
  cset_nan(z0); cset_nan(z1);
  if(f==NULL)                { FUNC_ERROR_ARG1("icvec_func",f); }
  else if(func_is(f,"nan"))  { cset_nan(z0); cset_nan(z1); }
  else if(func_is(f,"inf"))  { cset_inf(z0,1,1); cset_inf(z1,1,1); }
  else if(func_is_zero(f))   { e+=icset_dd(z0,z1,0,0); }
  else if(func_is_one(f))    { e+=icset_dd(z0,z1,1,0); }
  else if(func_is_bigint(f)) { e+=icset_bigint(z0,z1,func_bigint_p(f)); }
  else if(func_is_real(f))   { e+=ccopy_r(z0,func_real_p(f)); e+=ccopy_r(z1,func_real_p(f)); }
  else if(func_is_complex(f)){ e+=ccopy(z0,func_complex_p(f)); e+=ccopy(z1,func_complex_p(f)); }
  else if(func_is_var(f))    {
    icset_dd(z0,z1,1,0);
    for(i=0; i<func_var_size(f); i++){
      if(func_var_pow(f,i)!=0){
	if(0<=func_var_num(f,i) && func_var_num(f,i)<n){
	  e+=icpow_si(a0,a1,x0[func_var_num(f,i)],x1[func_var_num(f,i)],func_var_pow(f,i));
	}else{ cset_nan(a0); cset_nan(a1); }
	e+=icmul(z0,z1,z0,z1,a0,a1);
      }
    }
  }
  else if(func_is_add(f))  {
    icset_dd(z0,z1,0,0);
    for(i=0; i<func_asize(f); i++){
      e+=icvec_func(a0,a1,func_aget(f,i),n,x0,x1);
      e+=icadd(z0,z1,z0,z1,a0,a1);
    }
  }
  else if(func_is_mul(f))  {
    icset_dd(z0,z1,1,0);
    for(i=0; i<func_asize(f); i++){
      e+=icvec_func(a0,a1,func_aget(f,i),n,x0,x1);
      e+=icmul(z0,z1,z0,z1,a0,a1);
    }
  }
  else if(func_is(f,"sqrt")) { e+=icvec_func(a0,a1,func_aget(f,0),n,x0,x1); e+=icsqrt(z0,z1,a0,a1); }
  else if(func_is(f,"exp"))  { e+=icvec_func(a0,a1,func_aget(f,0),n,x0,x1); e+=icexp(z0,z1,a0,a1); }
  else if(func_is(f,"log"))  { e+=icvec_func(a0,a1,func_aget(f,0),n,x0,x1); e+=iclog(z0,z1,a0,a1); }
  else if(func_is(f,"sin"))  { e+=icvec_func(a0,a1,func_aget(f,0),n,x0,x1); e+=icsin(z0,z1,a0,a1); }
  else if(func_is(f,"cos"))  { e+=icvec_func(a0,a1,func_aget(f,0),n,x0,x1); e+=iccos(z0,z1,a0,a1); }
  else if(func_is(f,"tan"))  { e+=icvec_func(a0,a1,func_aget(f,0),n,x0,x1); e+=ictan(z0,z1,a0,a1); }
  else if(func_is(f,"asin")) { e+=icvec_func(a0,a1,func_aget(f,0),n,x0,x1); e+=icasin(z0,z1,a0,a1); }
  else if(func_is(f,"acos")) { e+=icvec_func(a0,a1,func_aget(f,0),n,x0,x1); e+=icacos(z0,z1,a0,a1); }
  else if(func_is(f,"atan")) { e+=icvec_func(a0,a1,func_aget(f,0),n,x0,x1); e+=icatan(z0,z1,a0,a1); }
  else if(func_is(f,"sinh")) { e+=icvec_func(a0,a1,func_aget(f,0),n,x0,x1); e+=icsinh(z0,z1,a0,a1); }
  else if(func_is(f,"cosh")) { e+=icvec_func(a0,a1,func_aget(f,0),n,x0,x1); e+=iccosh(z0,z1,a0,a1); }
  else if(func_is(f,"tanh")) { e+=icvec_func(a0,a1,func_aget(f,0),n,x0,x1); e+=ictanh(z0,z1,a0,a1); }
  else if(func_is(f,"asinh")){ e+=icvec_func(a0,a1,func_aget(f,0),n,x0,x1); e+=icasinh(z0,z1,a0,a1); }
  else if(func_is(f,"acosh")){ e+=icvec_func(a0,a1,func_aget(f,0),n,x0,x1); e+=icacosh(z0,z1,a0,a1); }
  else if(func_is(f,"atanh")){ e+=icvec_func(a0,a1,func_aget(f,0),n,x0,x1); e+=icatanh(z0,z1,a0,a1); }
  else if(func_is(f,"pow"))  {
    e+=icvec_func(a0,a1,func_aget(f,0),n,x0,x1);
    e+=icvec_func(b0,b1,func_aget(f,1),n,x0,x1);
    e+=cpow_c(z0,a0,b0); ERROR_AT;
  }
  if(func_has_power(f))      { e+=icpow_si(z0,z1,z0,z1,func_power(f)); }
  e+=iccopy(y0,y1,z0,z1);
  CF(a0); CF(a1); CF(b0); CF(b1); CF(z0); CF(z1);
  return e;
}

/**
 @brief ベクトル写像 [y0,y1]=f([x0,x1])
 */
int icvec_func_list(int m, cmulti **y0, cmulti **y1, func_t *f, int n, cmulti **x0, cmulti **x1)
{
  int i,e=0;
  if(func_is_list(f)){
    for(i=0; i<m && i<func_asize(f); i++){
      e+=icvec_func(y0[i],y1[i],func_aget(f,i),n,x0,x1);
    }
  }
  return e;
}

/** @} */

//EOF
