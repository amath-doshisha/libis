#include"is_rmulti.h"
#include"is_rmat.h"
#include"is_irmulti.h"
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
#define RF(X)        ((X)=rfree(X))
#define RVA(X,Y,N)   { X=rvec_allocate_prec(N,rvec_get_prec_max(N,Y)); }
#define RVF(X,N)     { X=rvec_free(N,X); }

/** @name 基本操作 */
/** @{ */

/**
 @brief rmulti型のベクトルの値を倍精度実数から設定.
*/
int irvec_set_d(int n, rmulti **y0, rmulti **y1, double *x)
{
  int i,e=0;
  for(i=0; i<n; i++){ e+=irset_d(y0[i],y1[i],x[i]); }
  return e;
}

/**
 @brief コピー [y0,y1]=[x0,x1]
*/
int irvec_copy(int n, rmulti **y0, rmulti **y1, rmulti **x0, rmulti **x1)
{
  int i,e=0;
  for(i=0; i<n; i++){ e+=ircopy(y0[i],y1[i],x0[i],x1[i]); }
  return e;
}

/**
 @brief 表示
 */
void irvec_print(int n, rmulti **x0, rmulti **x1, const char *name, const char *f, int digits)
{
  int i,k;
  char format[128];
  if(STR_EQ(f,"f") || STR_EQ(f,"F")){ k=3; }
  else if((strcmp(f,"e")==0) || (strcmp(f,"E")==0)){ k=7; }
  else if((strcmp(f,"g")==0) || (strcmp(f,"G")==0)){ k=6; }
  else{ k=3; }
  sprintf(format,"[%%%d.%dR%s, %%%d.%dR%s]\n",digits+k,digits,f,digits+k,digits,f);
  if(name!=NULL){ printf("%s\n",name); }
  if(x0==NULL || x1==NULL){ printf("NULL\n"); return; }
  for(i=0; i<n; i++){ mpfr_printf(format,x0[i],x1[i]); }
}


/**
 @brief 区間の中心 xc=(x1+x0)/2 と半径 xr=x1-x0
*/
int irvec_center_radius(int n, rmulti **xc, rmulti **xr, rmulti **x0, rmulti **x1)
{
  int e=0;
  mpfr_rnd_t mode;
  mode=get_round_mode();
  set_round_mode(MPFR_RNDU);  // up
  e+=rvec_sub(n,xc,x1,x0);    // xc=x1-x0
  e+=rvec_mul_d(n,xc,xc,0.5); // xc=(x1-x0)/2
  e+=rvec_add(n,xc,xc,x0);    // xc=(x1-x0)/2+x0
  e+=rvec_sub(n,xr,xc,x0);    // xr=xc-x0
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
int irvec_neg(int n, rmulti **y0, rmulti **y1, rmulti **x0, rmulti **x1)
{
  int i,e=0;
  for(i=0; i<n; i++){ e+=irneg(y0[i],y1[i],x0[i],x1[i]); }
  return e;
}

/**
 @brief 符号の正負 [y0,y1]=[-abs(x),abs(x)]
*/
int irvec_pm(int n, rmulti **y0, rmulti **y1, rmulti **x)
{
  int i,e=0;
  for(i=0; i<n; i++){ e+=irpm(y0[i],y1[i],x[i]); }
  return e;
}

/**
 @brief 区間の拡張 [z0,z1]=[x0,x1]+[-y,y]
*/
int irvec_add_pm(int n, rmulti **z0, rmulti **z1, rmulti **x0, rmulti **x1, rmulti **y)
{
  int i,e=0;
  for(i=0; i<n; i++){ e+=iradd_pm(z0[i],z1[i],x0[i],x1[i],y[i]); }
  return e;
}

/**
 @brief 区間の中心 [m-r,m+r]=[x0,x1]
 */
int irvec_mid(int n, rmulti **mid, rmulti **x0, rmulti **x1)
{
  int i,e=0;
  for(i=0; i<n; i++){ e+=irmid(mid[i],x0[i],x1[i]); }
  return e;
}

/**
 @brief 区間の半径 [m-r,m+r]=[x0,x1]
 */
int irvec_rad(int n, rmulti **rad, rmulti **x0, rmulti **x1)
{
  int i,e=0;
  for(i=0; i<n; i++){ e+=irrad(rad[i],x0[i],x1[i]); }
  return e;
}

/**
 @brief 区間の中心と半径 [m-r,m+r]=[x0,x1]
 */
int irvec_mr(int n, rmulti **mid, rmulti **rad, rmulti **x0, rmulti **x1)
{
  int i,e=0;
  for(i=0; i<n; i++){ e+=irmr(mid[i],rad[i],x0[i],x1[i]); }
  return e;
}

/**
 @brief 足し算 [z0,z1]=[x0,z1]+[y0,y1]
*/
int irvec_add(int n, rmulti **z0, rmulti **z1, rmulti **x0, rmulti **x1, rmulti **y0, rmulti **y1)
{
  int i,e=0;
  for(i=0; i<n; i++){ e+=iradd(z0[i],z1[i],x0[i],x1[i],y0[i],y1[i]); }
  return e;
}

/**
 @brief 足し算 [z0,z1]=[x0,z1]+[y0,y1]
*/
int irvec_add_r(int n, rmulti **z0, rmulti **z1, rmulti **x0, rmulti **x1, rmulti *y0, rmulti *y1)
{
  int i,e=0;
  for(i=0; i<n; i++){ e+=iradd(z0[i],z1[i],x0[i],x1[i],y0,y1); }
  return e;
}

/**
 @brief 引き算 [z0,z1]=[x0,z1]-[y0,y1]
*/
int irvec_sub(int n, rmulti **z0, rmulti **z1, rmulti **x0, rmulti **x1, rmulti **y0, rmulti **y1)
{
  int i,e=0;
  for(i=0; i<n; i++){ e+=irsub(z0[i],z1[i],x0[i],x1[i],y0[i],y1[i]); }
  return e;
}

/**
 @brief 掛け算 [z0,z1]=[x0,x1]*[y0,y1]
 */
int irvec_mul_r(int n, rmulti **z0, rmulti **z1, rmulti **x0, rmulti **x1, rmulti *y0, rmulti *y1)
{
  int i,e=0;
  for(i=0; i<n; i++){ e+=irmul(z0[i],z1[i],x0[i],x1[i],y0,y1); }
  return e;
}

/**
 @brief 割り算 [z0,z1]=[x0,x1]/[y0,y1]
 */
int irvec_div(int n, rmulti **z0, rmulti **z1, rmulti **x0, rmulti **x1, rmulti **y0, rmulti **y1)
{
  int i,e=0;
  for(i=0; i<n; i++){ e+=irdiv(z0[i],z1[i],x0[i],x1[i],y0[i],y1[i]); }
  return e;
}

/**
 @brief 平方の総和 [y0,y1]=sum([x0,x1].^2)
*/
int irvec_sum_pow2(rmulti *y0, rmulti *y1, int n, rmulti **x0, rmulti **x1)
{
  int i,e=0;
  e+=irset_d(y0,y1,0);
  for(i=0; i<n; i++){
    e+=iradd_mul(y0,y1,x0[i],x1[i],x0[i],x1[i]);
  }
  return e;
}

/** @} */

//////////////////////////////////////////////////////

/** @name 数学関数 */
/** @{ */

/**
 @brief 絶対値 [y0,y1]=abs([x0,x1])
*/
int irvec_abs(int n, rmulti **y0, rmulti **y1, rmulti **x0, rmulti **x1)
{
  int i,e=0;
  for(i=0; i<n; i++){
    e+=irabs(y0[i],y1[i],x0[i],x1[i]);
  }
  return e;
}

/**
 @brief 差の絶対値 [z0,z1]=abs([x0,x1]-[y0,y1])
*/
int irvec_abs_sub(int n, rmulti **z0, rmulti **z1, rmulti **x0, rmulti **x1, rmulti **y0, rmulti **y1)
{
  int i,e=0;
  for(i=0; i<n; i++){
    e+=irabs_sub(z0[i],z1[i],x0[i],x1[i],y0[i],y1[i]);
  }
  return e;
}

/**
 @brief 最大値 [y0,y1]=[x0,max(x1)]
*/
int irvec_umax(rmulti *y0, rmulti *y1, int n, rmulti **x0, rmulti **x1)
{
  int i,e=0;
  e+=ircopy(y0,y1,x0[0],x1[0]);      // y0=x0[0], y1=x1[0]
  for(i=1; i<n; i++){
    if(rgt(x1[i],y1)){               // x1[i]>y1
      e+=ircopy(y0,y1,x0[i],x1[i]);  // y0=x0[i], y1=x1[i]
    }
  }
  return e;
}

/**
 @brief 最大値 [y0,y1]=[max(x0),x1]
*/
int irvec_dmax(rmulti *y0, rmulti *y1, int n, rmulti **x0, rmulti **x1)
{
  int i,e=0;
  e+=ircopy(y0,y1,x0[0],x1[0]);      // y0=x0[0], y1=x1[0]
  for(i=1; i<n; i++){
    if(rgt(x0[i],y0)){               // x0[i]>y0
      e+=ircopy(y0,y1,x0[i],x1[i]);  // y0=x0[i], y1=x1[i]
    }
  }
  return e;
}

/**
 @brief 最小値 [y0,y1]=[x0,min(x1)]
*/
int irvec_umin(rmulti *y0, rmulti *y1, int n, rmulti **x0, rmulti **x1)
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
int irvec_dmin(rmulti *y0, rmulti *y1, int n, rmulti **x0, rmulti **x1)
{
  int i,e=0;
  e+=ircopy(y0,y1,x0[0],x1[0]);      // y0=x0[0], y1=x1[0]
  for(i=1; i<n; i++){
    if(rlt(x0[i],y0)){               // x0[i]<y0
      e+=ircopy(y0,y1,x0[i],x1[i]);  // y0=x0[i], y1=x1[i]
    }
  }
  return e;
}

/**
 @brief irmulti型のベクトルの要素の絶対値の最大値 value=max(abs(x))
*/
int irvec_umax_abs(rmulti *y0, rmulti *y1, int n, rmulti **x0, rmulti **x1)
{
  int prec,e=0;
  rmulti **ax0=NULL,**ax1=NULL;
  prec=MAX2(rget_prec(y0),rget_prec(y1));
  ax0=rvec_allocate_prec(n,prec);
  ax1=rvec_allocate_prec(n,prec);
  irvec_abs(n,ax0,ax1,x0,x1);
  irvec_umax(y0,y1,n,ax0,ax1);
  RVF(ax0,n); RVF(ax1,n);
  return e;
}

/**
 @brief irmulti型のベクトルの要素の絶対値の最小値 value=min(abs(x))
*/
int irvec_dmin_abs(rmulti *y0, rmulti *y1, int n, rmulti **x0, rmulti **x1)
{
  int prec,e=0;
  rmulti **ax0=NULL,**ax1=NULL;
  prec=MAX2(rget_prec(y0),rget_prec(y1));
  ax0=rvec_allocate_prec(n,prec);
  ax1=rvec_allocate_prec(n,prec);
  irvec_abs(n,ax0,ax1,x0,x1);
  irvec_dmin(y0,y1,n,ax0,ax1);
  RVF(ax0,n); RVF(ax1,n);
  return e;
}

//追加

/**
 @brief irmulti型のベクトルの要素の総和 value=sum(x)
 */
int irvec_sum(rmulti *y0, rmulti *y1, int n, rmulti **x0, rmulti **x1)
{
  int i,e=0;
  e+=irset_d(y0,y1,0);
  for(i=0; i<n; i++){ e+=iradd(y0,y1,y0,y1,x0[i],x1[i]); }
  return e;
}

//ここまで

/** @} */

//////////////////////////////////////////////////////////////////////////

/** @name 線形変換 */
/** @{ */

/**
 @brief 線形変換 [y0,y1]=[A0,A1]*[x0,x1]
*/
int irvec_lintr(int m, int n, rmulti **y0, rmulti **y1, rmulti **A0, int LDA0, rmulti **A1, int LDA1, rmulti **x0, rmulti **x1)
{
  int i,j,e=0;
  rmulti **z0=NULL,**z1=NULL;
  RVA(z0,y0,m); RVA(z1,y1,m);
  for(i=0; i<m; i++){
    e+=irset_d(z0[i],z1[i],0);
    for(j=0; j<n; j++){
      e+=iradd_mul(z0[i],z1[i],MAT(A0,i,j,LDA0),MAT(A1,i,j,LDA1),x0[j],x1[j]);
    }
  }
  e+=irvec_copy(m,y0,y1,z0,z1);
  RVF(z0,m); RVF(z1,m);
  return e;
}

/**
 @brief 線形変換の加算 [y0,y1]+=[A0,A1]*[x0,x1]
*/
int irvec_add_lintr(int m, int n, rmulti **y0, rmulti **y1, rmulti **A0, int LDA0, rmulti **A1, int LDA1, rmulti **x0, rmulti **x1)
{
  int i,j,e=0;
  rmulti **z0=NULL,**z1=NULL;
  RVA(z0,y0,m); RVA(z1,y1,m);
  for(i=0; i<m; i++){
    e+=ircopy(z0[i],z1[i],y0[i],y1[i]);
    for(j=0; j<n; j++){
      e+=iradd_mul(z0[i],z1[i],MAT(A0,i,j,LDA0),MAT(A1,i,j,LDA1),x0[j],x1[j]);
    }
  }
  e+=irvec_copy(m,y0,y1,z0,z1);
  RVF(z0,m); RVF(z1,m);
  return e;
}

/**
 @brief 線形変換の減算 [y0,y1]-=[A0,A1]*[x0,x1]
*/
int irvec_sub_lintr(int m, int n, rmulti **y0, rmulti **y1, rmulti **A0, int LDA0, rmulti **A1, int LDA1, rmulti **x0, rmulti **x1)
{
  int i,j,e=0;
  rmulti **z0=NULL,**z1=NULL;
  RVA(z0,y0,m); RVA(z1,y1,m);
  for(i=0; i<m; i++){
    e+=ircopy(z0[i],z1[i],y0[i],y1[i]);
    for(j=0; j<n; j++){
      e+=irsub_mul(z0[i],z1[i],MAT(A0,i,j,LDA0),MAT(A1,i,j,LDA1),x0[j],x1[j]);
    }
  }
  e+=irvec_copy(m,y0,y1,z0,z1);
  RVF(z0,m); RVF(z1,m);
  return e;
}

/**
 @brief 転置行列の線形変換 [y0,y1]=[A0,A1]'*[x0,x1]
*/
int irvec_lintr_t(int m, int n, rmulti **y0, rmulti **y1, rmulti **A0, int LDA0, rmulti **A1, int LDA1, rmulti **x0, rmulti **x1)
{
  int i,j,e=0;
  rmulti **z0=NULL,**z1=NULL;
  RVA(z0,y0,n); RVA(z1,y1,n);
  for(j=0; j<n; j++){
    e+=irset_d(z0[j],z1[j],0);
    for(i=0; i<m; i++){
      e+=iradd_mul(z0[j],z1[j],MAT(A0,i,j,LDA0),MAT(A1,i,j,LDA1),x0[i],x1[i]);
    }
  }
  e+=irvec_copy(n,y0,y1,z0,z1);
  RVF(z0,n); RVF(z1,n);
  return e;
}

/**
 @brief 転置行列の線形変換の加算 [y0,y1]+=[A0,A1]'*[x0,x1]
*/
int irvec_add_lintr_t(int m, int n, rmulti **y0, rmulti **y1, rmulti **A0, int LDA0, rmulti **A1, int LDA1, rmulti **x0, rmulti **x1)
{
  int i,j,e=0;
  rmulti **z0=NULL,**z1=NULL;
  RVA(z0,y0,n); RVA(z1,y1,n);
  for(j=0; j<n; j++){
    e+=ircopy(z0[j],z1[j],y0[j],y1[j]);
    for(i=0; i<m; i++){
      e+=iradd_mul(z0[j],z1[j],MAT(A0,i,j,LDA0),MAT(A1,i,j,LDA1),x0[i],x1[i]);
    }
  }
  e+=irvec_copy(n,y0,y1,z0,z1);
  RVF(z0,n); RVF(z1,n);
  return e;
}

/**
 @brief 転置行列の線形変換の減算 [y0,y1]-=[A0,A1]'*[x0,x1]
*/
int irvec_sub_lintr_t(int m, int n, rmulti **y0, rmulti **y1, rmulti **A0, int LDA0, rmulti **A1, int LDA1, rmulti **x0, rmulti **x1)
{
  int i,j,e=0;
  rmulti **z0=NULL,**z1=NULL;
  RVA(z0,y0,n); RVA(z1,y1,n);
  for(j=0; j<n; j++){
    e+=ircopy(z0[j],z1[j],y0[j],y1[j]);
    for(i=0; i<m; i++){
      e+=irsub_mul(z0[j],z1[j],MAT(A0,i,j,LDA0),MAT(A1,i,j,LDA1),x0[i],x1[i]);
    }
  }
  e+=irvec_copy(n,y0,y1,z0,z1);
  RVF(z0,n); RVF(z1,n);
  return e;
}

/** @} */

///////////////////////////////////////////////////////////////////////////////////////////////////////////

/** @name 写像 */
/** @{ */

/**
 @brief 写像 [y0,y1]=f([x0,x1])
*/
int irvec_func(rmulti *y0, rmulti *y1, func_t *f, int n, rmulti **x0, rmulti **x1)
{
  int i,e=0;
  rmulti *a0=NULL,*a1=NULL,*b0=NULL,*b1=NULL,*z0=NULL,*z1=NULL;
  RA2(a0,y0,y1); RA2(a1,y0,y1); RA2(b0,y0,y1); RA2(b1,y0,y1); RA2(z0,y0,y1); RA2(z1,y0,y1);
  rset_nan(z0); rset_nan(z1);  
  if(f==NULL)                { FUNC_ERROR_ARG1("irvec_func",f); }
  else if(func_is(f,"inf"))  { rset_inf(z0,1); rset_inf(z1,1); }
  else if(func_is_zero(f))   { e+=irset_d(z0,z1,0); }
  else if(func_is_one(f))    { e+=irset_d(z0,z1,1); }
  else if(func_is_bigint(f)) { e+=irset_bigint(z0,z1,func_bigint_p(f)); }
  else if(func_is_real(f))   { e+=ircopy(z0,z1,func_real_p(f),func_real_p(f)); }
  else if(func_is_var(f))    {
    e+=irset_d(z0,z1,1);
    for(i=0; i<func_var_size(f); i++){
      if(func_var_pow(f,i)!=0){
	if(0<=func_var_num(f,i) && func_var_num(f,i)<n){
	  e+=irpow_si(a0,a1,x0[func_var_num(f,i)],x1[func_var_num(f,i)],func_var_pow(f,i));
	}else{ rset_nan(a0); rset_nan(a1); }
	e+=irmul(z0,z1,z0,z1,a0,a1);
      }
    }
  }
  else if(func_is_add(f))    {
    e+=irset_d(z0,z1,0);
    for(i=0; i<func_asize(f); i++){
      e+=irvec_func(a0,a1,func_aget(f,i),n,x0,x1);
      e+=iradd(z0,z1,z0,z1,a0,a1);
    }
  }
  else if(func_is_mul(f))    {
    e+=irset_d(z0,z1,1);
    for(i=0; i<func_asize(f); i++){
      e+=irvec_func(a0,a1,func_aget(f,i),n,x0,x1);
      e+=irmul(z0,z1,z0,z1,a0,a1);
    }
  }
  else if(func_is(f,"sqrt")) { e+=irvec_func(a0,a1,func_aget(f,0),n,x0,x1); e+=irsqrt(z0,z1,a0,a1); }
  else if(func_is(f,"exp"))  { e+=irvec_func(a0,a1,func_aget(f,0),n,x0,x1); e+=irexp(z0,z1,a0,a1); }
  else if(func_is(f,"log"))  { e+=irvec_func(a0,a1,func_aget(f,0),n,x0,x1); e+=irlog(z0,z1,a0,a1); }
  else if(func_is(f,"sin"))  { e+=irvec_func(a0,a1,func_aget(f,0),n,x0,x1); e+=irsin(z0,z1,a0,a1); }
  else if(func_is(f,"cos"))  { e+=irvec_func(a0,a1,func_aget(f,0),n,x0,x1); e+=ircos(z0,z1,a0,a1); }
  else if(func_is(f,"tan"))  { e+=irvec_func(a0,a1,func_aget(f,0),n,x0,x1); e+=irtan(z0,z1,a0,a1); }
  else if(func_is(f,"asin")) { e+=irvec_func(a0,a1,func_aget(f,0),n,x0,x1); e+=irasin(z0,z1,a0,a1); }
  else if(func_is(f,"acos")) { e+=irvec_func(a0,a1,func_aget(f,0),n,x0,x1); e+=iracos(z0,z1,a0,a1); }
  else if(func_is(f,"atan")) { e+=irvec_func(a0,a1,func_aget(f,0),n,x0,x1); e+=iratan(z0,z1,a0,a1); }
  else if(func_is(f,"sinh")) { e+=irvec_func(a0,a1,func_aget(f,0),n,x0,x1); e+=irsinh(z0,z1,a0,a1); }
  else if(func_is(f,"cosh")) { e+=irvec_func(a0,a1,func_aget(f,0),n,x0,x1); e+=ircosh(z0,z1,a0,a1); }
  else if(func_is(f,"tanh")) { e+=irvec_func(a0,a1,func_aget(f,0),n,x0,x1); e+=irtanh(z0,z1,a0,a1); }
  else if(func_is(f,"asinh")){ e+=irvec_func(a0,a1,func_aget(f,0),n,x0,x1); e+=irasinh(z0,z1,a0,a1); }
  else if(func_is(f,"acosh")){ e+=irvec_func(a0,a1,func_aget(f,0),n,x0,x1); e+=iracosh(z0,z1,a0,a1); }
  else if(func_is(f,"atanh")){ e+=irvec_func(a0,a1,func_aget(f,0),n,x0,x1); e+=iratanh(z0,z1,a0,a1); }
  else if(func_is(f,"pow"))  {
    e+=irvec_func(a0,a1,func_aget(f,0),n,x0,x1);
    e+=irvec_func(b0,b1,func_aget(f,1),n,x0,x1);
    e+=rpow(z0,a0,b0); ERROR_AT;
  }
  if(func_has_power(f))      { e+=irpow_si(z0,z1,z0,z1,func_power(f)); }
  e+=ircopy(y0,y1,z0,z1);
  RF(a0); RF(a1); RF(b0); RF(b1); RF(z0); RF(z1);
  return e;
}

/**
 @brief ベクトル写像 [y0,y1]=f([x0,x1])
 */
int irvec_func_list(int m, rmulti **y0, rmulti **y1, func_t *f, int n, rmulti **x0, rmulti **x1)
{
  int i,e=0;
  for(i=0; i<m; i++){
    if(func_is_list(f) && func_is_list(f) && i<func_asize(f)){
      e+=irvec_func(y0[i],y1[i],func_aget(f,i),n,x0,x1);
    }else{
      rset_nan(y0[i]); rset_nan(y1[i]);
    }
  }
  return e;
}

/** @} */

//EOF
