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


/** @name irmulti型ベクトルの値の設定に関する関数 */
/** @{ */

/**
 @brief rmulti型のベクトルの値を倍精度実数から設定.
*/
void irvec_set_d(int n, rmulti **y0, rmulti **y1, double *x)
{
  int i;
  for(i=0; i<n; i++){ irset_d(y0[i],y1[i],x[i]); }
}

/**
 @brief コピー [y0,y1]=[x0,x1]
*/
void irvec_copy(int n, rmulti **y0, rmulti **y1, rmulti **x0, rmulti **x1)
{
  int i;
  for(i=0; i<n; i++){ ircopy(y0[i],y1[i],x0[i],x1[i]); }
}

/** @} */
/** @name irmulti型ベクトルの入出力に関する関数 */
/** @{ */

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


/** @} */
/** @name irmulti型ベクトルの型変換に関する関数 */
/** @{ */

/** @} */
/** @name irmulti型ベクトルに関する関数子 */
/** @{ */

/**
 @brief 区間の中心 xc=(x1+x0)/2 と半径 xr=x1-x0
*/
void irvec_center_radius(int n, rmulti **xc, rmulti **xr, rmulti **x0, rmulti **x1)
{
  mpfr_rnd_t mode;
  mode=get_round_mode();
  set_round_mode(MPFR_RNDU);  // up
  rvec_sub(n,xc,x1,x0);    // xc=x1-x0
  rvec_mul_d(n,xc,xc,0.5); // xc=(x1-x0)/2
  rvec_add(n,xc,xc,x0);    // xc=(x1-x0)/2+x0
  rvec_sub(n,xr,xc,x0);    // xr=xc-x0
  set_round_mode(mode);       // back
}

/**
 @brief 符号の反転 [y0,y1]=-[x0,x1]
*/
void irvec_neg(int n, rmulti **y0, rmulti **y1, rmulti **x0, rmulti **x1)
{
  int i;
  for(i=0; i<n; i++){ irneg(y0[i],y1[i],x0[i],x1[i]); }
}

/**
 @brief 符号の正負 [y0,y1]=[-abs(x),abs(x)]
*/
void irvec_pm(int n, rmulti **y0, rmulti **y1, rmulti **x)
{
  int i;
  for(i=0; i<n; i++){ irpm(y0[i],y1[i],x[i]); }
}

/**
 @brief 区間の拡張 [z0,z1]=[x0,x1]+[-y,y]
*/
void irvec_add_pm(int n, rmulti **z0, rmulti **z1, rmulti **x0, rmulti **x1, rmulti **y)
{
  int i;
  for(i=0; i<n; i++){ iradd_pm(z0[i],z1[i],x0[i],x1[i],y[i]); }
}

/**
 @brief 区間の中心 [m-r,m+r]=[x0,x1]
 */
void irvec_mid(int n, rmulti **mid, rmulti **x0, rmulti **x1)
{
  int i;
  for(i=0; i<n; i++){ irmid(mid[i],x0[i],x1[i]); }
}

/**
 @brief 区間の半径 [m-r,m+r]=[x0,x1]
 */
void irvec_rad(int n, rmulti **rad, rmulti **x0, rmulti **x1)
{
  int i;
  for(i=0; i<n; i++){ irrad(rad[i],x0[i],x1[i]); }
}

/**
 @brief 区間の中心と半径 [m-r,m+r]=[x0,x1]
 */
void irvec_mr(int n, rmulti **mid, rmulti **rad, rmulti **x0, rmulti **x1)
{
  int i;
  for(i=0; i<n; i++){ irmr(mid[i],rad[i],x0[i],x1[i]); }
}

/**
 @brief 足し算 [z0,z1]=[x0,z1]+[y0,y1]
*/
void irvec_add(int n, rmulti **z0, rmulti **z1, rmulti **x0, rmulti **x1, rmulti **y0, rmulti **y1)
{
  int i;
  for(i=0; i<n; i++){ iradd(z0[i],z1[i],x0[i],x1[i],y0[i],y1[i]); }
}

/**
 @brief 足し算 [z0,z1]=[x0,z1]+[y0,y1]
*/
void irvec_add_r(int n, rmulti **z0, rmulti **z1, rmulti **x0, rmulti **x1, rmulti *y0, rmulti *y1)
{
  int i;
  for(i=0; i<n; i++){ iradd(z0[i],z1[i],x0[i],x1[i],y0,y1); }
}

/**
 @brief 引き算 [z0,z1]=[x0,z1]-[y0,y1]
*/
void irvec_sub(int n, rmulti **z0, rmulti **z1, rmulti **x0, rmulti **x1, rmulti **y0, rmulti **y1)
{
  int i;
  for(i=0; i<n; i++){ irsub(z0[i],z1[i],x0[i],x1[i],y0[i],y1[i]); }
}

/**
 @brief 掛け算 [z0,z1]=[x0,x1]*[y0,y1]
 */
void irvec_mul_r(int n, rmulti **z0, rmulti **z1, rmulti **x0, rmulti **x1, rmulti *y0, rmulti *y1)
{
  int i;
  for(i=0; i<n; i++){ irmul(z0[i],z1[i],x0[i],x1[i],y0,y1); }
}

/**
 @brief 割り算 [z0,z1]=[x0,x1]/[y0,y1]
 */
void irvec_div(int n, rmulti **z0, rmulti **z1, rmulti **x0, rmulti **x1, rmulti **y0, rmulti **y1)
{
  int i;
  for(i=0; i<n; i++){ irdiv(z0[i],z1[i],x0[i],x1[i],y0[i],y1[i]); }
}

/**
 @brief 平方の総和 [y0,y1]=sum([x0,x1].^2)
*/
void irvec_sum_pow2(rmulti *y0, rmulti *y1, int n, rmulti **x0, rmulti **x1)
{
  int i;
  irset_d(y0,y1,0);
  for(i=0; i<n; i++){ iradd_mul(y0,y1,x0[i],x1[i],x0[i],x1[i]); }
}

/**
 @brief 絶対値 [y0,y1]=abs([x0,x1])
*/
void irvec_abs(int n, rmulti **y0, rmulti **y1, rmulti **x0, rmulti **x1)
{
  int i;
  for(i=0; i<n; i++){ irabs(y0[i],y1[i],x0[i],x1[i]); }
}

/**
 @brief 差の絶対値 [z0,z1]=abs([x0,x1]-[y0,y1])
*/
void irvec_abs_sub(int n, rmulti **z0, rmulti **z1, rmulti **x0, rmulti **x1, rmulti **y0, rmulti **y1)
{
  int i;
  for(i=0; i<n; i++){ irabs_sub(z0[i],z1[i],x0[i],x1[i],y0[i],y1[i]); }
}

/**
 @brief 最大値 [y0,y1]=[max(x0),max(x1)]
*/
void irvec_max(rmulti *y0, rmulti *y1, int n, rmulti **x0, rmulti **x1)
{
  int i;
  ircopy(y0,y1,x0[0],x1[0]);        // y0=x0[0], y1=x1[0]
  for(i=1; i<n; i++){
    if(rgt(x0[i],y0)){              // x0[i]>y0
      mpfr_set(y0,x0[i],MPFR_RNDD); // y0=x0[i]            
    }
    if(rgt(x1[i],y1)){              // x1[i]>y1
      mpfr_set(y1,x1[i],MPFR_RNDU); // y1=x1[i]              
    }
  }
}

/**
 @brief 最大値 [y0,y1]=[x0,max(x1)] 上限で比較
*/
void irvec_umax(rmulti *y0, rmulti *y1, int n, rmulti **x0, rmulti **x1)
{
  int i;
  ircopy(y0,y1,x0[0],x1[0]);      // y0=x0[0], y1=x1[0]
  for(i=1; i<n; i++){
    if(rgt(x1[i],y1)){            // x1[i]>y1
      ircopy(y0,y1,x0[i],x1[i]);  // y0=x0[i], y1=x1[i]
    }
  }
}

/**
 @brief 最大値 [y0,y1]=[max(x0),x1] 下限で比較
*/
void irvec_dmax(rmulti *y0, rmulti *y1, int n, rmulti **x0, rmulti **x1)
{
  int i;
  ircopy(y0,y1,x0[0],x1[0]);      // y0=x0[0], y1=x1[0]
  for(i=1; i<n; i++){
    if(rgt(x0[i],y0)){            // x0[i]>y0
      ircopy(y0,y1,x0[i],x1[i]);  // y0=x0[i], y1=x1[i]
    }
  }
}

/**
 @brief 最小値 [y0,y1]=[min(x0),min(x1)]
*/
void irvec_min(rmulti *y0, rmulti *y1, int n, rmulti **x0, rmulti **x1)
{
  int i;
  ircopy(y0,y1,x0[0],x1[0]);        // y0=x0[0], y1=x1[0]
  for(i=1; i<n; i++){
    if(rlt(x0[i],y0)){              // x0[i]<y0
      mpfr_set(y0,x0[i],MPFR_RNDD); // y0=x0[i]        
    }
    if(rlt(x1[i],y1)){              // x1[i]<y1
      mpfr_set(y1,x1[i],MPFR_RNDD); // y0=x0[i]       
    }
  }
}

/**
 @brief 最小値 [y0,y1]=[x0,min(x1)]
*/
void irvec_umin(rmulti *y0, rmulti *y1, int n, rmulti **x0, rmulti **x1)
{
  int i;
  ircopy(y0,y1,x0[0],x1[0]);     // y0=x0[0], y1=x1[0]
  for(i=1; i<n; i++){
    if(rlt(x1[i],y1)){           // x0[i]<y0
      ircopy(y0,y1,x0[i],x1[i]); // y0=x0[i], y1=x1[i]
    }
  }
}

/**
 @brief 最小値 [y0,y1]=[min(x0),x1]
*/
void irvec_dmin(rmulti *y0, rmulti *y1, int n, rmulti **x0, rmulti **x1)
{
  int i;
  ircopy(y0,y1,x0[0],x1[0]);     // y0=x0[0], y1=x1[0]
  for(i=1; i<n; i++){
    if(rlt(x0[i],y0)){           // x0[i]<y0
      ircopy(y0,y1,x0[i],x1[i]); // y0=x0[i], y1=x1[i]
    }
  }
}

/**
 @brief irmulti型のベクトルの要素の絶対値の最大値 value=max(abs(x))
*/
void irvec_umax_abs(rmulti *y0, rmulti *y1, int n, rmulti **x0, rmulti **x1)
{
  int prec;
  rmulti **ax0=NULL,**ax1=NULL;
  prec=MAX2(rget_prec(y0),rget_prec(y1));
  ax0=rvec_allocate_prec(n,prec);
  ax1=rvec_allocate_prec(n,prec);
  irvec_abs(n,ax0,ax1,x0,x1);
  irvec_umax(y0,y1,n,ax0,ax1);
  RVF(ax0,n); RVF(ax1,n);
}

/**
 @brief irmulti型のベクトルの要素の絶対値の最小値 value=min(abs(x))
*/
void irvec_dmin_abs(rmulti *y0, rmulti *y1, int n, rmulti **x0, rmulti **x1)
{
  int prec;
  rmulti **ax0=NULL,**ax1=NULL;
  prec=MAX2(rget_prec(y0),rget_prec(y1));
  ax0=rvec_allocate_prec(n,prec);
  ax1=rvec_allocate_prec(n,prec);
  irvec_abs(n,ax0,ax1,x0,x1);
  irvec_dmin(y0,y1,n,ax0,ax1);
  RVF(ax0,n); RVF(ax1,n);
}

/**
 @brief irmulti型のベクトルの要素の総和 value=sum(x)
 */
void irvec_sum(rmulti *y0, rmulti *y1, int n, rmulti **x0, rmulti **x1)
{
  int i;
  irset_d(y0,y1,0);
  for(i=0; i<n; i++){ iradd(y0,y1,y0,y1,x0[i],x1[i]); }
}

/**
 @brief 線形変換 [y0,y1]=[A0,A1]*[x0,x1]
*/
void irvec_lintr(int m, int n, rmulti **y0, rmulti **y1, rmulti **A0, int LDA0, rmulti **A1, int LDA1, rmulti **x0, rmulti **x1)
{
  int i,j;
  rmulti **z0=NULL,**z1=NULL;
  RVA(z0,y0,m); RVA(z1,y1,m);
  for(i=0; i<m; i++){
    irset_d(z0[i],z1[i],0);
    for(j=0; j<n; j++){
      iradd_mul(z0[i],z1[i],MAT(A0,i,j,LDA0),MAT(A1,i,j,LDA1),x0[j],x1[j]);
    }
  }
  irvec_copy(m,y0,y1,z0,z1);
  RVF(z0,m); RVF(z1,m);
}

/**
 @brief 線形変換の加算 [y0,y1]+=[A0,A1]*[x0,x1]
*/
void irvec_add_lintr(int m, int n, rmulti **y0, rmulti **y1, rmulti **A0, int LDA0, rmulti **A1, int LDA1, rmulti **x0, rmulti **x1)
{
  int i,j;
  rmulti **z0=NULL,**z1=NULL;
  RVA(z0,y0,m); RVA(z1,y1,m);
  for(i=0; i<m; i++){
    ircopy(z0[i],z1[i],y0[i],y1[i]);
    for(j=0; j<n; j++){
      iradd_mul(z0[i],z1[i],MAT(A0,i,j,LDA0),MAT(A1,i,j,LDA1),x0[j],x1[j]);
    }
  }
  irvec_copy(m,y0,y1,z0,z1);
  RVF(z0,m); RVF(z1,m);
}

/**
 @brief 線形変換の減算 [y0,y1]-=[A0,A1]*[x0,x1]
*/
void irvec_sub_lintr(int m, int n, rmulti **y0, rmulti **y1, rmulti **A0, int LDA0, rmulti **A1, int LDA1, rmulti **x0, rmulti **x1)
{
  int i,j;
  rmulti **z0=NULL,**z1=NULL;
  RVA(z0,y0,m); RVA(z1,y1,m);
  for(i=0; i<m; i++){
    ircopy(z0[i],z1[i],y0[i],y1[i]);
    for(j=0; j<n; j++){
      irsub_mul(z0[i],z1[i],MAT(A0,i,j,LDA0),MAT(A1,i,j,LDA1),x0[j],x1[j]);
    }
  }
  irvec_copy(m,y0,y1,z0,z1);
  RVF(z0,m); RVF(z1,m);
}

/**
 @brief 転置行列の線形変換 [y0,y1]=[A0,A1]'*[x0,x1]
*/
void irvec_lintr_t(int m, int n, rmulti **y0, rmulti **y1, rmulti **A0, int LDA0, rmulti **A1, int LDA1, rmulti **x0, rmulti **x1)
{
  int i,j;
  rmulti **z0=NULL,**z1=NULL;
  RVA(z0,y0,n); RVA(z1,y1,n);
  for(j=0; j<n; j++){
    irset_d(z0[j],z1[j],0);
    for(i=0; i<m; i++){
      iradd_mul(z0[j],z1[j],MAT(A0,i,j,LDA0),MAT(A1,i,j,LDA1),x0[i],x1[i]);
    }
  }
  irvec_copy(n,y0,y1,z0,z1);
  RVF(z0,n); RVF(z1,n);
}

/**
 @brief 転置行列の線形変換の加算 [y0,y1]+=[A0,A1]'*[x0,x1]
*/
void irvec_add_lintr_t(int m, int n, rmulti **y0, rmulti **y1, rmulti **A0, int LDA0, rmulti **A1, int LDA1, rmulti **x0, rmulti **x1)
{
  int i,j;
  rmulti **z0=NULL,**z1=NULL;
  RVA(z0,y0,n); RVA(z1,y1,n);
  for(j=0; j<n; j++){
    ircopy(z0[j],z1[j],y0[j],y1[j]);
    for(i=0; i<m; i++){
      iradd_mul(z0[j],z1[j],MAT(A0,i,j,LDA0),MAT(A1,i,j,LDA1),x0[i],x1[i]);
    }
  }
  irvec_copy(n,y0,y1,z0,z1);
  RVF(z0,n); RVF(z1,n);
}

/**
 @brief 転置行列の線形変換の減算 [y0,y1]-=[A0,A1]'*[x0,x1]
*/
void irvec_sub_lintr_t(int m, int n, rmulti **y0, rmulti **y1, rmulti **A0, int LDA0, rmulti **A1, int LDA1, rmulti **x0, rmulti **x1)
{
  int i,j;
  rmulti **z0=NULL,**z1=NULL;
  RVA(z0,y0,n); RVA(z1,y1,n);
  for(j=0; j<n; j++){
    ircopy(z0[j],z1[j],y0[j],y1[j]);
    for(i=0; i<m; i++){
      irsub_mul(z0[j],z1[j],MAT(A0,i,j,LDA0),MAT(A1,i,j,LDA1),x0[i],x1[i]);
    }
  }
  irvec_copy(n,y0,y1,z0,z1);
  RVF(z0,n); RVF(z1,n);
}

/**
 @brief 写像 [y0,y1]=f([x0,x1])
*/
void irvec_func(rmulti *y0, rmulti *y1, func_t *f, int n, rmulti **x0, rmulti **x1)
{
  int i;
  rmulti *a0=NULL,*a1=NULL,*b0=NULL,*b1=NULL,*z0=NULL,*z1=NULL;
  RA2(a0,y0,y1); RA2(a1,y0,y1); RA2(b0,y0,y1); RA2(b1,y0,y1); RA2(z0,y0,y1); RA2(z1,y0,y1);
  rset_nan(z0); rset_nan(z1);  
  if(f==NULL)                { FUNC_ERROR_ARG1("irvec_func",f); }
  else if(func_is(f,"inf"))  { rset_inf(z0,1); rset_inf(z1,1); }
  else if(func_is_zero(f))   { irset_d(z0,z1,0); }
  else if(func_is_one(f))    { irset_d(z0,z1,1); }
  else if(func_is_bigint(f)) { irset_bigint(z0,z1,func_bigint_p(f)); }
  else if(func_is_real(f))   { ircopy(z0,z1,func_real_p(f),func_real_p(f)); }
  else if(func_is_var(f))    {
    irset_d(z0,z1,1);
    for(i=0; i<func_var_size(f); i++){
      if(func_var_pow(f,i)!=0){
	if(0<=func_var_num(f,i) && func_var_num(f,i)<n){
	  irpow_si(a0,a1,x0[func_var_num(f,i)],x1[func_var_num(f,i)],func_var_pow(f,i));
	}else{ rset_nan(a0); rset_nan(a1); }
	irmul(z0,z1,z0,z1,a0,a1);
      }
    }
  }
  else if(func_is_add(f))    {
    irset_d(z0,z1,0);
    for(i=0; i<func_asize(f); i++){
      irvec_func(a0,a1,func_aget(f,i),n,x0,x1);
      iradd(z0,z1,z0,z1,a0,a1);
    }
  }
  else if(func_is_mul(f))    {
    irset_d(z0,z1,1);
    for(i=0; i<func_asize(f); i++){
      irvec_func(a0,a1,func_aget(f,i),n,x0,x1);
      irmul(z0,z1,z0,z1,a0,a1);
    }
  }
  else if(func_is(f,"sqrt")) { irvec_func(a0,a1,func_aget(f,0),n,x0,x1); irsqrt(z0,z1,a0,a1); }
  else if(func_is(f,"exp"))  { irvec_func(a0,a1,func_aget(f,0),n,x0,x1); irexp(z0,z1,a0,a1); }
  else if(func_is(f,"log"))  { irvec_func(a0,a1,func_aget(f,0),n,x0,x1); irlog(z0,z1,a0,a1); }
  else if(func_is(f,"sin"))  { irvec_func(a0,a1,func_aget(f,0),n,x0,x1); irsin(z0,z1,a0,a1); }
  else if(func_is(f,"cos"))  { irvec_func(a0,a1,func_aget(f,0),n,x0,x1); ircos(z0,z1,a0,a1); }
  else if(func_is(f,"tan"))  { irvec_func(a0,a1,func_aget(f,0),n,x0,x1); irtan(z0,z1,a0,a1); }
  else if(func_is(f,"asin")) { irvec_func(a0,a1,func_aget(f,0),n,x0,x1); irasin(z0,z1,a0,a1); }
  else if(func_is(f,"acos")) { irvec_func(a0,a1,func_aget(f,0),n,x0,x1); iracos(z0,z1,a0,a1); }
  else if(func_is(f,"atan")) { irvec_func(a0,a1,func_aget(f,0),n,x0,x1); iratan(z0,z1,a0,a1); }
  else if(func_is(f,"sinh")) { irvec_func(a0,a1,func_aget(f,0),n,x0,x1); irsinh(z0,z1,a0,a1); }
  else if(func_is(f,"cosh")) { irvec_func(a0,a1,func_aget(f,0),n,x0,x1); ircosh(z0,z1,a0,a1); }
  else if(func_is(f,"tanh")) { irvec_func(a0,a1,func_aget(f,0),n,x0,x1); irtanh(z0,z1,a0,a1); }
  else if(func_is(f,"asinh")){ irvec_func(a0,a1,func_aget(f,0),n,x0,x1); irasinh(z0,z1,a0,a1); }
  else if(func_is(f,"acosh")){ irvec_func(a0,a1,func_aget(f,0),n,x0,x1); iracosh(z0,z1,a0,a1); }
  else if(func_is(f,"atanh")){ irvec_func(a0,a1,func_aget(f,0),n,x0,x1); iratanh(z0,z1,a0,a1); }
  else if(func_is(f,"pow"))  {
    irvec_func(a0,a1,func_aget(f,0),n,x0,x1);
    irvec_func(b0,b1,func_aget(f,1),n,x0,x1);
    rpow(z0,a0,b0); ERROR_AT;
  }
  if(func_has_power(f))      { irpow_si(z0,z1,z0,z1,func_power(f)); }
  ircopy(y0,y1,z0,z1);
  RF(a0); RF(a1); RF(b0); RF(b1); RF(z0); RF(z1);
}

/**
 @brief ベクトル写像 [y0,y1]=f([x0,x1])
 */
void irvec_func_list(int m, rmulti **y0, rmulti **y1, func_t *f, int n, rmulti **x0, rmulti **x1)
{
  int i;
  for(i=0; i<m; i++){
    if(func_is_list(f) && func_is_list(f) && i<func_asize(f)){
      irvec_func(y0[i],y1[i],func_aget(f,i),n,x0,x1);
    }else{
      rset_nan(y0[i]); rset_nan(y1[i]);
    }
  }
}

/** @} */

//EOF
