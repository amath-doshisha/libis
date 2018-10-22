#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<stdarg.h>

#include"mt19937ar.h"
#include"is_macros.h"
#include"is_strings.h"
#include"is_dcomplex.h"
#include"is_rmulti.h"
#include"is_cmulti.h"
#include"is_irmulti.h"
#include"is_icmulti.h"
#include"is_svec.h"
#include"is_ivec.h"
#include"is_zvec.h"
#include"is_dvec.h"

/**
 @file  zvec.c
 @brief 倍精度複素数型dcomplexのベクトルに関する関数の定義
 */


/////////////////////////////////

#define FILE_NAME_LENGTH_MAX 100000

/////////////////////////////////

/** @name dcomplex型のベクトルの初期化 */
/** @{ */

/**
 @brief dcomplex型のメモリ割当
 */
dcomplex* zvec_allocate(int n)
{
  return (dcomplex*)malloc(sizeof(dcomplex)*n);
}

/**
 @brief dcomplex型のメモリ解放
 */
dcomplex* zvec_free(dcomplex *x)
{
  if(x==NULL) return NULL;
  free(x);
  return x=NULL;
}

/** @} */

/////////////////////////////////////////////////////

/** @name dcomplex型ベクトルの値の初期化 */
/** @{ */

/**
 @brief y=x
 */
void zvec_set_zvec(int n, dcomplex *y, dcomplex *x)
{
  int i;
  for(i=0; i<n; i++){ y[i]=x[i]; }
}

/**
 @brief y=x
 */
void zvec_set_dvec_dvec(int n, dcomplex *y, double *x_r, double *x_i)
{
  int i;
  for(i=0; i<n; i++){
    Z_R(y[i])=x_r[i];
    Z_I(y[i])=x_i[i];
  }
}

/**
 @brief x=nan(n,1)
 */
void zvec_set_nan(int n, dcomplex *x)
{
  int i;
  for(i=0; i<n; i++){ Z_SET(x[i],NAN,NAN); }
}

/**
 @brief x=inf(n,1)
 */
void zvec_set_inf(int n, dcomplex *x, int rsgn, int isgn)
{
  double zero=0;
  int i;
  for(i=0; i<n; i++){ Z_SET(x[i],rsgn/zero,isgn/zero); }
}

/**
 @brief x=zeros(n,1)
 */
void zvec_set_zeros(int n, dcomplex *x)
{
  int i;
  for(i=0; i<n; i++){ Z_SET(x[i],0,0); }
}

/**
 @brief x=ones(n,1)
 */
void zvec_set_ones(int n, dcomplex *x)
{
  int i;
  for(i=0; i<n; i++){ Z_SET(x[i],1,0); }
}

/**
 @brief x=ones(n,1)*a
 */
void zvec_set_all_z(int n, dcomplex *x, dcomplex a)
{
  int i;
  for(i=0; i<n; i++){ x[i]=a; }
}

/**
 @brief x=ones(n,1)*a
 */
void zvec_set_all_d(int n, dcomplex *x, double a)
{
  int i;
  for(i=0; i<n; i++){ Z_SET(x[i],a,0); }
}

/**
 @brief x=ones(n,1)*a
 */
void zvec_set_all_dd(int n, dcomplex *x, double real, double imag)
{
  int i;
  for(i=0; i<n; i++){ Z_SET(x[i],real,imag); }
}

/**
 @brief x=zeros(n,1); x[k]=1
 */
void zvec_set_unit(int n, dcomplex *x, int k)
{
  zvec_set_zeros(n,x);
  Z_SET(x[k],1,0);
}

/**
 @brief x[0]=0; x[1]=1; x[2]=2; ...; x[n-1]=n-1
 */
void zvec_set_grid(int n, dcomplex *x)
{
  int i;
  for(i=0; i<n; i++){ Z_SET(x[i],i,0); }
}

/**
 @brief x=rand(n,1)*a+b
 */
void zvec_set_rand(int n, dcomplex *x, double a, double b){
  double xr,xi;
  int i;
  for(i=0; i<n; i++){
    xr=genrand_int32();
    xi=genrand_int32();
    xr/=0xffffffff;
    xi/=0xffffffff;
    xr=xr*a+b;
    xi=xi*a+b;
    Z_SET(x[i],xr,xi);
  }
}

/** @} */

/////////////////////////////////////////////////////

/** @name dcomplex型ベクトルの型変換 */
/** @{ */

/**
 @brief char -> dcomplex
 */
void zvec_set_svec(int n, dcomplex *y, char **x)
{
  int i;
  cmulti *a=NULL;
  a=callocate();
  for(i=0; i<n; i++){
    cset_s(a,x[i]);
    y[i]=cget_z(a);
  }
  a=cmfree(a);
}

/**
 @brief dcomplex->char
 */
void zvec_get_svec(int n, char **y, dcomplex *x, char format, int digits)
{
  char f[1024];
  int i;
  sprintf(f,"%%-.%d%c%%+.%d%ci",digits,format,digits,format);
  for(i=0; i<n; i++){ y[i]=char_renew_sprintf(y[i],NULL,f,Z_R(x[i]),Z_I(x[i])); }
}

/**
 @brief double -> int
 */
void zvec_set_ivec(int n, dcomplex *y, int *x)
{
  int i;
  for(i=0; i<n; i++){ Z_SET(y[i],x[i],0); }
}

/**
 @brief dcomplex -> int
 */
void zvec_get_ivec(int n, int *y, dcomplex *x)
{
  int i;
  for(i=0; i<n; i++){ y[i]=Z_R(x[i]); }
}

/**
 @brief dcomplex -> int
 */
void ivec_set_zvec(int n, int *y, dcomplex *x)
{
  int i;
  for(i=0; i<n; i++){ y[i]=Z_R(x[i]); }
}

/**
 @brief double -> dcomplex
 */
void zvec_set_dvec(int n, dcomplex *y, double *x)
{
  int i;
  for(i=0; i<n; i++){ Z_SET(y[i],x[i],0); }
}

/**
 @brief dcomplex -> double
 */
void zvec_get_dvec(int n, double *y, dcomplex *x)
{
  int i;
  for(i=0; i<n; i++){ y[i]=Z_R(x[i]); }
}

/**
 @brief dcomplex -> double
 */
void dvec_set_zvec(int n, double *y, dcomplex *x)
{
  int i;
  for(i=0; i<n; i++){ y[i]=Z_R(x[i]); }
}

/** @} */

/////////////////////////////////////////////////////

/** @name dcomplex型ベクトルの要素の並び替え */
/** @{ */

/**
 @brief x <=> y
 */
void zvec_swap(int n, dcomplex *x, dcomplex *y)
{
  int i;
  dcomplex value;
  for(i=0; i<n; i++){
    value=x[i];
    x[i]=y[i];
    y[i]=value;
  }
}

/**
 @brief x[0]=x[n-1]; x[1]=x[n-2]; x[2]=x[n-3]; ...
 */
void zvec_reverse(int n, dcomplex *x)
{
  int i;
  dcomplex a;
  for(i=0; i<n/2; i++){
    a=x[i];
    x[i]=x[n-i-1];
    x[n-i-1]=a;
  }
}

/**
 @brief x[i] <=> x[j]
 */
void zvec_swap_at(dcomplex *x, int i, int j)
{
  dcomplex a;
  a=x[i];
  x[i]=x[j];
  x[j]=a;
}

/**
 @brief x[i] <=> x[I[i]]
 */
void zvec_swap_index(int n, dcomplex *x, int *I)
{
  dcomplex *y=NULL;
  y=zvec_allocate(n);
  zvec_copy_zvec_index(n,y,x,I);
  zvec_copy_zvec(n,x,y);
  y=zvec_free(y);
}

/**
 @brief sort x as x[0] <= x[1] <= x[2] <= ... <= x[n-1]
 @detail If I==NULL, then I is not ussed. If I!=NULL, then I is stored with sorted indexes.
 */
void zvec_sort(int n, dcomplex *x, int *I)
{
  if(I!=NULL){ ivec_set_grid(n,I); }
  zvec_quick_sort(n,x,I,0,n-1);
}

/**
 @brief Don't call this function directly!
 */
void zvec_quick_sort(int n, dcomplex *x, int *I, int left, int right)
{
  int i,last;
  if(left>=right) return;
  zvec_swap_at(x,left,(left+right)/2);
  if(I!=NULL) ivec_swap_at(I,left,(left+right)/2);
  last=left;
  for(i=left+1; i<=right; i++){
    if(lt_zz(x[i],x[left])){
      ++last;
      zvec_swap_at(x,last,i);
      if(I!=NULL) ivec_swap_at(I,last,i);
    }
  }
  zvec_swap_at(x,left,last);
  if(I!=NULL) ivec_swap_at(I,left,last);
  zvec_quick_sort(n,x,I,left,last-1);
  zvec_quick_sort(n,x,I,last+1,right);
}

/**
 @brief store the list of indexes, x is not destroyed
 */
void zvec_sort_index(int *I, int n, dcomplex *X)
{
  dcomplex *Y=NULL;
  Y=zvec_allocate(n);
  zvec_copy_zvec(n,Y,X);
  zvec_sort(n,Y,I);
  Y=zvec_free(Y);
}

/** @} */

/////////////////////////////////////////////////////

/** @name dcomplex型ベクトルの関数 y=f(x) */
/** @{ */

/**
 @brief y=x
 */
void zvec_copy_zvec(int n, dcomplex *y, dcomplex *x)
{
  int i;
  for(i=0; i<n; i++){ y[i]=x[i]; }
}

/**
 @brief Y[i]=X[I[i]]
 */
void zvec_copy_zvec_index(int n, dcomplex *Y, dcomplex *X, int *I)
{
  int i;
  for(i=0; i<n; i++) Y[i]=X[I[i]];
}

/**
 @brief y=-x
 */
void zvec_neg_zvec(int n, dcomplex *y, dcomplex *x)
{
  int i;
  for(i=0; i<n; i++){
    Z_R(y[i])=-Z_R(x[i]);
    Z_I(y[i])=-Z_I(x[i]);
  }
}

/**
 @brief y=real(x)
 */
void dvec_real_zvec(int n, double *y, dcomplex *x)
{
  int i;
  for(i=0; i<n; i++){ y[i]=Z_R(x[i]); }
}

/**
 @brief y=conj(x)
 */
void zvec_conj_zvec(int n, dcomplex *y, dcomplex *x)
{
  int i;
  for(i=0; i<n; i++){
    Z_R(y[i])= Z_R(x[i]);
    Z_I(y[i])=-Z_I(x[i]);
  }
}

/**
 @brief y=x/sqrt(x'*x)
 */
void zvec_normalize_zvec(int n, dcomplex *y, dcomplex *x)
{
  double norm=0;
  norm=dnorm2_zvec(n,x);
  zvec_mul_zvec_dscalar(n,y,x,(1.0/norm));
}

/**
 @brief y=x/sqrt(x'*x) where x[i]/|x[i]|=1, abs(x[i]) is max
 */
void zvec_normalize_sgn_zvec(int n, dcomplex *y, dcomplex *x)
{
  int k;
  double norm=0;
  dcomplex a;
  norm=dnorm2_zvec(n,x);
  dmax_abs_zvec_index(n,x,&k);
  a=znormalize_z(x[k]);
  Z_R(a)=Z_CONJ_R(a);
  Z_I(a)=Z_CONJ_I(a);
  Z_R(a)/=norm;
  Z_I(a)/=norm;
  zvec_mul_zvec_zscalar(n,y,x,a);
}

/**
 @brief y=sum(abs(x))
 */
double dnorm1_zvec(int n, dcomplex *x)
{
  int i;
  double value=0;
  for(i=0; i<n; i++){
    value+=Z_ABS(x[i]);
  }
  return value;
}

/**
 @brief y=sqrt(sum(abs(x).^2))
 */
double dnorm2_zvec(int n, dcomplex *x)
{
  int i;
  double value=0;
  for(i=0; i<n; i++){
    value+=Z_ABS2(x[i]);
  }
  value=sqrt(value);
  return value;
}

/**
 @brief y=max(abs(x))
 */
double dnorm_max_zvec(int n, dcomplex *x)
{
  int i;
  double value=0,foo;
  for(i=0; i<n; i++){
    foo=Z_ABS2(x[i]);
    if(foo>value) value=foo;
  }
  value=sqrt(value);
  return value;
}


/**
 @brief y=sum(abs(x).^2)
 */
double dsum_pow2_abs_zvec(int n, dcomplex *x)
{
  double value=0;
  int i;
  for(i=0; i<n; i++){
    value+=Z_ABS2(x[i]);
  }
  return value;
}


/**
 @brief y=max(abs(x))
 */
double dmax_abs_zvec(int n, dcomplex *x)
{
  int i;
  double value=0,foo;
  for(i=0; i<n; i++){
    foo=Z_ABS2(x[i]);
    if(foo>value){
      value=foo;
    }
  }
  value=sqrt(value);
  return value;
}

/**
 @brief y=min(abs(x))
 */
double dmin_abs_zvec(int n, dcomplex *x)
{
  int i;
  double value,foo;
  value=Z_ABS2(x[0]);
  for(i=1; i<n; i++){
    foo=Z_ABS2(x[i]);
    if(foo<value){
      value=foo;
    }
  }
  value=sqrt(value);
  return value;
}

/**
 @brief y=max(abs(x)), (*I)=k where k=max{ k | abs(x[k]) } unless I==NULL
 */
double dmax_abs_zvec_index(int n, dcomplex *x, int *I)
{
  int i;
  double value=0,foo;
  value=Z_ABS2(x[0]);
  if(I!=NULL) (*I)=0;
  for(i=1; i<n; i++){
    foo=Z_ABS2(x[i]);
    if(foo>value){
      value=foo;
      if(I!=NULL) (*I)=i;
    }
  }
  value=sqrt(value);
  return value;
}

/**
 @brief y=min(abs(x)), (*I)=k where k=min{ k | abs(x[k]) } unless I==NULL
 */
double dmin_abs_zvec_index(int n, dcomplex *x, int *I)
{
  int i;
  double value,foo;
  value=Z_ABS2(x[0]);
  if(I!=NULL) (*I)=0;
  for(i=1; i<n; i++){
    foo=Z_ABS2(x[i]);
    if(foo<value){
      value=foo;
      if(I!=NULL) (*I)=i;
    }
  }
  value=sqrt(value);
  return value;
}

/**
 @brief y=max(x)
 */
dcomplex zmax_zvec(int n, dcomplex *x)
{
  int i;
  dcomplex value;
  value=x[0];
  for(i=1; i<n; i++){
    if(gt_zz(x[i],value)) value=x[i];
  }
  return value;
}

/**
 @brief y=min(x)
 */
dcomplex zmin_zvec(int n, dcomplex *x)
{
  int i;
  dcomplex value=x[0];
  for(i=1; i<n; i++){
    if(lt_zz(x[i],value)) value=x[i];
  }
  return value;
}

/**
 @brief y=sum(x)
 */
dcomplex zsum_zvec(int n, dcomplex *x)
{
  dcomplex value;
  int i;
  Z_SET(value,0,0);
  for(i=0; i<n; i++){
    Z_ADD(value,x[i]);
  }
  return value;
}

/**
 @brief y=sum(x)/n
 */
dcomplex zvec_avarage(int n, dcomplex *x)
{
  dcomplex value;
  value=zsum_zvec(n,x);
  Z_R(value)/=n;
  Z_I(value)/=n;
  return value;
}

/**
 @brief dcomplex型ｍのベクトルの要素の平方根 y=sqrt(x)
 */
void zvec_sqrt_zvec(int n, dcomplex *y, dcomplex *x)
{
  int i;
  for(i=0; i<n; i++){ y[i]=zsqrt_z(x[i]); }
}

/**
 @brief dcomplex型ｍのベクトルの要素の平方根 y=sqrt(x)
 */
void zvec_sqrt_dvec(int n, dcomplex *y, double *x)
{
  int i;
  for(i=0; i<n; i++){ y[i]=zsqrt_d(x[i]); }
}

/**
 @brief dcomplex型ｍのベクトルの要素の対数 y=log(x)
 */
void zvec_log_zvec(int n, dcomplex *y, dcomplex *x)
{
  int i;
  for(i=0; i<n; i++){ y[i]=zlog_z(x[i]); }
}

/**
 @brief dcomplex型ｍのベクトルの要素の対数 y=log(x)
 */
void zvec_log_dvec(int n, dcomplex *y, double *x)
{
  int i;
  for(i=0; i<n; i++){ y[i]=zlog_d(x[i]); }
}

/**
 @brief dcomplex型ｍのベクトルの要素の対数 y=log10(x)
 */
void zvec_log10_zvec(int n, dcomplex *y, dcomplex *x)
{
  int i;
  for(i=0; i<n; i++){ y[i]=zlog10_z(x[i]); }
}

/**
 @brief dcomplex型ｍのベクトルの要素の対数 y=log10(x)
 */
void zvec_log10_dvec(int n, dcomplex *y, double *x)
{
  int i;
  for(i=0; i<n; i++){ y[i]=zlog10_d(x[i]); }
}

/** @} */

/////////////////////////////////////////////////////

/** @name dcomplex型ベクトルの関数 z=f(x,y) */
/** @{ */

/**
 @brief z=x+y
 */
void zvec_add_zvec_zvec(int n, dcomplex *z, dcomplex *x, dcomplex *y)
{
  int i;
  for(i=0; i<n; i++){ z[i]=zadd_zz(x[i],y[i]); }
}

/**
 @brief z=x+y
 */
void zvec_add_zvec_dvec(int n, dcomplex *z, dcomplex *x, double *y)
{
  int i;
  for(i=0; i<n; i++){ z[i]=zadd_zd(x[i],y[i]); }
}

/**
 @brief z=x+y
 */
void zvec_add_dvec_zvec(int n, dcomplex *z, double *x, dcomplex *y)
{
  int i;
  for(i=0; i<n; i++){ z[i]=zadd_dz(x[i],y[i]); }
}

/**
 @brief z=x+y
 */
void zvec_add_zvec_zscalar(int n, dcomplex *z, dcomplex *x, dcomplex y)
{
  int i;
  for(i=0; i<n; i++){ z[i]=zadd_zz(x[i],y); }
}

/**
 @brief z=x+y
 */
void zvec_add_zvec_dscalar(int n, dcomplex *z, dcomplex *x, double y)
{
  int i;
  for(i=0; i<n; i++){ z[i]=zadd_zd(x[i],y); }
}

/**
 @brief z=x+y
 */
void zvec_add_dvec_zscalar(int n, dcomplex *z, double *x, dcomplex y)
{
  int i;
  for(i=0; i<n; i++){ z[i]=zadd_dz(x[i],y); }
}

/**
 @brief z=x+y
 */
void zvec_add_zscalar_zvec(int n, dcomplex *z, dcomplex x, dcomplex *y)
{
  int i;
  for(i=0; i<n; i++){ z[i]=zadd_zz(x,y[i]); }
}

/**
 @brief z=x+y
 */
void zvec_add_zscalar_dvec(int n, dcomplex *z, dcomplex x, double *y)
{
  int i;
  for(i=0; i<n; i++){ z[i]=zadd_zd(x,y[i]); }
}

/**
 @brief z=x+y
 */
void zvec_add_dscalar_zvec(int n, dcomplex *z, double x, dcomplex *y)
{
  int i;
  for(i=0; i<n; i++){ z[i]=zadd_dz(x,y[i]); }
}

/**
 @brief z=x-y
 */
void zvec_sub_zvec_zvec(int n, dcomplex *z, dcomplex *x, dcomplex *y)
{
  int i;
  for(i=0; i<n; i++){ z[i]=zsub_zz(x[i],y[i]); }
}

/**
 @brief z=x-y
 */
void zvec_sub_zvec_dvec(int n, dcomplex *z, dcomplex *x, double *y)
{
  int i;
  for(i=0; i<n; i++){ z[i]=zsub_zd(x[i],y[i]); }
}

/**
 @brief z=x-y
 */
void zvec_sub_dvec_zvec(int n, dcomplex *z, double *x, dcomplex *y)
{
  int i;
  for(i=0; i<n; i++){ z[i]=zsub_dz(x[i],y[i]); }
}

/**
 @brief z=x-y
 */
void zvec_sub_zvec_zscalar(int n, dcomplex *z, dcomplex *x, dcomplex y)
{
  int i;
  for(i=0; i<n; i++){ z[i]=zsub_zz(x[i],y); }
}

/**
 @brief z=x-y
 */
void zvec_sub_zvec_dscalar(int n, dcomplex *z, dcomplex *x, double y)
{
  int i;
  for(i=0; i<n; i++){ z[i]=zsub_zd(x[i],y); }
}

/**
 @brief z=x-y
 */
void zvec_sub_dvec_zscalar(int n, dcomplex *z, double *x, dcomplex y)
{
  int i;
  for(i=0; i<n; i++){ z[i]=zsub_dz(x[i],y); }
}

/**
 @brief z=x-y
 */
void zvec_sub_zscalar_zvec(int n, dcomplex *z, dcomplex x, dcomplex *y)
{
  int i;
  for(i=0; i<n; i++){ z[i]=zsub_zz(x,y[i]); }
}

/**
 @brief z=x-y
 */
void zvec_sub_zscalar_dvec(int n, dcomplex *z, dcomplex x, double *y)
{
  int i;
  for(i=0; i<n; i++){ z[i]=zsub_zd(x,y[i]); }
}

/**
 @brief z=x-y
 */
void zvec_sub_dscalar_zvec(int n, dcomplex *z, double x, dcomplex *y)
{
  int i;
  for(i=0; i<n; i++){ z[i]=zsub_dz(x,y[i]); }
}

/**
 @brief z=x*y
 */
void zvec_mul_zvec_zvec(int n, dcomplex *z, dcomplex *x, dcomplex *y)
{
  int i;
  for(i=0; i<n; i++){ z[i]=zmul_zz(x[i],y[i]); }
}

/**
 @brief z=x*y
 */
void zvec_mul_zvec_dvec(int n, dcomplex *z, dcomplex *x, double *y)
{
  int i;
  for(i=0; i<n; i++){ z[i]=zmul_zd(x[i],y[i]); }
}

/**
 @brief z=x*y
 */
void zvec_mul_dvec_zvec(int n, dcomplex *z, double *x, dcomplex *y)
{
  int i;
  for(i=0; i<n; i++){ z[i]=zmul_dz(x[i],y[i]); }
}

/**
 @brief z=x*y
 */
void zvec_mul_zvec_zscalar(int n, dcomplex *z, dcomplex *x, dcomplex y)
{
  int i;
  for(i=0; i<n; i++){ z[i]=zmul_zz(x[i],y); }
}

/**
 @brief z=x*y
 */
void zvec_mul_zvec_dscalar(int n, dcomplex *z, dcomplex *x, double y)
{
  int i;
  for(i=0; i<n; i++){ z[i]=zmul_zd(x[i],y); }
}

/**
 @brief z=x*y
 */
void zvec_mul_dvec_zscalar(int n, dcomplex *z, double *x, dcomplex y)
{
  int i;
  for(i=0; i<n; i++){ z[i]=zmul_dz(x[i],y); }
}

/**
 @brief z=x*y
 */
void zvec_mul_zscalar_zvec(int n, dcomplex *z, dcomplex x, dcomplex *y)
{
  int i;
  for(i=0; i<n; i++){ z[i]=zmul_zz(x,y[i]); }
}

/**
 @brief z=x*y
 */
void zvec_mul_zscalar_dvec(int n, dcomplex *z, dcomplex x, double *y)
{
  int i;
  for(i=0; i<n; i++){ z[i]=zmul_zd(x,y[i]); }
}

/**
 @brief z=x*y
 */
void zvec_mul_dscalar_zvec(int n, dcomplex *z, double x, dcomplex *y)
{
  int i;
  for(i=0; i<n; i++){ z[i]=zmul_dz(x,y[i]); }
}

/**
 @brief z=x/y
 */
void zvec_div_zvec_zvec(int n, dcomplex *z, dcomplex *x, dcomplex *y)
{
  int i;
  for(i=0; i<n; i++){ z[i]=zdiv_zz(x[i],y[i]); }
}

/**
 @brief z=x/y
 */
void zvec_div_zvec_dvec(int n, dcomplex *z, dcomplex *x, double *y)
{
  int i;
  for(i=0; i<n; i++){ z[i]=zdiv_zd(x[i],y[i]); }
}

/**
 @brief z=x/y
 */
void zvec_div_dvec_zvec(int n, dcomplex *z, double *x, dcomplex *y)
{
  int i;
  for(i=0; i<n; i++){ z[i]=zdiv_dz(x[i],y[i]); }
}

/**
 @brief z=x/y
 */
void zvec_div_zvec_zscalar(int n, dcomplex *z, dcomplex *x, dcomplex y)
{
  int i;
  for(i=0; i<n; i++){ z[i]=zdiv_zz(x[i],y); }
}

/**
 @brief z=x/y
 */
void zvec_div_zvec_dscalar(int n, dcomplex *z, dcomplex *x, double y)
{
  int i;
  for(i=0; i<n; i++){ z[i]=zdiv_zd(x[i],y); }
}

/**
 @brief z=x/y
 */
void zvec_div_dvec_zscalar(int n, dcomplex *z, double *x, dcomplex y)
{
  int i;
  for(i=0; i<n; i++){ z[i]=zdiv_dz(x[i],y); }
}

/**
 @brief z=x/y
 */
void zvec_div_zscalar_zvec(int n, dcomplex *z, dcomplex x, dcomplex *y)
{
  int i;
  for(i=0; i<n; i++){ z[i]=zdiv_zz(x,y[i]); }
}

/**
 @brief z=x/y
 */
void zvec_div_zscalar_dvec(int n, dcomplex *z, dcomplex x, double *y)
{
  int i;
  for(i=0; i<n; i++){ z[i]=zdiv_zd(x,y[i]); }
}

/**
 @brief z=x/y
 */
void zvec_div_dscalar_zvec(int n, dcomplex *z, double x, dcomplex *y)
{
  int i;
  for(i=0; i<n; i++){ z[i]=zdiv_dz(x,y[i]); }
}

/**
 @brief z=x'*y
 */
dcomplex zinnprod_zvec_zvec(int n, dcomplex *x, dcomplex *y)
{
  int i;
  dcomplex value;
  Z_SET(value,0,0);
  for(i=0; i<n; i++){
    Z_ADD_DOT(value,x[i],y[i]);
  }
  return value;
}

/**
 @brief y=x^n
 */
dcomplex zpow_i(dcomplex x, int n)
{
  int i;
  dcomplex y={0,0},a;
  if(n==0){
    Z_SET(y,1,0);
  }else if(n>0){
    Z_LET(y,x);
    for(i=1 ;i<n ;i++){
      Z_SET_TIMES(a,y,x); Z_LET(y,a);
    }
  }else if(n<0){
    Z_LET(y,x);
    for(i=1; i<(-n); i++){
      Z_SET_TIMES(a,y,x); Z_LET(y,a);
    }
    Z_SET_INV(a,y); Z_LET(y,a);
  }else{
    // error
  }
  return y;
}

/**
 @brief z=abs(x'*y)/\sqrt(x'*x)/sqrt(y'*y)
 */
double ddcos_abs_zvec_zvec(int n, dcomplex *x, dcomplex *y)
{
  double norm_x,norm_y,dcos;
  dcomplex dot;
  norm_x=dnorm2_zvec(n,x);
  norm_y=dnorm2_zvec(n,y);
  dot=zinnprod_zvec_zvec(n,x,y);
  dcos=Z_ABS(dot);
  dcos/=norm_x*norm_y;
  return dcos;
}

/**
 @brief z=acos(abs(x'*y)/\sqrt(x'*x)/sqrt(y'*y))
 */
double dangle_zvec_zvec(int n, dcomplex *x, dcomplex *y)
{
  double theta;
  theta=ddcos_abs_zvec_zvec(n,x,y);
  if(theta>1) theta=1;
  theta=acos(theta);
  return theta;
}

/**
 @brief z=(180/PI)*acos(abs(x'*y)/\sqrt(x'*x)/sqrt(y'*y))
 */
double ddeg_angle_zvec_zvec(int n, dcomplex *x, dcomplex *y)
{
  return dangle_zvec_zvec(n,x,y)*(180.0/M_PI);
}

/**
 @brief z=|x[i]-y[i]| (0<=i<n)
 */
double dsum_abs_sub_zvec_zvec(int n, dcomplex *x, dcomplex *y)
{
  int i;
  double value;
  value=0.0;
  for(i=0; i<n; i++) value+=Z_DIST(x[i],y[i]);
  return value;
}

/**
 @brief z=max(|x-y|)
 */
double dmax_abs_sub_zvec_zvec(int n, dcomplex *x, dcomplex *y)
{
  int i;
  double value,a;
  value=Z_DIST2(x[0],y[0]);
  for(i=1; i<n; i++){
    a=Z_DIST2(x[i],y[i]);
    if(a>value) value=a;
  }
  value=sqrt(value);
  return value;
}

/**
 @brief z=abs(x-y)
 */
void dvec_abs_sub_zvec_zvec(int n, double *z, dcomplex *x, dcomplex *y)
{
  int i;
  for(i=0; i<n; i++){ z[i]=Z_DIST(x[i],y[i]); }
}

/**
 @brief z=abs(x-y)
 */
void dvec_abs_sub_zvec_zscalar(int n, double *z, dcomplex *x, dcomplex y)
{
  int i;
  for(i=0; i<n; i++){ z[i]=Z_DIST(x[i],y); }
}


/** @} */

/////////////////////////////////////////////////////

/** @name dcomplex型ベクトルの関数 y=y+f(x) */
/** @{ */

/**
 @brief y=y+a*x
 */
void zvec_add_scaled(int n, dcomplex *y, dcomplex a, dcomplex *x)
{
  int i;
  for(i=0; i<n; i++){
    Z_ADD_TIMES(y[i],a,x[i]);
  }
}

/**
 @biref y=y-a*x
 */
void zvec_sub_scaled(int n, dcomplex *y, dcomplex a, dcomplex *x)
{
  int i;
  for(i=0; i<n; i++){
    Z_SUB_TIMES(y[i],a,x[i]);
  }
}

/**
 @brief y=y-(x'*y)*x where x'*x=1
 */
void zvec_orthogonalize(int n, dcomplex *y, dcomplex *x)
{
  dcomplex dot;
  dot=zinnprod_zvec_zvec(n,x,y);
  Z_R(dot)=-Z_R(dot);
  Z_I(dot)=-Z_I(dot);
  zvec_add_scaled(n,y,dot,x);
}


/** @} */

/////////////////////////////////////////////////////

/** @name I/O */
/** @{ */


/**
 @brief print
 */
void zvec_print(int n, dcomplex *x, char *name, char format, int digits)
{
  char **s=NULL;
  if(x==NULL){
    if(name!=NULL){ printf("%s=NULL\n",name); }
    else          { printf("NULL\n"); }
    return;
  }
  s=svec_allocate(n);
  zvec_get_svec(n,s,x,format,digits);
  svec_print(n,s,name);
  s=svec_free(n,s);
  /*
   int i;
   char format[128];
   if(STR_EQ(f,"e")){ sprintf(format,"(%%%d.%d%s, %%%d.%d%s)\n",digits+7,digits,f,digits+7,digits,f); }
   else             { sprintf(format,"%%-.%d%s%%+.%d%si\n",digits,f,digits,f); }
   if(x==NULL){
   if(name!=NULL){ printf("%s=NULL\n",name); }
   else          { printf("NULL\n"); }
   return;
   }
   if(name!=NULL){ printf("%s\n",name); }
   if(x==NULL) return;
   for(i=0; i<n; i++){
   printf(format,Z_R(x[i]),Z_I(x[i]));
   }
   */
}

/**
 @brief save data as text format
 */
void zvec_save(int n, dcomplex *x, int offset, char* fmt, ...)
{
  int i;
  char fname[FILE_NAME_LENGTH_MAX+1];
  va_list argp;
  FILE *fid;
  va_start(argp,fmt); vsprintf(fname,fmt,argp);   // file name
  if((fid=fopen(fname,"w"))==0){ ERROR_AT; printf("Cant't open file: %s\n",fname); exit(0); } // open
  for(i=0; i<n; i++){ fprintf(fid,"%d\t%+.20e\t%+.20e\n",i+offset,Z_R(x[i]),Z_I(x[i])); } // write
  fclose(fid); // close
}

/**
 @brief save data of complex plane as text file
 */
void zvec_save_cplane(int n, dcomplex *x, char* fmt, ...)
{
  int i;
  char fname[FILE_NAME_LENGTH_MAX+1];
  va_list argp;
  FILE *fid;
  va_start(argp,fmt); vsprintf(fname,fmt,argp);   // file name
  if((fid=fopen(fname,"w"))==0){ ERROR_AT; printf("Cant't open file: %s\n",fname); exit(0); } // open
  for(i=0; i<n; i++){ fprintf(fid,"%+.20e\t%+.20e\n",Z_R(x[i]),Z_I(x[i])); } // write
  fclose(fid); // close
}

/**
 @brief save data of complex plane with label as text format
 */
void zvec_save_cplane_label(int n, dcomplex *x, int offset, int *label, char* fmt, ...)
{
  int i;
  char fname[FILE_NAME_LENGTH_MAX+1];
  va_list argp;
  FILE *fid;
  va_start(argp,fmt); vsprintf(fname,fmt,argp); // file name
  if((fid=fopen(fname,"w"))==0){ ERROR_AT; printf("Cant't open file: %s\n",fname); exit(0); } // open
  for(i=0; i<n; i++){ fprintf(fid,"%+.20e\t%+.20e\t%d\n",Z_R(x[i]),Z_I(x[i]),label[i]+offset); } // write
  fclose(fid); // close
}

/**
 @brief save data as binary format
 */
void zvec_bin_save(int n, dcomplex *x, char* fmt, ...)
{
  char fname[FILE_NAME_LENGTH_MAX+1];
  va_list argp;
  FILE *fid;
  va_start(argp,fmt); vsprintf(fname,fmt,argp); // file name
  if((fid=fopen(fname,"w"))==0){ ERROR_AT; printf("Cant't open file: %s\n",fname); exit(0); } // open
  fwrite("zvec",sizeof(char),strlen("zvec"),fid); // write header
  fwrite(&n,sizeof(int),1,fid); // write size
  fwrite(x,sizeof(dcomplex),n,fid); // write data
  fclose(fid); // close
}

/**
 @brief load binary formated data
 */
dcomplex *zvec_bin_load(int *n, char* fmt, ...)
{
  int l;
  size_t k;
  dcomplex *zx=NULL;
  double *dx=NULL;
  char fname[FILE_NAME_LENGTH_MAX+1],*buf=NULL;
  va_list argp;
  FILE *fid;
  // file name
  va_start(argp,fmt); vsprintf(fname,fmt,argp);
  // open file
  if((fid=fopen(fname,"r"))==0){ fclose(fid); zx=NULL; (*n)=0; }
  else{
    // read header
    l=strlen("zvec");
    buf=malloc(sizeof(char)*(l+1));
    k=fread(buf,sizeof(char),l,fid);
    if(k==(size_t)l && strncmp(buf,"zvec",l)==0){ /* zvec */
      // read size
      k=fread(n,sizeof(int),1,fid);
      if(k!=1 || (*n)<=0){ ERROR_AT; printf("Failed to load the size from the file '%s'.\n",fname); exit(0); }
      // allocate
      zx=zvec_allocate((*n));
      // read data
      k=fread(zx,sizeof(dcomplex),(*n),fid);
      if(k!=(size_t)(*n)){ ERROR_AT; printf("Failed to load the data from the file '%s'.\n",fname); exit(0); }
      // close
      fclose(fid);
    }else if(k==(size_t)l && strncmp(buf,"dvec",l)==0){ /* dvec */
      fclose(fid);               // close
      dx=dvec_bin_load(n,fname); // load
      zx=zvec_allocate((*n));    // allocate
      zvec_set_dvec((*n),zx,dx);    // copy
      dx=dvec_free(dx);          // free
    }else{ fclose(fid); zx=NULL; (*n)=0; }
  }
  // done
  free(buf);
  return zx;
}

/** @} */

/////////////////////////////////////////////////////

/** @name 作業中：cmat.cに移動する */
/** @{ */

// y=A*x
void zvec_mul_zmat_zvec(int m, int n, dcomplex *y, dcomplex *A, int LDA, dcomplex *x)
{
  int i,j;
  for(i=0; i<m; i++){
    Z_SET(y[i],0,0);
    for(j=0; j<n; j++){
      Z_ADD_TIMES(y[i],MAT(A,i,j,LDA),x[j]);
    }
  }
}

// y=y+A*x
void zvec_add_mul_zmat_zvec(int m, int n, dcomplex *y, dcomplex *A, int LDA, dcomplex *x)
{
  int i,j;
  for(i=0; i<m; i++){
    for(j=0; j<n; j++){
      Z_ADD_TIMES(y[i],MAT(A,i,j,LDA),x[j]);
    }
  }
}

// y=y-A*x
void zvec_sub_zmat_zvec(int m, int n, dcomplex *y, dcomplex *A, int LDA, dcomplex *x)
{
  int i,j;
  for(i=0; i<m; i++){
    for(j=0; j<n; j++){
      Z_SUB_TIMES(y[i],MAT(A,i,j,LDA),x[j]);
    }
  }
}

// y=A^T*x
void zvec_mul_zmat_t_zvec(int m, int n, dcomplex *y, dcomplex *A, int LDA, dcomplex *x)
{
  int i,j;
  for(j=0; j<n; j++){
    Z_SET(y[j],0,0);
    for(i=0; i<m; i++){
      Z_ADD_TIMES(y[j],MAT(A,i,j,LDA),x[i]);
    }
  }
}

// y=y+A^T*x
void zvec_add_mul_zmat_t_zvec(int m, int n, dcomplex *y, dcomplex *A, int LDA, dcomplex *x)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      Z_ADD_TIMES(y[j],MAT(A,i,j,LDA),x[i]);
    }
  }
}

// y=y-A^T*x
void zvec_sub_zmat_t_zvec(int m, int n, dcomplex *y, dcomplex *A, int LDA, dcomplex *x)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      Z_SUB_TIMES(y[j],MAT(A,i,j,LDA),x[i]);
    }
  }
}

// y=A'*x
void zvec_mul_zmat_ct_zvec(int m, int n, dcomplex *y, dcomplex *A, int LDA, dcomplex *x)
{
  int i,j;
  for(j=0; j<n; j++){
    Z_SET(y[j],0,0);
    for(i=0; i<m; i++){
      Z_ADD_DOT(y[j],MAT(A,i,j,LDA),x[i]);
    }
  }
}


// y=y+A'*x
void zvec_add_mul_zmat_ct_zvec(int m, int n, dcomplex *y, dcomplex *A, int LDA, dcomplex *x)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      Z_ADD_DOT(y[j],MAT(A,i,j,LDA),x[i]);
    }
  }
}

// y=y-A'*x
void zvec_sub_zmat_ct_zvec(int m, int n, dcomplex *y, dcomplex *A, int LDA, dcomplex *x)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      Z_SUB_DOT(y[j],MAT(A,i,j,LDA),x[i]);
    }
  }
}

/** @} */

//EOF
