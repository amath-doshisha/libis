#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<stdarg.h>
#include"is_macros.h"
#include"is_dcomplex.h"

/**
 @file   dcomplex.c
 @brief  倍精度複素数構造体 dcomplex に関する関数の定義
 */


/////////////////////////////////////////////

/** @name dcomplex型の設定に関する関数 */
/** @{ */

/**
 @brief z=r*exp(i*theta)
 */
void zget_polar(double *r, double *theta, dcomplex z)
{
  (*r)=dabs_z(z);
  (*theta)=atan2(Z_I(z),Z_R(z));
}

/**
 @brief z=r*exp(i*theta)
 */
dcomplex zset_polar(double r, double theta)
{
  dcomplex z;
  Z_R(z)=r*cos(theta);
  Z_I(z)=r*sin(theta);
  return z;
}

/** @} */

/////////////////////////////////////////////

/** @name dcomplex型に関する関数 y=f(x) */
/** @{ */


/**
 @brief y=-x
 */
dcomplex zneg_z(dcomplex x)
{
  dcomplex z;
  Z_R(z)=-Z_R(x);
  Z_I(z)=-Z_I(x);
  return z;
}

/**
 @brief y=-x
 */
dcomplex zneg_d(double x)
{
  dcomplex z;
  Z_R(z)=-x;
  Z_I(z)=0;
  return z;
}

/**
 @brief y=1/x
 */
dcomplex zinv_z(dcomplex x)
{
  dcomplex z;
  double den;
  den=Z_R(x)*Z_R(x)+Z_I(x)*Z_I(x);
  Z_R(z)= Z_R(x)/den;
  Z_I(z)=-Z_I(x)/den;
  return z;
}

/**
 @brief y=1/x
 */
dcomplex zinv_d(double x)
{
  dcomplex z;
  Z_R(z)=1/x;
  Z_I(z)=0;
  return z;
}

/**
 @brief y=abs(x)
 */
double dabs_z(dcomplex x)
{
  return sqrt(Z_R(x)*Z_R(x)+Z_I(x)*Z_I(x));
}

/**
 @brief y=abs(x)
 */
double dabs_d(double x)
{
  return fabs(x);
}

/**
 @brief y=abs(x)^2
 */
double dabs2_z(dcomplex x)
{
  return Z_R(x)*Z_R(x)+Z_I(x)*Z_I(x);
}

/**
 @brief y=abs(x)^2
 */
double dabs2_d(double x)
{
  return x*x;
}

/**
 @brief y=x/|x|
 */
dcomplex znormalize_z(dcomplex x)
{
  double a;
  dcomplex z;
  if(Z_R(x)==0 && Z_I(x)==0){
    Z_R(z)=0;
    Z_I(z)=0;
  }else{
    a=sqrt(Z_R(x)*Z_R(x)+Z_I(x)*Z_I(x));
    Z_R(z)=Z_R(x)/a;
    Z_I(z)=Z_I(x)/a;
  }
  return z;
}

/**
 @brief y=x/|x|
 */
dcomplex dnormalize_d(double x)
{
  dcomplex z;
  if(x>0){
    Z_R(z)=1;
    Z_I(z)=0;
  }else if(x<0){
    Z_R(z)=-1;
    Z_I(z)=0;
  }else{
    Z_R(z)=0;
    Z_I(z)=0;
  }
  return z;
}

/**
 @brief y=sqrt(x)
 */
double dsqrt_d(double x)
{
  return sqrt(x);
}


/**
 @brief y=sqrt(x)
 */
dcomplex zsqrt_d(double x)
{
  dcomplex y;
  if(x>=0){ Z_R(y)=sqrt(x); Z_I(y)=0; }
  else    { Z_R(y)=0;       Z_I(y)=sqrt(-x); }
  return y;
}

/**
 @brief y=sqrt(x)
 */
dcomplex zsqrt_z(dcomplex x)
{
  return zpow_zd(x,0.5);
}

/**
 @brief y=log(x)
 */
double dlog_d(double x)
{
  return log(x);
}

/**
 @brief y=log(x)
 */
dcomplex zlog_d(double x)
{
  dcomplex y;
  if(x==0)    { Z_R(y)=log(0);  Z_I(y)=log(0); }
  else if(x>0){ Z_R(y)=log(x);  Z_I(y)=0; }
  else        { Z_R(y)=log(-x); Z_I(y)=M_PI; }
  return y;
}

/**
 @brief y=log(x)
 */
dcomplex zlog_z(dcomplex x)
{
  dcomplex y;
  double r,theta;
  zget_polar(&r,&theta,x); // x=r*exp(i*theta)
  Z_R(y)=log(r);           // y=log(r)+i*theta
  Z_I(y)=theta;            // y=log(r)+i*theta
  return y;
}

/**
 @brief y=log10(x)
 */
double dlog10_d(double x)
{
  return log10(x);
}

/**
 @brief y=log10(x)
 */
dcomplex zlog10_d(double x)
{
  dcomplex y;
  if(x==0)    { Z_R(y)=log10(0);  Z_I(y)=log10(0); }
  else if(x>0){ Z_R(y)=log10(x);  Z_I(y)=0; }
  else        { Z_R(y)=log10(-x); Z_I(y)=M_PI/log(10); }
  return y;
}

/**
 @brief y=log10(x)
 */
dcomplex zlog10_z(dcomplex x)
{
  dcomplex y;
  double r,theta;
  zget_polar(&r,&theta,x); // x=r*exp(i*theta)
  Z_R(y)=log10(r);         // y=(log(r)+i*theta)/log(10)
  Z_I(y)=theta/log(10);    // y=(log(r)+i*theta)/log(10)
  return y;
}

/** @} */

/////////////////////////////////////////////

/** @name dcomplex型に関する関数 z=f(x,y) */
/** @{ */

/**
 @brief z=x+y
 */
dcomplex zadd_zz(dcomplex x, dcomplex y)
{
  dcomplex z;
  Z_R(z)=Z_R(x)+Z_R(y);
  Z_I(z)=Z_I(x)+Z_I(y);
  return z;
}

/**
 @brief z=x+y
 */
dcomplex zadd_zd(dcomplex x, double y)
{
  dcomplex z;
  Z_R(z)=Z_R(x)+y;
  Z_I(z)=Z_I(x);
  return z;
}

/**
 @brief z=x+y
 */
dcomplex zadd_dz(double x, dcomplex y)
{
  dcomplex z;
  Z_R(z)=x+Z_R(y);
  Z_I(z)=  Z_I(y);
  return z;
}

/**
 @brief z=x-y
 */
dcomplex zsub_zz(dcomplex x, dcomplex y)
{
  dcomplex z;
  Z_R(z)=Z_R(x)-Z_R(y);
  Z_I(z)=Z_I(x)-Z_I(y);
  return z;
}

/**
 @brief z=x-y
 */
dcomplex zsub_zd(dcomplex x, double y)
{
  dcomplex z;
  Z_R(z)=Z_R(x)-y;
  Z_I(z)=Z_I(x);
  return z;
}

/**
 @brief z=x-y
 */
dcomplex zsub_dz(double x, dcomplex y)
{
  dcomplex z;
  Z_R(z)=x-Z_R(y);
  Z_I(z)= -Z_I(y);
  return z;
}

/**
 @brief z=x*y
 */
dcomplex zmul_zz(dcomplex x, dcomplex y)
{
  dcomplex z;
  Z_R(z)=Z_R(x)*Z_R(y)-Z_I(x)*Z_I(y);
  Z_I(z)=Z_R(x)*Z_I(y)+Z_I(x)*Z_R(y);
  return z;
}

/**
 @brief z=x*y
 */
dcomplex zmul_zd(dcomplex x, double y)
{
  dcomplex z;
  Z_R(z)=Z_R(x)*y;
  Z_I(z)=Z_I(x)*y;
  return z;
}

/**
 @brief z=x*y
 */
dcomplex zmul_dz(double x, dcomplex y)
{
  dcomplex z;
  Z_R(z)=x*Z_R(y);
  Z_I(z)=x*Z_I(y);
  return z;
}

/**
 @brief z=conj(x)*y
 */
dcomplex zdot_zz(dcomplex x, dcomplex y)
{
  dcomplex z;
  Z_R(z)=Z_R(x)*Z_R(y)+Z_I(x)*Z_I(y);
  Z_I(z)=Z_R(x)*Z_I(y)-Z_I(x)*Z_R(y);
  return z;
}

/**
 @brief z=conj(x)*y
 */
dcomplex zdot_zd(dcomplex x, double y)
{
  dcomplex z;
  Z_R(z)= Z_R(x)*y;
  Z_I(z)=-Z_I(x)*y;
  return z;
}

/**
 @brief z=conj(x)*y
 */
dcomplex zdot_dz(double x, dcomplex y)
{
  dcomplex z;
  Z_R(z)=x*Z_R(y);
  Z_I(z)=x*Z_I(y);
  return z;
}

/**
 @brief z=x/y
 */
dcomplex zdiv_zz(dcomplex x, dcomplex y)
{
  dcomplex z;
  double den;
  den=Z_R(y)*Z_R(y)+Z_I(y)*Z_I(y);
  Z_R(z)=( Z_R(x)*Z_R(y)+Z_I(x)*Z_I(y))/den; // xr*yr+xi*yi
  Z_I(z)=(-Z_R(x)*Z_I(y)+Z_I(x)*Z_R(y))/den; // xi*yr-xr*yi
  return z;
}

/**
 @brief z=x/y
 */
dcomplex zdiv_zd(dcomplex x, double y)
{
  dcomplex z;
  Z_R(z)=Z_R(x)/y;
  Z_I(z)=Z_I(x)/y;
  return z;
}

/**
 @brief z=x/y
 */
dcomplex zdiv_dz(double x, dcomplex y)
{
  dcomplex z;
  double den;
  den=Z_R(y)*Z_R(y)+Z_I(y)*Z_I(y);
  Z_R(z)= x*Z_R(y)/den;
  Z_I(z)=-x*Z_I(y)/den;
  return z;
}

/**
 @brief z=x^y
 */
dcomplex zpow_zd(dcomplex x, double y)
{
  dcomplex z;
  double r,theta;
  zget_polar(&r,&theta,x); // x=r*exp(i*theta)
  r=pow(r,y);              // r=r^y
  theta=theta*y;           // theta=theta*y
  z=zset_polar(r,theta);   // z=r^y*cos(theta*y)+i*r^y*sin(theta*y)
  return z;
}

/** @} */

/////////////////////////////////////////////

/** @name dcomplex型に関する関数 x<=>y */
/** @{ */

static int __zcmp_order=0;

/**
 @brief dcomplex型の値の比較 x<=>y を実部，虚部で判定
 */
void zcmp_set_real_imag()
{
  __zcmp_order=0;
}

/**
 @brief dcomplex型の値の比較 x<=>y を絶対値，偏角で判定
 */
void zcmp_set_abs_arg()
{
  __zcmp_order=1;
}

/**
 @brief dcomplex型の値の比較 x<=>y の方法の取得
 */
int zcmp_get_type()
{
  return __zcmp_order;
}

/**
 @brief x > y
 */
int gt_zz(dcomplex x, dcomplex y)
{
  double abs_x,abs_y,arg_x,arg_y;
  if(zcmp_get_type()){
    abs_x=Z_ABS2(x);
    abs_y=Z_ABS2(y);
    if(abs_x>abs_y){ return 1; }
    if(abs_x<abs_y){ return 0; }
    arg_x=Z_ARG(x);
    arg_y=Z_ARG(x);
    if(arg_x>arg_y){ return 1; }
    return 0;
  }else{
    if(Z_R(x)>Z_R(y)){ return 1; }
    if(Z_R(x)<Z_R(y)){ return 0; }
    if(Z_I(x)>Z_I(y)){ return 1; }
    if(Z_I(x)<Z_I(y)){ return 0; }
    return 0;
  }
}

/**
 @brief x < y
 */
int lt_zz(dcomplex x, dcomplex y)
{
  double abs_x,abs_y,arg_x,arg_y;
  if(zcmp_get_type()){
    abs_x=Z_ABS2(x);
    abs_y=Z_ABS2(y);
    if(abs_x<abs_y){ return 1; }
    if(abs_x>abs_y){ return 0; }
    arg_x=Z_ARG(x);
    arg_y=Z_ARG(x);
    if(arg_x<arg_y){ return 1; }
    return 0;
  }else{
    if(Z_R(x)<Z_R(y)){ return 1; }
    if(Z_R(x)>Z_R(y)){ return 0; }
    if(Z_I(x)<Z_I(y)){ return 1; }
    if(Z_I(x)>Z_I(y)){ return 0; }
    return 0;
  }
}

/** @} */


//EOF
