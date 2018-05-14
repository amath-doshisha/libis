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
