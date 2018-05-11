#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<stdarg.h>
#include"is_macros.h"
#include"is_dcomplex.h"


///////////////////////////////////////////////////////////////
// dcomplex operators
///////////////////////////////////////////////////////////////

// y=-x
dcomplex zneg(dcomplex x)
{
  dcomplex z;
  Z_R(z)=Z_NEG_R(x);
  Z_I(z)=Z_NEG_I(x);
  return z;
}

/////////////////////////////////////////////

/**
 * @brief z=x+y
 */
dcomplex zadd_z(dcomplex x, dcomplex y)
{
  dcomplex z;
  Z_R(z)=Z_R(x)+Z_R(y);
  Z_I(z)=Z_I(x)+Z_I(y);
  return z;
}

/**
 * @brief z=x+y
 */
dcomplex zadd_d(dcomplex x, double y)
{
  dcomplex z;
  Z_R(z)=Z_R(x)+y;
  Z_I(z)=Z_I(x);
  return z;
}

/**
 * @brief z=x+y
 */
dcomplex dadd_z(double x, dcomplex y)
{
  dcomplex z;
  Z_R(z)=x+Z_R(y);
  Z_I(z)=  Z_I(y);
  return z;
}

/////////////////////////////////////////////

/**
 * @brief z=x-y
 */
dcomplex zsub_z(dcomplex x, dcomplex y)
{
  dcomplex z;
  Z_R(z)=Z_R(x)-Z_R(y);
  Z_I(z)=Z_I(x)-Z_I(y);
  return z;
}

/**
 * @brief z=x-y
 */
dcomplex zsub_d(dcomplex x, double y)
{
  dcomplex z;
  Z_R(z)=Z_R(x)-y;
  Z_I(z)=Z_I(x);
  return z;
}

/**
 * @brief z=x-y
 */
dcomplex dsub_z(double x, dcomplex y)
{
  dcomplex z;
  Z_R(z)=x-Z_R(y);
  Z_I(z)= -Z_I(y);
  return z;
}

/////////////////////////////////////////////

/**
 * @brief z=x*y
 */
dcomplex zmul_z(dcomplex x, dcomplex y)
{
  dcomplex z;
  Z_R(z)=Z_R(x)*Z_R(y)-Z_I(x)*Z_I(y);
  Z_I(z)=Z_R(x)*Z_I(y)+Z_I(x)*Z_R(y);
  return z;
}

/**
 * @brief z=x*y
 */
dcomplex zmul_d(dcomplex x, double y)
{
  dcomplex z;
  Z_R(z)=Z_R(x)*y;
  Z_I(z)=Z_I(x)*y;
  return z;
}

/**
 * @brief z=x*y
 */
dcomplex dmul_z(double x, dcomplex y)
{
  dcomplex z;
  Z_R(z)=x*Z_R(y);
  Z_I(z)=x*Z_I(y);
  return z;
}

/////////////////////////////////////////////

/**
 * @brief z=conj(x)*y
 */
dcomplex zdot_z(dcomplex x, dcomplex y)
{
  dcomplex z;
  Z_R(z)=Z_R(x)*Z_R(y)+Z_I(x)*Z_I(y);
  Z_I(z)=Z_R(x)*Z_I(y)-Z_I(x)*Z_R(y);
  return z;
}

/**
 * @brief z=conj(x)*y
 */
dcomplex zdot_d(dcomplex x, double y)
{
  dcomplex z;
  Z_R(z)= Z_R(x)*y;
  Z_I(z)=-Z_I(x)*y;
  return z;
}

/**
 * @brief z=conj(x)*y
 */
dcomplex ddot_z(double x, dcomplex y)
{
  dcomplex z;
  Z_R(z)=x*Z_R(y);
  Z_I(z)=x*Z_I(y);
  return z;
}


/////////////////////////////////////////////

/*
 * @brief z=x/y
 */
dcomplex zdiv_z(dcomplex x, dcomplex y)
{
  dcomplex z;
  double den;
  den=Z_DIV_DEN(y);
  Z_R(z)=Z_DIV_NUM_R(x,y);
  Z_I(z)=Z_DIV_NUM_I(x,y);
  Z_R(z)/=den;
  Z_I(z)/=den;
  return z;
}

/*
 * @brief z=x/y
 */
dcomplex zdiv_d(dcomplex x, double y)
{
  dcomplex z;
  Z_R(z)=Z_R(x)/y;
  Z_I(z)=Z_I(x)/y;
  return z;
}

/*
 * @brief z=x/y
 */
dcomplex ddiv_z(double x, dcomplex y)
{
  dcomplex z;
  Z_R(z)=x*Z_INV_R(y);
  Z_I(z)=x*Z_INV_I(y);
  return z;
}

/////////////////////////////////////////////////////////

// y=1/x
dcomplex zinv(dcomplex x)
{
  dcomplex z;
  double den;
  den=Z_INV_DEN(x);
  Z_R(z)=Z_INV_NUM_R(x);
  Z_I(z)=Z_INV_NUM_I(x);
  Z_R(z)/=den;
  Z_I(z)/=den;
  return z;
}

// y=x/|x|
dcomplex znormalize(dcomplex x)
{
  double a;
  dcomplex z;
  a=Z_ABS(x);
  Z_R(z)=Z_R(x)/a;
  Z_I(z)=Z_I(x)/a;
  return z;
}

// x>y
int zgt(dcomplex x, dcomplex y)
{
  if(Z_R(x)>Z_R(y)) return 1;
  if(Z_R(x)<Z_R(y)) return 0;
  if(Z_I(x)>Z_I(y)) return 1;
  if(Z_I(x)<Z_I(y)) return 0;
  return 0;
}

// x<y
int zlt(dcomplex x, dcomplex y)
{
  if(Z_R(x)<Z_R(y)) return 1;
  if(Z_R(x)>Z_R(y)) return 0;
  if(Z_I(x)<Z_I(y)) return 1;
  if(Z_I(x)>Z_I(y)) return 0;
  return 0;
  /*
  double abs_x,abs_y,arg_x,arg_y;
  abs_x=Z_ABS2(x);
  abs_y=Z_ABS2(y);
  if(abs_x<abs_y) return 1;
  if(abs_x>abs_y) return 0;
  arg_x=Z_ARG(x);
  arg_y=Z_ARG(x);
  if(arg_x<arg_y) return 1;
  return 0;
  */
}



//EOF
