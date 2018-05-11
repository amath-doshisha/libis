#ifndef IS_COMPLEX_H
#define IS_COMPLEX_H

/**
 @file   is_dcomplex.h
 @brief  倍精度複素数構造体に関する宣言と関数
 */

/**
 @brief 倍精度複素数構造体
 */
typedef struct {
  double r; /**< 虚部 */
  double i; /**< 実部 */
} dcomplex;

// functions of complex operators
// y=-x
dcomplex zneg(dcomplex x);
// y=1/x
dcomplex zinv(dcomplex x);
// y=x/|x|
dcomplex znormalize(dcomplex x);
// x>y
int zlt(dcomplex x, dcomplex y);
// x<y
int zgt(dcomplex x, dcomplex y);

/*
 * z=f(x,y)
 */
// z=x+y
dcomplex zadd_z(dcomplex x, dcomplex y);
dcomplex zadd_d(dcomplex x, double y);
dcomplex dadd_z(double x, dcomplex y);
// z=x-y
dcomplex zsub_z(dcomplex x, dcomplex y);
dcomplex zsub_d(dcomplex x, double y);
dcomplex dsub_z(double x, dcomplex y);
// z=x*y
dcomplex zmul_z(dcomplex x, dcomplex y);
dcomplex zmul_d(dcomplex x, double y);
dcomplex dmul_z(double x, dcomplex y);
// z=conj(x)*y
dcomplex zdot_z(dcomplex x, dcomplex y);
dcomplex zdot_d(dcomplex x, double y);
dcomplex ddot_z(double x, dcomplex y);
// z=x/y
dcomplex zdiv_z(dcomplex x, dcomplex y);
dcomplex zdiv_d(dcomplex x, double y);
dcomplex ddiv_z(double x, dcomplex y);



#endif
