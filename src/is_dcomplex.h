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
// z=x/y
dcomplex zdiv(dcomplex x, dcomplex y);
// y=1/x
dcomplex zinv(dcomplex x);
// y=x/|x|
dcomplex znormalize(dcomplex x);
// x>y
int zlt(dcomplex x, dcomplex y);
// x<y
int zgt(dcomplex x, dcomplex y);

#endif
