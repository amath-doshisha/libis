#ifndef IS_ICMAT_H
#define IS_ICMAT_H

#include<is_cmulti.h>

/**
 @file  is_icmat.h
 @brief 多倍長精度実数型cmultiの行列に関する関数の宣言
 */

/*
 * constructions and destruction
 */

/*
 * member variables
 */

/*
 * I/O
 */


/*
 * setting
 */
int icmat3_copy(int m, int n, int l, cmulti **B0, cmulti **B1, int LDB1, int LDB2, cmulti **A0, cmulti **A1, int LDA1, int LDA2);
//追加
int icmat3_copy_irmat3(int m, int n, int l, cmulti **B0, cmulti **B1, int LDB1, int LDB2, rmulti **A0, rmulti **A1, int LDA1, int LDA2);
int icmat3_copy_irmat3_irmat3(int m, int n, int l, cmulti **C0, cmulti **C1, int LDC1, int LDC2, rmulti **A0, rmulti **A1, int LDA1, int LDA2, rmulti **B0, rmulti **B1, int LDB1, int LDB2);
int icmat3_copy_ir_irmat3(int m, int n, int l, cmulti **C0, cmulti **C1, int LDC1, int LDC2, rmulti *a0, rmulti *a1, rmulti **B0, rmulti **B1, int LDB1, int LDB2);
int icmat3_copy_irmat3_ir(int m, int n, int l, cmulti **C0, cmulti **C1, int LDC1, int LDC2, rmulti **A0, rmulti **A1, int LDA1, int LDA2, rmulti *b0, rmulti *b1);
//ここまで
int icmat3_set_z(int m, int n, int l, cmulti **B0, cmulti **B1, int LDB1, int LDB2, dcomplex *A, int LDA1, int LDA2);
//追加
int icmat3_set_zz(int m, int n, int l, cmulti **C0, cmulti **C1, int LDC1, int LDC2, dcomplex *A, int LDA1, int LDA2, dcomplex *B, int LDB1, int LDB2);
//ここまで
int icmat3_set_d(int m, int n, int l, cmulti **B0, cmulti **B1, int LDB1, int LDB2, double *A, int LDA1, int LDA2);

/*
 * casting
 */

/*
 * rearange elements
 */

/*
 * operations with auto precision mode
 */


/*
 * oparations
 */
int icmat3_mid(int m, int n, int l, cmulti **mid, int LD1mid, int LD2mid, cmulti **A0, cmulti **A1, int LDA1, int LDA2);
int icmat3_rad(int m, int n, int l, cmulti **rad, int LD1rad, int LD2rad, cmulti **A0, cmulti **A1, int LDA1, int LDA2);
int icmat3_mr(int m, int n, int l, cmulti **mid, int LD1mid, int LD2mid, cmulti **rad, int LD1rad, int LD2rad, cmulti **A0, cmulti **A1, int LDA1, int LDA2);
//追加
int icmat3_conj(int m, int n, int l, cmulti **B0, cmulti **B1, int LDB1, int LDB2, cmulti **A0, cmulti **A1, int LDA1, int LDA2);
//ここまで
//編集済み
int icmat3_neg(int m, int n, int l, cmulti **B0, cmulti **B1, int LDB1, int LDB2, cmulti **A0, cmulti **A1, int LDA1, int LDA2);
//ここまで
int icmat3_add(int m, int n, int l, cmulti **Z0, cmulti **Z1, int LDZ1, int LDZ2, cmulti **X0, cmulti **X1, int LDX1, int LDX2, cmulti **Y0, cmulti **Y1, int LDY1, int LDY2);
int icmat3_add_rmat3(int m, int n, int l, cmulti **Z0, cmulti **Z1, int LDZ1, int LDZ2, cmulti **X0, cmulti **X1, int LDX1, int LDX2, rmulti **Y0, rmulti **Y1, int LDY1, int LDY2);
//編集済み
int icmat3_add_r(int m, int n, int l, cmulti **Z0, cmulti **Z1, int LDZ1, int LDZ2, cmulti **X0, cmulti **X1, int LDX1, int LDX2, rmulti *Y0, rmulti *Y1);
int icmat3_add_c(int m, int n, int l, cmulti **Z0, cmulti **Z1, int LDZ1, int LDZ2, cmulti **X0, cmulti **X1, int LDX1, int LDX2, cmulti *Y0, cmulti *Y1);
int irmat3_add_c(int m, int n, int l, cmulti **Z0, cmulti **Z1, int LDZ1, int LDZ2, rmulti **X0, rmulti **X1, int LDX1, int LDX2, cmulti *Y0, cmulti *Y1);
int icmat3_sub(int m, int n, int l, cmulti **Z0, cmulti **Z1, int LDZ1, int LDZ2, cmulti **X0, cmulti **X1, int LDX1, int LDX2, cmulti **Y0, cmulti **Y1, int LDY1, int LDY2);
int icmat3_sub_c1(int m, int n, int l, cmulti **Z0, cmulti **Z1, int LDZ1, int LDZ2, cmulti *X0, cmulti *X1, cmulti **Y0, cmulti **Y1, int LDY1, int LDY2);
int icmat3_sub_c2(int m, int n, int l, cmulti **Z0, cmulti **Z1, int LDZ1, int LDZ2, cmulti **X0, cmulti **X1, int LDX1, int LDX2, cmulti *Y0, cmulti *Y1);
int icmat3_sub_r1(int m, int n, int l, cmulti **Z0, cmulti **Z1, int LDZ1, int LDZ2, rmulti *X0, rmulti *X1, cmulti **Y0, cmulti **Y1, int LDY1, int LDY2);
int icmat3_sub_r2(int m, int n, int l, cmulti **Z0, cmulti **Z1, int LDZ1, int LDZ2, cmulti **X0, cmulti **X1, int LDX1, int LDX2, rmulti *Y0, rmulti *Y1);
int irmat3_sub_c1(int m, int n, int l, cmulti **Z0, cmulti **Z1, int LDZ1, int LDZ2, cmulti *X0, cmulti *X1, rmulti **Y0, rmulti **Y1, int LDY1, int LDY2);
int irmat3_sub_c2(int m, int n, int l, cmulti **Z0, cmulti **Z1, int LDZ1, int LDZ2, rmulti **X0, rmulti **X1, int LDX1, int LDX2, cmulti *Y0, cmulti *Y1);
//ここまで
int icmat3_sub_rmat3(int m, int n, int l, cmulti **Z0, cmulti **Z1, int LDZ1, int LDZ2, cmulti **X0, cmulti **X1, int LDX1, int LDX2, rmulti **Y0, rmulti **Y1, int LDY1, int LDY2);
int irmat3_sub_cmat3(int m, int n, int l, cmulti **Z0, cmulti **Z1, int LDZ1, int LDZ2, rmulti **X0, rmulti **X1, int LDX1, int LDX2, cmulti **Y0, cmulti **Y1, int LDY1, int LDY2);
int icmat3_mul(int m, int n, int l, cmulti **Z0, cmulti **Z1, int LDZ1, int LDZ2, cmulti **X0, cmulti **X1, int LDX1, int LDX2, cmulti **Y0, cmulti **Y1, int LDY1, int LDY2);
int icmat3_mul_rmat3(int m, int n, int l, cmulti **Z0, cmulti **Z1, int LDZ1, int LDZ2, cmulti **X0, cmulti **X1, int LDX1, int LDX2, rmulti **Y0, rmulti **Y1, int LDY1, int LDY2);
int icmat3_mul_c(int m, int n, int l, cmulti **Z0, cmulti **Z1, int LDZ1, int LDZ2, cmulti **X0, cmulti **X1, int LDX1, int LDX2, cmulti *Y0, cmulti *Y1);
int icmat3_mul_r(int m, int n, int l, cmulti **Z0, cmulti **Z1, int LDZ1, int LDZ2, cmulti **X0, cmulti **X1, int LDX1, int LDX2, rmulti *Y0, rmulti *Y1);
int irmat3_mul_c(int m, int n, int l, cmulti **Z0, cmulti **Z1, int LDZ1, int LDZ2, rmulti **X0, rmulti **X1, int LDX1, int LDX2, cmulti *Y0, cmulti *Y1);
int icmat3_div(int m, int n, int l, cmulti **Z0, cmulti **Z1, int LDZ1, int LDZ2, cmulti **X0, cmulti **X1, int LDX1, int LDX2, cmulti **Y0, cmulti **Y1, int LDY1, int LDY2);
int icmat3_div_rmat3(int m, int n, int l, cmulti **Z0, cmulti **Z1, int LDZ1, int LDZ2, cmulti **X0, cmulti **X1, int LDX1, int LDX2, rmulti **Y0, rmulti **Y1, int LDY1, int LDY2);
int irmat3_div_cmat3(int m, int n, int l, cmulti **Z0, cmulti **Z1, int LDZ1, int LDZ2, rmulti **X0, rmulti **X1, int LDX1, int LDX2, cmulti **Y0, cmulti **Y1, int LDY1, int LDY2);
int icmat3_div_c1(int m, int n, int l, cmulti **Z0, cmulti **Z1, int LDZ1, int LDZ2, cmulti *X0, cmulti *X1, cmulti **Y0, cmulti **Y1, int LDY1, int LDY2);
int icmat3_div_c2(int m, int n, int l, cmulti **Z0, cmulti **Z1, int LDZ1, int LDZ2, cmulti **X0, cmulti **X1, int LDX1, int LDX2, cmulti *Y0, cmulti *Y1);
int icmat3_div_r1(int m, int n, int l, cmulti **Z0, cmulti **Z1, int LDZ1, int LDZ2, rmulti *X0, rmulti *X1, cmulti **Y0, cmulti **Y1, int LDY1, int LDY2);
int icmat3_div_r2(int m, int n, int l, cmulti **Z0, cmulti **Z1, int LDZ1, int LDZ2, cmulti **X0, cmulti **X1, int LDX1, int LDX2, rmulti *Y0, rmulti *Y1);
int irmat3_div_c1(int m, int n, int l, cmulti **Z0, cmulti **Z1, int LDZ1, int LDZ2, cmulti *X0, cmulti *X1, rmulti **Y0, rmulti **Y1, int LDY1, int LDY2);
int irmat3_div_c2(int m, int n, int l, cmulti **Z0, cmulti **Z1, int LDZ1, int LDZ2, rmulti **X0, rmulti **X1, int LDX1, int LDX2, cmulti *Y0, cmulti *Y1);
int icmat3_abs(int m, int n, int l, rmulti **C0, rmulti **C1, int LDC1, int LDC2, cmulti **A0, cmulti **A1, int LDA1, int LDA2);
int icmat3_sum(int m, int n, int l, cmulti **B0, cmulti **B1, int LDB1, int LDB2, cmulti **A0, cmulti**A1, int LDA1, int LDA2);
int icmat3_max(int m, int n, int l, cmulti **B0, cmulti **B1, int LDB1, int LDB2, cmulti **A0, cmulti **A1, int LDA1, int LDA2);
int icmat3_umax(int m, int n, int l, cmulti **B0, cmulti **B1, int LDB1, int LDB2, cmulti **A0, cmulti **A1, int LDA1, int LDA2);
int icmat3_min(int m, int n, int l, cmulti **B0, cmulti **B1, int LDB1, int LDB2, cmulti **A0, cmulti **A1, int LDA1, int LDA2);
//ここまで

/*
 * compare
 */


#endif
