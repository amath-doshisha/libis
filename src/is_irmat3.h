#ifndef IS_IRMAT_H
#define IS_IRMAT_H

#include<is_rmulti.h>

/**
 @file  is_irmat.h
 @brief 多倍長精度実数型rmultiの行列に関する関数の宣言
 */

/*
 * setting
 */
void irmat3_copy(int m, int n, int l, rmulti **B0, rmulti **B1, int LDB1, int LDB2, rmulti **A0, rmulti **A1, int LDA1, int LDA2);
void irmat3_set_d(int m, int n, int l, rmulti **B0, rmulti **B1, int LDB1, int LDB2, double *A, int LDA1, int LDA2);
void irmat3_set_dd(int m, int n, int l, rmulti **C0, rmulti **C1, int LDC1, int LDC2, double *A, int LDA1, int LDA2, double *B, int LDB1, int LDB2);

/*
 * casting
 */

/*
 * oparations
 */
void irmat3_mid(int m, int n, int l, rmulti **mid, int LD1mid, int LD2mid, rmulti **A0, rmulti **A1, int LDA1, int LDA2);
void irmat3_rad(int m, int n, int l, rmulti **rad, int LD1rad, int LD2rad, rmulti **A0, rmulti **A1, int LDA1, int LDA2);
void irmat3_mr(int m, int n, int l, rmulti **mid, int LD1mid, int LD2mid, rmulti **rad, int LD1rad, int LD2rad, rmulti **A0, rmulti **A1, int LDA1, int LDA2);
void irmat3_neg(int m, int n, int l, rmulti **B0, rmulti **B1, int LDB1, int LDB2, rmulti **A0, rmulti **A1, int LDA1, int LDA2);
void irmat3_add(int m, int n, int l, rmulti **Z0, rmulti **Z1, int LDZ1, int LDZ2, rmulti **X0, rmulti **X1, int LDX1, int LDX2, rmulti **Y0, rmulti **Y1, int LDY1, int LDY2);
void irmat3_add_r(int m, int n, int l, rmulti **Z0, rmulti **Z1, int LDZ1, int LDZ2, rmulti **X0, rmulti **X1, int LDX1, int LDX2, rmulti *Y0, rmulti *Y1);
void irmat3_sub(int m, int n, int l, rmulti **Z0, rmulti **Z1, int LDZ1, int LDZ2, rmulti **X0, rmulti **X1, int LDX1, int LDX2, rmulti **Y0, rmulti **Y1, int LDY1, int LDY2);
void irmat3_sub_r1(int m, int n, int l, rmulti **Z0, rmulti **Z1, int LDZ1, int LDZ2, rmulti *X0, rmulti *X1, rmulti **Y0, rmulti **Y1, int LDY1, int LDY2);
void irmat3_sub_r2(int m, int n, int l, rmulti **Z0, rmulti **Z1, int LDZ1, int LDZ2, rmulti **X0, rmulti **X1, int LDX1, int LDX2, rmulti *Y0, rmulti *Y1);
void irmat3_mul(int m, int n, int l, rmulti **Z0, rmulti **Z1, int LDZ1, int LDZ2, rmulti **X0, rmulti **X1, int LDX1, int LDX2, rmulti **Y0, rmulti **Y1, int LDY1, int LDY2);
void irmat3_mul_r(int m, int n, int l, rmulti **Z0, rmulti **Z1, int LDZ1, int LDZ2, rmulti **X0, rmulti **X1, int LDX1, int LDX2, rmulti *Y0, rmulti *Y1);
void irmat3_div(int m, int n, int l, rmulti **Z0, rmulti **Z1, int LDZ1, int LDZ2, rmulti **X0, rmulti **X1, int LDX1, int LDX2, rmulti **Y0, rmulti **Y1, int LDY1, int LDY2);
void irmat3_div_r1(int m, int n, int l, rmulti **Z0, rmulti **Z1, int LDZ1, int LDZ2, rmulti *X0, rmulti *X1, rmulti **Y0, rmulti **Y1, int LDY1, int LDY2);
void irmat3_div_r2(int m, int n, int l, rmulti **Z0, rmulti **Z1, int LDZ1, int LDZ2, rmulti **X0, rmulti **X1, int LDX1, int LDX2, rmulti *Y0, rmulti *Y1);
void irmat3_abs(int m, int n, int l, rmulti **C0, rmulti **C1, int LDC1, int LDC2, rmulti **A0, rmulti **A1, int LDA1, int LDA2);
void irmat3_sqrt(int m, int n, int l, rmulti **C0, rmulti **C1, int LDC1, int LDC2, rmulti **A0, rmulti **A1, int LDA1, int LDA2);
void irmat3_log10(int m, int n, int l, rmulti **C0, rmulti **C1, int LDC1, int LDC2, rmulti **A0, rmulti **A1, int LDA1, int LDA2);
void irmat3_sum(int m, int n, int l, rmulti **B0, rmulti **B1, int LDB1, int LDB2, rmulti **A0, rmulti**A1, int LDA1, int LDA2);
void irmat3_max(int m, int n, int l, rmulti **B0, rmulti **B1, int LDB1, int LDB2, rmulti **A0, rmulti **A1, int LDA1, int LDA2);
void irmat3_umax(int m, int n, int l, rmulti **B0, rmulti **B1, int LDB1, int LDB2, rmulti **A0, rmulti **A1, int LDA1, int LDA2);
void irmat3_min(int m, int n, int l, rmulti **B0, rmulti **B1, int LDB1, int LDB2, rmulti **A0, rmulti **A1, int LDA1, int LDA2);

/*
 * compare
 */

#endif
