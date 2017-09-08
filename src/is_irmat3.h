#ifndef IS_IRMAT_H
#define IS_IRMAT_H

#include<is_rmulti.h>

/**
 @file  is_irmat.h
 @brief 多倍長精度実数型rmultiの行列に関する関数の宣言
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
int irmat3_copy(int m, int n, int l, rmulti **B0, rmulti **B1, int LDB1, int LDB2, rmulti **A0, rmulti **A1, int LDA1, int LDA2);
int irmat3_set_d(int m, int n, int l, rmulti **B0, rmulti **B1, int LDB1, int LDB2, double *A, int LDA1, int LDA2);

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
int irmat3_mid(int m, int n, int l, rmulti **mid, int LD1mid, int LD2mid, rmulti **A0, rmulti **A1, int LDA1, int LDA2);
int irmat3_rad(int m, int n, int l, rmulti **rad, int LD1rad, int LD2rad, rmulti **A0, rmulti **A1, int LDA1, int LDA2);
int irmat3_mr(int m, int n, int l, rmulti **mid, int LD1mid, int LD2mid, rmulti **rad, int LD1rad, int LD2rad, rmulti **A0, rmulti **A1, int LDA1, int LDA2);
int irmat3_add(int m, int n, int l, rmulti **Z0, rmulti **Z1, int LDZ1, int LDZ2, rmulti **X0, rmulti **X1, int LDX1, int LDX2, rmulti **Y0, rmulti **Y1, int LDY1, int LDY2);


/*
 * compare
 */

/*
 * mapping
 */

#endif
