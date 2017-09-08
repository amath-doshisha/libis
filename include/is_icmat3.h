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
int icmat3_set_z(int m, int n, int l, cmulti **B0, cmulti **B1, int LDB1, int LDB2, dcomplex *A, int LDA1, int LDA2);
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



/*
 * compare
 */


#endif
