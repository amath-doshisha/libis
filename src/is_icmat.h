#ifndef ISYS_ICMAT_H
#define ISYS_ICMAT_H

#include<is_cmulti.h>


/**
 @file  is_icmat.h
 @brief 多倍長精度実数型cmultiの機械区間演算の行列に関する関数の宣言.
*/

/*
 * setting
 */
// [B0,B1]=[A0,A1]
void icmat_copy(int m, int n, cmulti **B0, int LDB0, cmulti **B1, int LDB1, cmulti **A0, int LDA0, cmulti **A1, int LDA1);
// [B0,B1]=[A0,A1]
void icmat_copy_irmat(int m, int n, cmulti **B0, int LDB0, cmulti **B1, int LDB1, rmulti **A0, int LDA0, rmulti **A1, int LDA1);
// [B0,B1]=[A0,A1]
void icmat_copy_rmat(int m, int n, cmulti **B0, int LDB0, cmulti **B1, int LDB1, rmulti **A0, int LDA0);

/*
 * I/O
 */
void icmat_print(int m, int n, cmulti **A0, int LDA0, cmulti **A1, int LDA1, char *name, char format, int digits);

/*
 * casting
 */
void icmat_get_s(int m, int n, char **B, int LDB, cmulti **A0, int LDA0, cmulti **A1, int LDA1, char format, int digits);

/*
 * operations
 */
// [B0,B1]=[A0,A1]'
void icmat_copy_t(int m, int n, cmulti **B0, int LDB0, cmulti **B1, int LDB1, cmulti **A0, int LDA0, cmulti **A1, int LDA1);
// [B0,B1]=[A0,A1]'
void icmat_copy_ct(int m, int n, cmulti **B0, int LDB0, cmulti **B1, int LDB1, cmulti **A0, int LDA0, cmulti **A1, int LDA1);
// [m+r,m-r]=[A0,A1]
void icmat_mid(int m, int n, cmulti **mid, int LDmid, cmulti **A0, int LDA0, cmulti **A1, int LDA1);
// [m+r,m-r]=[A0,A1]
void icmat_rad(int m, int n, cmulti **rad, int LDrad, cmulti **A0, int LDA0, cmulti **A1, int LDA1);
// [m+r,m-r]=[A0,A1]
void icmat_mr(int m, int n, cmulti **mid, int LDmid, cmulti **rad, int LDrad, cmulti **A0, int LDA0, cmulti **A1, int LDA1);
// Ac=(A1+A0)/2, Ar=A1-A0
void icmat_center_radius(int m, int n, cmulti **Ac, int LDAc, cmulti **Ar, int LDAr, cmulti **A0, int LDA0, cmulti **A1, int LDA1);
// [A0,A1]=ones(m,n)*a
void icmat_set_all_d(int m, int n, cmulti **A0, int LDA0, cmulti **A1, int LDA1, double a);
void icmat_set_nan(int m, int n, cmulti **A0, int LDA0, cmulti **A1, int LDA1);
// [A0,A1]=[zeros(m,n),zeros(m,n)]
void icmat_set_zeros(int m, int n, cmulti **A0, int LDA0, cmulti **A1, int LDA1);
//[A0,A1]=eye(A0,A1)
void icmat_set_eye(int m, int n, cmulti **A0, cmulti **A1, int LDA);
// [C0,C1]=[A0,A1]*[B0,B1]
void icmat_prod(int l, int m, int n, cmulti **C0, int LDC0, cmulti **C1, int LDC1, cmulti **A0, int LDA0, cmulti **A1, int LDA1, cmulti **B0, int LDB0, cmulti **B1, int LDB1);
// [C0,C1]=[A0,A1]*[B0,B1]
void icmat_prod_r1(int l, int m, int n, cmulti **C0, int LDC0, cmulti **C1, int LDC1, rmulti **A0, int LDA0, rmulti **A1, int LDA1, cmulti **B0, int LDB0, cmulti **B1, int LDB1);
// [C0,C1]=[A0,A1]*[B0,B1]
void icmat_prod_r2(int l, int m, int n, cmulti **C0, int LDC0, cmulti **C1, int LDC1, cmulti **A0, int LDA0, cmulti **A1, int LDA1, rmulti **B0, int LDB0, rmulti **B1, int LDB1);
// [C0,C1]+=[A0,A1]*[B0,B1]
void icmat_add_prod(int l, int m, int n, cmulti **C0, int LDC0, cmulti **C1, int LDC1, cmulti **A0, int LDA0, cmulti **A1, int LDA1, cmulti **B0, int LDB0, cmulti **B1, int LDB1);
// [C0,C1]-=[A0,A1]*[B0,B1]
void icmat_sub_prod(int l, int m, int n, cmulti **C0, int LDC0, cmulti **C1, int LDC1, cmulti **A0, int LDA0, cmulti **A1, int LDA1, cmulti **B0, int LDB0, cmulti **B1, int LDB1);
// [B0,B1]=sum([A0,A1])
void icvec_sum_icmat(int m, int n, cmulti **B0, cmulti **B1, cmulti **A0, cmulti **A1, int LDA);
//[B0,B1]=inv(A0,A1);
void icmat_inv(int n, cmulti **B0, cmulti **B1, int LDB, cmulti **A0, cmulti **A1, int LDA);
// [B0,B1]=[max(A0),max(A1)]
void icvec_max_icmat(int m, int n, cmulti **B0, cmulti **B1, cmulti **A0, cmulti **A1, int LDA);
// [B0,B1]=[A0,max(A1)]
void icvec_umax_icmat(int m, int n, cmulti **B0, cmulti **B1, cmulti **A0, cmulti **A1, int LDA);
// [B0,B1]=[min(A0),min(A1)]
void icvec_min_icmat(int m, int n, cmulti **B0, cmulti **B1, cmulti **A0, cmulti **A1, int LDA);

#endif
