#ifndef ISYS_IRMAT_H
#define ISYS_IRMAT_H

#include<is_rmulti.h>
#include<is_rvec.h>
#include<is_func.h>


/**
 @file  is_irmat.h
 @brief 多倍長精度実数型rmultiの機械区間演算の行列に関する関数の宣言.
*/

/*
 * setting
 */
// [B0,B1]=A
void irmat_set_d(int m, int n, rmulti **B0, rmulti **B1, int LDB, const double *A, int LDA);
// [B0,B1]=[A0,A1]
void irmat_copy(int m, int n, rmulti **B0, int LDB0, rmulti **B1, int LDB1, rmulti **A0, int LDA0, rmulti **A1, int LDA1);

/*
 * I/O
 */
void irmat_print(int m, int n, rmulti **A0, int LDA0, rmulti **A1, int LDA1, char *name, char format, int digits);

/*
 * casting
 */
void irmat_get_s(int m, int n, char **B, int LDB, rmulti **A0, int LDA0, rmulti **A1, int LDA1, char format, int digits);

/*
 * operations
 */
// [B0,B1]=[A0,A1]'
void irmat_copy_t(int m, int n, rmulti **B0, int LDB0, rmulti **B1, int LDB1, rmulti **A0, int LDA0, rmulti **A1, int LDA1);
// Ac=(A1+A0)/2, Ar=A1-A0
void irmat_mid(int m, int n, rmulti **mid, int LDmid, rmulti **A0, int LDA0, rmulti **A1, int LDA1);
void irmat_rad(int m, int n, rmulti **rad, int LDrad, rmulti **A0, int LDA0, rmulti **A1, int LDA1);
void irmat_mr(int m, int n, rmulti **mid, int LDmid, rmulti **rad, int LDrad, rmulti **A0, int LDA0, rmulti **A1, int LDA1);
void irmat_center_radius(int m, int n, rmulti **Ac, int LDAc, rmulti **Ar, int LDAr, rmulti **A0, int LDA0, rmulti **A1, int LDA1);
//  [A0,A1]=[NaN(m,n),NaN(m,n)]
void irmat_set_nan(int m, int n, rmulti **A0, int LDA0, rmulti **A1, int LDA1);
// [A0,A1]=ones(m,n)*a
void irmat_set_all_d(int m, int n, rmulti **A0, int LDA0, rmulti **A1, int LDA1, double a);
// [A0,A1]=[zeros(m,n),zeros(m,n)]
void irmat_set_zeros(int m, int n, rmulti **A0, int LDA0, rmulti **A1, int LDA1);
// A=eye(m,n)
void irmat_set_eye(int m, int n, rmulti **A0, rmulti **A1, int LDA);
// [C0,C1]=[A0,A1]*[B0,B1]
void irmat_prod(int l, int m, int n, rmulti **C0, int LDC0, rmulti **C1, int LDC1, rmulti **A0, int LDA0, rmulti **A1, int LDA1, rmulti **B0, int LDB0, rmulti **B1, int LDB1);
// [C0,C1]=[C0,C1]+[A0,A1]*[B0,B1]
void irmat_add_prod(int l, int m, int n, rmulti **C0, int LDC0, rmulti **C1, int LDC1, rmulti **A0, int LDA0, rmulti **A1, int LDA1, rmulti **B0, int LDB0, rmulti **B1, int LDB1);
// [C0,C1]=[C0,C1]-[A0,A1]*[B0,B1]
void irmat_sub_prod(int l, int m, int n, rmulti **C0, int LDC0, rmulti **C1, int LDC1, rmulti **A0, int LDA0, rmulti **A1, int LDA1, rmulti **B0, int LDB0, rmulti **B1, int LDB1);
// [B0,B1]=inv(A0,A1)
void irmat_inv(int n, rmulti **B0, rmulti **B1, int LDB, rmulti **A0, rmulti **A1, int LDA);
// [B0,B1]=sum(A0,A1)
void irvec_sum_irmat(int m, int n, rmulti **B0, rmulti **B1, rmulti **A0, rmulti **A1, int LDA);
void irvec_max_irmat(int m, int n, rmulti **B0, rmulti **B1, rmulti **A0, rmulti **A1, int LDA);
// [B0,B1]=[A0,max(A1)]
void irvec_umax_irmat(int m, int n, rmulti **B0, rmulti **B1, rmulti **A0, rmulti **A1, int LDA);
void irvec_min_irmat(int m, int n, rmulti **B0, rmulti **B1, rmulti **A0, rmulti **A1, int LDA);
// [A0,A1]=F([x0,x1])
void irmat_func_list2(int m, int n, rmulti **A0, rmulti **A1, int LDA, func_t *f, int l, rmulti **x0, rmulti **x1);

#endif

