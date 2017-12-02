#ifndef IS_RMAT3_H
#define IS_RMAT3_H

#include<is_rmulti.h>
#include<is_func.h>

/**
 @file  is_rmat.h
 @brief 多倍長精度実数型rmultiの行列に関する関数の宣言
 */

/*
 * constructions and destruction
 */
rmulti **rmat3_allocate(int LDA1, int LDA2, int l);
rmulti **rmat3_allocate_prec(int LDA1, int LDA2, int l, int prec);
rmulti **rmat3_allocate_clone(int LDB1, int LDB2, int l, rmulti **B);
rmulti **rmat3_allocate_clone_index(int LDA1, int LDA2, int l, rmulti **B, int LDB1, int LDB2, int *I, int *J, int *K);
rmulti **rmat3_free(int LDA1, int LDA2, int l, rmulti **A);
int rmat3_round(int m, int n, int l, rmulti **A, int LDA1, int LDA2, int prec);
int rmat3_clone(int m, int n, int l, rmulti **B, int LDB1, int LDB2, rmulti **A, int LDA1, int LDA2);
int rmat3_clone_index(int m, int n, int l, rmulti **B, int LDB1, int LDB2, rmulti **A, int LDA1, int LDA2, int *I, int *J, int *K);
void rmat3_swap(int m, int n, int l, rmulti **A, int LDA1, int LDA2, rmulti **B, int LDB1, int LDB2);

/*
 * member variables
 */
void rmat3_get_prec(int m, int n, int l, int *P, int LDP1, int LDP2, rmulti **A, int LDA1, int LDA2);
int rmat3_get_prec_max(int m, int n, int l, rmulti **A, int LDA1, int LDA2);
void rmat3_get_exp(int m, int n, int l, int *P, int LDP1, int LDP2, rmulti **A, int LDA1, int LDA2);
void rmat3_get_sgn(int m, int n, int l, int *P, int LDP1, int LDP2, rmulti **A, int LDA1, int LDA2);

/*
 * I/O
 */
void rmat3_print(int m, int n, int l, rmulti **A, int LDA1, int LDA2, const char *name, const char *f, int digits);

/*
 * setting
 */
// B=rmulti(A)
int rmat3_set_d(int m, int n, int l, rmulti **B, int LDB1, int LDB2, double *A, int LDA1, int LDA2);
void rmat3_set_s(int m, int n, int l, rmulti **B, int LDB1, int LDB2, char **A, int LDA1, int LDA2);
// A=a*ones(m,n,l)
int rmat3_set_all_d(int m, int n, int l, rmulti **A, int LDA1, int LDA2, double a);
// A=zeros(m,n,l)
int rmat3_set_zeros(int m, int n, int l, rmulti **A, int LDA1, int LDA2);
// A=ones(m,n,l)
int rmat3_set_ones(int m, int n, int l, rmulti **A, int LDA1, int LDA2);
// A=nan(m,n,l)
int rmat3_set_nan(int m, int n, int l, rmulti **A, int LDA1, int LDA2);
// A=inf(m,n,l)
int rmat3_set_inf(int m, int n, int l, rmulti **A, int LDA1, int LDA2);


// A=a*rnad(m,n,l)+b
void rmat3_set_rand(int m, int n, int l, rmulti **A, int LDA1, int LDA2, double a, double b);

/*
 * casting
 */
// B=double(A)
void rmat3_get_d(int m, int n, int l, double *B, int LDB1, int LDB2, rmulti **A, int LDA1, int LDA2);

/*
 * rearange elements
 */

/*
 * operations with auto precision mode
 */
// B=A
int rmat3_copy(int m, int n, int l, rmulti **B, int LDB1, int LDB2, rmulti **A, int LDA1, int LDA2);
// B=A(I,J,K)
int rmat3_copy_index(int m, int n, int l, rmulti **B, int LDB1, int LDB2, rmulti **A, int LDA1, int LDA2, int *I, int *J, int *K);
// B(I,J,K)=A
int rmat3_index_copy(int m, int n, int l, rmulti **B, int LDB1, int LDB2, rmulti **A, int LDA1, int LDA2, int *I, int *J, int *K);
// B=-A
int rmat3_neg(int m, int n, int l, rmulti **B, int LDB1, int LDB2, rmulti **A, int LDA1, int LDA2);
// C=A+B
int rmat3_add(int m, int n, int l, rmulti **C, int LDC1, int LDC2, rmulti **A, int LDA1, int LDA2, rmulti **B, int LDB1, int LDB2);
// C=A+b
int rmat3_add_r(int m, int n, int l, rmulti **C, int LDC1, int LDC2, rmulti **A, int LDA1, int LDA2, rmulti *b);
// C=A-B
int rmat3_sub(int m, int n, int l, rmulti **C, int LDC1, int LDC2, rmulti **A, int LDA1, int LDA2, rmulti **B, int LDB1, int LDB2);
// C=a-B
int rmat3_sub_r1(int m, int n, int l, rmulti **C, int LDC1, int LDC2, rmulti *a, rmulti **B, int LDB1, int LDB2);
// C-A-b
int rmat3_sub_r2(int m, int n, int l, rmulti **C, int LDC1, int LDC2, rmulti **A, int LDA1, int LDA2, rmulti *b);
// C=A.*B
int rmat3_mul(int m, int n, int l, rmulti **C, int LDC1, int LDC2, rmulti **A, int LDA1, int LDA2, rmulti **B, int LDB1, int LDB2);
// C=A.*b
int rmat3_mul_r(int m, int n, int l, rmulti **C, int LDC1, int LDC2, rmulti **A, int LDA1, int LDA2, rmulti *b);

/*
 * oparations
 */
// C=abs(A)
int rmat3_abs(int m, int n, int l, rmulti **C, int LDC1, int LDC2, rmulti **A, int LDA1, int LDA2);
// B=max(A)
int rmat3_max(int m, int n, int l, rmulti **B, int LDB1, int LDB2, rmulti **A, int LDA1, int LDA2);
// B=min(A)
int rmat3_min(int m, int n, int l, rmulti **B, int LDB1, int LDB2, rmulti **A, int LDA1, int LDA2);
// C=A./B
int rmat3_div(int m, int n, int l, rmulti **C, int LDC1, int LDC2, rmulti **A, int LDA1, int LDA2, rmulti **B, int LDB1, int LDB2);
// C=a./B
int rmat3_div_r1(int m, int n, int l, rmulti **C, int LDC1, int LDC2, rmulti *a, rmulti **B, int LDB1, int LDB2);
// C=A./b
int rmat3_div_r2(int m, int n, int l, rmulti **C, int LDC1, int LDC2, rmulti **A, int LDA1, int LDA2, rmulti *b);
// C=A.^B
int rmat3_pow(int m, int n, int l, rmulti **C, int LDC1, int LDC2, rmulti **A, int LDA1, int LDA2, rmulti **B, int LDB1, int LDB2);
// C=a.^B
int rmat3_pow_r1(int m, int n, int l, rmulti **C, int LDC1, int LDC2, rmulti *a, rmulti **B, int LDB1, int LDB2);
// C=A.^b
int rmat3_pow_r2(int m, int n, int l, rmulti **C, int LDC1, int LDC2, rmulti **A, int LDA1, int LDA2, rmulti *b);

/*
 * compare
 */
// B=isnan(A)
void rmat3_isnan(int m, int n, int l, int *B, int LDB1, int LDB2, rmulti **A, int LDA1, int LDA2);
// B=isinf(A)
void rmat3_isinf(int m, int n, int l, int *B, int LDB1, int LDB2, rmulti **A, int LDA1, int LDA2);
// C=(A==B)
void rmat3_eq(int m, int n, int l, int *C, int LDC1, int LDC2, rmulti **A, int LDA1, int LDA2, rmulti **B, int LDB1, int LDB2);
// C=(A==b)
void rmat3_eq_r(int m, int n, int l, int *C, int LDC1, int LDC2, rmulti **A, int LDA1, int LDA2, rmulti *b);
// C=(A!=B)
void rmat3_ne(int m, int n, int l, int *C, int LDC1, int LDC2, rmulti **A, int LDA1, int LDA2, rmulti **B, int LDB1, int LDB2);
// C=(A!=b)
void rmat3_ne_r(int m, int n, int l, int *C, int LDC1, int LDC2, rmulti **A, int LDA1, int LDA2, rmulti *b);
// C=(A>=B)
void rmat3_ge(int m, int n, int l, int *C, int LDC1, int LDC2, rmulti **A, int LDA1, int LDA2, rmulti **B, int LDB1, int LDB2);
// C=(A>=b)
void rmat3_ge_r(int m, int n, int l, int *C, int LDC1, int LDC2, rmulti **A, int LDA1, int LDA2, rmulti *b);
// C=(A>B)
void rmat3_gt(int m, int n, int l, int *C, int LDC1, int LDC2, rmulti **A, int LDA1, int LDA2, rmulti **B, int LDB1, int LDB2);
// C=(A>b)
void rmat3_gt_r(int m, int n, int l, int *C, int LDC1, int LDC2, rmulti **A, int LDA1, int LDA2, rmulti *b);
// C=(A<=B)
void rmat3_le(int m, int n, int l, int *C, int LDC1, int LDC2, rmulti **A, int LDA1, int LDA2, rmulti **B, int LDB1, int LDB2);
// C=(A<=b)
void rmat3_le_r(int m, int n, int l, int *C, int LDC1, int LDC2, rmulti **A, int LDA1, int LDA2, rmulti *b);
// C=(A<B)
void rmat3_lt(int m, int n, int l, int *C, int LDC1, int LDC2, rmulti **A, int LDA1, int LDA2, rmulti **B, int LDB1, int LDB2);
// C=(A<b)
void rmat3_lt_r(int m, int n, int l, int *C, int LDC1, int LDC2, rmulti **A, int LDA1, int LDA2, rmulti *b);

/*
 * mapping
 */

#endif
