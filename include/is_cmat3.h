#ifndef IS_CMAT3_H
#define IS_CMAT3_H

#include<is_cmulti.h>
#include<is_func.h>

/**
 @file  is_cmat.h
 @brief 多倍長精度実数型cmultiの3次元配列に関する関数の宣言
 */

/*
 * constructions and destruction
 */
cmulti **cmat3_allocate(int LDA1, int LDA2, int l);
cmulti **cmat3_allocate_clone_index(int LDA1, int LDA2, int l, cmulti **B, int LDB1, int LDB2, int *I, int *J, int *K);
cmulti **cmat3_free(int LDA1, int LDA2, int l, cmulti **A);
int cmat3_clone(int m, int n, int l, cmulti **B, int LDB1, int LDB2, cmulti **A, int LDA1, int LDA2);
int cmat3_clone_index(int m, int n, int l, cmulti **B, int LDB1, int LDB2, cmulti **A, int LDA1, int LDA2, int *I, int *J, int *K);

/*
 * member variables
 */

/*
 * I/O
 */
void cmat3_print(int m, int n, int l, cmulti **A, int LDA1, int LDA2, const char *name, const char *f, int digits);
/*
 * setting
 */
// B=cmulti(A)
int cmat3_set_z(int m, int n, int l, cmulti **B, int LDB1, int LDB2, dcomplex *A, int LDA1, int LDA2);
// B=cmulti(A)
int cmat3_set_d(int m, int n, int l, cmulti **B, int LDB1, int LDB2, double *A, int LDA1, int LDA2);
// A=a*ones(m,n,l)
int cmat3_set_all_z(int m, int n, int l, cmulti **A, int LDA1, int LDA2, dcomplex a);
// A=a*ones(m,n,l)
int cmat3_set_all_dd(int m, int n, int l, cmulti **A, int LDA1, int LDA2, double a_r, double a_i);
// A=a*ones(m,n,l)
int cmat3_set_all_d(int m, int n, int l, cmulti **A, int LDA1, int LDA2, double a);
// A=zeros(m,n,l)
int cmat3_set_zeros(int m, int n, int l, cmulti **A, int LDA1, int LDA2);
// A=ones(m,n,l)
int cmat3_set_ones(int m, int n, int l, cmulti **A, int LDA1, int LDA2);
// A=a*rand(m,n,l)+b
void cmat3_set_rand(int m, int n, int l, cmulti **A, int LDA1, int LDA2, double a, double b);

/*
 * casting
 */
// B=dcomplexeA)
void cmat3_get_z(int m, int n, int l, dcomplex *B, int LDB1, int LDB2, cmulti **A, int LDA1, int LDA2);

/*
 * rearange elements
 */

/*
 * operations with auto precision mode
 */
// B=A
int cmat3_copy(int m, int n, int l, cmulti **B, int LDB1, int LDB2, cmulti **A, int LDA1, int LDA2);
// B=-A
int cmat3_neg(int m, int n, int l, cmulti **B, int LDB1, int LDB2, cmulti **A, int LDA1, int LDA2);
// C=A+B
int cmat3_add(int m, int n, int l, cmulti **C, int LDC1, int LDC2, cmulti **A, int LDA1, int LDA2, cmulti **B, int LDB1, int LDB2);
// C=A+B
int cmat3_add_rmat3(int m, int n, int l, cmulti **C, int LDC1, int LDC2, cmulti **A, int LDA1, int LDA2, rmulti **B, int LDB1, int LDB2);
// C=A+b
int cmat3_add_c(int m, int n, int l, cmulti **C, int LDC1, int LDC2, cmulti **A, int LDA1, int LDA2, cmulti *b);
// C=A+b
int cmat3_add_r(int m, int n, int l, cmulti **C, int LDC1, int LDC2, cmulti **A, int LDA1, int LDA2, rmulti *b);
// C=A+b
int rmat3_add_c(int m, int n, int l, cmulti **C, int LDC1, int LDC2, rmulti **A, int LDA1, int LDA2, cmulti *b);
// C=A.*B
int cmat3_mul(int m, int n, int l, cmulti **C, int LDC1, int LDC2, cmulti **A, int LDA1, int LDA2, cmulti **B, int LDB1, int LDB2);
// C=A.*B
int cmat3_mul_rmat3(int m, int n, int l, cmulti **C, int LDC1, int LDC2, cmulti **A, int LDA1, int LDA2, rmulti **B, int LDB1, int LDB2);
// C=A.*b
int cmat3_mul_c(int m, int n, int l, cmulti **C, int LDC1, int LDC2, cmulti **A, int LDA1, int LDA2, cmulti *b);
// C=A.*b
int cmat3_mul_r(int m, int n, int l, cmulti **C, int LDC1, int LDC2, cmulti **A, int LDA1, int LDA2, rmulti *b);
// C=A.*b
int rmat3_mul_c(int m, int n, int l, cmulti **C, int LDC1, int LDC2, rmulti **A, int LDA1, int LDA2, cmulti *b);

/*
 * oparations
 */
// C=A./B
int cmat3_div(int m, int n, int l, cmulti **C, int LDC1, int LDC2, cmulti **A, int LDA1, int LDA2, cmulti **B, int LDB1, int LDB2);
// C=A./B
int cmat3_div_rmat3(int m, int n, int l, cmulti **C, int LDC1, int LDC2, cmulti **A, int LDA1, int LDA2, rmulti **B, int LDB1, int LDB2);
// C=A./B
int rmat3_div_cmat3(int m, int n, int l, cmulti **C, int LDC1, int LDC2, rmulti **A, int LDA1, int LDA2, cmulti **B, int LDB1, int LDB2);
// C=A./b
int cmat3_div_c2(int m, int n, int l, cmulti **C, int LDC1, int LDC2, cmulti **A, int LDA1, int LDA2, cmulti *b);
// C=A./b
int cmat3_div_r2(int m, int n, int l, cmulti **C, int LDC1, int LDC2, cmulti **A, int LDA1, int LDA2, rmulti *b);
// C=A./b
int rmat3_div_c2(int m, int n, int l, cmulti **C, int LDC1, int LDC2, rmulti **A, int LDA1, int LDA2, cmulti *b);
// C=a./B
int cmat3_div_c1(int m, int n, int l, cmulti **C, int LDC1, int LDC2, cmulti *a, cmulti **B, int LDB1, int LDB2);
// C=a./B
int cmat3_div_r1(int m, int n, int l, cmulti **C, int LDC1, int LDC2, rmulti *a, cmulti **B, int LDB1, int LDB2);
// C=a./B
int rmat3_div_c1(int m, int n, int l, cmulti **C, int LDC1, int LDC2, cmulti *a, rmulti **B, int LDB1, int LDB2);
// C=A.^B
int cmat3_pow(int m, int n, int l, cmulti **C, int LDC1, int LDC2, cmulti **A, int LDA1, int LDA2, cmulti **B, int LDB1, int LDB2);
// C=A.^B
int cmat3_pow_rmat3(int m, int n, int l, cmulti **C, int LDC1, int LDC2, cmulti **A, int LDA1, int LDA2, rmulti **B, int LDB1, int LDB2);
// C=A.^B
int rmat3_pow_cmat3(int m, int n, int l, cmulti **C, int LDC1, int LDC2, rmulti **A, int LDA1, int LDA2, cmulti **B, int LDB1, int LDB2);
// C=A.^b
int cmat3_pow_c2(int m, int n, int l, cmulti **C, int LDC1, int LDC2, cmulti **A, int LDA1, int LDA2, cmulti *b);
// C=A.^b
int cmat3_pow_r2(int m, int n, int l, cmulti **C, int LDC1, int LDC2, cmulti **A, int LDA1, int LDA2, rmulti *b);
// C=A.^b
int rmat3_pow_c2(int m, int n, int l, cmulti **C, int LDC1, int LDC2, rmulti **A, int LDA1, int LDA2, cmulti *b);
// C=a.^B
int cmat3_pow_c1(int m, int n, int l, cmulti **C, int LDC1, int LDC2, cmulti *a, cmulti **B, int LDB1, int LDB2);
// C=a.^B
int cmat3_pow_r1(int m, int n, int l, cmulti **C, int LDC1, int LDC2, rmulti *a, cmulti **B, int LDB1, int LDB2);
// C=a.^B
int rmat3_pow_c1(int m, int n, int l, cmulti **C, int LDC1, int LDC2, cmulti *a, rmulti **B, int LDB1, int LDB2);

/*
 * compare
 */
// C=(A==B)
void cmat3_eq(int m, int n, int l, int *C, int LDC1, int LDC2, cmulti **A, int LDA1, int LDA2, cmulti **B, int LDB1, int LDB2);
// C=(A==B)
void cmat3_eq_rmat3(int m, int n, int l, int *C, int LDC1, int LDC2, cmulti **A, int LDA1, int LDA2, rmulti **B, int LDB1, int LDB2);
// C=(A==b)
void cmat3_eq_c(int m, int n, int l, int *C, int LDC1, int LDC2, cmulti **A, int LDA1, int LDA2, cmulti *b);
// C=(A==b)
void cmat3_eq_r(int m, int n, int l, int *C, int LDC1, int LDC2, cmulti **A, int LDA1, int LDA2, rmulti *b);
// C=(A==b)
void rmat3_eq_c(int m, int n, int l, int *C, int LDC1, int LDC2, rmulti **A, int LDA1, int LDA2, cmulti *b);
// C=(A!=B)
void cmat3_ne(int m, int n, int l, int *C, int LDC1, int LDC2, cmulti **A, int LDA1, int LDA2, cmulti **B, int LDB1, int LDB2);
// C=(A!=B)
void cmat3_ne_rmat3(int m, int n, int l, int *C, int LDC1, int LDC2, cmulti **A, int LDA1, int LDA2, rmulti **B, int LDB1, int LDB2);
// C=(A!=b)
void cmat3_ne_c(int m, int n, int l, int *C, int LDC1, int LDC2, cmulti **A, int LDA1, int LDA2, cmulti *b);
// C=(A!=b)
void cmat3_ne_r(int m, int n, int l, int *C, int LDC1, int LDC2, cmulti **A, int LDA1, int LDA2, rmulti *b);
// C=(A!=b)
void rmat3_ne_c(int m, int n, int l, int *C, int LDC1, int LDC2, rmulti **A, int LDA1, int LDA2, cmulti *b);
// C=(A>=B)
void cmat3_ge(int m, int n, int l, int *C, int LDC1, int LDC2, cmulti **A, int LDA1, int LDA2, cmulti **B, int LDB1, int LDB2);
// C=(A>=B)
void cmat3_ge_rmat3(int m, int n, int l, int *C, int LDC1, int LDC2, cmulti **A, int LDA1, int LDA2, rmulti **B, int LDB1, int LDB2);
// C=(A>=b)
void cmat3_ge_c(int m, int n, int l, int *C, int LDC1, int LDC2, cmulti **A, int LDA1, int LDA2, cmulti *b);
// C=(A>=b)
void cmat3_ge_r(int m, int n, int l, int *C, int LDC1, int LDC2, cmulti **A, int LDA1, int LDA2, rmulti *b);
// C=(A>=b)
void rmat3_ge_c(int m, int n, int l, int *C, int LDC1, int LDC2, rmulti **A, int LDA1, int LDA2, cmulti *b);
// C=(A>B)
void cmat3_gt(int m, int n, int l, int *C, int LDC1, int LDC2, cmulti **A, int LDA1, int LDA2, cmulti **B, int LDB1, int LDB2);
// C=(A>B)
void cmat3_gt_rmat3(int m, int n, int l, int *C, int LDC1, int LDC2, cmulti **A, int LDA1, int LDA2, rmulti **B, int LDB1, int LDB2);
// C=(A>b)
void cmat3_gt_c(int m, int n, int l, int *C, int LDC1, int LDC2, cmulti **A, int LDA1, int LDA2, cmulti *b);
// C=(A>b)
void cmat3_gt_r(int m, int n, int l, int *C, int LDC1, int LDC2, cmulti **A, int LDA1, int LDA2, rmulti *b);
// C=(A>b)
void rmat3_gt_c(int m, int n, int l, int *C, int LDC1, int LDC2, rmulti **A, int LDA1, int LDA2, cmulti *b);
// C=(A<=B)
void cmat3_le(int m, int n, int l, int *C, int LDC1, int LDC2, cmulti **A, int LDA1, int LDA2, cmulti **B, int LDB1, int LDB2);
// C=(A<=B)
void cmat3_le_rmat3(int m, int n, int l, int *C, int LDC1, int LDC2, cmulti **A, int LDA1, int LDA2, rmulti **B, int LDB1, int LDB2);
// C=(A<=b)
void cmat3_le_c(int m, int n, int l, int *C, int LDC1, int LDC2, cmulti **A, int LDA1, int LDA2, cmulti *b);
// C=(A<=b)
void cmat3_le_r(int m, int n, int l, int *C, int LDC1, int LDC2, cmulti **A, int LDA1, int LDA2, rmulti *b);
// C=(A<=b)
void rmat3_le_c(int m, int n, int l, int *C, int LDC1, int LDC2, rmulti **A, int LDA1, int LDA2, cmulti *b);
// C=(A<B)
void cmat3_lt(int m, int n, int l, int *C, int LDC1, int LDC2, cmulti **A, int LDA1, int LDA2, cmulti **B, int LDB1, int LDB2);
// C=(A<B)
void cmat3_lt_rmat3(int m, int n, int l, int *C, int LDC1, int LDC2, cmulti **A, int LDA1, int LDA2, rmulti **B, int LDB1, int LDB2);
// C=(A<b)
void cmat3_lt_c(int m, int n, int l, int *C, int LDC1, int LDC2, cmulti **A, int LDA1, int LDA2, cmulti *b);
// C=(A<b)
void cmat3_lt_r(int m, int n, int l, int *C, int LDC1, int LDC2, cmulti **A, int LDA1, int LDA2, rmulti *b);
// C=(A<b)
void rmat3_lt_c(int m, int n, int l, int *C, int LDC1, int LDC2, rmulti **A, int LDA1, int LDA2, cmulti *b);


/*
 * mapping
 */


#endif
