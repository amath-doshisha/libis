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
rmulti **rmat3_allocate_clone_index(int LDB1, int LDB2, int l, rmulti **B, int *I, int *J, int *K);
rmulti **rmat3_free(int LDA1, int LDA2, int l, rmulti **A);
int rmat3_round(int m, int n, int l, rmulti **A, int LDA1, int LDA2, int prec);
int rmat3_clone(int m, int n, int l, rmulti **B, int LDB1, int LDB2, rmulti **A, int LDA1, int LDA2);
void rmat3_swap(int m, int n, int l, rmulti **A, int LDA1, int LDA2, rmulti **B, int LDB1, int LDB2);

/*
 * member variables
 */

/*
 * I/O
 */
void rmat3_print(int m, int n, int l, rmulti **A, int LDA1, int LDA2, const char *name, const char *f, int digits);

/*
 * setting
 */
int rmat3_set_all_d(int m, int n, int l, rmulti **A, int LDA1, int LDA2, double a);
int rmat3_set_zeros(int m, int n, int l, rmulti **A, int LDA1, int LDA2);
int rmat3_set_ones(int m, int n, int l, rmulti **A, int LDA1, int LDA2);
void rmat3_set_rand(int m, int n, int l, rmulti **A, int LDA1, int LDA2, double a, double b);

/*
 * casting
 */

/*
 * rearange elements
 */

/*
 * operations with auto precision mode
 */
// B=A
int rmat3_copy(int m, int n, int l, rmulti **B, int LDB1, int LDB2, rmulti **A, int LDA1, int LDA2);
// B=-A
int rmat3_neg(int m, int n, int l, rmulti **B, int LDB1, int LDB2, rmulti **A, int LDA1, int LDA2);
// C=A+B
int rmat3_add(int m, int n, int l, rmulti **C, int LDC1, int LDC2, rmulti **A, int LDA1, int LDA2, rmulti **B, int LDB1, int LDB2);
// C=A+b
int rmat3_add_r(int m, int n, int l, rmulti **C, int LDC1, int LDC2, rmulti **A, int LDA1, int LDA2, rmulti *b);



/*
 * oparations
 */

/*
 * compare
 */

/*
 * mapping
 */

#endif
