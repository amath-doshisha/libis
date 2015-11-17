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
cmulti **cmat3_free(int LDA1, int LDA2, int l, cmulti **A);

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

/*
 * compare
 */

/*
 * mapping
 */


#endif
