#ifndef IS_RHSNBRG_H
#define IS_RHSNBRG_H

#include<is_rmulti.h>

/** @brief 密行列AをHessenberg型行列Bに相似変換する. */
void rhsnbrg_simtr(int n, rmulti **B, int LDB, rmulti **A, int LDA);

#endif
