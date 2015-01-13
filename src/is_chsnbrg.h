#ifndef IS_CHSNBRG_H
#define IS_CHSNBRG_H

#include<is_cmulti.h>

/** @brief 密行列AをHessenberg型行列Bに相似変換する. */
void chsnbrg_simtr(int n, cmulti **B, int LDB, cmulti **A, int LDA);

/** @brief 多倍長Hessenbergに相似変換する(mステップで終了) */
void chsnbrg_simtr_step(int n, cmulti **A, int LDA, int step);

#endif
