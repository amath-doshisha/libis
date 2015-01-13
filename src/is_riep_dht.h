#ifndef ISYS_IREIG_H
#define ISYS_IREIG_H

#include<is_rmulti.h>

// QEより行列生成
void riep_dhToda_matgen(int m, int M, rmulti **A, int LDA, rmulti ***Q, rmulti ***E, int debug);

// dhtodaの逆固有値問題(qd)
void riep_dhToda_TN(int m, int M, rmulti **A, int LDA, rmulti **Q, int LDQ, rmulti **E, rmulti **lambda, rmulti **c, int debug);

// dhtodaの逆固有値問題(dqd)
void riep_dhToda_TN_dqd(int m, int M, rmulti **A, int LDA, rmulti **Q, int LDQ, rmulti **E, rmulti **lambda, rmulti **c, int debug);

// dhtodaの逆固有値問題(タウ関数)
void riep_dhToda_TN_tau(int m, int M, rmulti **A, int LDA, rmulti **Q, int LDQ, rmulti **E, rmulti **lambda, rmulti **c, int debug);



#endif




