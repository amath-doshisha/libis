#ifndef ISYS_IS_RIEP_DHT_H
#define ISYS_IS_RIEP_DHT_H

#include<is_rmulti.h>


// QEより行列Aを生成
void riep_dhToda_QE_to_A(int m, int M, rmulti **A, int LDA, rmulti **Q, int LDQ, rmulti **E, int debug);
void riep_EXTdhToda_QE_to_A(int m, int N, int M, rmulti **A, int LDA, rmulti **Q, int LDQ, rmulti **E, int LDE, int debug);

// 指定固有値よりdhTodaによりTN行列を生成．漸化式は方程式．
void riep_dhToda_TN(int m, int M, rmulti **A, int LDA, rmulti **Q, int LDQ, rmulti **E, rmulti **lambda, rmulti **c, int debug);
void riep_EXTdhToda_TN(int m, int N, int M, rmulti **A, int LDA, rmulti **Q, int LDQ, rmulti **E, int LDE, rmulti **lambda, rmulti **c, int debug);
  
#endif




