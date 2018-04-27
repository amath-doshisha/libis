#ifndef ISYS_REIG_II_H
#define ISYS_REIG_II_H

#include<is_rmulti.h>

void reig_ii(int n, rmulti **X, int LDX, rmulti **A, int LDA, rmulti **lambda, int debug);
void reig_ii_1pair(int n, rmulti **X, rmulti **A, int LDA, rmulti *lambda, int debug);

#endif
