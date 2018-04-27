#ifndef ISYS_CEIG_II_H
#define ISYS_CEIG_II_H

#include<is_cmulti.h>

void ceig_ii(int n, cmulti **X, int LDX, cmulti **A, int LDA, cmulti **lambda, int debug);
void ceig_ii_1pair(int n, cmulti **X, cmulti **A, int LDA, cmulti *lambda, int debug);

#endif
