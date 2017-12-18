#ifndef IS_ICSOLVE_H
#define IS_ICSOLVE_H

#include<is_cmulti.h>

int icsolve(int n, int NRHS, cmulti **B0, cmulti **B1, int LDB, cmulti **A0, cmulti **A1, int LDA, int *info);
int icsolve_gauss_sweeper(int n, int NRHS, cmulti **B0, cmulti **B1, int LDB, cmulti **A0, cmulti **A1, int LDA, int *info);

#endif
