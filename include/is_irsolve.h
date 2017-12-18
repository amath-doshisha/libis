#ifndef IS_IRSOLVE_H
#define IS_IRSOLVE_H

#include<is_rmulti.h>

int irsolve(int n, int NRHS, rmulti **B0, rmulti **B1, int LDB, rmulti **A0, rmulti **A1, int LDA, int *info);


int irsolve_lu(int n, int NRHS, rmulti **B0, rmulti **B1, int LDB, rmulti **A0, rmulti **A1, int LDA, int *info);
int irsolve_lu_decomp(int n, rmulti **A0, rmulti **A1, int LDA, int *p0, int *P1, int *info);
int irsolve_lu_backsubs(int n, int NRHS, rmulti **B0, rmulti **B1, int LDB, rmulti **A0, rmulti **A1, int LDA, int *p0, int *p1);
int irsolve_gauss_sweeper(int n, int NRHS, rmulti **B0, rmulti **B1, int LDB, rmulti **A0, rmulti **A1, int LDA, int *info);

#endif
