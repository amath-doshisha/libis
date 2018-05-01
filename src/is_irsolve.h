#ifndef IS_IRSOLVE_H
#define IS_IRSOLVE_H

#include<is_rmulti.h>

void irsolve(int n, int NRHS, rmulti **B0, rmulti **B1, int LDB, rmulti **A0, rmulti **A1, int LDA, int *info);
void irsolve_lu(int n, int NRHS, rmulti **B0, rmulti **B1, int LDB, rmulti **A0, rmulti **A1, int LDA, int *info);
void irsolve_lu_decomp(int n, rmulti **A0, rmulti **A1, int LDA, int *p0, int *P1, int *info);
void irsolve_lu_backsubs(int n, int NRHS, rmulti **B0, rmulti **B1, int LDB, rmulti **A0, rmulti **A1, int LDA, int *p0, int *p1);
void irsolve_gauss_sweeper(int n, int NRHS, rmulti **B0, rmulti **B1, int LDB, rmulti **A0, rmulti **A1, int LDA, int *info);

#endif
