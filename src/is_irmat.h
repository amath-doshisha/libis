#ifndef ISYS_IRMAT_H
#define ISYS_IRMAT_H

#include<is_rmulti.h>
#include<is_rvec.h>
#include<is_func.h>

// [B0,B1]=A
int irmat_set_d(int m, int n, rmulti **B0, rmulti **B1, int LDB, const double *A, int LDA);
// print
void irmat_print(int m, int n, rmulti **A0, int LDA0, rmulti **A1, int LDA1, const char *name, const char *f, int digits);
// [B0,B1]=[A0,A1]
int irmat_copy(int m, int n, rmulti **B0, int LDB0, rmulti **B1, int LDB1, rmulti **A0, int LDA0, rmulti **A1, int LDA1);
// [B0,B1]=[A0,A1]'
int irmat_copy_t(int m, int n, rmulti **B0, int LDB0, rmulti **B1, int LDB1, rmulti **A0, int LDA0, rmulti **A1, int LDA1);
// Ac=(A1+A0)/2, Ar=A1-A0
int irmat_mid(int m, int n, rmulti **mid, int LDmid, rmulti **A0, int LDA0, rmulti **A1, int LDA1);
int irmat_rad(int m, int n, rmulti **rad, int LDrad, rmulti **A0, int LDA0, rmulti **A1, int LDA1);
int irmat_mr(int m, int n, rmulti **mid, int LDmid, rmulti **rad, int LDrad, rmulti **A0, int LDA0, rmulti **A1, int LDA1);
int irmat_center_radius(int m, int n, rmulti **Ac, int LDAc, rmulti **Ar, int LDAr, rmulti **A0, int LDA0, rmulti **A1, int LDA1);
// [A0,A1]=ones(m,n)*a
int irmat_set_all_d(int m, int n, rmulti **A0, int LDA0, rmulti **A1, int LDA1, double a);
// [A0,A1]=[zeros(m,n),zeros(m,n)]
int irmat_set_zeros(int m, int n, rmulti **A0, int LDA0, rmulti **A1, int LDA1);
// [C0,C1]=[A0,A1]*[B0,B1]
int irmat_prod(int l, int m, int n, rmulti **C0, int LDC0, rmulti **C1, int LDC1, rmulti **A0, int LDA0, rmulti **A1, int LDA1, rmulti **B0, int LDB0, rmulti **B1, int LDB1);
// [C0,C1]=[C0,C1]+[A0,A1]*[B0,B1]
int irmat_add_prod(int l, int m, int n, rmulti **C0, int LDC0, rmulti **C1, int LDC1, rmulti **A0, int LDA0, rmulti **A1, int LDA1, rmulti **B0, int LDB0, rmulti **B1, int LDB1);
// [C0,C1]=[C0,C1]-[A0,A1]*[B0,B1]
int irmat_sub_prod(int l, int m, int n, rmulti **C0, int LDC0, rmulti **C1, int LDC1, rmulti **A0, int LDA0, rmulti **A1, int LDA1, rmulti **B0, int LDB0, rmulti **B1, int LDB1);
// [A0,A1]=F([x0,x1])
int irmat_func_list2(int m, int n, rmulti **A0, rmulti **A1, int LDA, func_t *f, int l, rmulti **x0, rmulti **x1);

#endif

