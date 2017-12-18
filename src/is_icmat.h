#ifndef ISYS_ICMAT_H
#define ISYS_ICMAT_H

#include<is_cmulti.h>

int icmat_copy(int m, int n, cmulti **B0, int LDB0, cmulti **B1, int LDB1, cmulti **A0, int LDA0, cmulti **A1, int LDA1);          // [B0,B1]=[A0,A1]
//追加
int icmat_copy_irmat(int m, int n, cmulti **B0, int LDB0, cmulti **B1, int LDB1, rmulti **A0, int LDA0, rmulti **A1, int LDA1);    // [B0,B1]=[A0,A1]
int icmat_copy_rmat(int m, int n, cmulti **B0, int LDB0, cmulti **B1, int LDB1, rmulti **A0, int LDA0);                            // [B0,B1]=[A0,A1]
//ここまで
int icmat_copy_t(int m, int n, cmulti **B0, int LDB0, cmulti **B1, int LDB1, cmulti **A0, int LDA0, cmulti **A1, int LDA1);        // [B0,B1]=[A0,A1]'
int icmat_copy_ct(int m, int n, cmulti **B0, int LDB0, cmulti **B1, int LDB1, cmulti **A0, int LDA0, cmulti **A1, int LDA1);       // [B0,B1]=[A0,A1]'
void icmat_print(int m, int n, cmulti **A0, int LDA0, cmulti **A1, int LDA1, const char *name, const char *f, int digits);         // output
int icmat_mid(int m, int n, cmulti **mid, int LDmid, cmulti **A0, int LDA0, cmulti **A1, int LDA1); // [m+r,m-r]=[A0,A1]
int icmat_rad(int m, int n, cmulti **rad, int LDrad, cmulti **A0, int LDA0, cmulti **A1, int LDA1); // [m+r,m-r]=[A0,A1]
int icmat_mr(int m, int n, cmulti **mid, int LDmid, cmulti **rad, int LDrad, cmulti **A0, int LDA0, cmulti **A1, int LDA1); // [m+r,m-r]=[A0,A1]
int icmat_center_radius(int m, int n, cmulti **Ac, int LDAc, cmulti **Ar, int LDAr, cmulti **A0, int LDA0, cmulti **A1, int LDA1); // Ac=(A1+A0)/2, Ar=A1-A0
int icmat_set_all_d(int m, int n, cmulti **A0, int LDA0, cmulti **A1, int LDA1, double a);                                         // [A0,A1]=ones(m,n)*a
//追加
int icmat_set_nan(int m, int n, cmulti **A0, int LDA0, cmulti **A1, int LDA1);
//ここまで
int icmat_set_zeros(int m, int n, cmulti **A0, int LDA0, cmulti **A1, int LDA1);                                                   // [A0,A1]=[zeros(m,n),zeros(m,n)]
int icmat_prod(int l, int m, int n, cmulti **C0, int LDC0, cmulti **C1, int LDC1, cmulti **A0, int LDA0, cmulti **A1, int LDA1, cmulti **B0, int LDB0, cmulti **B1, int LDB1);     // [C0,C1]=[A0,A1]*[B0,B1]
//編集
int icmat_prod_r1(int l, int m, int n, cmulti **C0, int LDC0, cmulti **C1, int LDC1, rmulti **A0, int LDA0, rmulti **A1, int LDA1, cmulti **B0, int LDB0, cmulti **B1, int LDB1); // [C0,C1]=[A0,A1]*[B0,B1]
int icmat_prod_r2(int l, int m, int n, cmulti **C0, int LDC0, cmulti **C1, int LDC1, cmulti **A0, int LDA0, cmulti **A1, int LDA1, rmulti **B0, int LDB0, rmulti **B1, int LDB1); // [C0,C1]=[A0,A1]*[B0,B1]
//ここまで
int icmat_add_prod(int l, int m, int n, cmulti **C0, int LDC0, cmulti **C1, int LDC1, cmulti **A0, int LDA0, cmulti **A1, int LDA1, cmulti **B0, int LDB0, cmulti **B1, int LDB1); // [C0,C1]+=[A0,A1]*[B0,B1]
int icmat_sub_prod(int l, int m, int n, cmulti **C0, int LDC0, cmulti **C1, int LDC1, cmulti **A0, int LDA0, cmulti **A1, int LDA1, cmulti **B0, int LDB0, cmulti **B1, int LDB1); // [C0,C1]-=[A0,A1]*[B0,B1]
//追加
int icvec_sum_icmat(int m, int n, cmulti **B0, cmulti **B1, cmulti **A0, cmulti **A1, int LDA);  //[B0,B1]=sum([A0,A1])
//ここまで
int icmat_func_list2(int m, int n, cmulti **A0, cmulti **A1, int LDA, func_t *f, int l, cmulti **x0, cmulti **x1); // A=F(x)

#endif
