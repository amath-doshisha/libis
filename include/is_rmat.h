#ifndef IS_RMAT_H
#define IS_RMAT_H

#include<is_rmulti.h>
#include<is_func.h>

/**
 @file  is_rmat.h
 @brief 多倍長精度実数型rmultiの行列に関する関数の宣言
 */


/*
 * constructions and destruction
 */
rmulti **rmat_allocate(int LDA, int n);
rmulti **rmat_allocate_prec(int LDA, int n, int prec);
rmulti **rmat_allocate_clone(int LDB, int n, rmulti **B);
rmulti **rmat_free(int LDA, int n, rmulti **A);
void rmat_round(int m, int n, rmulti **A, int LDA, int prec);
void rmat_clone(int m, int n, rmulti **B, int LDB, rmulti **A, int LDA);   // B=A
void rmat_clone_t(int m, int n, rmulti **B, int LDB, rmulti **A, int LDA); // B=A'
void rmat_clone_index(int m, int n, rmulti **B, int LDB, rmulti **A, int LDA, const int *I);
void rmat_swap(int m, int n, rmulti **A, int LDA, rmulti **B, int LDB);   // A<->B

/*
 * member variables
 */
void rmat_get_prec(int m, int n, int *P, int LDP, rmulti **A, int LDA);
int rmat_get_prec_max(int m, int n, rmulti **A, int LDA);
void rmat_get_exp(int m, int n, int *P, int LDP, rmulti **A, int LDA);
int rmat_get_exp_max(int m, int n, rmulti **A, int LDA);
int rmat_is_number(int m, int n, rmulti **A, int LDA);
int rmat_has_nan(int m, int n, rmulti **A, int LDA);

/*
 * I/O
 */
void rmat_print(int m, int n, rmulti **A, int LDA, const char *name, const char *f, int digits);
void rmat_save(int m, int n, rmulti **A, int LDA, int digits, char* fmt, ...);
void rmat_load(int m, int n, rmulti **A, int LDA, char* fmt, ...);
void rmat_bin_save(int m, int n, rmulti **A, int LDA, char* fmt, ...);
rmulti **rmat_bin_load(int *m, int *n, char* fmt, ...);

/*
 * setting
 */
void rmat_set_nan(int m, int n, rmulti **A, int LDA);
void rmat_set_s(int m, int n, rmulti **B, int LDB, char **A, int LDA);
void rmat_set_d(int m, int n, rmulti **B, int LDB, const double *A, int LDA);
void rmat_set_z(int m, int n, rmulti **B, int LDB, const dcomplex *A, int LDA);
void rmat_set_all_d(int m, int n, rmulti **A, int LDA, double a);
void rmat_set_zeros(int m, int n, rmulti **A, int LDA);
void rmat_set_ones(int m, int n, rmulti **A, int LDA);
void rmat_set_eye(int m, int n, rmulti **A, int LDA);
void rmat_set_rand(int m, int n, rmulti **A, int LDA, double a, double b);
void rmat_set_diag_r(int m, int n, rmulti **A, int LDA, rmulti *a, int offset);
void rmat_set_diag_d(int m, int n, rmulti **A, int LDA, double a, int offset);

/*
 * casting
 */
void rmat_get_d(int m, int n, double *B, int LDB, rmulti **A, int LDA);
void rmat_get_z(int m, int n, dcomplex *B, int LDB, rmulti **A, int LDA);


/*
 * rearange elements
 */
void rmat_cols_swap(int m, int n, rmulti **A, int LDA, int k, int l);
void rmat_cols_swap_index(int m, int n, rmulti **A, int LDA, const int *I);
void rmat_rows_swap(int m, int n, rmulti **A, int LDA, int k, int l);

/*
 * oparations
 */
void rmat_copy(int m, int n, rmulti **B, int LDB, rmulti **A, int LDA);                                 // B=A
void rmat_copy_t(int m, int n, rmulti **B, int LDB, rmulti **A, int LDA);                               // B=A'
void rmat_copy_col_index(int m, int n, rmulti **B, int LDB, rmulti **A, int LDA, const int *I);        // B(:,j)=A(:,I(j))
void rmat_copy_index(int m, int n, rmulti **B, int LDB, rmulti **A, int LDA, int *I, int *J);           // B=A(I,J)
void rmat_index_copy(int m, int n, rmulti **B, int LDB, rmulti **A, int LDA, int *I, int *J);           // B(I,J)=A
void rmat_neg(int m, int n, rmulti **B, int LDB, rmulti **A, int LDA);                                  // B=-A
void rmat_abs(int m, int n, rmulti **B, int LDB, rmulti **A, int LDA);                                  // B=abs(A)
void rmat_add(int m, int n, rmulti **C, int LDC, rmulti **A, int LDA, rmulti **B, int LDB);             // C=A+B
void rmat_add_r(int m, int n, rmulti **C, int LDC, rmulti **A, int LDA, rmulti *b);                     // C=A+b
void rmat_sub(int m, int n, rmulti **C, int LDC, rmulti **A, int LDA, rmulti **B, int LDB);             // C=A-B
void rmat_sub_r1(int m, int n, rmulti **C, int LDC, rmulti *a, rmulti **B, int LDB);                    // C=a-B
void rmat_sub_r2(int m, int n, rmulti **C, int LDC, rmulti **A, int LDA, rmulti *b);                    // C=A-b
void rmat_mul(int m, int n, rmulti **C, int LDC, rmulti **A, int LDA, rmulti **B, int LDB);             // C=A.*B
void rmat_mul_r(int m, int n, rmulti **C, int LDC, rmulti **A, int LDA, rmulti *b);                     // C=A*b
void rmat_mul_d(int m, int n, rmulti **C, int LDC, rmulti **A, int LDA, double b);                      // C=A*b
void rmat_prod(int l, int m, int n, rmulti **C, int LDC, rmulti **A, int LDA, rmulti **B, int LDB);     // C=A*B
void rmat_add_prod(int l, int m, int n, rmulti **C, int LDC, rmulti **A, int LDA, rmulti **B, int LDB); // C+=A*B
void rmat_sub_prod(int l, int m, int n, rmulti **C, int LDC, rmulti **A, int LDA, rmulti **B, int LDB); // C-=A*B
void rmat_div(int m, int n, rmulti **C, int LDC, rmulti **A, int LDA, rmulti **B, int LDB);             // C=A./B
void rmat_div_r1(int m, int n, rmulti **C, int LDC, rmulti *a, rmulti **B, int LDB);                    // C=a./B
void rmat_div_r2(int m, int n, rmulti **C, int LDC, rmulti **A, int LDA, rmulti *b);                    // C=A./b
void rmat_div_d(int m, int n, rmulti **C, int LDC, rmulti **A, int LDA, double b);                      // C=A./b
void rmat_rank1op(int m, int n, rmulti **B, int LDB, rmulti **A, int LDA, rmulti *a, rmulti **x, rmulti **y); // B=A+a*x*y'
void rmat_diag_copy_rvec(int m, int n, rmulti **A, int LDA, rmulti **a);                                // A=diag(a)
void rmat_diag_copy_r(int m, int n, rmulti **A, int LDA, rmulti *a);                                    // A=diag(a)
void rmat_diag_copy_d(int m, int n, rmulti **A, int LDA, double a);                                     // A=diag(a)
void rmat_diag_add_rvec(int m, int n, rmulti **B, int LDB, rmulti **A, int LDA, rmulti **a);            // B=A+diag(a)
void rmat_diag_add_r(int m, int n, rmulti **B, int LDB, rmulti **A, int LDA, rmulti *a);                // B=A+diag(a)
void rmat_diag_sub_rvec(int m, int n, rmulti **B, int LDB, rmulti **A, int LDA, rmulti **a);            // B=A-diag(a)
void rmat_diag_sub_r(int m, int n, rmulti **B, int LDB, rmulti **A, int LDA, rmulti *a);                // B=A-diag(a)
void rvec_sum_rmat(int m, int n, rmulti **B, rmulti **A, int LDA);                                      // B=sum(A)
void rvec_max_rmat(int m, int n, rmulti **B, rmulti **A, int LDA);                                      // B=max(A)
void rvec_min_rmat(int m, int n, rmulti **B, rmulti **A, int LDA);                                      // B=min(A)
void rmat_max_abs(rmulti *value, int m, int n, rmulti **A, int LDA);                                    // value=max(abs(x))
void rmat_min_abs(rmulti *value, int m, int n, rmulti **A, int LDA);                                    // value=min(abs(x))
void rmat_cols_max_abs_sub(rmulti **x, int m, int n, rmulti **A, int LDA, rmulti **B, int LDB); // x(j)=max(abs(A(:,j)-B(:,j)))
void rmat_cols_max_abs_sub_rvec(rmulti **y, int m, int n, rmulti **A, int LDA, rmulti **x); // y(j)=max(abs(A(:,j)-x))
void rmat_inv(int n, rmulti **B, int LDB, rmulti **A, int LDA);                             // B=inv(A)
void rmat_power(int n, rmulti **B, int LDB, rmulti **A, int LDA, int p);                    // B=A^p
void rmat_cols_normalize(int m, int n, rmulti **B, int LDB, rmulti **A, int LDA);
void rmat_cols_normalize_sgn(int m, int n, rmulti **B, int LDB, rmulti **A, int LDA);

/*
 * compare
 */
int rmat_cmp(int m, int n, rmulti **A, int LDA, int k, int l, rmulti **B, int LDB);

/*
 * mapping
 */
void rmat_func_list2(int m, int n, rmulti **A, int LDA, func_t *f, int l, rmulti **x); // A=F(x)

#endif
