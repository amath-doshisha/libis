#ifndef IS_CMAT_H
#define IS_CMAT_H

#include<is_cmulti.h>
#include<is_func.h>

/**
 @file  is_cmat.h
 @brief 多倍長精度実数型cmultiの行列に関する関数の宣言
 */

/*
 * constructions and destruction
 */
cmulti **cmat_allocate(int LDA, int n);
cmulti **cmat_allocate_prec(int LDA, int n, int prec);
cmulti **cmat_allocate_clone(int LDB, int n, cmulti **B);
cmulti **cmat_allocate_clone_r(int LDB, int n, rmulti **B);
cmulti **cmat_free(int LDA, int n, cmulti **A);
void cmat_round(int m, int n, cmulti **A, int LDA, int prec);
void cmat_clone(int m, int n, cmulti **B, int LDB, cmulti **A, int LDA);
void cmat_clone_r(int m, int n, cmulti **B, int LDB, rmulti **A, int LDA);
void cmat_clone_t(int m, int n, cmulti **B, int LDB, cmulti **A, int LDA);
void cmat_swap(int m, int n, cmulti **A, int LDA, cmulti **B, int LDB);

/*
 * casting
 */
void cmat_get_s(int m, int n, char **B, int LDB, cmulti **A, int LDA, char format, int digits);

/*
 * member variables
 */
void cmat_get_prec(int m, int n, int *P, int LDP, cmulti **A, int LDA);
int cmat_get_prec_max(int m, int n, cmulti **A, int LDA);
void cmat_get_exp(int m, int n, int *P, int LDP, cmulti **A, int LDA);
int cmat_get_exp_max(int m, int n, cmulti **A, int LDA);
int cmat_is_number(int m, int n, cmulti **A, int LDA);
int cmat_has_nan(int m, int n, cmulti **A, int LDA);
int cmat_is_real(int m, int n, cmulti **A, int LDA);

/*
 * I/O
 */
void cmat_print(int m, int n, cmulti **A, int LDA, char *name, char format, int digits);
void cmat_load(int m, int n, cmulti **A, int LDA, char* fmt, ...);
void cmat_save(int m, int n, cmulti **A, int LDA, int digits, char* fmt, ...);
void cmat_bin_save(int m, int n, cmulti **A, int LDA, char* fmt, ...);
cmulti **cmat_bin_load(int *m, int *n, char* fmt, ...);

/*
 * setting
 */
void cmat_set_nan(int m, int n, cmulti **A, int LDA);
void cmat_set_z(int m, int n, cmulti **B, int LDB, const dcomplex *A, int LDA);
void cmat_set_d(int m, int n, cmulti **B, int LDB, const double *A, int LDA);
void cmat_set_dd(int m, int n, cmulti **B, int LDB, const double *Ar, int LDAr, const double *Ai, int LDAi);
void cmat_set_all_z(int m, int n, cmulti **A, int LDA, dcomplex a);
void cmat_set_all_dd(int m, int n, cmulti **A, int LDA, double a_r, double a_i);
void cmat_set_all_d(int m, int n, cmulti **A, int LDA, double a);
void cmat_set_zeros(int m, int n, cmulti **A, int LDA);
void cmat_set_ones(int m, int n, cmulti **A, int LDA);
void cmat_set_eye(int m, int n, cmulti **A, int LDA);
void cmat_set_rand(int m, int n, cmulti **A, int LDA, double a, double b);
void cmat_set_diag_c(int m, int n, cmulti **A, int LDA, cmulti *a, int offset);
void cmat_set_diag_d(int m, int n, cmulti **A, int LDA, double a, int offset);

/*
 * casting
 */
void cmat_get_z(int m, int n, dcomplex *B, int LDB, cmulti **A, int LDA);
void cmat_get_d(int m, int n, double *B, int LDB, cmulti **A, int LDA);

/*
 * rearange elements
 */
void cmat_cols_swap(int m, int n, cmulti **A, int LDA, int k, int l);
void cmat_cols_swap_index(int m, int n, cmulti **A, int LDA, const int *I);
void cmat_swap_rows(int m, int n, cmulti **A, int LDA, int k, int l);
void cmat_cols_rotate_right(int m, int n, cmulti **A, int LDA); // rotate columns
void cmat_cols_rotate_left(int m, int n, cmulti **A, int LDA);  // rotate columns

/*
 * oparations
 */
void cmat_copy(int m, int n, cmulti **B, int LDB, cmulti **A, int LDA);                            // B=A
void cmat_copy_rmat(int m, int n, cmulti **B, int LDB, rmulti **A, int LDA);                       // B=A
void cmat_copy_col_index(int m, int n, cmulti **B, int LDB, cmulti **A, int LDA, const int *I);   // B=A(I)
void cmat_copy_index(int m, int n, cmulti **B, int LDB, cmulti **A, int LDA, int *I, int *J);      // B=A(I,J)
void cmat_index_copy(int m, int n, cmulti **B, int LDB, cmulti **A, int LDA, int *I, int *J);      // B(I,J)=A
void cmat_index_copy_rmat(int m, int n, cmulti **B, int LDB, rmulti **A, int LDA, int *I, int *J); // B(I,J)=A
void cmat_copy_t(int m, int n, cmulti **B, int LDB, cmulti **A, int LDA);                         // B=A^T
void cmat_copy_rmat_t(int m, int n, cmulti **B, int LDB, rmulti **A, int LDA);                    // B=A'
void cmat_copy_ct(int m, int n, cmulti **B, int LDB, cmulti **A, int LDA);                        // B=A'
void cmat_real(int m, int n, rmulti **B, int LDB, cmulti **A, int LDA);                           // B=real(A)
void cmat_real_clone(int m, int n, rmulti **B, int LDB, cmulti **A, int LDA);                     // B=real(A)
void cmat_imag(int m, int n, rmulti **B, int LDB, cmulti **A, int LDA);                           // B=imag(A)
void cmat_imag_clone(int m, int n, rmulti **B, int LDB, cmulti **A, int LDA);                     // B=imag(A)
void cmat_conj(int m, int n, cmulti **B, int LDB, cmulti **A, int LDA);                           // B=conj(A)
void cmat_neg(int m, int n, cmulti **B, int LDB, cmulti **A, int LDA);                            // B=-A
void cmat_abs2(int m, int n, rmulti **B, int LDB, cmulti **A, int LDA);                           // B=abs(A).^2
void cmat_add(int m, int n, cmulti **C, int LDC, cmulti **A, int LDA, cmulti **B, int LDB);       // C=A+B
void cmat_add_rmat(int m, int n, cmulti **C, int LDC, cmulti **A, int LDA, rmulti **B, int LDB);  // C=A+B
void cmat_add_c(int m, int n, cmulti **C, int LDC, cmulti **A, int LDA, cmulti *b);               // C=A+b
void cmat_sub(int m, int n, cmulti **C, int LDC, cmulti **A, int LDA, cmulti **B, int LDB);       // C=A-B
void cmat_sub_rmat1(int m, int n, cmulti **C, int LDC, rmulti **A, int LDA, cmulti **B, int LDB); // C=A-B
void cmat_sub_rmat2(int m, int n, cmulti **C, int LDC, cmulti **A, int LDA, rmulti **B, int LDB); // C=A-B
void cmat_sub_c1(int m, int n, cmulti **C, int LDC, cmulti *a, cmulti **B, int LDB);              // C=a-B
void cmat_sub_c2(int m, int n, cmulti **C, int LDC, cmulti **A, int LDA, cmulti *b);              // C=A-b
void cmat_mul(int m, int n, cmulti **C, int LDC, cmulti **A, int LDA, cmulti **B, int LDB);       // C=A.*B
void cmat_mul_c(int m, int n, cmulti **C, int LDC, cmulti **A, int LDA, cmulti *b);               // C=A*b
void cmat_mul_r(int m, int n, cmulti **C, int LDC, cmulti **A, int LDA, rmulti *b);               // C=A*b
void cmat_mul_z(int m, int n, cmulti **C, int LDC, cmulti **A, int LDA, dcomplex b);              // C=A*b
void cmat_mul_d(int m, int n, cmulti **C, int LDC, cmulti **A, int LDA, double b);                // C=A*b
void cmat_prod(int l, int m, int n, cmulti **C, int LDC, cmulti **A, int LDA, cmulti **B, int LDB);     // C=A*B
void cmat_prod_r1(int l, int m, int n, cmulti **C, int LDC, rmulti **A, int LDA, cmulti **B, int LDB);  // C=A*B
void cmat_prod_r2(int l, int m, int n, cmulti **C, int LDC, cmulti **A, int LDA, rmulti **B, int LDB);  // C=A*B
void cmat_add_prod(int l, int m, int n, cmulti **C, int LDC, cmulti **A, int LDA, cmulti **B, int LDB); // C+=A*B
void cmat_sub_prod(int l, int m, int n, cmulti **C, int LDC, cmulti **A, int LDA, cmulti **B, int LDB); // C-=A*
void cmat_div(int m, int n, cmulti **C, int LDC, cmulti **A, int LDA, cmulti **B, int LDB);       // C=A./*B
void cmat_div_c1(int m, int n, cmulti **C, int LDC, cmulti *a, cmulti **B, int LDB);              // C=a./*B
void cmat_div_c2(int m, int n, cmulti **C, int LDC, cmulti **A, int LDA, cmulti *b);              // C=A./*b
void cmat_rank1op(int m, int n, cmulti **B, int LDB, cmulti **A, int LDA, cmulti *a, cmulti **x, cmulti **y); // B=A+a*x*y'
void cmat_diag_copy_cvec(int m, int n, cmulti **A, int LDA, cmulti **a);                       // A=diag(a)
void cmat_diag_copy_rvec(int m, int n, cmulti **A, int LDA, rmulti **a);                       // A=diag(a)
void cmat_diag_copy_c(int m, int n, cmulti **A, int LDA, cmulti *a);                           // A=diag(a)
void cmat_diag_copy_r(int m, int n, cmulti **A, int LDA, rmulti *a);                           // A=diag(a)
void cmat_diag_add_cvec(int m, int n, cmulti **B, int LDB, cmulti **A, int LDA, cmulti **a);   // B=A+diag(a).
void cmat_diag_add_rvec(int m, int n, cmulti **B, int LDB, cmulti **A, int LDA, rmulti **a);   // B=A+diag(a).
void cmat_diag_add_c(int m, int n, cmulti **B, int LDB, cmulti **A, int LDA, cmulti *a);       // B=A+diag(a).
void cmat_diag_add_r(int m, int n, cmulti **B, int LDB, cmulti **A, int LDA, rmulti *a);       // B=A+diag(a).
void cmat_diag_sub_cvec(int m, int n, cmulti **B, int LDB, cmulti **A, int LDA, cmulti **a);   // B=A-diag(a).
void cmat_diag_sub_rvec(int m, int n, cmulti **B, int LDB, cmulti **A, int LDA, rmulti **a);   // B=A-diag(a).
void cmat_diag_sub_c(int m, int n, cmulti **B, int LDB, cmulti **A, int LDA, cmulti *a);       // B=A-diag(a).
void cmat_diag_sub_r(int m, int n, cmulti **B, int LDB, cmulti **A, int LDA, rmulti *a);       // B=A-diag(a).
void cmat_inv(int n, cmulti **B, int LDB, cmulti **A, int LDA);                             // B=inv(A)
void cmat_power(int n, cmulti **B, int LDB, cmulti **A, int LDA, int p);                    // B=A^p
void cmat_abs(int m, int n, rmulti **B, int LDB, cmulti **A, int LDA);                      // B=abs(A)
void cmat_arg(int m, int n, rmulti **B, int LDB, cmulti **A, int LDA);                      // B=arg(A)
void cvec_sum_cmat(int m, int n, cmulti **B, cmulti **A, int LDA);                          // B=sum(A)
void cvec_max_cmat(int m, int n, cmulti **B, cmulti **A, int LDA);                          // B=max(A)
void cvec_min_cmat(int m, int n, cmulti **B, cmulti **A, int LDA);                          // B=min(A)
void cmat_max_abs(rmulti *value, int m, int n, cmulti **A, int LDA);
void cmat_max_absc(rmulti *value, int m, int n, cmulti **A, int LDA);                       // value=max(abs(x))
void cmat_min_abs(rmulti *value, int m, int n, cmulti **A, int LDA);                        // value=min(abs(x))
void cmat_cols_normalize(int m, int n, cmulti **B, int LDB, cmulti **A, int LDA);
void cmat_cols_normalize_sgn(int m, int n, cmulti **B, int LDB, cmulti **A, int LDA);
void cmat_cols_max_abs_sub(rmulti **x, int m, int n, cmulti **A, int LDA, cmulti **B, int LDB); // x(j)=max(abs(A(:,j)-B(:,j))) for j=1,2,..,n
void cmat_cols_max_abs_sub_cvec(rmulti **y, int m, int n, cmulti **A, int LDA, cmulti **x);     // y(j)=max(abs(A(:,j)-x)) for j=1,2,..,n

/*
 * compare
 */
int cmat_cmp(int m, int n, cmulti **A, int LDA, int k, int l, cmulti **B, int LDB);

/*
 * mapping
 */
void cmat_func_list2(int m, int n, cmulti **A, int LDA, func_t *f, int l, cmulti **x);

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void cmat_angle_deg_list(rmulti **angle, int m, int n, cmulti **A, int LDA);

#endif
