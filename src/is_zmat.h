#ifndef IS_ZMAT_H
#define IS_ZMAT_H

#include<is_dcomplex.h>
/*
 * allocation
 */
dcomplex* zmat_allocate(int m, int n);
dcomplex* zmat_free(dcomplex *A);

/*
 * initialization
 */
// A=zeros(m,n)
void zmat_set_zeros(int m, int n, dcomplex *A, int LDA);
// A=ones(m,n)
void zmat_set_ones(int m, int n, dcomplex *A, int LDA);
// A=ones(m,n)*a
void zmat_set_all_z(int m, int n, dcomplex *A, int LDA, dcomplex a);
// A=ones(m,n)*a
void zmat_set_all_d(int m, int n, dcomplex *A, int LDA, double a);
// A=eye(m,n); % identiy matrix
void zmat_set_eye(int m, int n, dcomplex *A, int LDA);
// A=rand(m,n)*a+b
void zmat_set_rand(int m, int n, dcomplex *A, int LDA, double a, double b);

/*
 * casting
 */
// y=char(x)
void zmat_get_s(int m, int n, char **B, int LDB, dcomplex *A, int LDA, char format, int digits);

/*
 * convert its elements
 */
// A(:,k) <-> A(:,l)
void zmat_swap_columns(int m, int n, dcomplex *A, int LDA, int k, int l);
// A(k,:) <-> A(l,:)
void zmat_swap_rows(int m, int n, dcomplex *A, int LDA, int k, int l);

/*
 * convert itself to itself
 */
// B=A
void zmat_copy(int m, int n, dcomplex *B, int LDB, dcomplex *A, int LDA);
// B=A
void zmat_copy_d(int m, int n, dcomplex *B, int LDB, double *A, int LDA);
// B=A^T
void zmat_copy_t(int m, int n, dcomplex *B, int LDB, dcomplex *A, int LDA);
// B=A'
void zmat_copy_ct(int m, int n, dcomplex *B, int LDB, dcomplex *A, int LDA);
// B(:,j)=A(:,I(j)) for 0<=j<n
void zmat_copy_col_index(int m, int n, dcomplex *B, int LDB, dcomplex *A, int LDA, int *I);
// A(:,j) <-> A(:,I[j]) for 0<=j<n
void zmat_swap_index(int m, int n, dcomplex *A, int LDA, int *I);
// A(k,k+offset)=a for k=0,1,..
void zmat_diag_set_scalar(int m, int n, dcomplex *A, int LDA, dcomplex a, int offset);
void zmat_diag_set_scalar_d(int m, int n, dcomplex *A, int LDA, double a, int offset);
// A(k,k)+=a for k=0,1,..,n-1
void zmat_diag_add_scalar(int n, dcomplex *A, int LDA, dcomplex a);
// A(k,k)-=a for k=0,1,..,n-1
void zmat_diag_sub_scalar(int n, dcomplex *A, int LDA, dcomplex a);
// B=B+A
void zmat_add(int m, int n, dcomplex *B, int LDB, dcomplex *A, int LDA);
// B=B-A
void zmat_sub(int m, int n, dcomplex *B, int LDB, dcomplex *A, int LDA);
// C=A*B
void zmat_prod(int l, int m, int n, dcomplex *C, int LDC, dcomplex *A, int LDA, dcomplex *B, int LDB);
// C=C+A*B
void zmat_add_prod(int l, int m, int n, dcomplex *C, int LDC, dcomplex *A, int LDA, dcomplex *B, int LDB);
// C=C-A*B
void zmat_sub_prod(int l, int m, int n, dcomplex *C, int LDC, dcomplex *A, int LDA, dcomplex *B, int LDB);
// A=A+a*x*y'
void zmat_rank1op(int m, int n, dcomplex *A, int LDA, dcomplex a, dcomplex *x, dcomplex *y);
// normalize all columns
void zmat_normalize(int m, int n, dcomplex *A, int LDA);
void zmat_normalize_sgn(int m, int n, dcomplex *A, int LDA);


/*
 * convert two matrices to real vector
 */
// x(j)=max(abs(A(:,j)-B(:,j))) for j=1,2,..,n
void zmat_sub_norm_max(double *x, int m, int n, dcomplex *A, int LDA, dcomplex *B, int LDB);

/*
 * convert the matrix to another real matrix
 */
// B=real(A)
void zmat_real(int m, int n, double *B, int LDB, dcomplex *A, int LDA);

/*
 * convert itself with vector to another real vector
 */
// y(j)=max(abs(A(:,j)-x)) for j=1,2,..,n
void zmat_dist_norm_max(double *y, int m, int n, dcomplex *A, int LDA, dcomplex *x);

/*
 * output and input
 */
void zmat_print(int m, int n, dcomplex *A, int LDA, char *name, char format, int digits);
void zmat_save(int m, int n, dcomplex *A, int LDA, char* fmt, ...);
void zmat_bin_save(int m, int n, dcomplex *A, int LDA, char* fmt, ...);
dcomplex *zmat_bin_load(int *m, int *n, char* fmt, ...);

#endif
