#ifndef IS_DMAT_H
#define IS_DMAT_H

/*
 * allocation
 */
double* dmat_allocate(int m, int n);
double* dmat_free(double *A);


/*
 * initialization
 */
// A=zeros(m,n)
void dmat_zeros(int m, int n, double *A, int LDA);
// A=ones(m,n)
void dmat_ones(int m, int n, double *A, int LDA);
// A=ones(m,n)*a
void dmat_set(int m, int n, double *A, int LDA, double a);
// A=eye(m,n); % identiy matrix
void dmat_eye(int m, int n, double *A, int LDA);
// A=rand(m,n)*a+b
void dmat_rand(int m, int n, double *A, int LDA, double a, double b);

/*
 * convert its elements
 */
// A(:,k) <-> A(:,l)
void dmat_swap_columns(int m, int n, double *A, int LDA, int k, int l);
// A(k,:) <-> A(l,:)
void dmat_swap_rows(int m, int n, double *A, int LDA, int k, int l);

/*
 * convert itself to itself
 */
// B=A
void dmat_copy(int m, int n, double *B, int LDB, const double *A, int LDA);
// B=A'
void dmat_copy_t(int m, int n, double *B, int LDB, const double *A, int LDA);
// B(:,j)=A(:,I(j)) for 0<=j<n
void dmat_copy_col_index(int m, int n, double *B, int LDB, const double *A, int LDA, const int *I);
// A(:,j) <-> A(:,I[j]) for 0<=j<n
void dmat_swap_index(int m, int n, double *A, int LDA, const int *I);
// A(k,k+offset)=a for k=0,1,..
void dmat_diag_set_scalar(int m, int n, double *A, int LDA, double a, int offset);
// A(k,k)+=a for k=0,1,..,n-1
void dmat_diag_add_scalar(int n, double *A, int LDA, double a);
// A(k,k)-=a for k=0,1,..,n-1
void dmat_diag_sub_scalar(int n, double *A, int LDA, double a);
// B=B+A
void dmat_add(int m, int n, double *B, int LDB, const double *A, int LDA);
// B=B-A
void dmat_sub(int m, int n, double *B, int LDB, const double *A, int LDA);
// C=A*B
void dmat_prod(int l, int m, int n, double *C, int LDC, const double *A, int LDA, const double *B, int LDB);
// C=C+A*B
void dmat_add_prod(int l, int m, int n, double *C, int LDC, const double *A, int LDA, const double *B, int LDB);
// C=C-A*B
void dmat_sub_prod(int l, int m, int n, double *C, int LDC, const double *A, int LDA, const double *B, int LDB);
// A=A+a*x*y'
void dmat_rank1op(int m, int n, double *A, int LDA, double a, const double *x, const double *y);
// normalize all columns
void dmat_normalize(int m, int n, double *A, int LDA);
void dmat_normalize_sgn(int m, int n, double *A, int LDA);


/*
 * convert two matrices to vector
 */
// x(j)=max(abs(A(:,j)-B(:,j))) for j=1,2,..,n
void dmat_sub_norm_max(double *x, int m, int n, const double *A, int LDA, const double *B, int LDB);

/*
 * convert itself with vector to another vector
 */
// y(j)=max(abs(A(:,j)-x)) for j=1,2,..,n
void dmat_dist_norm_max(double *y, int m, int n, const double *A, int LDA, const double *x);

/*
 * output and input
 */
void dmat_print(int m, int n, const double *A, int LDA, const char *name, const char *f, int digits);
void dmat_save(int m, int n, double *A, int LDA, char* fmt, ...);
void dmat_bin_save(int m, int n, double *A, int LDA, char* fmt, ...);
double *dmat_bin_load(int *m, int *n, char* fmt, ...);

#endif
