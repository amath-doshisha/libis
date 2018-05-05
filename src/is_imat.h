#ifndef IS_IMAT_H
#define IS_IMAT_H

/*
 * allocation
 */
int* imat_allocate(int m, int n);
int* imat_free(int *A);

/*
 * initialization
 */
// A=zeros(m,n)
void imat_zeros(int m, int n, int *A, int LDA);
// A=ones(m,n)
void imat_ones(int m, int n, int *A, int LDA);
// A=ones(m,n)*a
void imat_set(int m, int n, int *A, int LDA, int a);
// A=eye(m,n); % identiy matrix
void imat_eye(int m, int n, int *A, int LDA);
// A=rand(m,n)*a+b
void imat_rand(int m, int n, int *A, int LDA, int a, int b);

/*
 * casting
 */
// y=char(x)
void imat_get_s(int m, int n, char **B, int LDB, int *A, int LDA);

/*
 * operators
 */
// B+=A
void imat_add(int m, int n, int *B, int LDB, const int *A, int LDA);
// B-=A
void imat_sub(int m, int n, int *B, int LDB, const int *A, int LDA);
// B+=a
void imat_add_scalar(int m, int n, int *B, int LDB, int a);
// B*=a
void imat_scale(int m, int n, int *B, int LDB, int a);


/*
 * convert its elements
 */
// A(:,k) <-> A(:,l)
void imat_swap_columns(int m, int n, int *A, int LDA, int k, int l);
// A(k,:) <-> A(l,:)
void imat_swap_rows(int m, int n, int *A, int LDA, int k, int l);


/*
 * convert itself to itself
 */
// B=A
void imat_copy(int m, int n, int *B, int LDB, const int *A, int LDA);
// B(:,j)=A(:,I(j)) for 0<=j<n
void imat_copy_col_index(int m, int n, int *B, int LDB, const int *A, int LDA, const int *I);
// A(:,j) <-> A(:,I[j]) for 0<=j<n
void imat_swap_index(int m, int n, int *A, int LDA, const int *I);

/*
 * get informations
 */
int imat_row_count_if(int m, int n, int *A, int LDA, int i, int value);
int imat_column_count_if(int m, int n, int *A, int LDA, int j, int value);

/*
 * output and input
 */
void imat_print(int m, int n, int *A, int LDA, char *name);
void imat_save(int m, int n, int *A, int LDA, char* fmt, ...);
int* imat_load_allocate(int *pm, int *pn, char* fmt, ...);


#endif
