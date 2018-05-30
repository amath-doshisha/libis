#ifndef IS_IVEC_H
#define IS_IVEC_H

/*
 * allocation
 */
int* ivec_allocate(int n);
int* ivec_allocate_clone(int n, int *y);
int* ivec_free(int *x);

/*
 * initialization
 */
// X=zeros(n,1)
void ivec_set_zeros(int n, int *X);
// X=ones(n,1)
void ivec_set_ones(int n, int *X);
// X=ones(n,1)*a
void ivec_set_all(int n, int *X, int a);
// X[0]=0; X[1]=1; X[2]=2; ...; X[n-1]=n-1
void ivec_set_grid(int n, int *X);
// x=rand(n,1) under [a0,a1]
void ivec_set_rand(int n, int *x, int a0, int a1);
// y=x
void ivec_set_s(int n, int *y, char **x);

/*
 * casting
 */
// y=char(x)
void ivec_get_s(int n, char **y, int *x);

/*
 * convert itself to another vector
 */
// Y=X
void ivec_copy(int n, int *Y, int *X);
// Y[i]=X[I[i]], 0<=i<n
void ivec_copy_index(int n, int *Y, int *X, int *I);

/*
 * convert two vectors to themselves
 */
// x <=> y
void ivec_swap(int n, int *x, int *y);


/*
 * convert itself
 */
// y+=x
void ivec_add(int n, int *y, int *x);
void ivec_add_scalar(int n, int *y, int x);
// y-=x
void ivec_sub(int n, int *y, int *x);
// x*=a
void ivec_scale(int n, int *x, int a);

/*
 * convert its elements
 */
// X[0]=X[n-1]; X[1]=X[n-2]; X[2]=X[n-3]; ...
void ivec_reverse(int n, int *X);
// X[I[i]] => X[i]
void ivec_relocate(int n, int *X, int *I);
// X[i] <=> X[j]
void ivec_swap_at(int *X, int i, int j);
// x[i] <=> x[I[i]]
void ivec_swap_index(int n, int *x, int *I);
// sort X as X[0] <= X[1] <= X[2] <= ... <= X[n-1]
// if I==NULL, then I is not ussed.
// if I!=NULL, then I is stored with sorted indexes
void ivec_sort(int n, int *X, int *I);
// Don't call this function directly!
void ivec_quick_sort(int n, int *X, int *I, int left, int right);
// store the list of indexes, x is not destroyed
void ivec_sort_index(int *I, int n, int *X);


/*
 * convert to itself to scalar
 */
// max(X)
int ivec_max(int n, int *X);
// max(X)
int ivec_min(int n, int *X);
// sum(x)
int ivec_sum(int n, int *x);
// sum(x)
double ivec_average(int n, int *x);

/*
 * convert to boolean
 */
// X==Y
int ivec_eq(int n, int *X, int *Y);
// X!=Y
int ivec_ne(int n, int *X, int *Y);
// X<=>Y
int ivec_cmp(int n, int *X, int m, int *Y);
// x[i]>y for all i
int ivec_all_gt(int n, int *x, int y);
// x[i]>=y for all i
int ivec_all_ge(int n, int *x, int y);
// x[i]<y for all i
int ivec_all_lt(int n, int *x, int y);
// x[i]<=y for all i
int ivec_all_le(int n, int *x, int y);


/*
 * index operations
 */
void ivec_copy_index(int n, int *Y, int *X, int *I);
int ivec_count_if(int n, int *X, int value);
int ivec_count_if_not(int n, int *X, int value);
int* ivec_allocate_index_if(int n, int *X, int value, int *pm);
int* ivec_allocate_index_if_not(int n, int *X, int value, int *pm);
int ivec_first_index_if(int n, int *X, int value);
int ivec_first_index_if_not(int n, int *X, int value);


/*
 * input and output
 */
void ivec_put(int n, int *x, char *sep);
void ivec_print(int n, int *x, char *name);
void ivec_save(int n, int *x, int offset, char *fmt, ...);

#endif
