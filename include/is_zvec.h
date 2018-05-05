#ifndef IS_ZVEC_H
#define IS_ZVEC_H

#include<is_dcomplex.h>

/*
 * allocation
 */
dcomplex* zvec_allocate(int n);
dcomplex* zvec_free(dcomplex *x);

/*
 * initialization
 */
// x=nan(n,1)
void zvec_set_nan(int n, dcomplex *x);
// x=inf(n,1)
void zvec_set_inf(int n, dcomplex *x, int rsgn, int isgn);
// x=zeros(n,1)
void zvec_set_zeros(int n, dcomplex *x);
// x=ones(n,1)
void zvec_set_ones(int n, dcomplex *x);
// x=ones(n,1)*a
void zvec_set_all(int n, dcomplex *x, dcomplex a);
// x=ones(n,1)*a
void zvec_set_all_d(int n, dcomplex *x, double a);
// x=ones(n,1)*a
void zvec_set_all_dd(int n, dcomplex *x, double a_r, double a_i);
// x=zeros(n,1); x[k]=1
void zvec_set_unit(int n, dcomplex *x, int k);
// x[0]=0; x[1]=1; x[2]=2; ...; x[n-1]=n-1
void zvec_set_grid(int n, dcomplex *x);
// x=rand(n,1)*a+b
void zvec_set_rand(int n, dcomplex *x, double a, double b);
// y=x
void zvec_set_s(int n, dcomplex *y, char **x);

/*
 * casting
 */
// y=int(x)
void zvec_get_si(int n, int *y, dcomplex *x);
// y=char(x)
void zvec_get_s(int n, char **y, dcomplex *x, char format, int digits);


/*
 * convert its elements
 */
// x[0]=x[n-1]; x[1]=x[n-2]; x[2]=x[n-3]; ...
void zvec_reverse(int n, dcomplex *x);
// x[i] <=> x[j]
void zvec_swap_at(dcomplex *x, int i, int j);
// x[i] <=> x[I[i]]
void zvec_swap_index(int n, dcomplex *x, const int *I);
// sort x as x[0] <= x[1] <= x[2] <= ... <= x[n-1]
// if I==NULL, then I is not ussed.
// if I!=NULL, then I is stored with sorted indexes
void zvec_sort(int n, dcomplex *x, int *I);
// Don't call this function directly!
void zvec_quick_sort(int n, dcomplex *x, int *I, int left, int right);
// store the list of indexes, x is not destroyed
void zvec_sort_index(int *I, int n, dcomplex *X);

/*
 * convert itself to itself
 */
// y=x
void zvec_copy(int n, dcomplex *y, const dcomplex *x);
// y=x
void zvec_copy_d(int n, dcomplex *y, const double *x);
// y=x
void zvec_copy_dd(int n, dcomplex *y, const double *x_r, const double *x_i);
// Y[i]=X[I[i]], 0<=i<n
void zvec_copy_index(int n, dcomplex *Y, const dcomplex *X, const int *I);
// x=x+a
void zvec_add_scalar(int n, dcomplex *x, dcomplex a);
// x=x-a
void zvec_sub_scalar(int n, dcomplex *x, dcomplex a);
// y=y+x
void zvec_add(int n, dcomplex *y, const dcomplex *x);
// y=y-x
void zvec_sub(int n, dcomplex *y, const dcomplex *x);
// y=y+a*x
void zvec_add_scaled(int n, dcomplex *y, dcomplex a, const dcomplex *x);
// y=y-a*x
void zvec_sub_scaled(int n, dcomplex *y, dcomplex a, const dcomplex *x);
// y=A*x
void zvec_lintr(int m, int n, dcomplex *y, const dcomplex *A, int LDA, const dcomplex *x);
// y=y+A*x
void zvec_add_lintr(int m, int n, dcomplex *y, const dcomplex *A, int LDA, const dcomplex *x);
// y=y-A*x
void zvec_sub_lintr(int m, int n, dcomplex *y, const dcomplex *A, int LDA, const dcomplex *x);
// y=A^T*x
void zvec_lintr_t(int m, int n, dcomplex *y, const dcomplex *A, int LDA, const dcomplex *x);
// y=y+A^T*x
void zvec_add_lintr_t(int m, int n, dcomplex *y, const dcomplex *A, int LDA, const dcomplex *x);
// y=y-A^T*x
void zvec_sub_lintr_t(int m, int n, dcomplex *y, const dcomplex *A, int LDA, const dcomplex *x);
// y=A'*x
void zvec_lintr_ct(int m, int n, dcomplex *y, const dcomplex *A, int LDA, const dcomplex *x);
// y=y+A'*x
void zvec_add_lintr_ct(int m, int n, dcomplex *y, const dcomplex *A, int LDA, const dcomplex *x);
// y=y-A'*x
void zvec_sub_lintr_ct(int m, int n, dcomplex *y, const dcomplex *A, int LDA, const dcomplex *x);
// x=-x
void zvec_neg(int n, dcomplex *x);
// x=a*x
void zvec_scale(int n, dcomplex *x, dcomplex a);
// x=a*x
void zvec_scale_d(int n, dcomplex *x, double a);
// x=(-a)*x
void zvec_scale_neg(int n, dcomplex *x, dcomplex a);
// x=x/a
void zvec_scale_div(int n, dcomplex *x, dcomplex a);
// x=x/sqrt(x'*x)
void zvec_normalize(int n, dcomplex *x);
// x=x/sqrt(x'*x) where x[i]/|x[i]|=1, abs(x[i]) is max
void zvec_normalize_sgn(int n, dcomplex *x);
// y=y-(x'*y)*x where x'*x=1
void zvec_orthogonalize(int n, dcomplex *y, const dcomplex *x);

/*
 * convert two vectors to themselves
 */
// x <=> y
void zvec_swap(int n, dcomplex *x, dcomplex *y);

/*
 * convert two vectors to complex scalar
 */
// x'*y
dcomplex zvec_dot(int n, const dcomplex *x, const dcomplex *y);

/*
 * convert two vectors to real scalar
 */
// abs(x'*y)/\sqrt(x'*x)/sqrt(y'*y)
double zvec_dcos_abs(int n, const dcomplex *x, const dcomplex *y);
// acos(abs(x'*y)/\sqrt(x'*x)/sqrt(y'*y))
double zvec_angle(int n, const dcomplex *x, const dcomplex *y);
// (180/PI)*acos(abs(x'*y)/\sqrt(x'*x)/sqrt(y'*y))
double zvec_angle_deg(int n, const dcomplex *x, const dcomplex *y);
// value=sum(abs(x-y))
double zvec_dist_norm1(int n, const dcomplex *x, const dcomplex *y);
// max(|x-y|)
double zvec_dist_norm_max(int n, const dcomplex *x, const dcomplex *y);

/*
 * convert the vector to another real vector
 */
// y=real(x)
void zvec_real(int n, double *y, dcomplex *x);

/*
 * convert two vectors to another real vector
 */
// z=abs(x-y)
void zvec_abs_sub(int n, double *z, const dcomplex *x, const dcomplex *y);
void zvec_abs_sub_scalar(int n, double *z, const dcomplex *x, dcomplex y);

/*
 * convert itself to real scalar
 */
// sum(abs(x))
double zvec_norm1(int n, const dcomplex *x);
// sqrt(sum(abs(x).^2))
double zvec_norm2(int n, const dcomplex *x);
// max(abs(x))
double zvec_norm_max(int n, const dcomplex *x);
// sum(abs(x).^2)
double zvec_sum_pow2_abs(int n, const dcomplex *x);
// max(abs(x))
double zvec_max_abs(int n, const dcomplex *x);
// min(abs(x))
double zvec_min_abs(int n, const dcomplex *x);
// max(abs(x))
// (*I)=k where k=max{ k | abs(x[k]) } unless I==NULL
double zvec_max_abs_index(int n, const dcomplex *x, int *I);
// min(abs(x))
// (*I)=k where k=min{ k | abs(x[k]) } unless I==NULL
double zvec_min_abs_index(int n, const dcomplex *x, int *I);


/*
 * convert itself to complex scalar
 */
// max(x)
dcomplex zvec_max(int n, const dcomplex *x);
// min(x)
dcomplex zvec_min(int n, const dcomplex *x);
// sum(x)
dcomplex zvec_sum(int n, const dcomplex *x);
// sum(x)/n
dcomplex zvec_average(int n, const dcomplex *x);

/*
 * input and output
 */
void zvec_print(int n, dcomplex *x, char *name, char format, int digits);
void zvec_save(int n, dcomplex *x, int offset, char* fmt, ...);
void zvec_save_cplane(int n, dcomplex *x, char* fmt, ...);
void zvec_save_cplane_label(int n, dcomplex *x, int offset, int *label, char* fmt, ...);
void zvec_bin_save(int n, dcomplex *x, char* fmt, ...);
dcomplex *zvec_bin_load(int *n, char* fmt, ...);

#endif
