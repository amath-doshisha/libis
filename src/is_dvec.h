#ifndef IS_DVEC_H
#define IS_DVEC_H

/*
 * allocation
 */
double* dvec_allocate(int n);
double* dvec_free(double *x);
double *dvec_allocate_s(char *str, int *n);

/*
 * initialization
 */
// y=x
void dvec_set_d(int n, double *y, double *x);
// x=nan(n,1)
void dvec_set_nan(int n, double *x);
// x=inf(n,1)
void dvec_set_inf(int n, double *x, int sgn);
// x=zeros(n,1)
void dvec_set_zeros(int n, double *x);
// x=ones(n,1)
void dvec_set_ones(int n, double *x);
// x=ones(n,1)*a
void dvec_set_all(int n, double *x, double a);
// x=zeros(n,1); x[k]=1
void dvec_set_unit(int n, double *x, int k);
// x[0]=0; x[1]=1; x[2]=2; ...; x[n-1]=n-1
void dvec_set_grid(int n, double *x);
// x=rand(n,1)*a+b
void dvec_set_rand(int n, double *x, double a, double b);

/*
 * casting
 */
// int <-> double
void dvec_get_si(int n, int *y, double *x);
void dvec_set_si(int n, double *y, int *x);
// char <-> double
void dvec_get_s(int n, char **y, double *x, char format, int digits);
void dvec_set_s(int n, double *y, char **x);

/*
 * convert its elements
 */
// x[0]=x[n-1]; x[1]=x[n-2]; x[2]=x[n-3]; ...
void dvec_reverse(int n, double *x);
// x[i] <=> x[j]
void dvec_swap_at(double *x, int i, int j);
// x[i] <=> x[I[i]]
void dvec_swap_index(int n, double *x, int *I);
// sort x as x[0] <= x[1] <= x[2] <= ... <= x[n-1]
// if I==NULL, then I is not ussed.
// if I!=NULL, then I is stored with sorted indexes
void dvec_sort(int n, double *x, int *I);
// Don't call this function directly!
void dvec_quick_sort(int n, double *x, int *I, int left, int right);
// store the list of indexes, x is not destroyed
void dvec_sort_index(int *I, int n, double *X);

/*
 * y=f(x)
 */
// y=x
void dvec_copy(int n, double *y, double *x);
// Y[i]=X[I[i]], 0<=i<n
void dvec_copy_index(int n, double *Y, double *X, int *I);


/*
 * z=f(x,y)
 */
// z=x+y
void dvec_add_dvec(int n, double *z, double *x, double *y);
void dvec_add_dscalar(int n, double *z, double *x, double y);
void dscalar_add_dvec(int n, double *z, double x, double *y);
// z=x-y
void dvec_sub_dvec(int n, double *z, double *x, double *y);
void dvec_sub_dscalar(int n, double *z, double *x, double y);
void dscalar_sub_dvec(int n, double *z, double x, double *y);
// z=x*y
void dvec_mul_dvec(int n, double *z, double *x, double *y);
void dvec_mul_dscalar(int n, double *z, double *x, double y);
void dscalar_mul_dvec(int n, double *z, double x, double *y);
// z=x/y
void dvec_div_dvec(int n, double *z, double *x, double *y);
void dvec_div_dscalar(int n, double *z, double *x, double y);
void dscalar_div_dvec(int n, double *z, double x, double *y);


/*
 * others
 */


// y=y+a*x
void dvec_add_scaled(int n, double *y, double a, double *x);
// y=y-a*x
void dvec_sub_scaled(int n, double *y, double a, double *x);
// y=A*x
void dvec_lintr(int m, int n, double *y, double *A, int LDA, double *x);
// y=y+A*x
void dvec_add_lintr(int m, int n, double *y, double *A, int LDA, double *x);
// y=y-A*x
void dvec_sub_lintr(int m, int n, double *y, double *A, int LDA, double *x);
// y=A'*x
void dvec_lintr_t(int m, int n, double *y, double *A, int LDA, double *x);
// y=y+A'*x
void dvec_add_lintr_t(int m, int n, double *y, double *A, int LDA, double *x);
// y=y-A'*x
void dvec_sub_lintr_t(int m, int n, double *y, double *A, int LDA, double *x);
// x=-x
void dvec_neg(int n, double *x);
// x=a*x
void dvec_scale(int n, double *x, double a);
// x=x/sqrt(x'*x)
void dvec_normalize(int n, double *x);
// x=x/sqrt(x'*x) where x[i]>0, abs(x[i]) is max
void dvec_normalize_sgn(int n, double *x);
// y=y-(x'*y)*x where x'*x=1
void dvec_orthogonalize(int n, double *y, double *x);
// z=x.^y
void dvec_pow(int n, double *z, double *x, double *y);
double pow_si(double x, int n);
// z=x^y
void dvec_pow_si(int n, double *z, double *x, long y);
void dvec_pow_d(int n, double *z, double *x, double y);
// y=x.^2
void dvec_pow2(int n, double *y, double *x);
// y=log2(abs(x))
void dvec_log2_abs(int n, double *y, double *x);

/*
 * convert two vectors to themselves
 */
// x <=> y
void dvec_swap(int n, double *x, double *y);

/*
 * convert two vectors to scalar
 */
// x'*y
double dvec_dot(int n, double *x, double *y);
// abs(x'*y)/\sqrt(x'*x)/sqrt(y'*y)
double dvec_dcos_abs(int n, double *x, double *y);
// acos(abs(x'*y)/\sqrt(x'*x)/sqrt(y'*y))
double dvec_angle(int n, double *x, double *y);
// (180/PI)*acos(abs(x'*y)/\sqrt(x'*x)/sqrt(y'*y))
double dvec_angle_deg(int n, double *x, double *y);
// value=sum(abs(x-y))
double dvec_dist_norm1(int n, double *x, double *y);
// max(abs(x-y))
double dvec_dist_norm_max(int n, double *x, double *y);

/*
 * convert two vectors to another vector
 */
// z=abs(x-y)
void dvec_abs_sub(int n, double *z, double *x, double *y);
void dvec_abs_sub_scalar(int n, double *z, double *x, double y);

/*
 * convert itself to scalar
 */
// sum(abs(x))
double dvec_norm1(int n, double *x);
// sqrt(sum(x.^2))
double dvec_norm2(int n, double *x);
// max(abs(x))
double dvec_norm_max(int n, double *x);
// sum(x)
double dvec_sum(int n, double *x);
// sum(x)/n
double dvec_average(int n, double *x);
// sum(x.^2)
double dvec_sum_pow2(int n, double *x);
// max(abs(x))
double dvec_max_abs(int n, double *x);
// min(abs(x))
double dvec_min_abs(int n, double *x);
// max(abs(x))
// (*I)=k where k=max{ k | abs(x[k]) } unless I==NULL
double dvec_max_abs_index(int n, double *x, int *I);
// min(abs(x))
// (*I)=k where k=min{ k | abs(x[k]) } unless I==NULL
double dvec_min_abs_index(int n, double *x, int *I);
// max(x)
double dvec_max(int n, double *x);
// min(x)
double dvec_min(int n, double *x);


/*
 * input and output
 */
void dvec_print(int n, double *x, char *name, char format, int digits);
void dvec_load(int n, double *x, char* fmt, ...);
void dvec_save(int n, double *x, int offset, char* fmt, ...);
void dvec_save2(int n, double *x, double *y, char* fmt, ...);
void dvec_save_log2_abs(int n, double *x, int offset, char* fmt, ...);
void dvec_bin_save(int n, double *x, char* fmt, ...);
double *dvec_bin_load(int *n, char* fmt, ...);

#endif
