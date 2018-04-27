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
// x=zeros(n,1)
void dvec_zeros(int n, double *x);
// x=ones(n,1)
void dvec_ones(int n, double *x);
// x=ones(n,1)*a
void dvec_set(int n, double *x, double a);
// x=zeros(n,1); x[k]=1
void dvec_unit(int n, double *x, int k);
// x[0]=0; x[1]=1; x[2]=2; ...; x[n-1]=n-1
void dvec_set_grid(int n, double *x);
// x=rand(n,1)*a+b
void dvec_rand(int n, double *x, double a, double b);

/*
 * convert its elements
 */
// x[0]=x[n-1]; x[1]=x[n-2]; x[2]=x[n-3]; ...
void dvec_reverse(int n, double *x);
// x[i] <=> x[j]
void dvec_swap_at(double *x, int i, int j);
// x[i] <=> x[I[i]]
void dvec_swap_index(int n, double *x, const int *I);
// sort x as x[0] <= x[1] <= x[2] <= ... <= x[n-1]
// if I==NULL, then I is not ussed.
// if I!=NULL, then I is stored with sorted indexes
void dvec_sort(int n, double *x, int *I);
// Don't call this function directly!
void dvec_quick_sort(int n, double *x, int *I, int left, int right);
// store the list of indexes, x is not destroyed
void dvec_sort_index(int *I, int n, double *X);

/*
 * convert itself to itself
 */
void dvec_add_d(int n, double *z, const double *x, double y);  // z=x+y
void dvec_div_d2(int n, double *z, const double *x, double y); // z=x/y

// y=x
void dvec_copy(int n, double *y, const double *x);
// Y[i]=X[I[i]], 0<=i<n
void dvec_copy_index(int n, double *Y, const double *X, const int *I);
// x=x+a
void dvec_add_scalar(int n, double *x, double a);
// x=x-a
void dvec_sub_scalar(int n, double *x, double a);
// y=y+x
void dvec_add(int n, double *y, const double *x);
// y=y-x
void dvec_sub(int n, double *y, const double *x);
// y=y+a*x
void dvec_add_scaled(int n, double *y, double a, const double *x);
// y=y-a*x
void dvec_sub_scaled(int n, double *y, double a, const double *x);
// y=A*x
void dvec_lintr(int m, int n, double *y, const double *A, int LDA, const double *x);
// y=y+A*x
void dvec_add_lintr(int m, int n, double *y, const double *A, int LDA, const double *x);
// y=y-A*x
void dvec_sub_lintr(int m, int n, double *y, const double *A, int LDA, const double *x);
// y=A'*x
void dvec_lintr_t(int m, int n, double *y, const double *A, int LDA, const double *x);
// y=y+A'*x
void dvec_add_lintr_t(int m, int n, double *y, const double *A, int LDA, const double *x);
// y=y-A'*x
void dvec_sub_lintr_t(int m, int n, double *y, const double *A, int LDA, const double *x);
// x=-x
void dvec_neg(int n, double *x);
// x=a*x
void dvec_scale(int n, double *x, double a);
// x=x/sqrt(x'*x)
void dvec_normalize(int n, double *x);
// x=x/sqrt(x'*x) where x[i]>0, abs(x[i]) is max
void dvec_normalize_sgn(int n, double *x);
// y=y-(x'*y)*x where x'*x=1
void dvec_orthogonalize(int n, double *y, const double *x);
// z=x.^y
void dvec_pow(int n, double *z, double *x, double *y);
double pow_si(double x, int n);
// z=x^y
void dvec_pow_si(int n, double *z, double *x, long y);
void dvec_pow_d(int n, double *z, double *x, double y);


// y=x.^2
void dvec_pow2(int n, double *y, const double *x);
// y=log2(abs(x))
void dvec_log2_abs(int n, double *y, const double *x);

/*
 * convert two vectors to themselves
 */
// x <=> y
void dvec_swap(int n, double *x, double *y);

/*
 * convert two vectors to scalar
 */
// x'*y
double dvec_dot(int n, const double *x, const double *y);
// abs(x'*y)/\sqrt(x'*x)/sqrt(y'*y)
double dvec_dcos_abs(int n, const double *x, const double *y);
// acos(abs(x'*y)/\sqrt(x'*x)/sqrt(y'*y))
double dvec_angle(int n, const double *x, const double *y);
// (180/PI)*acos(abs(x'*y)/\sqrt(x'*x)/sqrt(y'*y))
double dvec_angle_deg(int n, const double *x, const double *y);
// value=sum(abs(x-y))
double dvec_dist_norm1(int n, const double *x, const double *y);
// max(abs(x-y))
double dvec_dist_norm_max(int n, const double *x, const double *y);

/*
 * convert two vectors to another vector
 */
// z=abs(x-y)
void dvec_abs_sub(int n, double *z, const double *x, const double *y);
void dvec_abs_sub_scalar(int n, double *z, const double *x, double y);

/*
 * convert itself to scalar
 */
// sum(abs(x))
double dvec_norm1(int n, const double *x);
// sqrt(sum(x.^2))
double dvec_norm2(int n, const double *x);
// max(abs(x))
double dvec_norm_max(int n, const double *x);
// sum(x)
double dvec_sum(int n, const double *x);
// sum(x)/n
double dvec_average(int n, const double *x);
// sum(x.^2)
double dvec_sum_pow2(int n, const double *x);
// max(abs(x))
double dvec_max_abs(int n, const double *x);
// min(abs(x))
double dvec_min_abs(int n, const double *x);
// max(abs(x))
// (*I)=k where k=max{ k | abs(x[k]) } unless I==NULL
double dvec_max_abs_index(int n, const double *x, int *I);
// min(abs(x))
// (*I)=k where k=min{ k | abs(x[k]) } unless I==NULL
double dvec_min_abs_index(int n, const double *x, int *I);
// max(x)
double dvec_max(int n, const double *x);
// min(x)
double dvec_min(int n, const double *x);


/*
 * input and output
 */
void dvec_print(int n, const double *x, const char *name, const char *format, int digits);
void dvec_load(int n, double *x, char* fmt, ...);
void dvec_save(int n, double *x, int offset, char* fmt, ...);
void dvec_save2(int n, double *x, double *y, char* fmt, ...);
void dvec_save_log2_abs(int n, double *x, int offset, char* fmt, ...);
void dvec_bin_save(int n, double *x, char* fmt, ...);
double *dvec_bin_load(int *n, char* fmt, ...);

#endif
