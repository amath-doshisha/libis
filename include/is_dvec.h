#ifndef IS_DVEC_H
#define IS_DVEC_H

/**
 @file  is_rvec.h
 @brief 倍精度型doubleのベクトルに関する関数の宣言
 */

/*
 * allocation
 */
double *dvec_allocate(int n);
double *dvec_allocate_clone_dvec(int n, double *y);
double* dvec_free(double *x);

/*
 * I/O
 */
void dvec_print(int n, double *x, char *name, char format, int digits);
void dvec_load(int n, double *x, char* fmt, ...);
void dvec_save(int n, double *x, int offset, char* fmt, ...);
void dvec_save2(int n, double *x, double *y, char* fmt, ...);
void dvec_save_log2_abs(int n, double *x, int offset, char* fmt, ...);
void dvec_bin_save(int n, double *x, char* fmt, ...);
double *dvec_bin_load(int *n, char* fmt, ...);

/*
 * setting
 */
void dvec_set_nan(int n, double *x);                      // x=nan(n,1)
void dvec_set_inf(int n, double *x, int sgn);             // x=sgn*inf(n,1)
void dvec_set_zeros(int n, double *x);                    // x=zeros(n,1)
void dvec_set_ones(int n, double *x);                     // x=ones(n,1)
void dvec_set_all_d(int n, double *x, double a);          // x=ones(n,1)*a
void dvec_set_unit(int n, double *x, int k);              // x=zeros(n,1); x[k]=1
void dvec_set_grid(int n, double *x);                     // x=0:(n-1)
void dvec_set_rand(int n, double *x, double a, double b); // x=rand(n,1)*a+b

/*
 * casting
 */
void dvec_set_dvec(int n, double *y, double *x); // double <- double
void dvec_set_ivec(int n, double *y, int *x);    // double <- int
void dvec_get_ivec(int n, int *y, double *x);    // double -> int
void ivec_set_dvec(int n, int *y, double *x);    // double -> int
void dvec_set_svec(int n, double *y, char **x);                          // double <- char
void dvec_get_svec(int n, char **y, double *x, char format, int digits); // double -> char

/*
 * convert itself
 */
void dvec_swap(int n, double *x, double *y);    // x <=> y
void dvec_reverse(int n, double *x);            // x[0]=x[n-1]; x[1]=x[n-2]; x[2]=x[n-3]; ...
void dvec_swap_at(double *x, int i, int j);     // x[i] <=> x[j]
void dvec_swap_index(int n, double *x, int *I); // x[i] <=> x[I[i]]
void dvec_sort(int n, double *x, int *I);       // sort x as x[0] <= x[1] <= x[2] <= ... <= x[n-1]. If I==NULL, then I is not ussed. If I!=NULL, then I is stored with sorted indexes.
void dvec_quick_sort(int n, double *x, int *I, int left, int right); // Don't call this function directly!
void dvec_sort_index(int *I, int n, double *X); // store the list of indexes, x is not destroyed

/*
 * y=f(x)
 */
void dvec_copy_dvec(int n, double *y, double *x);                // y=x
void dvec_copy_dvec_index(int n, double *y, double *x, int *I);  // y=x(I)
void dvec_index_copy_dvec(int n, double *y, int *I, double *x);  // y(I)=x
void dvec_abs_dvec(int n, double *y, double *x);                 // y=abs(x)
void dvec_neg_dvec(int n, double *y, double *x);                 // y=-x
void dvec_pow2_dvec(int n, double *y, double *x);                // y=x.^2
void dvec_log2_abs_dvec(int n, double *y, double *x);            // y=log2(abs(x))
void dvec_normalize_dvec(int n, double *y, double *x);           // y=x/sqrt(x'*x)
void dvec_normalize_sgn_dvec(int n, double *y, double *x);       // y=x/sqrt(x'*x) where x[i]>0, abs(x[i]) is max
double dsum_dvec(int n, double *x);                              // y=sum(x)
double daverage_dvec(int n, double *x);                          // y=sum(x)/n
double dsum_abs_dvec(int n, double *x);                          // y=sum(abs(x))
double dsum_pow2_abs_dvec(int n, double *x);                     // y=sum(abs(x)^2)
double dnorm1_dvec(int n, double *x);                            // y=sum(abs(x))
double dnorm2_dvec(int n, double *x);                            // y=sqrt(sum(abs(x)^2))
double dnorm_max_dvec(int n, double *x);                         // y=max(abs(x))
double dmax_dvec(int n, double *x);                              // y=max(x)
double dmin_dvec(int n, double *x);                              // y=min(x)
double dmax_abs_dvec(int n, double *x);                          // y=max(abs(x))
double dmin_abs_dvec(int n, double *x);                          // y=min(abs(x))
double dmax_abs_dvec_index(int n, double *x, int *I);            // y=max(abs(x)), (*I)=k where k=max{ k | abs(x[k]) } unless I==NULL
double dmin_abs_dvec_index(int n, double *x, int *I);            // y=min(abs(x)), (*I)=k where k=min{ k | abs(x[k]) } unless I==NULL

/*
 * z=f(x,y)
 */
// z=x+y
void dvec_add_dvec_dvec(int n, double *z, double *x, double *y);
void dvec_add_dvec_dscalar(int n, double *z, double *x, double y);
void dvec_add_dscalar_dvec(int n, double *z, double x, double *y);
// z=x-y
void dvec_sub_dvec_dvec(int n, double *z, double *x, double *y);
void dvec_sub_dvec_dscalar(int n, double *z, double *x, double y);
void dvec_sub_dscalar_dvec(int n, double *z, double x, double *y);
// z=x*y
void dvec_mul_dvec_dvec(int n, double *z, double *x, double *y);
void dvec_mul_dvec_dscalar(int n, double *z, double *x, double y);
void dvec_mul_dscalar_dvec(int n, double *z, double x, double *y);
// z=x/y
void dvec_div_dvec_dvec(int n, double *z, double *x, double *y);
void dvec_div_dvec_dscalar(int n, double *z, double *x, double y);
void dvec_div_dscalar_dvec(int n, double *z, double x, double *y);
// z=x^y
double pow_i(double x, int n);
void dvec_pow_dvec_dvec(int n, double *z, double *x, double *y);
void dvec_pow_dvec_dscalar(int n, double *z, double *x, double y);
void dvec_pow_dvec_iscalar(int n, double *z, double *x, int y);
// z=abs(x-y)
void dvec_abs_sub_dvec_dvec(int n, double *z, double *x, double *y);
void dvec_abs_sub_dvec_dscalar(int n, double *z, double *x, double y);
// z=sum(abs(x-y))
double dsum_abs_sub_dvec_dvec(int n, double *x, double *y);
// z=max(abs(x-y))
double dmax_abs_sub_dvec_dvec(int n, double *x, double *y);
// z=max(abs(x-y)^2)
double dmax_pow2_abs_sub_dvec_dvec(int n, double *x, double *y);
// z=x'*y
double dinnprod_dvec_dvec(int n, double *x, double *y);
// z=(x'*y)/\sqrt(x'*x)/sqrt(y'*y)
double ddcos_dvec_dvec(int n, double *x, double *y);
// z=abs(x'*y)/\sqrt(x'*x)/sqrt(y'*y)
double ddcos_abs_dvec_dvec(int n, double *x, double *y);
// others
// z=acos(abs(x'*y)/\sqrt(x'*x)/sqrt(y'*y))
double dangle_dvec_dvec(int n, double *x, double *y);
// z=(180/PI)*acos(abs(x'*y)/\sqrt(x'*x)/sqrt(y'*y))
double ddeg_angle_dvec_dvec(int n, double *x, double *y);
// z=max(abs(x)/abs(y))
double dmax_div_abs_dvec_dvec(int n, double *x, double *y);

/*
 * y=y+f(x)
 */
void dvec_add_scaled(int n, double *y, double a, double *x);                   // y=y+a*x
void dvec_sub_scaled(int n, double *y, double a, double *x);                   // y=y-a*x
void dvec_orthogonalize(int n, double *y, double *x);                          // y=y-(x'*y)*x where x'*x=1


#endif
