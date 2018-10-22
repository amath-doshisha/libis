#ifndef IS_RVEC_H
#define IS_RVEC_H

#include<is_dcomplex.h>
#include<is_rmulti.h>

/**
 @file  is_rvec.h
 @brief 多倍長精度実数型rmultiのベクトルに関する関数の宣言
 */

/*
 * allocation
 */
rmulti **rvec_allocate(int n);
rmulti **rvec_allocate_prec(int n, int prec);
rmulti **rvec_allocate_clone_rvec(int n, rmulti **y);
rmulti **rvec_allocate_clone_zvec(int n, dcomplex *y);
rmulti **rvec_allocate_clone_dvec(int n, double *y);
rmulti **rvec_free(int n, rmulti **x);

/*
 * I/O
 */
void rvec_print(int n, rmulti **x, char *name, char format, int digits);
void rvec_print_bin(int n, rmulti **x, const char *name, int digits);
void rvec_print_prec(int n, rmulti **x, const char *name, const char *f, int digits);
void rvec_print_exp(int n, rmulti **x, const char *name);
void rvec_save(int n, rmulti **x, int offset, int digits, char* fmt, ...);
void rvec_save_log2_abs(int n, rmulti **x,int offset, int digits, char* fmt, ...);
void rvec_save_itrmap(int n, rmulti **x, int digits, char* fmt, ...);
void rvec_load(int n, rmulti **x, char* fmt, ...);
void rvec_bin_save(int n, rmulti **x, char* fmt, ...);
rmulti **rvec_bin_load(int *n, char* fmt, ...);

/*
 * setting
 */
void rvec_set_nan(int n, rmulti **x);                      // x=nan(n,1)
void rvec_set_inf(int n, rmulti **x, int sgn);             // x=sgn*inf(n,1)
void rvec_set_zeros(int n, rmulti **x);                    // x=zeros(n,1)
void rvec_set_ones(int n, rmulti **x);                     // x=ones(n,1)
void rvec_set_all_r(int n, rmulti **x, rmulti *a);         // x=ones(n,1)*a
void rvec_set_all_z(int n, rmulti **x, dcomplex a);        // x=ones(n,1)*a
void rvec_set_all_d(int n, rmulti **x, double a);          // x=ones(n,1)*a
void rvec_set_unit(int n, rmulti **x, int k);              // x=zeros(n,1); x[k]=1
void rvec_set_grid(int n, rmulti **x);                     // x=0:(n-1)
void rvec_set_rand(int n, rmulti **x, double a, double b); // x=rand(n,1)*a+b

/*
 * casting
 */
void rvec_set_rvec(int n, rmulti **y, rmulti **x);  // rmulti <- rmulti
void rvec_set_zvec(int n, rmulti **y, dcomplex *x); // rmulti <- dcomplex
void rvec_get_zvec(int n, dcomplex *y, rmulti **x); // rmulti -> dcomplex
void zvec_set_rvec(int n, dcomplex *y, rmulti **x); // rmulti -> dcomplex
void rvec_set_dvec(int n, rmulti **y, double *x);   // rmulti <- double
void rvec_get_dvec(int n, double *y, rmulti **x);   // rmulti -> double
void dvec_set_rvec(int n, double *y, rmulti **x);   // rmulti -> double
void rvec_set_ivec(int n, rmulti **y, int *x);      // rmulti <- int
void rvec_get_ivec(int n, int *y, rmulti **x);      // rmulti -> int
void ivec_set_rvec(int n, int *y, rmulti **x);      // rmulti -> int
void rvec_set_svec(int n, rmulti **x, char **str);                        // rmulti <- char
void rvec_get_svec(int n, char **y, rmulti **x, char format, int digits); // rmulti -> char

/*
 * convert itself
 */
void rvec_swap(int n, rmulti **x, rmulti **y); // x <=> y
void rvec_reverse(int n, rmulti **x);
void rvec_swap_at(rmulti **x, int i, int j);
void rvec_swap_index(int n, rmulti **x, const int *I);
void rvec_sort(int n, rmulti **x, int *I);
void rvec_quick_sort(int n, rmulti **x, int *I, int left, int right);
void rvec_sort_index(int *I, int n, rmulti **X);

/*
 * member variables
 */
void rvec_get_prec(int n, int *p, rmulti **x);
int  rvec_get_prec_max(int n, rmulti **x);
void rvec_get_exp(int n, int *p, rmulti **x);
int  rvec_get_exp_max(int n, rmulti **x);
int  rvec_is_number(int n, rmulti **x);
int  rvec_has_nan(int n, rmulti **x);
int  rvec_has_negative(int n, rmulti **x);

/*
 * y=f(x)
 */
void rvec_round_rvec(int n, rmulti **x, int prec);                 // x=round(x)
void rvec_clone_rvec(int n, rmulti **y, rmulti **x);               // y=x
void rvec_clone_rvec_index(int n, rmulti **y, rmulti **x, int *I); // y=x(I)
void rvec_index_clone_rvec(int n, rmulti **y, rmulti **x, int *I); // y(I)=x
void rvec_copy_rvec(int n, rmulti **y, rmulti **x);                // y=x
void rvec_copy_rvec_index(int n, rmulti **y, rmulti **x, int *I);  // y=x(I)
void rvec_index_copy_rvec(int n, rmulti **y, int *I, rmulti **x);  // y(I)=x
void rvec_abs_rvec(int n, rmulti **y, rmulti **x);                 // y=abs(x)
void rvec_neg_rvec(int n, rmulti **y, rmulti **x);                 // y=-x
void rvec_pow2_rvec(int n, rmulti **y, rmulti **x);                // y=x^2
void rvec_log2_abs_rvec(int n, rmulti **y, rmulti **x);            // y=log2(abs(x))
void rvec_normalize_rvec(int n, rmulti **y, rmulti **x);           // y=x/sqrt(x'*x)
void rvec_normalize_sgn_rvec(int n, rmulti **y, rmulti **x);       // y=x/sqrt(x'*x)
void rsum_rvec(rmulti *y, int n, rmulti **x);                      // y=sum(x)
void raverage_rvec(rmulti *y, int n, rmulti **x);                  // y=sum(x)/n
void rsum_abs_rvec(rmulti *y, int n, rmulti **x);                  // y=sum(abs(x))
void rsum_pow2_abs_rvec(rmulti *y, int n, rmulti **x);             // y=sum(abs(x)^2)
void rnorm1_rvec(rmulti *y, int n, rmulti **x);                    // y=sum(abs(x))
void rnorm2_rvec(rmulti *y, int n, rmulti **x);                    // y=sqrt(sum(abs(x)^2))
void rnorm_max_rvec(rmulti *y, int n, rmulti **x);                 // y=max(abs(x))
void rmax_rvec(rmulti *y, int n, rmulti **x);                      // y=max(x)
void rmin_rvec(rmulti *y, int n, rmulti **x);                      // y=min(x)
void rmax_abs_rvec(rmulti *y, int n, rmulti **x);                  // y=max(abs(x))
void rmin_abs_rvec(rmulti *y, int n, rmulti **x);                  // y=min(abs(x))
void rmax_abs_rvec_index(rmulti *y, int n, rmulti **x, int *I);    // y=max(abs(x))
void rmin_abs_rvec_index(rmulti *y, int n, rmulti **x, int *I);    // y=min(abs(x))
void rvec_mul_2exp_rvec(int n, rmulti **y, rmulti **x, int p);     // y=x*2^p
void rvec_div_2exp_rvec(int n, rmulti **y, rmulti **x, int p);     // y=x/2^p
void rvec_sqrt_rvec(int n, rmulti **y, rmulti **x);                // y=sqrt(x)
void rvec_log_rvec(int n, rmulti **y, rmulti **x);                 // y=log(x)
void rvec_log10_rvec(int n, rmulti **y, rmulti **x);               // y=log10(x)

/*
 * z=f(x,y)
 */
// z=x+y
void rvec_add_rvec_rvec(int n, rmulti **z, rmulti **x, rmulti **y);
void rvec_add_rvec_dvec(int n, rmulti **z, rmulti **x, double *y);
void rvec_add_dvec_rvec(int n, rmulti **z, double *x, rmulti **y);
void rvec_add_rvec_rscalar(int n, rmulti **z, rmulti **x, rmulti *y);
void rvec_add_rvec_dscalar(int n, rmulti **z, rmulti **x, double y);
void rvec_add_dvec_rscalar(int n, rmulti **z, double *x, rmulti *y);
void rvec_add_rscalar_rvec(int n, rmulti **z, rmulti *x, rmulti **y);
void rvec_add_rscalar_dvec(int n, rmulti **z, rmulti *x, double *y);
void rvec_add_dscalar_rvec(int n, rmulti **z, double x, rmulti **y);
// z=x-y
void rvec_sub_rvec_rvec(int n, rmulti **z, rmulti **x, rmulti **y);
void rvec_sub_rvec_dvec(int n, rmulti **z, rmulti **x, double *y);
void rvec_sub_dvec_rvec(int n, rmulti **z, double *x, rmulti **y);
void rvec_sub_rvec_rscalar(int n, rmulti **z, rmulti **x, rmulti *y);
void rvec_sub_rvec_dscalar(int n, rmulti **z, rmulti **x, double y);
void rvec_sub_dvec_rscalar(int n, rmulti **z, double *x, rmulti *y);
void rvec_sub_rscalar_rvec(int n, rmulti **z, rmulti *x, rmulti **y);
void rvec_sub_rscalar_dvec(int n, rmulti **z, rmulti *x, double *y);
void rvec_sub_dscalar_rvec(int n, rmulti **z, double x, rmulti **y);
// z=x*y
void rvec_mul_rvec_rvec(int n, rmulti **z, rmulti **x, rmulti **y);
void rvec_mul_rvec_dvec(int n, rmulti **z, rmulti **x, double *y);
void rvec_mul_dvec_rvec(int n, rmulti **z, double *x, rmulti **y);
void rvec_mul_rvec_rscalar(int n, rmulti **z, rmulti **x, rmulti *y);
void rvec_mul_rvec_dscalar(int n, rmulti **z, rmulti **x, double y);
void rvec_mul_dvec_rscalar(int n, rmulti **z, double *x, rmulti *y);
void rvec_mul_rscalar_rvec(int n, rmulti **z, rmulti *x, rmulti **y);
void rvec_mul_rscalar_dvec(int n, rmulti **z, rmulti *x, double *y);
void rvec_mul_dscalar_rvec(int n, rmulti **z, double x, rmulti **y);
// z=x/y
void rvec_div_rvec_rvec(int n, rmulti **z, rmulti **x, rmulti **y);
void rvec_div_rvec_dvec(int n, rmulti **z, rmulti **x, double *y);
void rvec_div_dvec_rvec(int n, rmulti **z, double *x, rmulti **y);
void rvec_div_rvec_rscalar(int n, rmulti **z, rmulti **x, rmulti *y);
void rvec_div_rvec_dscalar(int n, rmulti **z, rmulti **x, double y);
void rvec_div_dvec_rscalar(int n, rmulti **z, double *x, rmulti *y);
void rvec_div_rscalar_rvec(int n, rmulti **z, rmulti *x, rmulti **y);
void rvec_div_rscalar_dvec(int n, rmulti **z, rmulti *x, double *y);
void rvec_div_dscalar_rvec(int n, rmulti **z, double x, rmulti **y);
// z=x^y
void rvec_pow_rvec_rvec(int n, rmulti **z, rmulti **x, rmulti **y);
void rvec_pow_rvec_dvec(int n, rmulti **z, rmulti **x, double *y);
void rvec_pow_rvec_rscalar(int n, rmulti **z, rmulti **x, rmulti *y);
void rvec_pow_rvec_dscalar(int n, rmulti **z, rmulti **x, double y);
void rvec_pow_rvec_iscalar(int n, rmulti **z, rmulti **x, int y);
// z=abs(x-y)
void rvec_abs_sub_rvec_rvec(int n, rmulti **z, rmulti **x, rmulti **y);
void rvec_abs_sub_rvec_rscalar(int n, rmulti **z, rmulti **x, rmulti *y);
// z=sum(abs(x-y))
void rsum_abs_sub_rvec_rvec(rmulti *value, int n, rmulti **x, rmulti **y);
// z=max(abs(x-y))
void rmax_abs_sub_rvec_rvec(rmulti *value, int n, rmulti **x, rmulti **y);
// z=max(abs(x-y)^2)
void rmax_pow2_abs_sub_rvec_rvec(rmulti *z, int n, rmulti **x, rmulti **y);
// z=x'*y
void rinnprod_rvec_rvec(rmulti *value, int n, rmulti **x, rmulti **y);
// z=(x'*y)/sqrt(x'*x)/sqrt(y'*y)
void rdcos_rvec_rvec(rmulti *z, int n, rmulti **x, rmulti **y);
// z=abs(x'*y)/sqrt(x'*x)/sqrt(y'*y)
void rdcos_abs_rvec_rvec(rmulti *z, int n, rmulti **x, rmulti **y);
// z=acos(abs(x'*y)/sqrt(x'*x)/sqrt(y'*y))
void rangle_rvec_rvec(rmulti *theta, int n, rmulti **x, rmulti **y);
// z=(180/PI)*acos(abs(x'*y)/sqrt(x'*x)/sqrt(y'*y))
void rdeg_angle_rvec_rvec(rmulti *theta, int n, rmulti **x, rmulti **y);
// z=max(abs(x)/abs(y))
void rmax_div_abs_rvec_rvec(rmulti *z, int n, rmulti **x, rmulti **y);

/*
 * y=y+f(x)
 */
// y=y-(x'*y)*x where x'*x=1
void rvec_orthogonalize(int n, rmulti **y, rmulti **x);

/*
 * z=z+f(x,y)
 */
// z=z+x*y
void rvec_add_mul_rvec_rvec(int n, rmulti **z, rmulti **x, rmulti **y);
void rvec_add_mul_rvec_rscalar(int n, rmulti **z, rmulti **x, rmulti *y);
void rvec_add_mul_rvec_dscalar(int n, rmulti **z, rmulti **x, double y);
// z=z-x*y
void rvec_sub_mul_rvec_rvec(int n, rmulti **z, rmulti **x, rmulti **y);
void rvec_sub_mul_rvec_rscalar(int n, rmulti **z, rmulti **x, rmulti *y);
void rvec_sub_mul_rvec_dscalar(int n, rmulti **z, rmulti **x, double y);

/*
 * compare
 */
// X<=>Y
int rvec_cmp(int n, rmulti **x, int m, rmulti **y); // X<=>Y
int rvec_eq(int n, rmulti **x, rmulti **y);         // X==Y
int rvec_gt(int n, rmulti **x, rmulti **y);         // X>Y
int rvec_ge(int n, rmulti **x, rmulti **y);         // X>=Y
int rvec_lt(int n, rmulti **x, rmulti **y);         // X<Y
int rvec_le(int n, rmulti **x, rmulti **y);         // X<=Y
int rvec_eq_d(int n, rmulti **x, double a);         // x[i]==a 
int rvec_gt_d2(int n, rmulti **x, double a);        // x[i]>a 
int rvec_ge_d2(int n, rmulti **x, double a);        // x[i]>=a 
int rvec_lt_d2(int n, rmulti **x, double a);        // x[i]>a 
int rvec_le_d2(int n, rmulti **x, double a);        // x[i]<=a

/*
 * 作業中： rmat.cに移動する．
 */
// y=A*x
void rvec_mul_rmat_rvec(int m, int n, rmulti **y, rmulti **A, int LDA, rmulti **x);
void rvec_add_lintr(int m, int n, rmulti **y, rmulti **A, int LDA, rmulti **x);   // y+=A*x
void rvec_sub_lintr(int m, int n, rmulti **y, rmulti **A, int LDA, rmulti **x);   // y-=A*x
void rvec_lintr_t(int m, int n, rmulti **y, rmulti **A, int LDA, rmulti **x);     // y=A'*x
void rvec_add_lintr_t(int m, int n, rmulti **y, rmulti **A, int LDA, rmulti **x); // y+=A'*x
void rvec_sub_lintr_t(int m, int n, rmulti **y, rmulti **A, int LDA, rmulti **x); // y-=A'*x


#endif
