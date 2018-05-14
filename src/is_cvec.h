#ifndef IS_CVEC_H
#define IS_CVEC_H

#include<is_cmulti.h>

/**
 @file  is_cvec.h
 @brief 多倍長精度複素数型cmultiのベクトルに関する関数の宣言.
*/

/*
 * allocation
 */
cmulti **cvec_allocate(int n);
cmulti **cvec_allocate_prec(int n, int prec);
cmulti **cvec_allocate_clone_cvec(int n, cmulti **y);
cmulti **cvec_allocate_clone_rvec(int n, rmulti **y);
cmulti **cvec_allocate_clone_zvec(int n, dcomplex *y);
cmulti **cvec_allocate_clone_dvec(int n, double *y);
cmulti **cvec_free(int n, cmulti **x);

/*
 * I/O
 */
void cvec_print(int n, cmulti **x, char *name, char format, int digits);
void cvec_print_prec(int n, cmulti **x, const char *name, const  char *format, int digits);
void cvec_print_exp(int n, cmulti **x, const char *name);
void cvec_save(int n, cmulti **x, int offset, int digits, char* fmt, ...);
void cvec_save_cplane(int n, cmulti **x, int digits, char* fmt, ...);
void cvec_load(int n, cmulti **x, char* fmt, ...);
void cvec_bin_save(int n, cmulti **x, char* fmt, ...);
cmulti **cvec_bin_load(int *n, char* fmt, ...);


/*
 * setting
 */
void cvec_set_nan(int n, cmulti **x);                                    // x=nan(n,1)
void cvec_set_inf(int n, cmulti **x, int rsgn, int isgn);                // x=sgn*inf(n,1)
void cvec_set_zeros(int n, cmulti **x);                                  // x=zeros(n,1)
void cvec_set_ones(int n, cmulti **x);                                   // x=ones(n,1)
void cvec_set_all_c(int n, cmulti **x, cmulti *a);                       // x=ones(n,1)*a
void cvec_set_all_r(int n, cmulti **x, rmulti *a);                       // x=ones(n,1)*a
void cvec_set_all_z(int n, cmulti **x, dcomplex a);                      // x=ones(n,1)*a
void cvec_set_all_dd(int n, cmulti **x, double a_r, double a_i);         // x=ones(n,1)*a
void cvec_set_all_d(int n, cmulti **x, double a);                        // x=ones(n,1)*a
void cvec_set_unit(int n, cmulti **x, int k);                            // x=zeros(n,1); x[k]=1
void cvec_set_grid(int n, cmulti **x);                                   // x=0:(n-1)
void cvec_set_rand(int n, cmulti **x, double a, double b);               // x=rand(n,1)*a+b

/*
 * casting
 */
void cvec_set_cvec(int n, cmulti **y, cmulti **x);                      // cmulti <- cmulti
void cvec_set_rvec_rvec(int n, cmulti **y, rmulti **x_r, rmulti **x_i); // cmulti <- (rmulti,rmulti)
void cvec_set_rvec(int n, cmulti **y, rmulti **x);                      // cmulti <- rmulti
void cvec_get_rvec(int n, rmulti **y, cmulti **x);                      // cmulti -> rmulti
void rvec_set_cvec(int n, rmulti **y, cmulti **x);                      // cmulti -> rmulti
void cvec_set_zvec(int n, cmulti **y, dcomplex *x);                     // cmulti <- dcomplex
void cvec_get_zvec(int n, dcomplex *y, cmulti **x);                     // cmulti -> dcomplex
void zvec_set_cvec(int n, dcomplex *y, cmulti **x);                     // cmulti -> dcomplex
void cvec_set_dvec_dvec(int n, cmulti **y, double *x_r, double *x_i);   // cmulti <- (double,double)
void cvec_set_dvec(int n, cmulti **y, double *x);                       // cmulti <- double
void cvec_get_dvec(int n, double *y, cmulti **x);                       // cmulti -> double
void dvec_set_cvec(int n, double *y, cmulti **x);                       // cmulti -> double
void cvec_set_ivec(int n, cmulti **y, int *x);                          // cmulti <- int
void cvec_get_ivec(int n, int *y, cmulti **x);                          // cmulti -> int
void ivec_set_cvec(int n, int *y, cmulti **x);                          // cmulti -> int
void cvec_set_svec(int n, cmulti **x, char **str);                        // cmulti <- char
void cvec_get_svec(int n, char **y, cmulti **x, char format, int digits); // cmulti -> char

/*
 * member variables
 */
int cvec_get_prec_max(int n, cmulti **x);
void cvec_get_prec(int n, int *p, cmulti **x);
int cvec_get_exp_max(int n, cmulti **x);
void cvec_get_exp(int n, int *p, cmulti **x);
int cvec_is_number(int n, cmulti **x);
int cvec_has_nan(int n, cmulti **x);
int cvec_is_real(int n, cmulti **x);

/*
 * rearange elements
 */
void cvec_swap(int n, cmulti **x, cmulti **y);
void cvec_reverse(int n, cmulti **x);
void cvec_swap_at(cmulti **x, int i, int j);
void cvec_swap_index(int n, cmulti **x, const int *I);
void cvec_sort(int n, cmulti **x, int *I);
void cvec_quick_sort(int n, cmulti **x, int *I, int left, int right);
void cvec_sort_index(int *I, int n, cmulti **X);

/*
 * y=f(x)
 */
void cvec_round(int n, cmulti **x, int prec);
void cvec_clone(int n, cmulti **y, cmulti **x);
void cvec_clone_r(int n, cmulti **y, rmulti **x);
void cvec_clone_rr(int n, cmulti **y, rmulti **x_r, rmulti **x_i);
void cvec_clone_index(int n, cmulti **y, cmulti **x, int *I); // y=x(I)
void cvec_index_clone(int n, cmulti **y, cmulti **x, int *I); // y(I)=x
void cvec_copy(int n, cmulti **y, cmulti **x);                     // y=x
void cvec_copy_index(int n, cmulti **y, cmulti **x, int *I);       // y=x(I)
void cvec_index_copy(int n, cmulti **y, int *I, cmulti **x);       // y(I)=x
void cvec_index_copy_rvec(int n, cmulti **y, int *I, rmulti **x);  // y(I)=x
void cvec_neg(int n, cmulti **y, cmulti **x);                      // y=-x
void cvec_real(int n, rmulti **y, cmulti **x);                     // y=real(x)
void cvec_real_clone(int n, rmulti **y, cmulti **x);               // y=real(x)
void cvec_imag(int n, rmulti **y, cmulti **x);                     // y=imag(x)
void cvec_imag_clone(int n, rmulti **y, cmulti **x);               // y=imag(x)
void cvec_conj(int n, cmulti **y, cmulti **x);                     // y=conj(x)
void cvec_absc(int n, cmulti **y, cmulti **x);                     // y=abs(real(x)+i*imag(x))
void cvec_abs(int n, rmulti **y, cmulti **x);                      // y=abs(x)
void cvec_sum_abs(rmulti *value, int n, cmulti **x);               // value=sum(abs(x))
void cvec_max_abs(rmulti *value, int n, cmulti **x);               // value=max(abs(x))
void cvec_max_abs_index(rmulti *value, int n, cmulti **x, int *I); // value=max(abs(x))
void cvec_min_abs(rmulti *value, int n, cmulti **x);               // value=min(abs(x))
void cvec_min_abs_index(rmulti *value, int n, cmulti **x, int *I); // value=min(abs(x))
void cvec_max_absc(rmulti *value, int n, cmulti **x);              // value=max(abs(x))
void cvec_arg(int n, rmulti **y, cmulti **x);                      // y=arg(x)
void cvec_sum(cmulti *value, int n, cmulti **x);                   // value=sum(x)
void cvec_sum_abs2(rmulti *value, int n, cmulti **x);              // value=sum(abs(x).^2)
void cvec_max(cmulti *value, int n, cmulti **x);                   // value=max(x)
void cvec_min(cmulti *value, int n, cmulti **x);                   // value=min(x)
void cvec_log2_abs(int n, rmulti **y, cmulti **x);                 // y=log2(abs(x))
void cvec_normalize(int n, cmulti **y, cmulti **x);                // y=x/sqrt(x'*x)
void cvec_normalize_sgn(int n, cmulti **y, cmulti **x);            // y=x/sqrt(x'*x)
void cvec_norm2(rmulti *value, int n, cmulti **x);                 // value=sqrt(sum(abs(x).^2))
void cvec_average(cmulti *value, int n, cmulti **x);               // value=sum(x)/n

/*
 * z=f(x,y)
 */
// z=x+y
void cvec_add_cvec(int n, cmulti **z, cmulti **x, cmulti **y);
void cvec_add_rvec(int n, cmulti **z, cmulti **x, rmulti **y);
void rvec_add_cvec(int n, cmulti **z, rmulti **x, cmulti **y);
void cvec_add_zvec(int n, cmulti **z, cmulti **x, dcomplex *y);
void zvec_add_cvec(int n, cmulti **z, dcomplex *x, cmulti **y);
void cvec_add_dvec(int n, cmulti **z, cmulti **x, double *y);
void dvec_add_cvec(int n, cmulti **z, double *x, cmulti **y);
void rvec_add_zvec(int n, cmulti **z, rmulti **x, dcomplex *y);
void zvec_add_rvec(int n, cmulti **z, dcomplex *x, rmulti **y);
void cvec_add_cscalar(int n, cmulti **z, cmulti **x, cmulti *y);
void cvec_add_rscalar(int n, cmulti **z, cmulti **x, rmulti *y);
void rvec_add_cscalar(int n, cmulti **z, rmulti **x, cmulti *y);
void cvec_add_zscalar(int n, cmulti **z, cmulti **x, dcomplex y);
void zvec_add_cscalar(int n, cmulti **z, dcomplex *x, cmulti *y);
void cvec_add_dscalar(int n, cmulti **z, cmulti **x, double y);
void dvec_add_cscalar(int n, cmulti **z, double *x, cmulti *y);
void rvec_add_zscalar(int n, cmulti **z, rmulti **x, dcomplex y);
void zvec_add_rscalar(int n, cmulti **z, dcomplex *x, rmulti *y);
void cscalar_add_cvec(int n, cmulti **z, cmulti *x, cmulti **y);
void cscalar_add_rvec(int n, cmulti **z, cmulti *x, rmulti **y);
void rscalar_add_cvec(int n, cmulti **z, rmulti *x, cmulti **y);
void cscalar_add_zvec(int n, cmulti **z, cmulti *x, dcomplex *y);
void zscalar_add_cvec(int n, cmulti **z, dcomplex x, cmulti **y);
void cscalar_add_dvec(int n, cmulti **z, cmulti *x, double *y);
void dscalar_add_cvec(int n, cmulti **z, double x, cmulti **y);
void rscalar_add_zvec(int n, cmulti **z, rmulti *x, dcomplex *y);
void zscalar_add_rvec(int n, cmulti **z, dcomplex x, rmulti **y);
// z=x-y
void cvec_sub_cvec(int n, cmulti **z, cmulti **x, cmulti **y);
void cvec_sub_rvec(int n, cmulti **z, cmulti **x, rmulti **y);
void rvec_sub_cvec(int n, cmulti **z, rmulti **x, cmulti **y);
void cvec_sub_zvec(int n, cmulti **z, cmulti **x, dcomplex *y);
void zvec_sub_cvec(int n, cmulti **z, dcomplex *x, cmulti **y);
void cvec_sub_dvec(int n, cmulti **z, cmulti **x, double *y);
void dvec_sub_cvec(int n, cmulti **z, double *x, cmulti **y);
void rvec_sub_zvec(int n, cmulti **z, rmulti **x, dcomplex *y);
void zvec_sub_rvec(int n, cmulti **z, dcomplex *x, rmulti **y);
void cvec_sub_cscalar(int n, cmulti **z, cmulti **x, cmulti *y);
void cvec_sub_rscalar(int n, cmulti **z, cmulti **x, rmulti *y);
void rvec_sub_cscalar(int n, cmulti **z, rmulti **x, cmulti *y);
void cvec_sub_zscalar(int n, cmulti **z, cmulti **x, dcomplex y);
void zvec_sub_cscalar(int n, cmulti **z, dcomplex *x, cmulti *y);
void cvec_sub_dscalar(int n, cmulti **z, cmulti **x, double y);
void dvec_sub_cscalar(int n, cmulti **z, double *x, cmulti *y);
void rvec_sub_zscalar(int n, cmulti **z, rmulti **x, dcomplex y);
void zvec_sub_rscalar(int n, cmulti **z, dcomplex *x, rmulti *y);
void cscalar_sub_cvec(int n, cmulti **z, cmulti *x, cmulti **y);
void cscalar_sub_rvec(int n, cmulti **z, cmulti *x, rmulti **y);
void rscalar_sub_cvec(int n, cmulti **z, rmulti *x, cmulti **y);
void cscalar_sub_zvec(int n, cmulti **z, cmulti *x, dcomplex *y);
void zscalar_sub_cvec(int n, cmulti **z, dcomplex x, cmulti **y);
void cscalar_sub_dvec(int n, cmulti **z, cmulti *x, double *y);
void dscalar_sub_cvec(int n, cmulti **z, double x, cmulti **y);
void rscalar_sub_zvec(int n, cmulti **z, rmulti *x, dcomplex *y);
void zscalar_sub_rvec(int n, cmulti **z, dcomplex x, rmulti **y);
// z=x.*y
void cvec_mul_cvec(int n, cmulti **z, cmulti **x, cmulti **y);
void cvec_mul_rvec(int n, cmulti **z, cmulti **x, rmulti **y);
void rvec_mul_cvec(int n, cmulti **z, rmulti **x, cmulti **y);
void cvec_mul_zvec(int n, cmulti **z, cmulti **x, dcomplex *y);
void zvec_mul_cvec(int n, cmulti **z, dcomplex *x, cmulti **y);
void cvec_mul_dvec(int n, cmulti **z, cmulti **x, double *y);
void dvec_mul_cvec(int n, cmulti **z, double *x, cmulti **y);
void rvec_mul_zvec(int n, cmulti **z, rmulti **x, dcomplex *y);
void zvec_mul_rvec(int n, cmulti **z, dcomplex *x, rmulti **y);
void cvec_mul_cscalar(int n, cmulti **z, cmulti **x, cmulti *y);
void cvec_mul_rscalar(int n, cmulti **z, cmulti **x, rmulti *y);
void rvec_mul_cscalar(int n, cmulti **z, rmulti **x, cmulti *y);
void cvec_mul_zscalar(int n, cmulti **z, cmulti **x, dcomplex y);
void zvec_mul_cscalar(int n, cmulti **z, dcomplex *x, cmulti *y);
void cvec_mul_dscalar(int n, cmulti **z, cmulti **x, double y);
void dvec_mul_cscalar(int n, cmulti **z, double *x, cmulti *y);
void rvec_mul_zscalar(int n, cmulti **z, rmulti **x, dcomplex y);
void zvec_mul_rscalar(int n, cmulti **z, dcomplex *x, rmulti *y);
void cscalar_mul_cvec(int n, cmulti **z, cmulti *x, cmulti **y);
void cscalar_mul_rvec(int n, cmulti **z, cmulti *x, rmulti **y);
void rscalar_mul_cvec(int n, cmulti **z, rmulti *x, cmulti **y);
void cscalar_mul_zvec(int n, cmulti **z, cmulti *x, dcomplex *y);
void zscalar_mul_cvec(int n, cmulti **z, dcomplex x, cmulti **y);
void cscalar_mul_dvec(int n, cmulti **z, cmulti *x, double *y);
void dscalar_mul_cvec(int n, cmulti **z, double x, cmulti **y);
void rscalar_mul_zvec(int n, cmulti **z, rmulti *x, dcomplex *y);
void zscalar_mul_rvec(int n, cmulti **z, dcomplex x, rmulti **y);
// z=x./y
void cvec_div_cvec(int n, cmulti **z, cmulti **x, cmulti **y);
void cvec_div_rvec(int n, cmulti **z, cmulti **x, rmulti **y);
void rvec_div_cvec(int n, cmulti **z, rmulti **x, cmulti **y);
void cvec_div_zvec(int n, cmulti **z, cmulti **x, dcomplex *y);
void zvec_div_cvec(int n, cmulti **z, dcomplex *x, cmulti **y);
void cvec_div_dvec(int n, cmulti **z, cmulti **x, double *y);
void dvec_div_cvec(int n, cmulti **z, double *x, cmulti **y);
void rvec_div_zvec(int n, cmulti **z, rmulti **x, dcomplex *y);
void zvec_div_rvec(int n, cmulti **z, dcomplex *x, rmulti **y);
void cvec_div_cscalar(int n, cmulti **z, cmulti **x, cmulti *y);
void cvec_div_rscalar(int n, cmulti **z, cmulti **x, rmulti *y);
void rvec_div_cscalar(int n, cmulti **z, rmulti **x, cmulti *y);
void cvec_div_zscalar(int n, cmulti **z, cmulti **x, dcomplex y);
void zvec_div_cscalar(int n, cmulti **z, dcomplex *x, cmulti *y);
void cvec_div_dscalar(int n, cmulti **z, cmulti **x, double y);
void dvec_div_cscalar(int n, cmulti **z, double *x, cmulti *y);
void rvec_div_zscalar(int n, cmulti **z, rmulti **x, dcomplex y);
void zvec_div_rscalar(int n, cmulti **z, dcomplex *x, rmulti *y);
void cscalar_div_cvec(int n, cmulti **z, cmulti *x, cmulti **y);
void cscalar_div_rvec(int n, cmulti **z, cmulti *x, rmulti **y);
void rscalar_div_cvec(int n, cmulti **z, rmulti *x, cmulti **y);
void cscalar_div_zvec(int n, cmulti **z, cmulti *x, dcomplex *y);
void zscalar_div_cvec(int n, cmulti **z, dcomplex x, cmulti **y);
void cscalar_div_dvec(int n, cmulti **z, cmulti *x, double *y);
void dscalar_div_cvec(int n, cmulti **z, double x, cmulti **y);
void rscalar_div_zvec(int n, cmulti **z, rmulti *x, dcomplex *y);
void zscalar_div_rvec(int n, cmulti **z, dcomplex x, rmulti **y);
// others
void cvec_mul_2exp(int n, cmulti **y, cmulti **x, int pr, int pi); // y=x*2^p
void cvec_div_2exp(int n, cmulti **y, cmulti **x, int pr, int pi); // y=x/2^p
void cvec_dot(int n, cmulti **z, cmulti **x, cmulti **y);          // z=conj(x).*y
void cvec_dot_c1(int n, cmulti **z, cmulti *x, cmulti **y);        // z=conj(x)*y
void cvec_dot_c2(int n, cmulti **z, cmulti **x, cmulti *y);        // z=conj(x)*y
void cvec_dot_z1(int n, cmulti **z, dcomplex x, cmulti **y);       // z=conj(x)*y
void cvec_dot_z2(int n, cmulti **z, cmulti **x, dcomplex y);       // z=conj(x)*y
void cvec_pow(int n, cmulti **z, cmulti **x, cmulti **y);          // z=x.^y
void cvec_pow_si(int n, cmulti **z, cmulti **x, int y);            // z=x^y
void cvec_pow_r(int n, cmulti **z, cmulti **x, rmulti *y);         // z=x^y
void cvec_abs_sub(int n, rmulti **z, cmulti **x, cmulti **y);         // z=abs(x-y)
void cvec_abs_sub_c(int n, rmulti **z, cmulti **x, cmulti *y);        // z=abs(x-y)
void cvec_abs_sub_r(int n, rmulti **z, cmulti **x, rmulti *y);        // z=abs(x-y)
void cvec_sum_abs_sub(rmulti *value, int n, cmulti **x, cmulti **y);  // value=sum(abs(x-y))
void cvec_max_abs_sub(rmulti *value, int n, cmulti **x, cmulti **y);  // value=max(abs(x-y))
void cvec_dcos(cmulti *value, int n, cmulti **x, cmulti **y);         // value=x'*y/\sqrt(x'*x)/sqrt(y'*y)
void cvec_abs_dcos(rmulti *value, int n, cmulti **x, cmulti **y);     // value=abs(x'*y)/\sqrt(x'*x)/sqrt(y'*y)
void cvec_angle(rmulti *theta, int n, cmulti **x, cmulti **y);        // theta=acos(abs(x'*y)/\sqrt(x'*x)/sqrt(y'*y))
void cvec_angle_deg(rmulti *theta, int n, cmulti **x, cmulti **y);    // theta=(180/PI)*acos(abs(x'*y)/\sqrt(x'*x)/sqrt(y'*y))
void cvec_sum_dot(cmulti *value, int n, cmulti **x, cmulti **y);   // value=sum(conj(x).*y)
void cvec_max_abs2_sub(rmulti *value, int n, cmulti **x, cmulti **y); // value=max(abs(x-y).^2)


/*
 * y=y+f(x)
 */
void cvec_orthogonalize(int n, cmulti **y, cmulti **x); // y=y-(x'*y)*x where x'*x=1


/*
 * z=z+f(x,y)
 */
// z=z+x.*y
void cvec_add_mul_cvec_cvec(int n, cmulti **z, cmulti **x, cmulti **y);
void cvec_add_mul_cvec_cscalar(int n, cmulti **z, cmulti **x, cmulti *y);
void cvec_add_mul_cvec_rscalar(int n, cmulti **z, cmulti **x, rmulti *y);
void cvec_add_mul_cvec_zscalar(int n, cmulti **z, cmulti **x, dcomplex y);
void cvec_add_mul_cvec_dscalar(int n, cmulti **z, cmulti **x, double y);
// z=z-x.*y
void cvec_sub_mul_cvec_cvec(int n, cmulti **z, cmulti **x, cmulti **y);
void cvec_sub_mul_cvec_cscalar(int n, cmulti **z, cmulti **x, cmulti *y);
void cvec_sub_mul_cvec_zscalar(int n, cmulti **z, cmulti **x, dcomplex y);
void cvec_sub_mul_cvec_dscalar(int n, cmulti **z, cmulti **x, double y);
// z=z+conj(x).*y
void cvec_add_dot_cvec_cvec(int n, cmulti **z, cmulti **x, cmulti **y);
void cvec_add_dot_cvec_cscalar(int n, cmulti **z, cmulti **x, cmulti *y);
void cvec_add_dot_cvec_zscalar(int n, cmulti **z, cmulti **x, dcomplex y);
void cvec_add_dot_cscalar_cvec(int n, cmulti **z, cmulti *x, cmulti **y);
void cvec_add_dot_zscalar_cvec(int n, cmulti **z, dcomplex x, cmulti **y);
// z=z-conj(x).*y
void cvec_sub_dot_cvec_cvec(int n, cmulti **z, cmulti **x, cmulti **y);
void cvec_sub_dot_cvec_cscalar(int n, cmulti **z, cmulti **x, cmulti *y);
void cvec_sub_dot_cvec_zscalar(int n, cmulti **z, cmulti **x, dcomplex y);
void cvec_sub_dot_cscalar_cvec(int n, cmulti **z, cmulti *x, cmulti **y);
void cvec_sub_dot_zscalar_cvec(int n, cmulti **z, dcomplex x, cmulti **y);

/*
 * compare
 */
int cvec_cmp(int n, cmulti **x, int m, cmulti **y); // X<=>Y
int cvec_eq(int n, cmulti **x, cmulti **y);         // X==Y
int cvec_gt(int n, cmulti **x, cmulti **y);         // X>Y
int cvec_ge(int n, cmulti **x, cmulti **y);         // X>=Y
int cvec_lt(int n, cmulti **x, cmulti **y);         // X<Y
int cvec_le(int n, cmulti **x, cmulti **y);         // X<=Y
int cvec_eqc(int n, cmulti **x, cmulti **y);        // X==Y
int cvec_gtc(int n, cmulti **x, cmulti **y);        // X>Y
int cvec_gec(int n, cmulti **x, cmulti **y);        // X>=Y
int cvec_ltc(int n, cmulti **x, cmulti **y);        // X<Y
int cvec_lec(int n, cmulti **x, cmulti **y);        // X<=Y

/*
 * 作業中：cmat.cに移動する
 */
void cvec_lintr(int m, int n, cmulti **y, cmulti **A, int LDA, cmulti **x);        // y=A*x
void cvec_add_lintr(int m, int n, cmulti **y, cmulti **A, int LDA, cmulti **x);    // y+=A*x
void cvec_sub_lintr(int m, int n, cmulti **y, cmulti **A, int LDA, cmulti **x);    // y-=A*x
void cvec_lintr_t(int m, int n, cmulti **y, cmulti **A, int LDA, cmulti **x);      // y=A^T*x
void cvec_add_lintr_t(int m, int n, cmulti **y, cmulti **A, int LDA, cmulti **x);  // y+=A^T*x
void cvec_sub_lintr_t(int m, int n, cmulti **y, cmulti **A, int LDA, cmulti **x);  // y-=A^T*x
void cvec_lintr_ct(int m, int n, cmulti **y, cmulti **A, int LDA, cmulti **x);     // y=A'*x
void cvec_add_lintr_ct(int m, int n, cmulti **y, cmulti **A, int LDA, cmulti **x); // y+=A'*x
void cvec_sub_lintr_ct(int m, int n, cmulti **y, cmulti **A, int LDA, cmulti **x); // y-=A'*x


#endif
