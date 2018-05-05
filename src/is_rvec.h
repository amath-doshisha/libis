#ifndef IS_RVEC_H
#define IS_RVEC_H

#include<is_rmulti.h>
#include<is_func.h>

/**
 @file  is_rvec.h
 @brief 多倍長精度実数型rmultiのベクトルに関する関数の宣言
 */


/*
 * constructions and destruction
 */
rmulti **rvec_allocate(int n);
rmulti **rvec_allocate_prec(int n, int prec);
rmulti **rvec_allocate_clone(int n, rmulti **y);
rmulti **rvec_free(int n, rmulti **x);
void rvec_round(int n, rmulti **x, int prec);
void rvec_clone(int n, rmulti **y, rmulti **x); // y=x
void rvec_clone_index(int n, rmulti **y, rmulti **x, int *I); // y=x(I)
void rvec_index_clone(int n, rmulti **y, rmulti **x, int *I); // y(I)=x
void rvec_swap(int n, rmulti **x, rmulti **y); // x <=> y

/*
 * member variables
 */
void rvec_get_prec(int n, int *p, rmulti **x);
int rvec_get_prec_max(int n, rmulti **x);
void rvec_get_exp(int n, int *p, rmulti **x);
int rvec_get_exp_max(int n, rmulti **x);
int rvec_is_number(int n, rmulti **x);
int rvec_has_nan(int n, rmulti **x);

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
void rvec_set_nan(int n, rmulti **x);
void rvec_set_inf(int n, rmulti **x, int sgn);
void rvec_set_s(int n, rmulti **x, char **str);
void rvec_set_si(int n, rmulti **y, int *x);
void rvec_set_d(int n, rmulti **y, double *x);
void rvec_set_z(int n, rmulti **y, dcomplex *x);
void rvec_set_all_d(int n, rmulti **x, double a);
void rvec_set_zeros(int n, rmulti **x);
void rvec_set_ones(int n, rmulti **x);
void rvec_set_unit(int n, rmulti **x, int k);
void rvec_set_grid(int n, rmulti **x);
void rvec_set_rand(int n, rmulti **x, double a, double b);

/*
 * casting
 */
void rvec_get_si(int n, int *y, rmulti **x);     // y=int(x)
void rvec_get_d(int n, double *y, rmulti **x);   // y=double(x)
void rvec_get_z(int n, dcomplex *y, rmulti **x); // y=dcomplex(x)
void rvec_get_s(int n, char **y, rmulti **x, char format, int digits); // y=char(x)

/*
 * rearange elements
 */
void rvec_reverse(int n, rmulti **x);
void rvec_swap_at(rmulti **x, int i, int j);
void rvec_swap_index(int n, rmulti **x, const int *I);
void rvec_sort(int n, rmulti **x, int *I);
void rvec_quick_sort(int n, rmulti **x, int *I, int left, int right);
void rvec_sort_index(int *I, int n, rmulti **X);

/*
 * operations
 */
void rvec_copy(int n, rmulti **y, rmulti **x);                       // y=x
void rvec_copy_index(int n, rmulti **y, rmulti **x, const int *I);   // y=x(I)
void rvec_index_copy(int n, rmulti **y, rmulti **x, int *I);         // y(I)=x
void rvec_mul_2exp(int n, rmulti **y, rmulti **x, int p);            // y=x*2^p
void rvec_div_2exp(int n, rmulti **y, rmulti **x, int p);            // y=x/2^p
void rvec_neg(int n, rmulti **y, rmulti **x);                        // y=-x
void rvec_abs(int n, rmulti **y, rmulti **x);                        // y=abs(x)
void rvec_add(int n, rmulti **z, rmulti **x, rmulti **y);            // z=x+y
void rvec_add_r(int n, rmulti **z, rmulti **x, rmulti *y);           // z=x+y
void rvec_add_d(int n, rmulti **z, rmulti **x, double y);            // z=x+y
void rvec_sub(int n, rmulti **z, rmulti **x, rmulti **y);            // z=x-y
void rvec_sub_r1(int n, rmulti **z, rmulti *x, rmulti **y);          // z=x-y
void rvec_sub_r2(int n, rmulti **z, rmulti **x, rmulti *y);          // z=x-y
void rvec_mul(int n, rmulti **z, rmulti **x, rmulti **y);            // z=x*y
void rvec_mul_r(int n, rmulti **z, rmulti **x, rmulti *y);           // z=x*y
void rvec_mul_d(int n, rmulti **z, rmulti **x, double y);            // z=x*y
void rvec_add_mul(int n, rmulti **z, rmulti **x, rmulti **y);        // z+=x*y
void rvec_add_mul_r(int n, rmulti **z, rmulti **x, rmulti *y);       // z+=x*y
void rvec_add_mul_d(int n, rmulti **z, rmulti **x, double y);        // z+=x*y
void rvec_sub_mul(int n, rmulti **z, rmulti **x, rmulti **y);        // z-=x*y
void rvec_sub_mul_r(int n, rmulti **z, rmulti **x, rmulti *y);       // z-=x*y
void rvec_sub_mul_d(int n, rmulti **z, rmulti **x, double y);        // z-=x*y
void rvec_div(int n, rmulti **z, rmulti **x, rmulti **y);            // z=x./y
void rvec_div_r1(int n, rmulti **z, rmulti *x, rmulti **y);          // z=x./y
void rvec_div_d1(int n, rmulti **z, double x, rmulti **y);           // z=x./y
void rvec_div_r2(int n, rmulti **z, rmulti **x, rmulti *y);          // z=x/y
void rvec_div_d2(int n, rmulti **z, rmulti **x, double y);           // z=x/y
void rvec_pow_ui(int n, rmulti **z, rmulti **x, ulong y);            // z=x^y
void rvec_pow(int n, rmulti **z, rmulti **x, rmulti **y);            // z=x.^y
void rvec_pow_si(int n, rmulti **z, rmulti **x, long y);             // z=x^y
void rvec_pow_r(int n, rmulti **z, rmulti **x, rmulti *y);           // z=x^y
void rvec_abs_sub(int n, rmulti **z, rmulti **x, rmulti **y);        // z=abs(x-y)
void rvec_abs_sub_r(int n, rmulti **z, rmulti **x, rmulti *y);       // z=abs(x-y)
void rvec_sum(rmulti *value, int n, rmulti **x);                     // value=sum(x)
void rvec_sum_abs(rmulti *value, int n, rmulti **x);                 // value=sum(abs(x))
void rvec_sum_pow2(rmulti *value, int n, rmulti **x);                // value=sum(x.^2)
void rvec_sum_mul(rmulti *value, int n, rmulti **x, rmulti **y);     // value=sum(x.*y)
void rvec_sum_abs_sub(rmulti *value, int n, rmulti **x, rmulti **y); // value=sum(abs(x-y))
void rvec_max(rmulti *value, int n, rmulti **x);                     // value=max(x)
void rvec_max_abs(rmulti *value, int n, rmulti **x);                 // value=max(abs(x))
void rvec_max_abs_index(rmulti *value, int n, rmulti **x, int *I);   // value=max(abs(x))
void rvec_max_abs_sub(rmulti *value, int n, rmulti **x, rmulti **y); // value=max(abs(x-y))
void rvec_min(rmulti *value, int n, rmulti **x);                     // value=min(x)
void rvec_min_abs(rmulti *value, int n, rmulti **x);                 // value=min(abs(x))
void rvec_min_abs_index(rmulti *value, int n, rmulti **x, int *I);   // value=min(abs(x))
void rvec_log2_abs(int n, rmulti **y, rmulti **x);                   // y=log2(abs(x))
void rvec_normalize(int n, rmulti **y, rmulti **x);                  // y=x/sqrt(x'*x)
void rvec_normalize_sgn(int n, rmulti **y, rmulti **x);              // y=x/sqrt(x'*x)
void rvec_orthogonalize(int n, rmulti **y, rmulti **x);              // y-=(x'*y)*x where x'*x=1
void rvec_norm2(rmulti *value, int n, rmulti **x);                   // value=sqrt(sum(x.^2))
void rvec_average(rmulti *value, int n, rmulti **x);                 // value=sum(x)/n
void rvec_max_div_abs(rmulti *value, int n, rmulti **x, rmulti **y); // value=max(abs(x)./abs(y))
void rvec_dcos(rmulti *value, int n, rmulti **x, rmulti **y);        // value=(x'*y)/sqrt(x'*x)/sqrt(y'*y)
void rvec_abs_dcos(rmulti *value, int n, rmulti **x, rmulti **y);    // value=abs(x'*y)/sqrt(x'*x)/sqrt(y'*y)
void rvec_angle(rmulti *theta, int n, rmulti **x, rmulti **y);       // theta=acos(abs(x'*y)/sqrt(x'*x)/sqrt(y'*y))
void rvec_angle_deg(rmulti *theta, int n, rmulti **x, rmulti **y);   // theta=(180/PI)*acos(abs(x'*y)/sqrt(x'*x)/sqrt(y'*y))
void rvec_lintr(int m, int n, rmulti **y, rmulti **A, int LDA, rmulti **x);       // y=A*x
void rvec_add_lintr(int m, int n, rmulti **y, rmulti **A, int LDA, rmulti **x);   // y+=A*x
void rvec_sub_lintr(int m, int n, rmulti **y, rmulti **A, int LDA, rmulti **x);   // y-=A*x
void rvec_lintr_t(int m, int n, rmulti **y, rmulti **A, int LDA, rmulti **x);     // y=A'*x
void rvec_add_lintr_t(int m, int n, rmulti **y, rmulti **A, int LDA, rmulti **x); // y+=A'*x
void rvec_sub_lintr_t(int m, int n, rmulti **y, rmulti **A, int LDA, rmulti **x); // y-=A'*x

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
 * mapping
 */
// y=f(x)
void rvec_func(rmulti *y, func_t *f, int n, rmulti **x);
void rvec_func_list(int m, rmulti **y, func_t *f, int n, rmulti **x);
void r1_func(rmulti *y, func_t *f, rmulti *x0);
void r2_func(rmulti *y, func_t *f, rmulti *x0, rmulti *x1);
void r3_func(rmulti *y, func_t *f, rmulti *x0, rmulti *x1, rmulti *x2);

#endif
