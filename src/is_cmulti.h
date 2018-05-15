#ifndef IS_CMULTI_H
#define IS_CMULTI_H

#include<is_macros.h>
#include<is_dcomplex.h>
#include<is_rmulti.h>

/**
 @file  is_cmulti.h
 @brief 多倍長精度複素数型cmulti
 */

/**
 @brief 多倍長精度複素数構造体cmulti
 */
typedef struct {
  rmulti *r; /**< 実部 */
  rmulti *i; /**< 虚部 */
} cmulti_struct;

/**
 @brief 多倍長精度複素数型
 */
typedef cmulti_struct cmulti;

/**
 @brief cmulit型の実部
 */
#define C_R(X) ((X)->r)

/**
 @brief cmulit型の虚部
 */
#define C_I(X) ((X)->i)

/*
 * constructions and destruction
 */
cmulti *callocate(void);
cmulti *callocate_prec(int prec);
cmulti *callocate_clone(cmulti *y);
cmulti *callocate_clone_r(rmulti *y);
cmulti *cfree(cmulti *x);
void cround(cmulti *x, int prec);
void cclone_c(cmulti *y, cmulti *x);
void cclone_r(cmulti *y, rmulti *x);
void cclone_rr(cmulti *y, rmulti *x_r, rmulti *x_i);
void cswap(cmulti *x, cmulti *y);

/*
 * I/O
 */
void cbin_save(cmulti *x, FILE *fid);
cmulti *cbin_load(FILE *fid);

/*
 * member variables
 */
int cget_prec(cmulti *x); // max(x_r.prec,y_i.prec)
int cget_exp(cmulti *x);  // max(x_r.exp,y_i.exp)

/*
 * query
 */
int cis_nan(cmulti *x);              // x==nan
int cis_inf(cmulti *x);              // x==inf
int cis_zero(cmulti *x);             // x==0
int cis_number(cmulti *x);           // number?
int cis_real(cmulti *x);             // real number?
int cis_pure_imaginary(cmulti *x);   // pure imaginary number?
int cis_imaginary(cmulti *x);        // imaginary number?

/*
 * setting
 */
void cset_nan(cmulti *x);                             // x=nan
void cset_inf(cmulti *x, int sgn_r, int sgn_i);       // x=sgn_r*inf+i*sng_i*inf
void cset_zero(cmulti *x);                            // x=0
void cset_one(cmulti *x);                             // x=1
void cset_one_neg(cmulti *x);                         // x=-1
void cset_rand(cmulti *x);                            // x=rand()
void cset_polar(cmulti *z, rmulti *r, rmulti *theta); // z=r*exp(i*theta)
void cget_polar(rmulti *r, rmulti *theta, cmulti *z); // z=r*exp(i*theta)

/*
 * casting
 */
void     cset_c(cmulti *y, cmulti *x);                 // cmulti <- cmulti
void     cset_rr(cmulti *y, rmulti *x_r, rmulti *x_i); // cmulti <- (rmulti,rmulti)
void     cset_z(cmulti *y, dcomplex x);                // cmulti <- dcomplex
dcomplex cget_z(cmulti *x);                            // cmulti -> dcomplex
void     cset_dd(cmulti *y, double xr, double xi);     // cmulti <- (double,double)
void     cset_r(cmulti *y, rmulti *x);                 // cmulti <- rmulti
void     cget_r(rmulti *y, rmulti *x);                 // cmulti -> rmulti
void     cset_d(cmulti *y, double x);                  // cmulti <- double
double   cget_d(cmulti *x);                            // cmulti -> double
void     cset_i(cmulti *y, int x);                     // cmulti <- int
int      cget_i(cmulti *x);                            // cmulti -> int
void     cset_s(cmulti *y, char *x);                   // cmulti <- char
void     cget_s(char *y, cmulti *x);                   // cmulti -> char

/*
 * y=f(x)
 */
// y=conj(x)
void cconj_c(cmulti *y, cmulti *x);
void cconj_r(cmulti *y, rmulti *x);
void cconj_z(cmulti *y, dcomplex x);
void cconj_d(cmulti *y, double x);
// y=real(x)
void rreal_c(rmulti *y, cmulti *x);
void rreal_z(rmulti *y, dcomplex x);
// y=imag(x)
void rimag_c(rmulti *y, cmulti *x);
void rimag_z(rmulti *y, dcomplex x);
// y=abs(x)
void rabs_c_ws(rmulti *y, cmulti *x, int n_rws, rmulti **rws); // n_rws>=1
void rabs_c(rmulti *y, cmulti *x);
// y=abs(x)^2
void rabs2_c_ws(rmulti *y, cmulti *x, int n_rws, rmulti **rws); // n_rws>=1
void rabs2_c(rmulti *y, cmulti *x);
// y=abs(x_r)+i*abs(x_i)
void cabsc_c(cmulti *y, cmulti *x);
// y=max(abs(x_r),abs(x_i))
void rmax_absc_c(rmulti *y, cmulti *x);
// y=-x
void cneg_c(cmulti *y, cmulti *x);
void cneg_r(cmulti *y, rmulti *x);
void cneg_z(cmulti *y, dcomplex x);
void cneg_d(cmulti *y, double x);
// y=1/x
void cinv_c_ws(cmulti *y, cmulti *x, int n_rws, rmulti **rws, int n_cws, cmulti **cws); // n_rws>=2,n_cws>=1
void cinv_c(cmulti *y, cmulti *x);
void cinv_r(cmulti *y, rmulti *x);
void cinv_z(cmulti *y, dcomplex x);
void cinv_d(cmulti *y, double x);
// y=x/abs(x)
void cnormalize_c(cmulti *y, cmulti *x);
void cnormalize_r(cmulti *y, rmulti *x);
void cnormalize_z(cmulti *y, dcomplex x);
void cnormalize_d(cmulti *y, double x);
// theta=arg(z)
void rarg_c(rmulti *theta, cmulti *z);
void rarg_z(rmulti *theta, dcomplex z);
// y=floor(x)
void cfloor_c(cmulti *y, cmulti *x);
void cfloor_r(cmulti *y, rmulti *x);
void cfloor_z(cmulti *y, dcomplex x);
void cfloor_d(cmulti *y, double x);
// y=ceil(x)
void cceil_c(cmulti *y, cmulti *x);
void cceil_r(cmulti *y, rmulti *x);
void cceil_z(cmulti *y, dcomplex x);
void cceil_d(cmulti *y, double x);
// y=trunc(x)
void ctrunc_c(cmulti *y, cmulti *x);
void ctrunc_r(cmulti *y, rmulti *x);
void ctrunc_z(cmulti *y, dcomplex x);
void ctrunc_d(cmulti *y, double x);
// y=sqrt(x)
void csqrt_c(cmulti *y, cmulti *x);
void csqrt_r(cmulti *y, rmulti *x);
void csqrt_z(cmulti *y, dcomplex x);
void csqrt_d(cmulti *y, double x);
// y=exp(x)
void cexp_c(cmulti *y, cmulti *x);
void cexp_r(cmulti *y, rmulti *x);
void cexp_z(cmulti *y, dcomplex x);
void cexp_d(cmulti *y, double x);
// y=exp2(x)
void cexp2_c(cmulti *y, cmulti *x);
void cexp2_r(cmulti *y, rmulti *x);
void cexp2_z(cmulti *y, dcomplex x);
void cexp2_d(cmulti *y, double x);
// y=exp10(x)
void cexp10_c(cmulti *y, cmulti *x);
void cexp10_r(cmulti *y, rmulti *x);
void cexp10_z(cmulti *y, dcomplex x);
void cexp10_d(cmulti *y, double x);
// y=log(x)
void clog_c(cmulti *y, cmulti *x);
void clog_r(cmulti *y, rmulti *x);
void clog_z(cmulti *y, dcomplex x);
void clog_d(cmulti *y, double x);
// y=log2(x)
void clog2_c(cmulti *y, cmulti *x);
void clog2_r(cmulti *y, rmulti *x);
void clog2_z(cmulti *y, dcomplex x);
void clog2_d(cmulti *y, double x);
// y=log10(x)
void clog10_c(cmulti *y, cmulti *x);
void clog10_r(cmulti *y, rmulti *x);
void clog10_z(cmulti *y, dcomplex x);
void clog10_d(cmulti *y, double x);
// y=sin(x)
void csin_c(cmulti *y, cmulti *x);
void csin_r(cmulti *y, rmulti *x);
void csin_z(cmulti *y, dcomplex x);
void csin_d(cmulti *y, double x);
// y=cos(x)
void ccos_c(cmulti *y, cmulti *x);
void ccos_r(cmulti *y, rmulti *x);
void ccos_z(cmulti *y, dcomplex x);
void ccos_d(cmulti *y, double x);
// y=tan(x)
void ctan_c(cmulti *y, cmulti *x);
void ctan_r(cmulti *y, rmulti *x);
void ctan_z(cmulti *y, dcomplex x);
void ctan_d(cmulti *y, double x);
// y=asin(x)
void casin_r(cmulti *y, rmulti *x);
void casin_r(cmulti *y, rmulti *x);
void casin_z(cmulti *y, dcomplex x);
void casin_d(cmulti *y, double x);
// y=asin(x)
void casin_c(cmulti *y, cmulti *x);
void casin_r(cmulti *y, rmulti *x);
void casin_z(cmulti *y, dcomplex x);
void casin_d(cmulti *y, double x);
// y=cos(x)
void cacos_c(cmulti *y, cmulti *x);
void cacos_r(cmulti *y, rmulti *x);
void cacos_z(cmulti *y, dcomplex x);
void cacos_d(cmulti *y, double x);
// y=tan(x)
void catan_c(cmulti *y, cmulti *x);
void catan_r(cmulti *y, rmulti *x);
void catan_z(cmulti *y, dcomplex x);
void catan_d(cmulti *y, double x);
// y=sinh(x)
void csinh_c(cmulti *y, cmulti *x);
void csinh_r(cmulti *y, rmulti *x);
void csinh_z(cmulti *y, dcomplex x);
void csinh_d(cmulti *y, double x);
// y=cosh(x)
void ccosh_c(cmulti *y, cmulti *x);
void ccosh_r(cmulti *y, rmulti *x);
void ccosh_z(cmulti *y, dcomplex x);
void ccosh_d(cmulti *y, double x);
// y=tanh(x)
void ctanh_c(cmulti *y, cmulti *x);
void ctanh_r(cmulti *y, rmulti *x);
void ctanh_z(cmulti *y, dcomplex x);
void ctanh_d(cmulti *y, double x);
// y=asinh(x)
void casinh_c(cmulti *y, cmulti *x);
void casinh_r(cmulti *y, rmulti *x);
void casinh_z(cmulti *y, dcomplex x);
void casinh_d(cmulti *y, double x);
// y=acosh(x)
void cacosh_c(cmulti *y, cmulti *x);
void cacosh_r(cmulti *y, rmulti *x);
void cacosh_z(cmulti *y, dcomplex x);
void cacosh_d(cmulti *y, double x);
// y=atanh(x)
void catanh_c(cmulti *y, cmulti *x);
void catanh_r(cmulti *y, rmulti *x);
void catanh_z(cmulti *y, dcomplex x);
void catanh_d(cmulti *y, double x);
// z=x^n
void cpow_c(cmulti *y, cmulti *x, int n);
void cpow_r(cmulti *y, rmulti *x, int n);
void cpow_z(cmulti *y, dcomplex x, int n);
void cpow_d(cmulti *y, double x, int n);
// y=x*2^n
void cmul_2exp(cmulti *y, cmulti *x, int nr, int ni);
// y=x/2^n
void cdiv_2exp(cmulti *y, cmulti *x, int nr, int ni);

/*
 * z=f(x,y)
 */
// z=x+y
void cadd_cc(cmulti *z, cmulti *x, cmulti *y);
void cadd_cz(cmulti *z, cmulti *x, dcomplex y);
void cadd_zc(cmulti *z, dcomplex x, cmulti *y);
void cadd_cr(cmulti *z, cmulti *x, rmulti *y);
void cadd_rc(cmulti *z, rmulti *x, cmulti *y);
void cadd_cd(cmulti *z, cmulti *x, double y);
void cadd_dc(cmulti *z, double x, cmulti *y);
void cadd_rz(cmulti *z, rmulti *x, dcomplex y);
void cadd_zr(cmulti *z, dcomplex x, rmulti *y);
// z=x-y
void csub_cc(cmulti *z, cmulti *x, cmulti *y);
void csub_cz(cmulti *z, cmulti *x, dcomplex y);
void csub_zc(cmulti *z, dcomplex x, cmulti *y);
void csub_cr(cmulti *z, cmulti *x, rmulti *y);
void csub_rc(cmulti *z, rmulti *x, cmulti *y);
void csub_cd(cmulti *z, cmulti *x, double y);
void csub_dc(cmulti *z, double x, cmulti *y);
void csub_rz(cmulti *z, rmulti *x, dcomplex y);
void csub_zr(cmulti *z, dcomplex x, rmulti *y);
// z=x*y
void cmul_cc_ws(cmulti *z, cmulti *x, cmulti *y, int n_rws, rmulti **rws, int n_cws, cmulti **cws);  // n_rws>=1,n_cws>=1
void cmul_cc(cmulti *z, cmulti *x, cmulti *y);
void cmul_cz(cmulti *z, cmulti *x, dcomplex y);
void cmul_zc(cmulti *z, dcomplex x, cmulti *y);
void cmul_cr(cmulti *z, cmulti *x, rmulti *y);
void cmul_rc(cmulti *z, rmulti *x, cmulti *y);
void cmul_cd(cmulti *z, cmulti *x, double y);
void cmul_dc(cmulti *z, double x, cmulti *y);
void cmul_rz(cmulti *z, rmulti *x, dcomplex y);
void cmul_zr(cmulti *z, dcomplex x, rmulti *y);
// z=conj(x)*y
void cdot_cc(cmulti *z, cmulti *x, cmulti *y);
void cdot_cz(cmulti *z, cmulti *x, dcomplex y);
void cdot_zc(cmulti *z, dcomplex x, cmulti *y);
void cdot_cr(cmulti *z, cmulti *x, rmulti *y);
void cdot_rc(cmulti *z, rmulti *x, cmulti *y);
void cdot_cd(cmulti *z, cmulti *x, double y);
void cdot_dc(cmulti *z, double x, cmulti *y);
void cdot_rz(cmulti *z, rmulti *x, dcomplex y);
void cdot_zr(cmulti *z, dcomplex x, rmulti *y);
// z=x/y
void cdiv_cc(cmulti *z, cmulti *x, cmulti *y);
void cdiv_cz(cmulti *z, cmulti *x, dcomplex y);
void cdiv_zc(cmulti *z, dcomplex x, cmulti *y);
void cdiv_cr(cmulti *z, cmulti *x, rmulti *y);
void cdiv_rc(cmulti *z, rmulti *x, cmulti *y);
void cdiv_cd(cmulti *z, cmulti *x, double y);
void cdiv_dc(cmulti *z, double x, cmulti *y);
void cdiv_rz(cmulti *z, rmulti *x, dcomplex y);
void cdiv_zr(cmulti *z, dcomplex x, rmulti *y);
// z=x^y
void cpow_cc(cmulti *z, cmulti *x, cmulti *y);
void cpow_cz(cmulti *z, cmulti *x, dcomplex y);
void cpow_zc(cmulti *z, dcomplex x, cmulti *y);
void cpow_cr(cmulti *z, cmulti *x, rmulti *y);
void cpow_rc(cmulti *z, rmulti *x, cmulti *y);
void cpow_cd(cmulti *z, cmulti *x, double y);
void cpow_dc(cmulti *z, double x, cmulti *y);
void cpow_rz(cmulti *z, rmulti *x, dcomplex y);
void cpow_zr(cmulti *z, dcomplex x, rmulti *y);
// z=abs(x-y)
void rabs_sub_cc(rmulti *z, cmulti *x, cmulti *y);
void rabs_sub_cz(rmulti *z, cmulti *x, dcomplex y);
void rabs_sub_zc(rmulti *z, dcomplex x, cmulti *y);
void rabs_sub_cr(rmulti *z, cmulti *x, rmulti *y);
void rabs_sub_rc(rmulti *z, rmulti *x, cmulti *y);
void rabs_sub_cd(rmulti *z, cmulti *x, double y);
void rabs_sub_dc(rmulti *z, double x, cmulti *y);
void rabs_sub_rz(rmulti *z, rmulti *x, dcomplex y);
void rabs_sub_zr(rmulti *z, dcomplex x, rmulti *y);
// z=abs(x)+abs(y)
void radd_abs_abs_cc(rmulti *z, cmulti *x, cmulti *y);
void radd_abs_abs_cz(rmulti *z, cmulti *x, dcomplex y);
void radd_abs_abs_zc(rmulti *z, dcomplex x, cmulti *y);
void radd_abs_abs_cr(rmulti *z, cmulti *x, rmulti *y);
void radd_abs_abs_rc(rmulti *z, rmulti *x, cmulti *y);
void radd_abs_abs_cd(rmulti *z, cmulti *x, double y);
void radd_abs_abs_dc(rmulti *z, double x, cmulti *y);
void radd_abs_abs_rz(rmulti *z, rmulti *x, dcomplex y);
void radd_abs_abs_zr(rmulti *z, dcomplex x, rmulti *y);
// z=10^(floor(log10(abs(x))-y))
void cexp10_floor_log10_abs_sub(cmulti *z, cmulti *x, double y);
// z=2^(floor(log2(abs(x))-y))
void cexp2_floor_log2_abs_sub(cmulti *z, cmulti *x, double y);

/*
 * y=y+f(x)
 */
// y=y+abs(x)
void radd_abs_c(rmulti *y, cmulti *x);
void radd_abs_z(rmulti *y, dcomplex x);
// y=y+abs(x)^2
void radd_abs2_c(rmulti *y, cmulti *x);
void radd_abs2_z(rmulti *y, dcomplex x);

/*
 * z=z+f(x,y)
 */
// z=z+x*y
void cadd_mul_cc(cmulti *z, cmulti *x, cmulti *y);
void cadd_mul_cr(cmulti *z, cmulti *x, rmulti *y);
void cadd_mul_rc(cmulti *z, rmulti *x, cmulti *y);
void cadd_mul_cz(cmulti *z, cmulti *x, dcomplex y);
void cadd_mul_zc(cmulti *z, dcomplex x, cmulti *y);
void cadd_mul_cd(cmulti *z, cmulti *x, double y);
void cadd_mul_dc(cmulti *z, double x, cmulti *y);
void cadd_mul_rz(cmulti *z, rmulti *x, dcomplex y);
void cadd_mul_zr(cmulti *z, dcomplex x, rmulti *y);
// z=z-x*y
void csub_mul_cc_ws(cmulti *z, cmulti *x, cmulti *y, int n_rws, rmulti **rws, int n_cws, cmulti **cws);  // n_rws>=1,n_cws>=2
void csub_mul_cc(cmulti *z, cmulti *x, cmulti *y);
void csub_mul_cr(cmulti *z, cmulti *x, rmulti *y);
void csub_mul_rc(cmulti *z, rmulti *x, cmulti *y);
void csub_mul_cz(cmulti *z, cmulti *x, dcomplex y);
void csub_mul_zc(cmulti *z, dcomplex x, cmulti *y);
void csub_mul_cd(cmulti *z, cmulti *x, double y);
void csub_mul_dc(cmulti *z, double x, cmulti *y);
void csub_mul_rz(cmulti *z, rmulti *x, dcomplex y);
void csub_mul_zr(cmulti *z, dcomplex x, rmulti *y);
// z=z+conj(x)*y
void cadd_dot_cc(cmulti *z, cmulti *x, cmulti *y);
void cadd_dot_cr(cmulti *z, cmulti *x, rmulti *y);
void cadd_dot_rc(cmulti *z, rmulti *x, cmulti *y);
void cadd_dot_cz(cmulti *z, cmulti *x, dcomplex y);
void cadd_dot_zc(cmulti *z, dcomplex x, cmulti *y);
void cadd_dot_cd(cmulti *z, cmulti *x, double y);
void cadd_dot_dc(cmulti *z, double x, cmulti *y);
void cadd_dot_rz(cmulti *z, rmulti *x, dcomplex y);
void cadd_dot_zr(cmulti *z, dcomplex x, rmulti *y);
// z=z-conj(x)*y
void csub_dot_cc(cmulti *z, cmulti *x, cmulti *y);
void csub_dot_cr(cmulti *z, cmulti *x, rmulti *y);
void csub_dot_rc(cmulti *z, rmulti *x, cmulti *y);
void csub_dot_cz(cmulti *z, cmulti *x, dcomplex y);
void csub_dot_zc(cmulti *z, dcomplex x, cmulti *y);
void csub_dot_cd(cmulti *z, cmulti *x, double y);
void csub_dot_dc(cmulti *z, double x, cmulti *y);
void csub_dot_rz(cmulti *z, rmulti *x, dcomplex y);
void csub_dot_zr(cmulti *z, dcomplex x, rmulti *y);

/*
 * x<=>y
 */
void ccmp_set_real_imag(void);
void ccmp_set_abs_arg(void);
int ccmp_get_type(void);
// x<=>y
int cmp_cc(cmulti *x, cmulti *y);
int cmp_cz(cmulti *x, dcomplex y);
int cmp_zc(dcomplex x, cmulti *y);
int cmp_cr(cmulti *x, rmulti *y);
int cmp_rc(rmulti *x, cmulti *y);
int cmp_cd(cmulti *x, double y);
int cmp_dc(double x, cmulti *y);
int cmp_rz(rmulti *x, dcomplex y);
int cmp_zr(dcomplex x, rmulti *y);
// x==y
int eq_cc(cmulti *x, cmulti *y);
int eq_cz(cmulti *x, dcomplex y);
int eq_zc(dcomplex x, cmulti *y);
int eq_cr(cmulti *x, rmulti *y);
int eq_rc(rmulti *x, cmulti *y);
int eq_cd(cmulti *x, double y);
int eq_dc(double x, cmulti *y);
int eq_rz(rmulti *x, dcomplex y);
int eq_zr(dcomplex x, rmulti *y);
// x!=y
int ne_cc(cmulti *x, cmulti *y);
int ne_cz(cmulti *x, dcomplex y);
int ne_zc(dcomplex x, cmulti *y);
int ne_cr(cmulti *x, rmulti *y);
int ne_rc(rmulti *x, cmulti *y);
int ne_cd(cmulti *x, double y);
int ne_dc(double x, cmulti *y);
int ne_rz(rmulti *x, dcomplex y);
int ne_zr(dcomplex x, rmulti *y);
// x>y
int gt_cc(cmulti *x, cmulti *y);
int gt_rc(rmulti *x, cmulti *y);
int gt_cr(cmulti *x, rmulti *y);
int gt_zc(dcomplex x, cmulti *y);
int gt_cz(cmulti *x, dcomplex y);
int gt_dc(double x, cmulti *y);
int gt_cd(cmulti *x, double y);
int gt_rz(rmulti *x, dcomplex y);
int gt_zr(dcomplex x, rmulti *y);
// x>=y
int ge_cc(cmulti *x, cmulti *y);
int ge_rc(rmulti *x, cmulti *y);
int ge_cr(cmulti *x, rmulti *y);
int ge_zc(dcomplex x, cmulti *y);
int ge_cz(cmulti *x, dcomplex y);
int ge_dc(double x, cmulti *y);
int ge_cd(cmulti *x, double y);
int ge_rz(rmulti *x, dcomplex y);
int ge_zr(dcomplex x, rmulti *y);
// x<y
int lt_cc(cmulti *x, cmulti *y);
int lt_rc(rmulti *x, cmulti *y);
int lt_cr(cmulti *x, rmulti *y);
int lt_zc(dcomplex x, cmulti *y);
int lt_cz(cmulti *x, dcomplex y);
int lt_dc(double x, cmulti *y);
int lt_cd(cmulti *x, double y);
int lt_rz(rmulti *x, dcomplex y);
int lt_zr(dcomplex x, rmulti *y);
// x<=y
int le_cc(cmulti *x, cmulti *y);
int le_rc(rmulti *x, cmulti *y);
int le_cr(cmulti *x, rmulti *y);
int le_zc(dcomplex x, cmulti *y);
int le_cz(cmulti *x, dcomplex y);
int le_dc(double x, cmulti *y);
int le_cd(cmulti *x, double y);
int le_rz(rmulti *x, dcomplex y);
int le_zr(dcomplex x, rmulti *y);
// others
int cmp_abs_cc(cmulti *x, cmulti *y); // abs(x)<=>abs(y)
int eq_abs_cc(cmulti *x, cmulti *y);  // abs(x)==abs(y)
int gt_abs_cc(cmulti *x, cmulti *y);  // abs(x)>abs(y)
int ge_abs_cc(cmulti *x, cmulti *y);  // abs(x)>=abs(y)
int lt_abs_cc(cmulti *x, cmulti *y);  // abs(x)<abs(y)
int le_abs_cc(cmulti *x, cmulti *y);  // abs(x)<=abs(y)
int ceq_cc(cmulti *x, cmulti *y);     // x.r==y.r && x.i==y.i
int cgt_cc(cmulti *x, cmulti *y);     // x.r> y.r && x.i> y.i
int cge_cc(cmulti *x, cmulti *y);     // x.r>=y.r && x.i>=y.i
int clt_cc(cmulti *x, cmulti *y);     // x.r< y.r && x.i< y.i
int cle_cc(cmulti *x, cmulti *y);     // x.r<=y.r && x.i<=y.i

#endif

