#ifndef IS_RMULTI_H
#define IS_RMULTI_H

/**
 @file  is_rmulti.h
 @brief 多倍長精度実数型rmulti
 */

#include<is_macros.h>
#include<mpfr.h>

/**
 @brief 多倍長精度複素数構造体rmulti
 @details 実体は MPFR(http://www.mpfr.org/)の__mpfr_struct 構造体．
 */
typedef __mpfr_struct rmulti;

/*
 * precision
 */
void set_default_prec(int prec);
int get_default_prec(void);

/*
 * rouding mode
 */
mpfr_rnd_t get_round_mode(void);
void set_round_mode(mpfr_rnd_t default_round);

/*
 * allocation
 */
rmulti *rallocate(void);
rmulti *rallocate_prec(int prec);
rmulti *rallocate_clone(rmulti *y);
rmulti *rmfree(rmulti *x);
rmulti *rmepsilon(int prec);         // epsilon=1/2^prec
void rround(rmulti *x, int prec);    // y=round(x)
void rclone_r(rmulti *y, rmulti *x); // y=x
void rswap(rmulti *x, rmulti *y);    // x<->y

/*
 * I/O
 */
void rprint(rmulti *x, const char *name, const char *f, int digits);
void rput(rmulti *x);
void rput_info(rmulti *x, char *name);
void rbin_save(rmulti *x, FILE *fid);
rmulti *rbin_load(FILE *fid);

/*
 * member variables
 */
int rget_prec(rmulti *x);
int rget_exp(rmulti *x);
int rget_sgn(rmulti *x);
int rget_sgn_not0(rmulti *x);
int rget_size(rmulti *x);
int rget_length(rmulti *x);

/*
 * query
 */
int ris_nan(rmulti *x);       // x==nan
int ris_inf(rmulti *x);       // x==inf
int ris_zero(rmulti *x);      // x==0
int ris_positive(rmulti *x);  // x > 0
int ris_negative(rmulti *x);  // x < 0
int ris_number(rmulti *x);    // number?
int ris_integer(rmulti *x);   // integer?

/*
 * setting
 */
void rset_nan(rmulti *x);           // x=nan
void rset_inf(rmulti *x, int sgn);  // x=sgn*inf
void rset_zero(rmulti *x);          // x=0
void rset_one(rmulti *x);           // x=1
void rset_one_neg(rmulti *x);       // x=-1
void rset_rand(rmulti *x);          // x=rand()

/*
 * casting
 */
void   rset_r(rmulti *y, rmulti *x); // rmulti <- rmulti
void   rset_d(rmulti *y, double x);  // rmulti <- double
double rget_d(rmulti *x);            // rmulti -> double
void   rset_i(rmulti *y, int x);     // rmulti <- int
int    rget_i(rmulti *x);            // rmulti -> int
void   rset_s(rmulti *y, char *x);   // rmulti <- char
void   rget_s(char *y, rmulti *x);   // rmulti -> char

/*
 * y=f(x)
 */
// y=conj(x)
void rconj_r(rmulti *y, rmulti *x);
void rconj_d(rmulti *y, double x);
// y=real(x)
void rreal_r(rmulti *y, rmulti *x);
void rreal_d(rmulti *y, double x);
// y=imag(x)
void rimag_r(rmulti *y, rmulti *x);
void rimag_d(rmulti *y, double x);
// y=abs(x)
void rabs_r(rmulti *y, rmulti *x);
// y=abs(x)^2
void rabs2_r(rmulti *y, rmulti *x);
// y=abs(x)+i*abs(0)
void rabsc_r(rmulti *y, rmulti *x);
// y=max(abs(x),abs(0))
void rmax_absc_r(rmulti *y, rmulti *x);
// y=-x
void rneg_r(rmulti *y, rmulti *x);
void rneg_d(rmulti *y, double x);
// y=1/x
void rinv_r(rmulti *y, rmulti *x);
void rinv_d(rmulti *y, double x);
// y=x/abs(x)
void rnormalize_r(rmulti *y, rmulti *x);
void rnormalize_d(rmulti *y, double x);
// theta=arg(z)
void rarg_r(rmulti *theta, rmulti *z);
void rarg_d(rmulti *theta, double z);
// y=floor(x)
void rfloor_r(rmulti *y, rmulti *x);
void rfloor_d(rmulti *y, double x);
// y=ceil(x)
void rceil_r(rmulti *y, rmulti *x);
void rceil_d(rmulti *y, double x);
// y=trunc(x)
void rtrunc_r(rmulti *y, rmulti *x);
void rtrunc_d(rmulti *y, double x);
// y=sqrt(x)
void rsqrt_r(rmulti *y, rmulti *x);
void rsqrt_d(rmulti *y, double x);
// y=exp(x)
void rexp_r(rmulti *y, rmulti *x);
void rexp_d(rmulti *y, double x);
// y=exp2(x)
void rexp2_r(rmulti *y, rmulti *x);
void rexp2_d(rmulti *y, double x);
// y=exp10(x)
void rexp10_r(rmulti *y, rmulti *x);
void rexp10_d(rmulti *y, double x);
// y=log(x)
void rlog_r(rmulti *y, rmulti *x);
void rlog_d(rmulti *y, double x);
// y=log2(x)
void rlog2_r(rmulti *y, rmulti *x);
void rlog2_d(rmulti *y, double x);
// y=log10(x)
void rlog10_r(rmulti *y, rmulti *x);
void rlog10_d(rmulti *y, double x);
// y=sin(x)
void rsin_r(rmulti *y, rmulti *x);
void rsin_d(rmulti *y, double x);
// y=cos(x)
void rcos_r(rmulti *y, rmulti *x);
void rcos_d(rmulti *y, double x);
// y=tan(x)
void rtan_r(rmulti *y, rmulti *x);
void rtan_d(rmulti *y, double x);
// y=sec(x)
void rsec_r(rmulti *y, rmulti *x);
void rsec_d(rmulti *y, double x);
// y=csc(x)
void rcsc_r(rmulti *y, rmulti *x);
void rcsc_d(rmulti *y, double x);
// y=cot(x)
void rcot_r(rmulti *y, rmulti *x);
void rcot_d(rmulti *y, double x);
// y=asin(x)
void rasin_r(rmulti *y, rmulti *x);
void rasin_d(rmulti *y, double x);
// y=acos(x)
void racos_r(rmulti *y, rmulti *x);
void racos_d(rmulti *y, double x);
// y=atan(x)
void ratan_r(rmulti *y, rmulti *x);
void ratan_d(rmulti *y, double x);
// y=sinh(x)
void rsinh_r(rmulti *y, rmulti *x);
void rsinh_d(rmulti *y, double x);
// y=cosh(x)
void rcosh_r(rmulti *y, rmulti *x);
void rcosh_d(rmulti *y, double x);
// y=tanh(x)
void rtanh_r(rmulti *y, rmulti *x);
void rtanh_d(rmulti *y, double x);
// y=asinh(x)
void rasinh_r(rmulti *y, rmulti *x);
void rasinh_d(rmulti *y, double x);
// y=acosh(x)
void racosh_r(rmulti *y, rmulti *x);
void racosh_d(rmulti *y, double x);
// y=atanh(x)
void ratanh_r(rmulti *y, rmulti *x);
void ratanh_d(rmulti *y, double x);
// z=x^n
void rpow_r(rmulti *z, rmulti *x, int n);
void rpow_d(rmulti *z, double x, int n);
// y=x*2^n
void rmul_2exp(rmulti *y, rmulti *x, int n);
// y=x/2^n
void rdiv_2exp(rmulti *y, rmulti *x, int n);

/*
 * z=f(x,y)
 */
// z=x+y
int radd_rr_exact(rmulti *z, rmulti *x, rmulti *y);
void radd_rr(rmulti *z, rmulti *x, rmulti *y);
void radd_rd(rmulti *z, rmulti *x, double y);
void radd_dr(rmulti *z, double x, rmulti *y);
// z=x-y
int rsub_rr_exact(rmulti *z, rmulti *x, rmulti *y);
void rsub_rr(rmulti *z, rmulti *x, rmulti *y);
void rsub_rd(rmulti *z, rmulti *x, double y);
void rsub_dr(rmulti *z, double x, rmulti *y);
// z=x*y
int rmul_rr_exact(rmulti *z, rmulti *x, rmulti *y);
void rmul_rr(rmulti *z, rmulti *x, rmulti *y);
void rmul_rd(rmulti *z, rmulti *x, double y);
void rmul_dr(rmulti *z, double x, rmulti *y);
// z=x/y
int rdiv_rr_rouding_check(rmulti *z, rmulti *x, rmulti *y);
void rdiv_rr(rmulti *z, rmulti *x, rmulti *y);
void rdiv_rd(rmulti *z, rmulti *x, double y);
void rdiv_dr(rmulti *z, double x, rmulti *y);
// z=x^y
void rpow_rr(rmulti *z, rmulti *x, rmulti *y);
void rpow_rd(rmulti *z, rmulti *x, double y);
void rpow_dr(rmulti *z, double x, rmulti *y);
// z=atan(x/y)
void ratan2_rr(rmulti *z, rmulti *x, rmulti *y);
// z=abs(x-y)
void rabs_sub_rr(rmulti *z, rmulti *x, rmulti *y);
void rabs_sub_rd(rmulti *z, rmulti *x, double y);
void rabs_sub_dr(rmulti *z, double x, rmulti *y);
void rabs_sub_dd(rmulti *z, double x, double y);
// z=abs(x)+abs(y)
void radd_abs_abs_rr(rmulti *z, rmulti *x, rmulti *y);
void radd_abs_abs_rd(rmulti *z, rmulti *x, double y);
void radd_abs_abs_dr(rmulti *z, double x, rmulti *y);
void radd_abs_abs_dd(rmulti *z, double x, double y);
// z=max2(x,y)
void rmax2_up(rmulti *z, rmulti *x, rmulti *y);
// z=min2(x,y)
void rmin2_down(rmulti *z, rmulti *x, rmulti *y);
// z=10^(floor(log10(abs(x))-y))
void rexp10_floor_log10_abs_sub(rmulti *z, rmulti *x, double y);
// z=2^(floor(log2(abs(x))-y))
void rexp2_floor_log2_abs_sub(rmulti *z, rmulti *x, double y);

/*
 * y=y+f(x)
 */
// y=y+abs(x)
void radd_abs_r(rmulti *y, rmulti *x);
void radd_abs_d(rmulti *y, double x);
// y=y+abs(x)^2
void radd_abs2_r(rmulti *y, rmulti *x);
void radd_abs2_d(rmulti *y, double x);

/*
 * z=z+f(x,y)
 */
// z=z+x*y
int radd_mul_rr_exact(rmulti *z, rmulti *x, rmulti *y);
void radd_mul_rr_ws(rmulti *z, rmulti *x, rmulti *y, int n_rws, rmulti **rws);  // n_rws>=1
void radd_mul_rr(rmulti *z, rmulti *x, rmulti *y);
void radd_mul_rd(rmulti *z, rmulti *x, double y);
void radd_mul_dr(rmulti *z, double x, rmulti *y);
// z=z-x*y
int rsub_mul_rr_exact(rmulti *z, rmulti *x, rmulti *y);
void rsub_mul_rr_ws(rmulti *z, rmulti *x, rmulti *y, int n_rws, rmulti **rws);  // n_rws>=1
void rsub_mul_rr(rmulti *z, rmulti *x, rmulti *y);
void rsub_mul_rd(rmulti *z, rmulti *x, double y);
void rsub_mul_dd(rmulti *z, double x, double y);

/*
 * w=f(x,y,z)
 */
void radd_abs_abs_abs(rmulti *z, rmulti *x0, rmulti *x1, rmulti *x2); // z=abs(x0)+abs(x1)+abs(x2)

/*
 * x <=> y
 */
int cmp_rr(rmulti *x, rmulti *y);  // x<=>y
int cmp_rd(rmulti *x, double y);   // x<=>y
int eq_rr(rmulti *x, rmulti *y);   // x==y
int eq_rd(rmulti *x, double y);    // x==y
int ne_rr(rmulti *x, rmulti *y);   // x!=y
int ne_rd(rmulti *x, double y);    // x!=y
int gt_rr(rmulti *x, rmulti *y);   // x>y
int gt_dr(double x, rmulti *y);    // x>y
int gt_rd(rmulti *x, double y);    // x>y
int ge_rr(rmulti *x, rmulti *y);   // x>=y
int ge_dr(double x, rmulti *y);    // x>=y
int ge_rd(rmulti *x, double y);    // x>=y
int lt_rr(rmulti *x, rmulti *y);   // x<y
int lt_dr(double x, rmulti *y);    // x<y
int lt_rd(rmulti *x, double y);    // x<y
int le_rr(rmulti *x, rmulti *y);   // x<=y
int le_dr(double x, rmulti *y);    // x<=y
int le_rd(rmulti *x, double y);    // x<=y
int cmp_abs_rr(rmulti *x, rmulti *y); // abs(x)<=>abs(y)
int eq_abs_rr(rmulti *x, rmulti *y);  // abs(x)==abs(y)
int gt_abs_rr(rmulti *x, rmulti *y);  // abs(x)> abs(y)
int ge_abs_rr(rmulti *x, rmulti *y);  // abs(x)>=abs(y)
int lt_abs_rr(rmulti *x, rmulti *y);  // abs(x)< xabs(y)
int le_abs_rr(rmulti *x, rmulti *y);  // abs(x)<=abs(y)
int cmp_dist_rrr(rmulti *x, rmulti *y, rmulti *z); // abs(x-y)<=>z
int cmp_dist_rdd(rmulti *x, double y, double z);   // abs(x-y)<=>z
int eq_dist_rrr(rmulti *x, rmulti *y, rmulti *z);  // abs(x-y)==z
int eq_dist_rdd(rmulti *x, double y, double z);    // abs(x-y)==z
int ge_dist_rrr(rmulti *x, rmulti *y, rmulti *z);  // abs(x-y)>=z
int ge_dist_rdd(rmulti *x, double y, double z);    // abs(x-y)>=z
int gt_dist_rrr(rmulti *x, rmulti *y, rmulti *z);  // abs(x-y)> z
int gt_dist_rdd(rmulti *x, double y, double z);    // abs(x-y)> z
int le_dist_rrr(rmulti *x, rmulti *y, rmulti *z);  // abs(x-y)<=z
int le_dist_rdd(rmulti *x, double y, double z);    // abs(x-y)<=z
int lt_dist_rrr(rmulti *x, rmulti *y, rmulti *z);  // abs(x-y)< z
int lt_dist_rdd(rmulti *x, double y, double z);    // abs(x-y)< z

#endif
