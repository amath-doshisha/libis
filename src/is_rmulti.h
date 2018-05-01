#ifndef IS_RMULTI_H
#define IS_RMULTI_H

/**
 @file  is_rmulti.h
 @brief 多倍長精度実数型rmultiの宣言と関数．
 */

#include<is_macros.h>
#include<mpfr.h>

/**
 @brief   多倍長精度実数構造体の宣言．
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
 * constructions and destruction
 */
rmulti *rallocate(void);
rmulti *rallocate_prec(int prec);
rmulti *rallocate_clone(rmulti *y);
rmulti *rfree(rmulti *x);
void rround(rmulti *x, int prec);
void rclone(rmulti *y, rmulti *x); // y=x
void rcopy(rmulti *y, rmulti *x);  // y=x
void rswap(rmulti *x, rmulti *y);  // x<->y
rmulti *rmepsilon(int prec);

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
 * type query
 */
int ris_nan(rmulti *x);
int ris_inf(rmulti *x);
int ris_number(rmulti *x);
int ris_zero(rmulti *x);
int ris_positive(rmulti *x);
int ris_negative(rmulti *x);
int ris_integer(rmulti *x);

/*
 * I/O
 */
void rprint(rmulti *x, const char *name, const char *f, int digits);
void rput(rmulti *x);
void rput_info(rmulti *x, char *name);
void rbin_save(rmulti *x, FILE *fid);
rmulti *rbin_load(FILE *fid);

/*
 * setting
 */
void rset_s(rmulti *x, const char *value);
void rset_d(rmulti *x, double value);
void rset_ui(rmulti *x, ulong value);
void rset_si(rmulti *x, long value);
void rset_inf(rmulti *x, int sgn);
void rset_nan(rmulti *x);
void rset_zero(rmulti *x);
void rset_one(rmulti *x);
void rset_one_neg(rmulti *x);
void rset_rand(rmulti *x);

/*
 * casting
 */
double rget_d(rmulti *x); // y=(double)x
ulong rget_ui(rmulti *x); // y=(unsigned long int)x
long rget_si(rmulti *x);  // y=(signed long int)x

/*
 * operatior for one argument
 */
void rneg(rmulti *y, rmulti *x);     // y=-x
void rinv(rmulti *y, rmulti *x);     // y=1/x
void rinv_d(rmulti *y, double x);    // y=1/x
void rfloor(rmulti *y, rmulti *x);   // y=floor(x)
void rceil(rmulti *y, rmulti *x);    // y=ceil(x)
void rtrunc(rmulti *y, rmulti *x);   // y=trunc(x)
void rabs(rmulti *y, rmulti *x);     // y=abs(x)
void radd_abs(rmulti *y, rmulti *x); // y+=abs(x)
void rsqrt(rmulti *y, rmulti *x);    // y=sqrt(x)
void rexp(rmulti *y, rmulti *x);     // y=exp(x)
void rexp2(rmulti *y, rmulti *x);    // y=exp2(x)
void rexp10(rmulti *y, rmulti *x);   // y=exp10(x)
void rlog(rmulti *y, rmulti *x);     // y=log(x)
void rlog2(rmulti *y, rmulti *x);    // y=log2(x)
void rlog10(rmulti *y, rmulti *x);   // y=log10(x)
void rsin(rmulti *y, rmulti *x);     // y=sin(x)
void rcos(rmulti *y, rmulti *x);     // y=cos(x)
void rtan(rmulti *y, rmulti *x);     // y=tan(x)
void rsec(rmulti *y, rmulti *x);     // y=sec(x)
void rcsc(rmulti *y, rmulti *x);     // y=csc(x)
void rcot(rmulti *y, rmulti *x);     // y=cot(x)
void rasin(rmulti *y, rmulti *x);    // y=asin(x)
void racos(rmulti *y, rmulti *x);    // y=acos(x)
void ratan(rmulti *y, rmulti *x);    // y=atan(x)
void rsinh(rmulti *y, rmulti *x);    // y=sinh(x)
void rcosh(rmulti *y, rmulti *x);    // y=cosh(x)
void rtanh(rmulti *y, rmulti *x);    // y=tanh(x)
void rasinh(rmulti *y, rmulti *x);   // y=asinh(x)
void racosh(rmulti *y, rmulti *x);   // y=acosh(x)
void ratanh(rmulti *y, rmulti *x);   // y=atanh(x)

/*
 * operatior for two arguments
 */
void radd(rmulti *z, rmulti *x, rmulti *y);        // z=x+y
void radd_d(rmulti *z, rmulti *x, double y);       // z=x+y
void radd_ui(rmulti *z, rmulti *x, ulong y);       // z=x+y
void radd_si(rmulti *z, rmulti *x, long int y);    // z=x+y
void rsub(rmulti *z, rmulti *x, rmulti *y);        // z=x-y
void rsub_d1(rmulti *z, double x, rmulti *y);      // z=x-y
void rsub_d2(rmulti *z, rmulti *x, double y);      // z=x-y
void rsub_ui1(rmulti *z, ulong x, rmulti *y);      // z=x-y
void rsub_ui2(rmulti *z, rmulti *x, ulong y);      // z=x-y
void rsub_si1(rmulti *z, long x, rmulti *y);       // z=x-y
void rsub_si2(rmulti *z, rmulti *x, long y);       // z=x-y
void rmul(rmulti *z, rmulti *x, rmulti *y);        // z=x*y
void rmul_d(rmulti *z, rmulti *x, double y);       // z=x*y
void rmul_ui(rmulti *z, rmulti *x, ulong y);       // z=x*y
void rmul_si(rmulti *z, rmulti *x, long y);        // z=x*y
void radd_mul(rmulti *z, rmulti *x, rmulti *y);    // z+=x*y
void radd_mul_ws(rmulti *z, rmulti *x, rmulti *y, int n_rws, rmulti **rws);  // n_rws>=1
void radd_mul_d(rmulti *z, rmulti *x, double y);   // z+=x*y
void rsub_mul(rmulti *z, rmulti *x, rmulti *y);    // z-=x*y
void rsub_mul_ws(rmulti *z, rmulti *x, rmulti *y, int n_rws, rmulti **rws);  // n_rws>=1
void rsub_mul_d(rmulti *z, rmulti *x, double y);   // z-=x*y
void rdiv(rmulti *z, rmulti *x, rmulti *y);        // z=x/y
void rdiv_d1(rmulti *z, double x, rmulti *y);      // z=x/y
void rdiv_d2(rmulti *z, rmulti *x, double y);      // z=x/y
void rdiv_ui2(rmulti *z, rmulti *x, ulong y);      // z=x/y
void rdiv_ui1(rmulti *z, ulong x, rmulti *y);      // z=x/y
void rdiv_si2(rmulti *z, rmulti *x, long y);       // z=x/y
void rdiv_si1(rmulti *z, long x, rmulti *y);       // z=x/y
void radd_abs_abs(rmulti *z, rmulti *x, rmulti *y); // z=abs(x)+abs(y)
void rabs_sub(rmulti *z, rmulti *x, rmulti *y);    // z=abs(x-y)
void rabs_sub_d(rmulti *z, rmulti *x, double y);   // z=abs(x-y)
void rdiv_abs(rmulti *z, rmulti *x, rmulti *y);    // z=x/abs(y)
void rpow(rmulti *z, rmulti *x, rmulti *y);        // z=x^y
void rpow_d2(rmulti *z, rmulti *x, double y);      // z=x^y
void rpow_ui(rmulti *z, rmulti *x, ulong y);       // z=x^y
void rpow_si(rmulti *z, rmulti *x, long y);        // z=x^y
void rmul_2exp(rmulti *y, rmulti *x, int n);       // y=x*2^n
void rdiv_2exp(rmulti *y, rmulti *x, int n);       // y=x/2^n
void ratan2(rmulti *z, rmulti *x, rmulti *y);      // z=atan(x/y)
void rmax2_up(rmulti *z, rmulti *x, rmulti *y);    // z=max2(x,y)
void rmin2_down(rmulti *z, rmulti *x, rmulti *y);  // z=min2(x,y)
void rexp10_floor_log10_abs_sub(rmulti *z, rmulti *x, double y); // z=10^(floor(log10(abs(x))-y))
void rexp2_floor_log2_abs_sub(rmulti *z, rmulti *x, double y);   // z=2^(floor(log2(abs(x))-y))
int radd_exact(rmulti *z, rmulti *x, rmulti *y);         // z=x+y
int rsub_exact(rmulti *z, rmulti *x, rmulti *y);         // z=x-y
int rmul_exact(rmulti *z, rmulti *x, rmulti *y);         // z=x*y
int radd_mul_exact(rmulti *z, rmulti *x, rmulti *y);     // z+=x*y
int rsub_mul_exact(rmulti *z, rmulti *x, rmulti *y);     // z-=x*y
int rdiv_rouding_check(rmulti *z, rmulti *x, rmulti *y); // z=x/y

/*
 * operatior for tree arguments
 */
void radd_abs_abs_abs(rmulti *z, rmulti *x0, rmulti *x1, rmulti *x2); // z=abs(x0)+abs(x1)+abs(x2)

/*
 * comparisions
 */
int rcmp(rmulti *x, rmulti *y);  // x<=>y
int rcmp_ui(rmulti *x, ulong y); // x<=>y
int rcmp_si(rmulti *x, long y);  // x<=>y
int rcmp_d(rmulti *x, double y); // x<=>y
int req(rmulti *x, rmulti *y);   // x==y
int req_d(rmulti *x, double y);  // x==y
int req_ui(rmulti *x, ulong y);  // x==y
int req_si(rmulti *x, long y);   // x==y
int rne(rmulti *x, rmulti *y);   // x!=y
int rne_d(rmulti *x, double y);  // x!=y
int rne_ui(rmulti *x, ulong y);  // x!=y
int rne_si(rmulti *x, long y);   // x!=y
int rgt(rmulti *x, rmulti *y);   // x>y
int rgt_d1(double x, rmulti *y); // x>y
int rgt_d2(rmulti *x, double y); // x>y
int rgt_ui1(ulong x, rmulti *y); // x>y
int rgt_ui2(rmulti *x, ulong y); // x>y
int rgt_si1(long x, rmulti *y);  // x>y
int rgt_si2(rmulti *x, long y);  // x>y
int rge(rmulti *x, rmulti *y);   // x>=y
int rge_d1(double x, rmulti *y); // x>=y
int rge_d2(rmulti *x, double y); // x>=y
int rge_ui1(ulong x, rmulti *y); // x>=y
int rge_ui2(rmulti *x, ulong y); // x>=y
int rge_si1(long x, rmulti *y);  // x>=y
int rge_si2(rmulti *x, long y);  // x>=y
int rlt(rmulti *x, rmulti *y);   // x<y
int rlt_d1(double x, rmulti *y); // x<y
int rlt_d2(rmulti *x, double y); // x<y
int rlt_ui1(ulong x, rmulti *y); // x<y
int rlt_ui2(rmulti *x, ulong y); // x<y
int rlt_si1(long x, rmulti *y);  // x<y
int rlt_si2(rmulti *x, long y);  // x<y
int rle(rmulti *x, rmulti *y);   // x<=y
int rle_d1(double x, rmulti *y); // x<=y
int rle_d2(rmulti *x, double y); // x<=y
int rle_ui1(ulong x, rmulti *y); // x<=y
int rle_ui2(rmulti *x, ulong y); // x<=y
int rle_si1(long x, rmulti *y);  // x<=y
int rle_si2(rmulti *x, long y);  // x<=y
int rabs_cmp(rmulti *x, rmulti *y); // abs(x)<=>abs(y)
int rabs_eq(rmulti *x, rmulti *y);  // abs(x)==abs(y)
int rabs_gt(rmulti *x, rmulti *y);  // abs(x)> abs(y)
int rabs_ge(rmulti *x, rmulti *y);  // abs(x)>=abs(y)
int rabs_lt(rmulti *x, rmulti *y);  // abs(x)< xabs(y)
int rabs_le(rmulti *x, rmulti *y);  // abs(x)<=abs(y)
int rdist_cmp(rmulti *x, rmulti *y, rmulti *z); // abs(x-y)<=>z
int rdist_eq(rmulti *x, rmulti *y, rmulti *z);  // abs(x-y)==z
int rdist_ge(rmulti *x, rmulti *y, rmulti *z);  // abs(x-y)>=z
int rdist_gt(rmulti *x, rmulti *y, rmulti *z);  // abs(x-y)> z
int rdist_le(rmulti *x, rmulti *y, rmulti *z);  // abs(x-y)<=z
int rdist_lt(rmulti *x, rmulti *y, rmulti *z);  // abs(x-y)< z
int rdist_cmp_d(rmulti *x, double y, double z); // abs(x-y)<=>z
int rdist_eq_d(rmulti *x, double y, double z);  // abs(x-y)==z
int rdist_ge_d(rmulti *x, double y, double z);  // abs(x-y)>=z
int rdist_gt_d(rmulti *x, double y, double z);  // abs(x-y)> z
int rdist_le_d(rmulti *x, double y, double z);  // abs(x-y)<=z
int rdist_lt_d(rmulti *x, double y, double z);  // abs(x-y)< z

#endif
