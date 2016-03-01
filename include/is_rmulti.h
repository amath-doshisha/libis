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
void set_auto_prec_enabled(void);
void set_auto_prec_disabled(void);
void set_auto_prec_mode(int mode);
int get_auto_prec_mode(void);
int __get_next_prec(int prec);

/*
 * constructions and destruction
 */
rmulti *rallocate(void);
rmulti *rallocate_prec(int prec);
rmulti *rallocate_clone(rmulti *y);
rmulti *rfree(rmulti *x);
int rround(rmulti *x, int prec);
int rclone(rmulti *y, rmulti *x); // y=x
void rswap(rmulti *x, rmulti *y); // x<->y
rmulti *rmepsilon(int prec);

/*
 * rouding mode
 */
mpfr_rnd_t get_round_mode(void);
void set_round_mode(mpfr_rnd_t default_round);

/*
 * member variables
 */
int rget_prec(rmulti *x);
int rget_exp(rmulti *x);
int rget_sgn(rmulti *x);
int rget_sgn_not0(rmulti *x);
int rget_size(rmulti *x);
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
int rset_d(rmulti *x, double value);
int rset_ui(rmulti *x, ulong value);
int rset_si(rmulti *x, long value);
void rset_inf(rmulti *x, int sgn);
void rset_nan(rmulti *x);
int rset_zero(rmulti *x);
int rset_one(rmulti *x);
int rset_one_neg(rmulti *x);
void rset_rand(rmulti *x);

/*
 * casting
 */
double rget_d(rmulti *x);
ulong rget_ui(rmulti *value);
long rget_si(rmulti *value);

/*
 * operations with auto precision mode
 */
int rcopy(rmulti *y, rmulti *x);                // y=x
int rmul_2exp(rmulti *y, rmulti *x, int n);     // y=x*2^n
int rdiv_2exp(rmulti *y, rmulti *x, int n);     // y=x/2^n
int rneg(rmulti *y, rmulti *x);                 // y=-x
int rabs(rmulti *y, rmulti *x);                 // y=abs(x)
int radd(rmulti *z, rmulti *x, rmulti *y);      // z=x+y
int radd_d(rmulti *z, rmulti *x, double y);     // z=x+y
int radd_ui(rmulti *z, rmulti *x, ulong y);     // z=x+y
int radd_si(rmulti *z, rmulti *x, long int y);  // z=x+y
int rsub(rmulti *z, rmulti *x, rmulti *y);      // z=x-y
int rsub_d1(rmulti *z, double x, rmulti *y);    // z=x-y
int rsub_d2(rmulti *z, rmulti *x, double y);    // z=x-y
int rsub_ui1(rmulti *z, ulong x, rmulti *y);    // z=x-y
int rsub_ui2(rmulti *z, rmulti *x, ulong y);    // z=x-y
int rsub_si1(rmulti *z, long x, rmulti *y);     // z=x-y
int rsub_si2(rmulti *z, rmulti *x, long y);     // z=x-y
int rmul(rmulti *z, rmulti *x, rmulti *y);      // z=x*y
int rmul_d(rmulti *z, rmulti *x, double y);     // z=x*y
int rmul_ui(rmulti *z, rmulti *x, ulong y);     // z=x*y
int rmul_si(rmulti *z, rmulti *x, long y);      // z=x*y
int radd_mul(rmulti *z, rmulti *x, rmulti *y);  // z+=x*y
int radd_mul_d(rmulti *z, rmulti *x, double y); // z+=x*y
int rsub_mul(rmulti *z, rmulti *x, rmulti *y);  // z-=x*y
int rsub_mul_d(rmulti *z, rmulti *x, double y); // z-=x*y
int radd_abs(rmulti *y, rmulti *x);             // y+=abs(x)
int rabs_sub(rmulti *z, rmulti *x, rmulti *y);  // z=abs(x-y)
int rabs_sub_d(rmulti *z, rmulti *x, double y); // z=abs(x-y)
int rdiv_abs(rmulti *z, rmulti *x, rmulti *y);  // z=x/abs(y)
int rpow_ui(rmulti *z, rmulti *x, ulong y);     // z=x^y
int rsum_abs_r2(rmulti *z, rmulti *x, rmulti *y); // z=abs(x)+abs(y)
int rsum_abs_r3(rmulti *z, rmulti *x0, rmulti *x1, rmulti *x2); // z=abs(x0)+abs(x1)+abs(x2)

/*
 * operations
 */
int rfloor(rmulti *y, rmulti *x);            // y=floor(x)
int rceil(rmulti *y, rmulti *x);             // y=ceil(x)
int rtrunc(rmulti *y, rmulti *x);            // y=trunc(x)
int rinv(rmulti *z, rmulti *x);              // z=1/x
int rinv_d(rmulti *z, double x);             // z=1/x
int rdiv(rmulti *z, rmulti *x, rmulti *y);   // z=x/y
int rdiv_d1(rmulti *z, double x, rmulti *y); // z=x/y
int rdiv_d2(rmulti *z, rmulti *x, double y); // z=x/y
int rdiv_ui2(rmulti *z, rmulti *x, ulong y); // z=x/y
int rdiv_ui1(rmulti *z, ulong x, rmulti *y); // z=x/y
int rdiv_si2(rmulti *z, rmulti *x, long y);  // z=x/y
int rdiv_si1(rmulti *z, long x, rmulti *y);  // z=x/y
int rpow_si(rmulti *z, rmulti *x, long y);   // z=x^y
int rpow(rmulti *z, rmulti *x, rmulti *y);   // z=x^y
int rpow_d2(rmulti *z, rmulti *x, double y); // z=x^y
int rsqrt(rmulti *y, rmulti *x);             // y=sqrt(x)
int rexp(rmulti *y, rmulti *x);              // y=exp(x)
int rexp2(rmulti *y, rmulti *x);             // y=exp2(x)
int rexp10(rmulti *y, rmulti *x);            // y=exp10(x)
int rlog(rmulti *y, rmulti *x);              // y=log(x)
int rlog2(rmulti *y, rmulti *x);             // y=log2(x)
int rlog10(rmulti *y, rmulti *x);            // y=log10(x)
int rsin(rmulti *y, rmulti *x);              // y=sin(x)
int rcos(rmulti *y, rmulti *x);              // y=cos(x)
int rtan(rmulti *y, rmulti *x);              // y=tan(x)
int rsec(rmulti *y, rmulti *x);              // y=sec(x)
int rcsc(rmulti *y, rmulti *x);              // y=csc(x)
int rcot(rmulti *y, rmulti *x);              // y=cot(x)
int rasin(rmulti *y, rmulti *x);             // y=asin(x)
int racos(rmulti *y, rmulti *x);             // y=acos(x)
int ratan(rmulti *y, rmulti *x);             // y=atan(x)
int ratan2(rmulti *z, rmulti *x, rmulti *y); // z=atan(x/y)
int rsinh(rmulti *y, rmulti *x);             // y=sinh(x)
int rcosh(rmulti *y, rmulti *x);             // y=cosh(x)
int rtanh(rmulti *y, rmulti *x);             // y=tanh(x)
int rasinh(rmulti *y, rmulti *x);            // y=asinh(x)
int racosh(rmulti *y, rmulti *x);            // y=acosh(x)
int ratanh(rmulti *y, rmulti *x);            // y=atanh(x)

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
