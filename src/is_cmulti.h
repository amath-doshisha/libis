#ifndef IS_CMULTI_H
#define IS_CMULTI_H

#include<is_macros.h>
#include<is_dcomplex.h>
#include<is_rmulti.h>

/**
 @file  is_cmulti.h
 @brief 多倍長精度複素数型cmultiに関する宣言
 */

/**
 @brief 多倍長精度複素数構造体
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
cmulti *callocate_prec2(int prec_r, int prec_i);
cmulti *callocate_clone(cmulti *y);
cmulti *callocate_clone_r(rmulti *y);
cmulti *cfree(cmulti *x);
void cround(cmulti *x, int prec);
void cclone(cmulti *y, cmulti *x);
void cclone_r(cmulti *y, rmulti *x);
void cclone_rr(cmulti *y, rmulti *x_r, rmulti *x_i);
void ccopy(cmulti *y, cmulti *x);                     // y=x
void ccopy_r(cmulti *y, rmulti *x_r);                 // y=x
void ccopy_rr(cmulti *y, rmulti *x_r, rmulti *x_i);   // y=x
void cswap(cmulti *x, cmulti *y); 

/*
 * member variables
 */
int cget_prec(cmulti *x);
int cget_exp(cmulti *x); // max exponent


/*
 * type query
 */
int cis_nan(cmulti *x); // nan?
int cis_inf(cmulti *x); // inf?
int cis_number(cmulti *x); // number?
int cis_zero(cmulti *x); // zeror?
int cis_real(cmulti *x);
int cis_pure_imaginary(cmulti *x);
int cis_imaginary(cmulti *x);

/*
 * I/O
 */
void cbin_save(cmulti *x, FILE *fid);
cmulti *cbin_load(FILE *fid);

/*
 * setting
 */
void cset_s(cmulti *x, const char *str_real);
void cset_ss(cmulti *x, const char *str_real, const char *str_imag);
void cset_script(cmulti *x, const char *str);
void cset_z(cmulti *x, dcomplex value);
void cset_dd(cmulti *x, double value_r, double value_i);
void cset_d(cmulti *x, double value_r);
void cset_ui(cmulti *x, ulong real);
void cset_si(cmulti *x, long real);
void cset_inf(cmulti *x, int sgn_r, int sgn_i);
void cset_nan(cmulti *x);
void cset_zero(cmulti *x);
void cset_one(cmulti *x);
void cset_one_neg(cmulti *x);
void cset_rand(cmulti *x);

/*
 * casting
 */
dcomplex cget_z(cmulti *x);
double cget_d(cmulti *x);
ulong cget_ui(cmulti *x);
long cget_si(cmulti *x);

/*
 * operatior for one argument
 */
void cneg(cmulti *y, cmulti *x);                      // y=-x
void cinv(cmulti *y, cmulti *x);                      // y=1/x
void cinv_ws(cmulti *z, cmulti *x, int n_rws, rmulti **rws, int n_cws, cmulti **cws); // n_rws>=2,n_cws>=1
void cconj(cmulti *y, cmulti *x);                     // y=conj(x)
void cconj_clone(cmulti *y, cmulti *x);               // y=conj(x)
void cabsv(rmulti *y, cmulti *x);                     // y=abs(x)
void cabsv_ws(rmulti *y, cmulti *x, int n_rws, rmulti **rws); // n_rws>=1
void cadd_abs(rmulti *y, cmulti *x);                  // y+=abs(x)
void cabs2(rmulti *y, cmulti *x);                     // y=abs(x)^2
void cabs2_ws(rmulti *y, cmulti *x, int n_rws, rmulti **rws); // n_rws>=1
void cadd_abs2(rmulti *y, cmulti *x);                 // y+=abs(x)^2
void cabsc(cmulti *y, cmulti *x);                     // y=abs(x.r)+i*abs(r.i)
void cmax_absc(rmulti *y, cmulti *x);                 // y=max(abs(x.r),abs(r.i))
void cnormalize(cmulti *y, cmulti *x);                // y=x/abs(x)
void cargument(rmulti *theta, cmulti *z);             // theta=arg(z)
void csqrt_r(cmulti *y, rmulti *x);                   // y=sqrt(x)
void csqrt_c(cmulti *y, cmulti *x);                   // y=sqrt(x)
void cexp_c(cmulti *y, cmulti *x);                    // y=exp(x)
void clog_r(cmulti *y, rmulti *x);                    // y=log(x)
void clog_c(cmulti *y, cmulti *x);                    // y=log(x)
void csin_c(cmulti *y, cmulti *x);                    // y=sin(x)
void ccos_c(cmulti *y, cmulti *x);                    // y=cos(x)
void ctan_c(cmulti *y, cmulti *x);                    // y=tan(x)
void casin_r(cmulti *y, rmulti *x);                   // y=asin(x)
void casin_c(cmulti *y, cmulti *x);                   // y=asin(x)
void cacos_r(cmulti *y, rmulti *x);                   // y=cos(x)
void cacos_c(cmulti *y, cmulti *x);                   // y=cos(x)
void catan_r(cmulti *y, rmulti *x);                   // y=tan(x)
void catan_c(cmulti *y, cmulti *x);                   // y=tan(x)
void csinh_c(cmulti *y, cmulti *x);                   // y=sinh(x)
void ccosh_c(cmulti *y, cmulti *x);                   // y=cosh(x)
void ctanh_c(cmulti *y, cmulti *x);                   // y=tanh(x)
void casinh_r(cmulti *y, rmulti *x);                  // y=asinh(x)
void casinh_c(cmulti *y, cmulti *x);                  // y=asinh(x)
void cacosh_r(cmulti *y, rmulti *x);                  // y=acosh(x)
void cacosh_c(cmulti *y, cmulti *x);                  // y=acosh(x)
void catanh_r(cmulti *y, rmulti *x);                  // y=atanh(x)
void catanh_c(cmulti *y, cmulti *x);                  // y=atanh(x)

/*
 * operatior for two arguments
 */
void cadd(cmulti *z, cmulti *x, cmulti *y);           // z=x+y
void cadd_z(cmulti *z, cmulti *x, dcomplex y);        // z=x+y
void cadd_r(cmulti *z, cmulti *x, rmulti *y);         // z=x+y
void cadd_d(cmulti *z, cmulti *x, double y);          // z=x+y
void cadd_ui(cmulti *z, cmulti *x, ulong y);          // z=x+y
void cadd_si(cmulti *z, cmulti *x, long y);           // z=x+y
void csub(cmulti *z, cmulti *x, cmulti *y);           // z=x-y
void csub_z1(cmulti *z, dcomplex x, cmulti *y);       // z=x-y
void csub_z2(cmulti *z, cmulti *x, dcomplex y);       // z=x-y
void csub_r1(cmulti *z, rmulti *x, cmulti *y);        // z=x-y
void csub_r2(cmulti *z, cmulti *x, rmulti *y);        // z=x-y
void csub_d1(cmulti *z, double x, cmulti *y);         // z=x-y
void csub_d2(cmulti *z, cmulti *x, double y);         // z=x-y
void csub_ui1(cmulti *z, ulong x, cmulti *y);         // z=x-y
void csub_ui2(cmulti *z, cmulti *x, ulong y);         // z=x-y
void csub_si1(cmulti *z, long x, cmulti *y);          // z=x-y
void csub_si2(cmulti *z, cmulti *x, long y);          // z=x-y
void cmul(cmulti *z, cmulti *x, cmulti *y);           // z=x*y
void cmul_ws(cmulti *z, cmulti *x, cmulti *y, int n_rws, rmulti **rws, int n_cws, cmulti **cws);  // n_rws>=1,n_cws>=1
void cmul_z(cmulti *z, cmulti *x, dcomplex y);        // z=x*y
void cmul_r(cmulti *z, cmulti *x, rmulti *y);         // z=x*y
void cmul_d(cmulti *z, cmulti *x, double y);          // z=x*y
void cmul_ui(cmulti *z, cmulti *x, ulong y);          // z=x*y
void cmul_si(cmulti *z, cmulti *x, long y);           // z=x*y
void cdot(cmulti *z, cmulti *x, cmulti *y);           // z=conj(x)*y
void cdot_z1(cmulti *z, dcomplex x, cmulti *y);       // z=conj(x)*y
void cdot_z2(cmulti *z, cmulti *x, dcomplex y);       // z=conj(x)*y
void cadd_mul(cmulti *z, cmulti *x, cmulti *y);       // z+=x*y
void cadd_mul_r(cmulti *z, cmulti *x, rmulti *y);     // z+=x*y
void cadd_mul_z(cmulti *z, cmulti *x, dcomplex y);    // z+=x*y
void cadd_mul_d(cmulti *z, cmulti *x, double y);      // z+=x*y
void csub_mul(cmulti *z, cmulti *x, cmulti *y);       // z-=x*y
void csub_mul_ws(cmulti *z, cmulti *x, cmulti *y, int n_rws, rmulti **rws, int n_cws, cmulti **cws);  // n_rws>=1,n_cws>=2
void csub_mul_r(cmulti *z, cmulti *x, rmulti *y);     // z-=x*y
void csub_mul_z(cmulti *z, cmulti *x, dcomplex y);    // z-=x*y
void csub_mul_d(cmulti *z, cmulti *x, double y);      // z-=x*y
void cadd_dot(cmulti *z, cmulti *x, cmulti *y);       // z+=conj(x)*y
void cadd_dot_z1(cmulti *z, dcomplex x, cmulti *y);   // z+=conj(x)*y
void cadd_dot_z2(cmulti *z, cmulti *x, dcomplex y);   // z+=conj(x)*y
void csub_dot(cmulti *z, cmulti *x, cmulti *y);       // z-=conj(x)*y
void csub_dot_z1(cmulti *z, dcomplex x, cmulti *y);   // z-=conj(x)*y
void csub_dot_z2(cmulti *z, cmulti *x, dcomplex y);   // z-=conj(x)*y
void cdiv(cmulti *z, cmulti *x, cmulti *y);           // z=x/y
void cdiv_z1(cmulti *z, dcomplex x, cmulti *y);       // z=x/y
void cdiv_z2(cmulti *z, cmulti *x, dcomplex y);       // z=x/y
void cdiv_r1(cmulti *z, rmulti *x, cmulti *y);        // z=x/y
void cdiv_r2(cmulti *z, cmulti *x, rmulti *y);        // z=x/y
void cdiv_d1(cmulti *z, double x, cmulti *y);         // z=x/y
void cdiv_d2(cmulti *z, cmulti *x, double y);         // z=x/y
void cdiv_ui1(cmulti *z, ulong x, cmulti *y);         // z=x/y
void cdiv_ui2(cmulti *z, cmulti *x, ulong y);         // z=x/y
void cdiv_si1(cmulti *z, long x, cmulti *y);          // z=x/y
void cdiv_si2(cmulti *z, cmulti *x, long y);          // z=x/y
void cabs_sub(rmulti *z, cmulti *x, cmulti *y);       // z=abs(x-y)
void cabs_sub_r(rmulti *z, cmulti *x, rmulti *y);     // z=abs(x-y)
void cdiv_abs(rmulti *z, rmulti *x, cmulti *y);       // z=x/abs(y)
void cpow_c(cmulti *z, cmulti *x, cmulti *y);         // z=x^y
void cpow_r1(cmulti *z, rmulti *x, cmulti *y);        // z=x^y
void cpow_r2(cmulti *z, cmulti *x, rmulti *y);        // z=x^y
void cpow_d2(cmulti *z, cmulti *x, double y);         // z=x^y
void cpow_si(cmulti *y, cmulti *x, long n);           // y=x^n
void cpow_ui(cmulti *y, cmulti *x, ulong n);          // y=x^n
void cmul_2exp(cmulti *y, cmulti *x, int nr, int ni); // y=x*2^n
void cdiv_2exp(cmulti *y, cmulti *x, int nr, int ni); // y=x/2^n
void cexp10_floor_log10_abs_sub(cmulti *z, cmulti *x, double y); // z=10^(floor(log10(abs(x))-y))
void cexp2_floor_log2_abs_sub(cmulti *z, cmulti *x, double y);   // z=2^(floor(log2(abs(x))-y))
void cget_polar(rmulti *r, rmulti *theta, cmulti *z); // z=r*exp(i*theta)
void cset_polar(cmulti *z, rmulti *r, rmulti *theta); // z=r*exp(i*theta)

/*
 * comparisions
 */
void ccmp_set_real_imag();
void ccmp_set_abs_arg();
int ccmp_get_type();
int ccmp(cmulti *x, cmulti *y);       // x<=>y
int ccmp_z(cmulti *x, dcomplex y);    // x<=>y
int ccmp_r1(rmulti *x, cmulti *y);    // x<=>y
int ccmp_r2(cmulti *x, rmulti *y);    // x<=>y
int ccmp_d(cmulti *x, double y);      // x<=>y
int ccmp_ui(cmulti *x, ulong y);      // x<=>y
int ccmp_si(cmulti *x, long y);       // x<=>y
int ceq(cmulti *x, cmulti *y);        // x==y
int ceq_z(cmulti *x, dcomplex y);     // x==y
int ceq_r(cmulti *x, rmulti *y);      // x==y
int ceq_d(cmulti *x, double y);       // x==y
int ceq_ui(cmulti *x, ulong y);       // x==y
int ceq_si(cmulti *x, long y);        // x==y
int cne(cmulti *x, cmulti *y);        // x!=y
int cne_z(cmulti *x, dcomplex y);     // x!=y
int cne_r(cmulti *x, rmulti *y);      // x!=y
int cne_d(cmulti *x, double y);       // x!=y
int cne_ui(cmulti *x, ulong y);       // x!=y
int cne_si(cmulti *x, long y);        // x!=y
int cgt(cmulti *x, cmulti *y);        // x>y
int cgt_r1(rmulti *x, cmulti *y);     // x>y
int cgt_r2(cmulti *x, rmulti *y);     // x>y
int cgt_z1(dcomplex x, cmulti *y);    // x>y
int cgt_z2(cmulti *x, dcomplex y);    // x>y
int cgt_d1(double x, cmulti *y);      // x>y
int cgt_d2(cmulti *x, double y);      // x>y
int cgt_ui1(ulong x, cmulti *y);      // x>y
int cgt_ui2(cmulti *x, ulong y);      // x>y
int cgt_si1(long x, cmulti *y);       // x>y
int cgt_si2(cmulti *x, long y);       // x>y
int cge(cmulti *x, cmulti *y);        // x>=y
int cge_r1(rmulti *x, cmulti *y);     // x>=y
int cge_r2(cmulti *x, rmulti *y);     // x>=y
int cge_z1(dcomplex x, cmulti *y);    // x>=y
int cge_z2(cmulti *x, dcomplex y);    // x>=y
int cge_d1(double x, cmulti *y);      // x>=y
int cge_d2(cmulti *x, double y);      // x>=y
int cge_ui1(ulong x, cmulti *y);      // x>=y
int cge_ui2(cmulti *x, ulong y);      // x>=y
int cge_si1(long x, cmulti *y);       // x>=y
int cge_si2(cmulti *x, long y);       // x>=y
int clt(cmulti *x, cmulti *y);        // x<y
int clt_r1(rmulti *x, cmulti *y);     // x<y
int clt_r2(cmulti *x, rmulti *y);     // x<y
int clt_z1(dcomplex x, cmulti *y);    // x<y
int clt_z2(cmulti *x, dcomplex y);    // x<y
int clt_d1(double x, cmulti *y);      // x<y
int clt_d2(cmulti *x, double y);      // x<y
int clt_ui1(ulong x, cmulti *y);      // x<y
int clt_ui2(cmulti *x, ulong y);      // x<y
int clt_si1(long x, cmulti *y);       // x<y
int clt_si2(cmulti *x, long y);       // x<y
int cle(cmulti *x, cmulti *y);        // x<=y
int cle_r1(rmulti *x, cmulti *y);     // x<=y
int cle_r2(cmulti *x, rmulti *y);     // x<=y
int cle_z1(dcomplex x, cmulti *y);    // x<=y
int cle_z2(cmulti *x, dcomplex y);    // x<=y
int cle_d1(double x, cmulti *y);      // x<=y
int cle_d2(cmulti *x, double y);      // x<=y
int cle_ui1(ulong x, cmulti *y);      // x<=y
int cle_ui2(cmulti *x, ulong y);      // x<=y
int cle_si1(long x, cmulti *y);       // x<=y
int cle_si2(cmulti *x, long y);       // x<=y
int cabs_cmp(cmulti *x, cmulti *y);   // abs(x)<=>abs(y)
int cabs_eq(cmulti *x, cmulti *y);    // abs(x)==abs(y)
int cabs_gt(cmulti *x, cmulti *y);    // abs(x)>abs(y)
int cabs_ge(cmulti *x, cmulti *y);    // abs(x)>=abs(y)
int cabs_lt(cmulti *x, cmulti *y);    // abs(x)<abs(y)
int cabs_le(cmulti *x, cmulti *y);    // abs(x)<=abs(y)
int ceqc(cmulti *x, cmulti *y);       // x.r==y.r && x.i==y.i
int cgtc(cmulti *x, cmulti *y);       // x.r> y.r && x.i> y.i
int cgec(cmulti *x, cmulti *y);       // x.r>=y.r && x.i>=y.i
int cltc(cmulti *x, cmulti *y);       // x.r< y.r && x.i< y.i
int clec(cmulti *x, cmulti *y);       // x.r<=y.r && x.i<=y.i

#endif

