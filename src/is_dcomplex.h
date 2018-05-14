#ifndef IS_COMPLEX_H
#define IS_COMPLEX_H

/**
 @file   is_dcomplex.h
 @brief  倍精度複素数構造体dcomplex
 */

/**
 @brief 倍精度複素数構造体dcomplex
 */
typedef struct {
  double r; /**< 虚部 */
  double i; /**< 実部 */
} dcomplex;


/*
 * query
 */
int dis_nan(double x);       // x==nan
int dis_inf(double x);       // x==inf
int dis_zero(double x);      // x==0
int dis_positive(double x);  // x > 0
int dis_negative(double x);  // x < 0
int dis_number(double x);    // number?
int dis_integer(double x);   // integer?
int zis_nan(dcomplex z);              // z==nan
int zis_inf(dcomplex z);              // z==inf
int zis_zero(dcomplex z);             // z==0
int zis_number(dcomplex z);           // number?
int zis_real(dcomplex z);             // real number?
int zis_pure_imaginary(dcomplex z);   // pure imaginary number?
int zis_imaginary(dcomplex z);        // imaginary number?

/*
 * setting
 */
dcomplex zset_polar(double r, double theta);           // z=r*exp(i*theta)
void zget_polar(double *r, double *theta, dcomplex z); // z=r*exp(i*theta)

/*
 * y=f(x)
 */
// y=conj(x)
dcomplex zconj_z(dcomplex x);
dcomplex zconj_d(double x);
double   dconj_d(double x);
// y=real(x)
double dreal_z(dcomplex x);
double dreal_d(double x);
// y=imag(x)
double dimag_z(dcomplex x);
double dimag_d(double x);
// y=abs(x)
double dabs_z(dcomplex x);
double dabs_d(double x);
// y=abs(x)^2
double dabs2_z(dcomplex x);
double dabs2_d(double x);
// y=abs(x_r)+i*abs(x_i)
dcomplex zabsc_z(dcomplex x);
double   dabsc_d(double x);
// y=max(abs(x_r),abs(x_i))
double dmax_absc_z(dcomplex x);
// y=-x
dcomplex zneg_z(dcomplex x);
dcomplex zneg_d(double x);
double   dneg_d(double x);
// y=1/x
dcomplex zinv_z(dcomplex x);
dcomplex zinv_d(double x);
double   dinv_d(double x);
// y=x/|x|
dcomplex znormalize_z(dcomplex x);
dcomplex dnormalize_d(double x);
// theta=arg(z)
double darg_z(dcomplex z);
double darg_d(double z);
// y=floor(x)
dcomplex zfloor_z(dcomplex x);
dcomplex zfloor_d(double x);
double dfloor_d(double x);
// y=ceil(x)
dcomplex zceil_z(dcomplex x);
dcomplex zceil_d(double x);
double dceil_d(double x);
// y=trunc(x)
dcomplex ztrunc_z(dcomplex x);
double dtrunc_d(double x);
// y=sqrt(x)
dcomplex zsqrt_z(dcomplex x);
dcomplex zsqrt_d(double x);
double dsqrt_d(double x);
// y=exp(x)
dcomplex zexp_z(dcomplex x);
dcomplex zexp_d(double x);
double dexp_d(double x);
// y=exp2(x)
dcomplex zexp2_z(dcomplex x);
dcomplex zexp2_d(double x);
double dexp2_d(double x);
// y=exp10(x)
dcomplex zexp10_z(dcomplex x);
dcomplex zexp10_d(double x);
double dexp10_d(double x);
// y=log(x)
dcomplex zlog_z(dcomplex x);
dcomplex zlog_d(double x);
double dlog_d(double x);
// y=log2(x)
dcomplex zlog2_z(dcomplex x);
dcomplex zlog2_d(double x);
double dlog2_d(double x);
// y=log10(x)
dcomplex zlog10_z(dcomplex x);
dcomplex zlog10_d(double x);
double dlog10_d(double x);
// y=sin(x)
dcomplex zsin_z(dcomplex x);
dcomplex zsin_d(double x);
double dsin_d(double x);
// y=cos(x)
dcomplex zcos_z(dcomplex x);
dcomplex zcos_d(double x);
double dcos_d(double x);
// y=tan(x)
dcomplex ztan_z(dcomplex x);
dcomplex ztan_d(double x);
double dtan_d(double x);
// y=sec(x)
dcomplex zsec_z(dcomplex x);
dcomplex zsec_d(double x);
double dsec_d(double x);
// y=csc(x)
dcomplex zcsc_z(dcomplex x);
dcomplex zcsc_d(double x);
double dcsc_d(double x);
// y=cot(x)
dcomplex zcot_z(dcomplex x);
dcomplex zcot_d(double x);
double dcot_d(double x);
// y=asin(x)
dcomplex zasin_z(dcomplex x);
dcomplex zasin_d(double x);
double dasin_d(double x);
// y=acos(x)
dcomplex zacos_z(dcomplex x);
dcomplex zacos_d(double x);
double dacos_d(double x);
// y=atan(x)
dcomplex zatan_z(dcomplex x);
dcomplex zatan_d(double x);
double datan_d(double x);
// y=sinh(x)
dcomplex zsinh_z(dcomplex x);
dcomplex zsinh_d(double x);
double dsinh_d(double x);
// y=cosh(x)
dcomplex zcosh_z(dcomplex x);
dcomplex zcosh_d(double x);
double dcosh_d(double x);
// y=tanh(x)
dcomplex ztanh_z(dcomplex x);
dcomplex ztanh_d(double x);
double dtanh_d(double x);
// y=asinh(x)
dcomplex zasinh_z(dcomplex x);
dcomplex zasinh_d(double x);
double dasinh_d(double x);
// y=acosh(x)
dcomplex zacosh_z(dcomplex x);
dcomplex zacosh_d(double x);
double dacosh_d(double x);
// y=atanh(x)
dcomplex zatanh_z(dcomplex x);
dcomplex zatanh_d(double x);
double datanh_d(double x);
// z=x^n
dcomplex zpow_z(dcomplex x, int n);
double   dpow_d(double x, int n);

/*
 * z=f(x,y)
 */
// z=x+y
dcomplex zadd_zz(dcomplex x, dcomplex y);
dcomplex zadd_zd(dcomplex x, double y);
dcomplex zadd_dz(double x, dcomplex y);
double   dadd_dd(double x, double y);
// z=x-y
dcomplex zsub_zz(dcomplex x, dcomplex y);
dcomplex zsub_zd(dcomplex x, double y);
dcomplex zsub_dz(double x, dcomplex y);
double   dsub_dd(double x, double y);
// z=x*y
dcomplex zmul_zz(dcomplex x, dcomplex y);
dcomplex zmul_zd(dcomplex x, double y);
dcomplex zmul_dz(double x, dcomplex y);
double   dmul_dd(double x, double y);
// z=conj(x)*y
dcomplex zdot_zz(dcomplex x, dcomplex y);
dcomplex zdot_zd(dcomplex x, double y);
dcomplex zdot_dz(double x, dcomplex y);
double   ddot_dd(double x, double y);
// z=x/y
dcomplex zdiv_zz(dcomplex x, dcomplex y);
dcomplex zdiv_zd(dcomplex x, double y);
dcomplex zdiv_dz(double x, dcomplex y);
double   ddiv_dd(double x, double y);
// z=x^y
dcomplex zpow_zz(dcomplex x, dcomplex y);
dcomplex zpow_zd(dcomplex x, double y);
dcomplex zpow_dz(double x, dcomplex y);
double   dpow_dd(double x, double y);
// z=abs(x-y)
double dabs_sub_zz(dcomplex x, dcomplex y);
double dabs_sub_zd(dcomplex x, double y);
double dabs_sub_dz(double x, dcomplex y);
double dabs_sub_dd(double x, double y);
// z=abs(x)+abs(y)
double dadd_abs_abs_zz(dcomplex x, dcomplex y);
double dadd_abs_abs_zd(dcomplex x, double y);
double dadd_abs_abs_dz(double x, dcomplex y);
double dadd_abs_abs_dd(double x, double y);

/*
 * y=y+f(x)
 */
// y=y+abs(x)
double dadd_abs_z(double y, dcomplex x);
double dadd_abs_d(double y, double x);
// y=y+abs(x)^2
double dadd_abs2_z(double y, dcomplex x);
double dadd_abs2_d(double y, double x);

/*
 * z=z+f(x,y)
 */
// z=z+x*y
dcomplex zadd_mul_zz(dcomplex z, dcomplex x, dcomplex y);
dcomplex zadd_mul_zd(dcomplex z, dcomplex x, double y);
dcomplex zadd_mul_dz(dcomplex z, double x, dcomplex y);
dcomplex dadd_mul_dd(double z, dcomplex x, dcomplex y);
// z=z-x*y
dcomplex zsub_mul_zz(dcomplex z, dcomplex x, dcomplex y);
dcomplex zsub_mul_zd(dcomplex z, dcomplex x, double y);
dcomplex zsub_mul_dz(dcomplex z, double x, dcomplex y);
dcomplex dsub_mul_dd(double z, dcomplex x, dcomplex y);
// z=z+conj(x)*y
dcomplex zadd_dot_zz(dcomplex z, dcomplex x, dcomplex y);
dcomplex zadd_dot_zd(dcomplex z, dcomplex x, double y);
dcomplex zadd_dot_dz(dcomplex z, double x, dcomplex y);
dcomplex dadd_dot_dd(double z, dcomplex x, dcomplex y);
// z=z-conj(x)*y
dcomplex zsub_dot_zz(dcomplex z, dcomplex x, dcomplex y);
dcomplex zsub_dot_zd(dcomplex z, dcomplex x, double y);
dcomplex zsub_dot_dz(dcomplex z, double x, dcomplex y);
dcomplex dsub_dot_dd(double z, dcomplex x, dcomplex y);

/*
 * x<=>y
 */
void zcmp_set_real_imag(void); // dcomplex型の値の比較 x<=>y を実部，虚部で判定
void zcmp_set_abs_arg(void);   // dcomplex型の値の比較 x<=>y を絶対値，偏角で判定
int zcmp_get_type(void);       // dcomplex型の値の比較 x<=>y の方法の取得
int cmp_zz(dcomplex x, dcomplex y); // x<=>y
int cmp_zd(dcomplex x, double y);   // x<=>y
int cmp_dz(double x, dcomplex y);   // x<=>y
int cmp_dd(double x, double y);     // x<=>y
int eq_zz(dcomplex x, dcomplex y);  // x==y
int eq_zd(dcomplex x, double y);    // x==y
int eq_dz(double x, dcomplex y);    // x==y
int eq_dd(double x, double y);      // x==y
int ne_zz(dcomplex x, dcomplex y);  // x!=y
int ne_zd(dcomplex x, double y);    // x!=y
int ne_dz(double x, dcomplex y);    // x!=y
int ne_dd(double x, double y);      // x!=y
int gt_zz(dcomplex x, dcomplex y);  // x>y
int gt_zd(dcomplex x, double y);    // x>y
int gt_dz(double x, dcomplex y);    // x>y
int gt_dd(double x, double y);      // x>y
int ge_zz(dcomplex x, dcomplex y);  // x>=y
int ge_zd(dcomplex x, double y);    // x>=y
int ge_dz(double x, dcomplex y);    // x>=y
int ge_dd(double x, double y);      // x>=y
int lt_zz(dcomplex x, dcomplex y);  // x<y
int lt_zd(dcomplex x, double y);    // x<y
int lt_dz(double x, dcomplex y);    // x<y
int lt_dd(double x, double y);      // x<y
int le_zz(dcomplex x, dcomplex y);  // x<=y
int le_zd(dcomplex x, double y);    // x<=y
int le_dz(double x, dcomplex y);    // x<=y
int le_dd(double x, double y);      // x<=y

#endif
