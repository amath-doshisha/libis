#ifndef ISYS_ICMULTI_H
#define ISYS_ICMULTI_H

#include<is_cmulti.h>

/**
 @file  is_icmulti.h
 @brief 多倍長精度実数型cmultiの機械区間演算
*/

/*
 * casting
 */
void icset_c(cmulti *y0, cmulti *y1, cmulti *x0, cmulti *x1);                                // [cmulti,cmulti] <- [cmulti,cmulti]
void icset_rrrr(cmulti *y0, cmulti *y1, rmulti *x0r, rmulti *x0i, rmulti *x1r, rmulti *x1i); // [cmulti,cmulti] <- [(rmulti,rmulti),(rmulti,rmulti)]
void icset_z(cmulti *y0, cmulti *y1, dcomplex x0, dcomplex x1);                              // [cmulti,cmulti] <- [dcomplex,dcomplex]
void icset_r(cmulti *y0, cmulti *y1, rmulti *x0, rmulti *x1);                                // [cmulti,cmulti] <- [rmulti,rmulti]
void icset_d(cmulti *y0, cmulti *y1, double x0, double x1);                                  // [cmulti,cmulti] <- [double,double]
void icset_cr(cmulti *y0, cmulti *y1, cmulti *x0, rmulti *x1);                               // [cmulti,cmulti] <- [cmulti,rmulti]
void icset_rc(cmulti *y0, cmulti *y1, rmulti *x0, cmulti *x1);                               // [cmulti,cmulti] <- [rmulti,cmulti]
void icset_zd(cmulti *y0, cmulti *y1, dcomplex x0, double x1);                               // [cmulti,cmulti] <- [dcomplex,double]
void icset_dz(cmulti *y0, cmulti *y1, double x0, dcomplex x1);                               // [cmulti,cmulti] <- [double,dcomplex]
void icset_dddd(cmulti *y0, cmulti *y1, double x0r, double x0i, double x1r, double x1i);     // [cmulti,cmulti] <- [(double,double),(double,double)]
void icset_s(cmulti *x0, cmulti *x1, char *s);                                               // [cmulti,cmulti] <- char

/*
 * y=f(x)
 */
// y=conj(x)
void icconj_c(cmulti *y0, cmulti *y1, cmulti *x0, cmulti *x1);
void icconj_r(cmulti *y0, cmulti *y1, rmulti *x0, rmulti *x1);
void icconj_z(cmulti *y0, cmulti *y1, dcomplex x0, dcomplex x1);
void icconj_d(cmulti *y0, cmulti *y1, double x0, double x1);
// y=real(x)
void irreal_c(rmulti *y0, rmulti *y1, cmulti *x0, cmulti *x1);
void irreal_z(rmulti *y0, rmulti *y1, dcomplex x0, dcomplex x1);
// y=imag(x)
void irimag_c(rmulti *y0, rmulti *y1, cmulti *x0, cmulti *x1);
void irimag_z(rmulti *y0, rmulti *y1, dcomplex x0, dcomplex x1);
// y=abs(x)
void irabs_c_ws(rmulti *y0, rmulti *y1, cmulti *x0, cmulti *x1, int n_cws, cmulti **cws); // n_cws>=2
void irabs_c(rmulti *y0, rmulti *y1, cmulti *x0, cmulti *x1);
// y=abs(x)^2
void irabs2_c(rmulti *y0, rmulti *y1, cmulti *x0, cmulti *x1);
// y=abs(x.r)+i*abs(r.i)
void icabsc_c(cmulti *y0, cmulti *y1, cmulti *x0, cmulti *x1);
// y=max(abs(x_r),abs(x_i))
void irmax_absc_c(rmulti *y0, rmulti *y1, cmulti *x0, cmulti *x1);
// y=-x
void icneg_c(cmulti *y0, cmulti *y1, cmulti *x0, cmulti *x1);
void icneg_z(cmulti *y0, cmulti *y1, dcomplex x0, dcomplex x1);
void icneg_r(cmulti *y0, cmulti *y1, rmulti *x0, rmulti *x1);
void icneg_d(cmulti *y0, cmulti *y1, double x0, double x1);
// y=1/x
void icinv_c(cmulti *y0, cmulti *y1, cmulti *x0, cmulti *x1);
void icinv_z(cmulti *y0, cmulti *y1, dcomplex x0, dcomplex x1);
void icinv_r(cmulti *y0, cmulti *y1, rmulti *x0, rmulti *x1);
void icinv_d(cmulti *y0, cmulti *y1, double x0, double x1);
// y=x/abs(x)
void icnormalize_c(cmulti *y0, cmulti *y1, cmulti *x0, cmulti *x1);
void icnormalize_r(cmulti *y0, cmulti *y1, dcomplex x0, dcomplex x1);
void icnormalize_z(cmulti *y0, cmulti *y1, rmulti *x0, rmulti *x1);
void icnormalize_d(cmulti *y0, cmulti *y1, double x0, double x1);
// theta=arg(z)
void irarg_c(rmulti *theta0, rmulti *theta1, cmulti *z0, cmulti *z1);
void irarg_z(rmulti *theta0, rmulti *theta1, dcomplex z0, dcomplex z1);
// y=floor(x)
void icfloor_c(cmulti *y0, cmulti *y1, cmulti *x0, cmulti *x1);
void icfloor_r(cmulti *y0, cmulti *y1, dcomplex x0, dcomplex x1);
void icfloor_z(cmulti *y0, cmulti *y1, rmulti *x0, rmulti *x1);
void icfloor_d(cmulti *y0, cmulti *y1, double x0, double x1);
// y=ceil(x)
void icceil_c(cmulti *y0, cmulti *y1, cmulti *x0, cmulti *x1);
void icceil_r(cmulti *y0, cmulti *y1, dcomplex x0, dcomplex x1);
void icceil_z(cmulti *y0, cmulti *y1, rmulti *x0, rmulti *x1);
void icceil_d(cmulti *y0, cmulti *y1, double x0, double x1);
// y=trunc(x)
void ictrunc_c(cmulti *y0, cmulti *y1, cmulti *x0, cmulti *x1);
void ictrunc_r(cmulti *y0, cmulti *y1, dcomplex x0, dcomplex x1);
void ictrunc_z(cmulti *y0, cmulti *y1, rmulti *x0, rmulti *x1);
void ictrunc_d(cmulti *y0, cmulti *y1, double x0, double x1);
// z=x^n
void icpow_c(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, int n);
void icpow_r(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, int n);
void icpow_z(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, int n);
void icpow_d(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, int n);
// [y0,y1]=[-absc(x),absc(x)]
void icpm(cmulti *y0, cmulti *y1, cmulti *x);
// [m-r,m+r]=[x0,x1]
void icmid(cmulti *mid, cmulti *x0, cmulti *x1);
void icrad(cmulti *rad, cmulti *x0, cmulti *x1);
void icmr(cmulti *mid, cmulti *rad, cmulti *x0, cmulti *x1);

/*
 * z=f(x,y)
 */
// z=x+y
void icadd_cc(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, cmulti *y0, cmulti *y1);
void icadd_cr(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, rmulti *y0, rmulti *y1);
void icadd_rc(cmulti *z0, cmulti *z1, rmulti *x0, rmulti *x1, cmulti *y0, cmulti *y1);
void icadd_cz(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, dcomplex y0, dcomplex y1);
void icadd_zc(cmulti *z0, cmulti *z1, dcomplex x0, dcomplex x1, cmulti *y0, cmulti *y1);
void icadd_cd(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, double y0, double y1);
void icadd_dc(cmulti *z0, cmulti *z1, double x0, double x1, cmulti *y0, cmulti *y1);
void icadd_rz(cmulti *z0, cmulti *z1, rmulti *x0, rmulti *x1, dcomplex y0, dcomplex y1);
void icadd_zr(cmulti *z0, cmulti *z1, dcomplex x0, dcomplex x1, rmulti *y0, rmulti *y1);
// z=x-y
void icsub_cc_ws(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, cmulti *y0, cmulti *y1, int n_rws, rmulti **rws); // n_rws>=2
void icsub_cc(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, cmulti *y0, cmulti *y1);
void icsub_cr(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, rmulti *y0, rmulti *y1);
void icsub_rc(cmulti *z0, cmulti *z1, rmulti *x0, rmulti *x1, cmulti *y0, cmulti *y1);
void icsub_cz(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, dcomplex y0, dcomplex y1);
void icsub_zc(cmulti *z0, cmulti *z1, dcomplex x0, dcomplex x1, cmulti *y0, cmulti *y1);
void icsub_cd(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, double y0, double y1);
void icsub_dc(cmulti *z0, cmulti *z1, double x0, double x1, cmulti *y0, cmulti *y1);
void icsub_rz(cmulti *z0, cmulti *z1, rmulti *x0, rmulti *x1, dcomplex y0, dcomplex y1);
void icsub_zr(cmulti *z0, cmulti *z1, dcomplex x0, dcomplex x1, rmulti *y0, rmulti *y1);
// z=x*y
void icmul_cc_ws(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, cmulti *y0, cmulti *y1, int n_rws, rmulti **rws, int n_cws, cmulti **cws); // n_rws>=8,n_cws>=2
void icmul_cc(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, cmulti *y0, cmulti *y1);
void icmul_cr(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, rmulti *y0, rmulti *y1);
void icmul_rc(cmulti *z0, cmulti *z1, rmulti *x0, rmulti *x1, cmulti *y0, cmulti *y1);
void icmul_cz(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, dcomplex y0, dcomplex y1);
void icmul_zc(cmulti *z0, cmulti *z1, dcomplex x0, dcomplex x1, cmulti *y0, cmulti *y1);
void icmul_cd(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, double y0, double y1);
void icmul_dc(cmulti *z0, cmulti *z1, double x0, double x1, cmulti *y0, cmulti *y1);
void icmul_rz(cmulti *z0, cmulti *z1, rmulti *x0, rmulti *x1, dcomplex y0, dcomplex y1);
void icmul_zr(cmulti *z0, cmulti *z1, dcomplex x0, dcomplex x1, rmulti *y0, rmulti *y1);
// z=conj(x)*y
void icdot_cc_ws(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, cmulti *y0, cmulti *y1, int n_rws, rmulti **rws, int n_cws, cmulti **cws); // n_rws>=8,n_cws>=2
void icdot_cc(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, cmulti *y0, cmulti *y1);
void icdot_cr(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, rmulti *y0, rmulti *y1);
void icdot_rc(cmulti *z0, cmulti *z1, rmulti *x0, rmulti *x1, cmulti *y0, cmulti *y1);
void icdot_cz(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, dcomplex y0, dcomplex y1);
void icdot_zc(cmulti *z0, cmulti *z1, dcomplex x0, dcomplex x1, cmulti *y0, cmulti *y1);
void icdot_cd(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, double y0, double y1);
void icdot_dc(cmulti *z0, cmulti *z1, double x0, double x1, cmulti *y0, cmulti *y1);
void icdot_rz(cmulti *z0, cmulti *z1, rmulti *x0, rmulti *x1, dcomplex y0, dcomplex y1);
void icdot_zr(cmulti *z0, cmulti *z1, dcomplex x0, dcomplex x1, rmulti *y0, rmulti *y1);
// z=x/y
void icdiv_cc(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, cmulti *y0, cmulti *y1);
void icdiv_cr(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, rmulti *y0, rmulti *y1);
void icdiv_rc(cmulti *z0, cmulti *z1, rmulti *x0, rmulti *x1, cmulti *y0, cmulti *y1);
void icdiv_cz(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, dcomplex y0, dcomplex y1);
void icdiv_zc(cmulti *z0, cmulti *z1, dcomplex x0, dcomplex x1, cmulti *y0, cmulti *y1);
void icdiv_cd(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, double y0, double y1);
void icdiv_dc(cmulti *z0, cmulti *z1, double x0, double x1, cmulti *y0, cmulti *y1);
void icdiv_rz(cmulti *z0, cmulti *z1, rmulti *x0, rmulti *x1, dcomplex y0, dcomplex y1);
void icdiv_zr(cmulti *z0, cmulti *z1, dcomplex x0, dcomplex x1, rmulti *y0, rmulti *y1);
// z=abs(x-y)
void irabs_sub_cc(rmulti *z0, rmulti *z1, cmulti *x0, cmulti *x1, cmulti *y0, cmulti *y1);
void irabs_sub_cr(rmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, rmulti *y0, rmulti *y1);
void irabs_sub_rc(rmulti *z0, cmulti *z1, rmulti *x0, rmulti *x1, cmulti *y0, cmulti *y1);
void irabs_sub_cz(rmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, dcomplex y0, dcomplex y1);
void irabs_sub_zc(rmulti *z0, cmulti *z1, dcomplex x0, dcomplex x1, cmulti *y0, cmulti *y1);
void irabs_sub_cd(rmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, double y0, double y1);
void irabs_sub_dc(rmulti *z0, cmulti *z1, double x0, double x1, cmulti *y0, cmulti *y1);
void irabs_sub_rz(rmulti *z0, cmulti *z1, rmulti *x0, rmulti *x1, dcomplex y0, dcomplex y1);
void irabs_sub_zr(rmulti *z0, cmulti *z1, dcomplex x0, dcomplex x1, rmulti *y0, rmulti *y1);
// z=abs(x)+abs(y)
void iradd_abs_abs_cc(rmulti *z0, rmulti *z1, cmulti *x0, cmulti *x1, cmulti *y0, cmulti *y1);
void iradd_abs_abs_cr(rmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, rmulti *y0, rmulti *y1);
void iradd_abs_abs_rc(rmulti *z0, cmulti *z1, rmulti *x0, rmulti *x1, cmulti *y0, cmulti *y1);
void iradd_abs_abs_cz(rmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, dcomplex y0, dcomplex y1);
void iradd_abs_abs_zc(rmulti *z0, cmulti *z1, dcomplex x0, dcomplex x1, cmulti *y0, cmulti *y1);
void iradd_abs_abs_cd(rmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, double y0, double y1);
void iradd_abs_abs_dc(rmulti *z0, cmulti *z1, double x0, double x1, cmulti *y0, cmulti *y1);
void iradd_abs_abs_rz(rmulti *z0, cmulti *z1, rmulti *x0, rmulti *x1, dcomplex y0, dcomplex y1);
void iradd_abs_abs_zr(rmulti *z0, cmulti *z1, dcomplex x0, dcomplex x1, rmulti *y0, rmulti *y1);

/*
 * y=y+f(x)
 */
// y=y+abs(x)
void iradd_abs_c(rmulti *y0, rmulti *y1, cmulti *x0, cmulti *x1);
void iradd_abs_z(rmulti *z0, rmulti *z1, dcomplex x0, dcomplex x1);
// y=y+abs(x)^2
void iradd_abs2_c(rmulti *y0, rmulti *y1, cmulti *x0, cmulti *x1);
void iradd_abs2_z(rmulti *z0, rmulti *z1, dcomplex x0, dcomplex x1);
// [z]=[x]+[-absc(y),absc(y)]
void icadd_pm_cc(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, cmulti *y);

/*
 * z=z+f(x,y)
 */
// z=z+x*y
void icadd_mul_rr(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, cmulti *y0, cmulti *y1);
void icadd_mul_cr(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, rmulti *y0, rmulti *y1);
void icadd_mul_rc(cmulti *z0, cmulti *z1, rmulti *x0, rmulti *x1, cmulti *y0, cmulti *y1);
void icadd_mul_cz(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, dcomplex y0, dcomplex y1);
void icadd_mul_zc(cmulti *z0, cmulti *z1, dcomplex x0, dcomplex x1, cmulti *y0, cmulti *y1);
void icadd_mul_cd(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, double y0, double y1);
void icadd_mul_dc(cmulti *z0, cmulti *z1, double x0, double x1, cmulti *y0, cmulti *y1);
void icadd_mul_rz(cmulti *z0, cmulti *z1, rmulti *x0, rmulti *x1, dcomplex y0, dcomplex y1);
void icadd_mul_zr(cmulti *z0, cmulti *z1, dcomplex x0, dcomplex x1, rmulti *y0, rmulti *y1);
// z=z-x*y
void icsub_mul_cc_ws(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, cmulti *y0, cmulti *y1, int n_rws, rmulti **rws, int n_cws, cmulti **cws); // n_rws>=8,n_cws>=4
void icsub_mul_cc(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, cmulti *y0, cmulti *y1);
void icsub_mul_cr(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, rmulti *y0, rmulti *y1);
void icsub_mul_rc(cmulti *z0, cmulti *z1, rmulti *x0, rmulti *x1, cmulti *y0, cmulti *y1);
void icsub_mul_cz(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, dcomplex y0, dcomplex y1);
void icsub_mul_zc(cmulti *z0, cmulti *z1, dcomplex x0, dcomplex x1, cmulti *y0, cmulti *y1);
void icsub_mul_cd(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, double y0, double y1);
void icsub_mul_dc(cmulti *z0, cmulti *z1, double x0, double x1, cmulti *y0, cmulti *y1);
void icsub_mul_rz(cmulti *z0, cmulti *z1, rmulti *x0, rmulti *x1, dcomplex y0, dcomplex y1);
void icsub_mul_zr(cmulti *z0, cmulti *z1, dcomplex x0, dcomplex x1, rmulti *y0, rmulti *y1);
// z=z+conj(x)*y
void icadd_dot_cc(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, cmulti *y0, cmulti *y1);
void icadd_dot_cr(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, rmulti *y0, rmulti *y1);
void icadd_dot_rc(cmulti *z0, cmulti *z1, rmulti *x0, rmulti *x1, cmulti *y0, cmulti *y1);
void icadd_dot_cz(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, dcomplex y0, dcomplex y1);
void icadd_dot_zc(cmulti *z0, cmulti *z1, dcomplex x0, dcomplex x1, cmulti *y0, cmulti *y1);
void icadd_dot_cd(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, double y0, double y1);
void icadd_dot_dc(cmulti *z0, cmulti *z1, double x0, double x1, cmulti *y0, cmulti *y1);
void icadd_dot_rz(cmulti *z0, cmulti *z1, rmulti *x0, rmulti *x1, dcomplex y0, dcomplex y1);
void icadd_dot_zr(cmulti *z0, cmulti *z1, dcomplex x0, dcomplex x1, rmulti *y0, rmulti *y1);
// z=z-conj(x)*y
void icsub_dot_cc(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, cmulti *y0, cmulti *y1);
void icsub_dot_cr(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, rmulti *y0, rmulti *y1);
void icsub_dot_rc(cmulti *z0, cmulti *z1, rmulti *x0, rmulti *x1, cmulti *y0, cmulti *y1);
void icsub_dot_cz(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, dcomplex y0, dcomplex y1);
void icsub_dot_zc(cmulti *z0, cmulti *z1, dcomplex x0, dcomplex x1, cmulti *y0, cmulti *y1);
void icsub_dot_cd(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, double y0, double y1);
void icsub_dot_dc(cmulti *z0, cmulti *z1, double x0, double x1, cmulti *y0, cmulti *y1);
void icsub_dot_rz(cmulti *z0, cmulti *z1, rmulti *x0, rmulti *x1, dcomplex y0, dcomplex y1);
void icsub_dot_zr(cmulti *z0, cmulti *z1, dcomplex x0, dcomplex x1, rmulti *y0, rmulti *y1);

/*
 * x<=>y
 */
int ic_in_ic(cmulti *y0, cmulti *y1, cmulti *x0, cmulti *x1); // [y0,y1] in [x0,x1]
int ic_in_icpm(cmulti *y0, cmulti *y1, cmulti *x);            // [y0,y1] in [-abs(x),abs(x)]

#endif
