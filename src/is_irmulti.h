#ifndef ISYS_IRMULTI_H
#define ISYS_IRMULTI_H

#include<is_rmulti.h>

/**
 @file  is_irmulti.h
 @brief 多倍長精度実数型rmultiの機械区間演算
*/

/*
 * casting
 */
void irset_r(rmulti *y0, rmulti *y1, rmulti *x0, rmulti *x1);   // [rmulti,rmulti] <- [rmulti,rmulti]
void irset_d(rmulti *y0, rmulti *y1, double x0, double x1);     // [rmulti,rmulti] <- [double,double]
void irset_dr(rmulti *y0, rmulti *y1, double x0, rmulti *x1);   // [rmulti,rmulti] <- [double,rmulti]
void irset_rd(rmulti *y0, rmulti *y1, rmulti *x0, double x1);   // [rmulti,rmulti] <- [rmulti,double]
void irset_s(rmulti *x0, rmulti *x1, char *str);                // [rmulti,rmulti] <- char

/*
 * y=f(x)
 */
// y=conj(x)
void irconj_r(rmulti *y0, rmulti *y1, rmulti *x0, rmulti *x1);
void irconj_d(rmulti *y0, rmulti *y1, double x0, double x1);
// y=real(x)
void irreal_r(rmulti *y0, rmulti *y1, rmulti *x0, rmulti *x1);
void irreal_d(rmulti *y0, rmulti *y1, double x0, double x1);
// y=imag(x)
void irimag_r(rmulti *y0, rmulti *y1, rmulti *x0, rmulti *x1);
void irimag_d(rmulti *y0, rmulti *y1, double x0, double x1);
// y=abs(x)
void irabs_r(rmulti *y0, rmulti *y1, rmulti *x0, rmulti *x1);
// y=abs(x)^2
void irabs2_r(rmulti *y0, rmulti *y1, rmulti *x0, rmulti *x1);
// y=abs(x)+i*abs(0)
void irabsc_r(rmulti *y0, rmulti *y1, rmulti *x0, rmulti *x1);
// y=max(abs(x),abs(0))
void irmax_absc_r(rmulti *y0, rmulti *y1, rmulti *x0, rmulti *x1);
// y=-x
void irneg_r(rmulti *y0, rmulti *y1, rmulti *x0, rmulti *x1);
void irneg_d(rmulti *y0, rmulti *y1, double x0, double x1);
// y=1/x
void irinv_r(rmulti *y0, rmulti *y1, rmulti *x0, rmulti *x1);
void irinv_d(rmulti *y0, rmulti *y1, double x0, double x1);
// y=x/abs(x)
void irnormalize_r(rmulti *y0, rmulti *y1, rmulti *x0, rmulti *x1);
void irnormalize_d(rmulti *y0, rmulti *y1, double x0, double x1);
// theta=arg(z)
void irarg_r(rmulti *theta0, rmulti *theta1, rmulti *z0, rmulti *z1);
void irarg_d(rmulti *theta0, rmulti *theta1, double z0, double z1);
// y=floor(x)
void irfloor_r(rmulti *y0, rmulti *y1, rmulti *x0, rmulti *x1);
void irfloor_d(rmulti *y0, rmulti *y1, double x0, double x1);
// y=ceil(x)
void irceil_r(rmulti *y0, rmulti *y1, rmulti *x0, rmulti *x1);
void irceil_d(rmulti *y0, rmulti *y1, double x0, double x1);
// y=trunc(x)
void irtrunc_r(rmulti *y0, rmulti *y1, rmulti *x0, rmulti *x1);
void irtrunc_d(rmulti *y0, rmulti *y1, double x0, double x1);
// y=sqrt(x)
void irsqrt_r(rmulti *y0, rmulti *y1, rmulti *x0, rmulti *x1);
void irsqrt_d(rmulti *y0, rmulti *y1, double x0, double x1);
// z=x^n
void irpow_r(rmulti *z0, rmulti *z1, rmulti *x0, rmulti *x1, int n);
void irpow_d(rmulti *z0, rmulti *z1, double x0, double x1, int n);
// [y0,y1]=[-abs(x),abs(x)]
void irpm(rmulti *y0, rmulti *y1, rmulti *x);
// [m-r,m+r]=[x0,x1]
void irmid(rmulti *mid, rmulti *x0, rmulti *x1);
void irrad(rmulti *rad, rmulti *x0, rmulti *x1);
void irmr(rmulti *mid, rmulti *rad, rmulti *x0, rmulti *x1);

/*
 * z=f(x,y)
 */
// z=x+y
void iradd_rr(rmulti *z0, rmulti *z1, rmulti *x0, rmulti *x1, rmulti *y0, rmulti *y1);
void iradd_rd(rmulti *z0, rmulti *z1, rmulti *x0, rmulti *x1, double y0, double y1);
void iradd_dr(rmulti *z0, rmulti *z1, double x0, double x1, rmulti *y0, rmulti *y1);
void iradd_dd(rmulti *z0, rmulti *z1, double x0, double x1, double y0, double y1);
// z=x-y
void irsub_rr_ws(rmulti *z0, rmulti *z1, rmulti *x0, rmulti *x1, rmulti *y0, rmulti *y1, int n_rws, rmulti **rws); // n_rws>=2
void irsub_rr(rmulti *z0, rmulti *z1, rmulti *x0, rmulti *x1, rmulti *y0, rmulti *y1);
void irsub_rd(rmulti *z0, rmulti *z1, rmulti *x0, rmulti *x1, double y0, double y1);
void irsub_dr(rmulti *z0, rmulti *z1, double x0, double x1, rmulti *y0, rmulti *y1);
void irsub_dd(rmulti *z0, rmulti *z1, double x0, double x1, double y0, double y1);
// z=x*y
void irmul_rr_ws(rmulti *z0, rmulti *z1, rmulti *x0, rmulti *x1, rmulti *y0, rmulti *y1, int n_rws, rmulti **rws); // n_rws>=4
void irmul_rr(rmulti *z0, rmulti *z1, rmulti *x0, rmulti *x1, rmulti *y0, rmulti *y1);
void irmul_rd(rmulti *z0, rmulti *z1, rmulti *x0, rmulti *x1, double y0, double y1);
void irmul_dr(rmulti *z0, rmulti *z1, double x0, double x1, rmulti *y0, rmulti *y1);
void irmul_dd(rmulti *z0, rmulti *z1, double x0, double x1, double y0, double y1);
// z=x/y
void irdiv_rr(rmulti *z0, rmulti *z1, rmulti *x0, rmulti *x1, rmulti *y0, rmulti *y1);
void irdiv_rd(rmulti *z0, rmulti *z1, rmulti *x0, rmulti *x1, double y0, double y1);
void irdiv_dr(rmulti *z0, rmulti *z1, double x0, double x1, rmulti *y0, rmulti *y1);
void irdiv_dd(rmulti *z0, rmulti *z1, double x0, double x1, double y0, double y1);
// z=abs(x-y)
void irabs_sub_rr(rmulti *z0, rmulti *z1, rmulti *x0, rmulti *x1, rmulti *y0, rmulti *y1);
void irabs_sub_rd(rmulti *z0, rmulti *z1, rmulti *x0, rmulti *x1, double y0, double y1);
void irabs_sub_dr(rmulti *z0, rmulti *z1, double x0, double x1, rmulti *y0, rmulti *y1);
void irabs_sub_dd(rmulti *z0, rmulti *z1, double x0, double x1, double y0, double y1);
// z=abs(x)+abs(y)
void iradd_abs_abs_rr(rmulti *z0, rmulti *z1, rmulti *x0, rmulti *x1, rmulti *y0, rmulti *y1);
void iradd_abs_abs_rd(rmulti *z0, rmulti *z1, rmulti *x0, rmulti *x1, double y0, double y1);
void iradd_abs_abs_dr(rmulti *z0, rmulti *z1, double x0, double x1, rmulti *y0, rmulti *y1);
void iradd_abs_abs_dd(rmulti *z0, rmulti *z1, double x0, double x1, double y0, double y1);
// [z]=[x]+[-abs(y),abs(y)]
void iradd_pm_rr(rmulti *z0, rmulti *z1, rmulti *x0, rmulti *x1, rmulti *y);

/*
 * y=y+f(x)
 */
// y=y+abs(x)
void iradd_abs_r(rmulti *y0, rmulti *y1, rmulti *x0, rmulti *x1);
void iradd_abs_d(rmulti *y0, rmulti *y1, double x0, double x1);
// y=y+abs(x)^2
void iradd_abs2_r(rmulti *y0, rmulti *y1, rmulti *x0, rmulti *x1);
void iradd_abs2_d(rmulti *y0, rmulti *y1, double x0, double x1);

/*
 * z=z+f(x,y)
 */
// z=z+x*y
void iradd_mul_rr_ws(rmulti *z0, rmulti *z1, rmulti *x0, rmulti *x1, rmulti *y0, rmulti *y1, int n_rws, rmulti **rws);  // n_rws>=6
void iradd_mul_rr(rmulti *z0, rmulti *z1, rmulti *x0, rmulti *x1, rmulti *y0, rmulti *y1);
void iradd_mul_rd(rmulti *z0, rmulti *z1, rmulti *x0, rmulti *x1, double y0, double y1);
void iradd_mul_dr(rmulti *z0, rmulti *z1, double x0, double x1, rmulti *y0, rmulti *y1);
void iradd_mul_dd(rmulti *z0, rmulti *z1, double x0, double x1, double y0, double y1);
// z=z-x*y
void irsub_mul_rr_ws(rmulti *z0, rmulti *z1, rmulti *x0, rmulti *x1, rmulti *y0, rmulti *y1, int n_rws, rmulti **rws); // n_rws>=6
void irsub_mul_rr(rmulti *z0, rmulti *z1, rmulti *x0, rmulti *x1, rmulti *y0, rmulti *y1);
void irsub_mul_rd(rmulti *z0, rmulti *z1, rmulti *x0, rmulti *x1, double y0, double y1);
void irsub_mul_dr(rmulti *z0, rmulti *z1, double x0, double x1, rmulti *y0, rmulti *y1);
void irsub_mul_dd(rmulti *z0, rmulti *z1, double x0, double x1, double y0, double y1);

/*
 * x<=>y
 */
int ir_in_ir(rmulti *y0, rmulti *y1, rmulti *x0, rmulti *x1); // [y0,y1] in [x0,x1]
int ir_in_irpm(rmulti *y0, rmulti *y1, rmulti *x);            // [y0,y1] in [-abs(x),abs(x)]

#endif
