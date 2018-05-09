#ifndef ISYS_IRMULTI_H
#define ISYS_IRMULTI_H

#include<is_rmulti.h>
#include<is_rvec.h>
#include<is_func.h>

/**
 @file  is_irmulti.h
 @brief 多倍長精度実数型rmultiの機械区間演算に関する関数の宣言.
*/


/*
 * setting
 */
// [y0,y1]=[x0,x1]
void irset(rmulti *y0, rmulti *y1, rmulti *x0, rmulti *x1);
void irset_d(rmulti *y0, rmulti *y1, double x0, double x1);
void irset_dr(rmulti *y0, rmulti *y1, double x0, rmulti *x1);
void irset_rd(rmulti *y0, rmulti *y1, rmulti *x0, double x1);
// [y0,y1]='[x0,x1]'
void irset_s(rmulti *x0, rmulti *x1, char *str);
// [z0,z1]=x.num/x.den
void irset_bigint(rmulti *z0, rmulti *z1, bigint *x);

/*
 * operatior of one argument
 */
// [y0,y1]=[x0,x1]
void ircopy(rmulti *y0, rmulti *y1, rmulti *x0, rmulti *x1);
// [y0,y1]=-[x0,x1]
void irneg(rmulti *y0, rmulti *y1, rmulti *x0, rmulti *x1);
// [y0,y1]=[-abs(x),abs(x)]
void irpm(rmulti *y0, rmulti *y1, rmulti *x);
void irmid(rmulti *mid, rmulti *x0, rmulti *x1);                                      // [m-r,m+r]=[x0,x1]
void irrad(rmulti *rad, rmulti *x0, rmulti *x1);                                      // [m-r,m+r]=[x0,x1]
void irmr(rmulti *mid, rmulti *rad, rmulti *x0, rmulti *x1);                          // [m-r,m+r]=[x0,x1]
void irinv(rmulti *z0, rmulti *z1, rmulti *x0, rmulti *x1);                           // [y0,y1]=1/[x0,x1]
void irabs(rmulti *y0, rmulti *y1, rmulti *x0, rmulti *x1);                           // [y0,y1]=abs([x0,x1])
void irsqrt(rmulti *y0, rmulti *y1, rmulti *x0, rmulti *x1);                          // [y0,y1]=sqrt([x0,x1])
void irexp(rmulti *y0, rmulti *y1, rmulti *x0, rmulti *x1);                           // [y0,y1]=exp([x0,x1])
void irlog(rmulti *y0, rmulti *y1, rmulti *x0, rmulti *x1);                           // [y0,y1]=log([x0,x1])
void irlog10(rmulti *y0, rmulti *y1, rmulti *x0, rmulti *x1);                         // [y0,y1]=log10([x0,x1])   
void irsin(rmulti *y0, rmulti *y1, rmulti *x0, rmulti *x1);                           // [y0,y1]=sin([x0,x1])
void ircos(rmulti *y0, rmulti *y1, rmulti *x0, rmulti *x1);                           // [y0,y1]=cos([x0,x1])
void irtan(rmulti *y0, rmulti *y1, rmulti *x0, rmulti *x1);                           // [y0,y1]=tan([x0,x1])
void irasin(rmulti *y0, rmulti *y1, rmulti *x0, rmulti *x1);                          // [y0,y1]=asin([x0,x1])
void iracos(rmulti *y0, rmulti *y1, rmulti *x0, rmulti *x1);                          // [y0,y1]=acos([x0,x1])
void iratan(rmulti *y0, rmulti *y1, rmulti *x0, rmulti *x1);                          // [y0,y1]=atan([x0,x1])
void irsinh(rmulti *y0, rmulti *y1, rmulti *x0, rmulti *x1);                          // [y0,y1]=sinh([x0,x1])
void ircosh(rmulti *y0, rmulti *y1, rmulti *x0, rmulti *x1);                          // [y0,y1]=cosh([x0,x1])
void irtanh(rmulti *y0, rmulti *y1, rmulti *x0, rmulti *x1);                          // [y0,y1]=tanh([x0,x1])
void irasinh(rmulti *y0, rmulti *y1, rmulti *x0, rmulti *x1);                         // [y0,y1]=asinh([x0,x1])
void iracosh(rmulti *y0, rmulti *y1, rmulti *x0, rmulti *x1);                         // [y0,y1]=acosh([x0,x1])
void iratanh(rmulti *y0, rmulti *y1, rmulti *x0, rmulti *x1);                         // [y0,y1]=atanh([x0,x1])

/*
 * operatior of two arguments
 */
// [z0,z1]=[x0,x1]+[y0,y1]
void iradd(rmulti *z0, rmulti *z1, rmulti *x0, rmulti *x1, rmulti *y0, rmulti *y1);
void idadd_r(rmulti *z0, rmulti *z1, double x0, double x1, rmulti *y0, rmulti *y1);
void iradd_d(rmulti *z0, rmulti *z1, rmulti *x0, rmulti *x1, double y0, double y1);
void iradd_pm(rmulti *z0, rmulti *z1, rmulti *x0, rmulti *x1, rmulti *y);               // [z0,z1]=[x0,x1]+[-y,y]
void irsub(rmulti *z0, rmulti *z1, rmulti *x0, rmulti *x1, rmulti *y0, rmulti *y1);     // [z0,z1]=[x0,x1]-[y0,y1]
void irsub_ws(rmulti *z0, rmulti *z1, rmulti *x0, rmulti *x1, rmulti *y0, rmulti *y1, int n_rws, rmulti **rws); // n_rws>=2
void irsub_d2(rmulti *z0, rmulti *z1, rmulti *x0, rmulti *x1, double y);                // [z0,z1]=[x0,x1]-[y,y]
void irmul(rmulti *z0, rmulti *z1, rmulti *x0, rmulti *x1, rmulti *y0, rmulti *y1);     // [z0,z1]=[x0,x1]*[y0,y1]
void irmul_ws(rmulti *z0, rmulti *z1, rmulti *x0, rmulti *x1, rmulti *y0, rmulti *y1, int n_rws, rmulti **rws); // n_rws>=4
void irmul_d(rmulti *z0, rmulti *z1, rmulti *x0, rmulti *x1, double y);                 // [z0,z1]=[x0,x1]*y
void iradd_mul(rmulti *z0, rmulti *z1, rmulti *x0, rmulti *x1, rmulti *y0, rmulti *y1); // [z0,z1]=[z0,z1]+[x0,x1]*[y0,y1]
void iradd_mul_ws(rmulti *z0, rmulti *z1, rmulti *x0, rmulti *x1, rmulti *y0, rmulti *y1, int n_rws, rmulti **rws);  // n_rws>=6
void irsub_mul(rmulti *z0, rmulti *z1, rmulti *x0, rmulti *x1, rmulti *y0, rmulti *y1); // [z0,z1]=[z0,z1]-[x0,x1]*[y0,y1]
void irsub_mul_ws(rmulti *z0, rmulti *z1, rmulti *x0, rmulti *x1, rmulti *y0, rmulti *y1, int n_rws, rmulti **rws); // n_rws>=6
void irdiv(rmulti *z0, rmulti *z1, rmulti *x0, rmulti *x1, rmulti *y0, rmulti *y1);     // [z0,z1]=[x0,x1]/[y0,y1]
void irabs_sub(rmulti *z0, rmulti *z1, rmulti *x0, rmulti *x1, rmulti *y0, rmulti *y1); // [z0,z1]=abs([x0,x1]-[y0,y1])
void irdiv_abs(rmulti *z0, rmulti *z1, rmulti *x0, rmulti *x1, rmulti *y0, rmulti *y1); // [z0,z1]=[x0,x1]/abs([y0,y1])
void irdiv_abs_c(rmulti *z0, rmulti *z1, rmulti *x0, rmulti *x1, cmulti *y0, cmulti *y1); // [z0,z1]=[x0,x1]/abs([y0,y1])
void irpow_si(rmulti *y0, rmulti *y1, rmulti *x0, rmulti *x1, long n);                // [y0,y1]=[x0,x1]^n
void iratan2(rmulti *z0, rmulti *z1, rmulti *x0, rmulti *x1, rmulti *y0, rmulti *y1); // [z0,z1]=atan2([x0,x1]/[y0,y1])

/*
 * comparison
 */
int irin(rmulti *y0, rmulti *y1, rmulti *x0, rmulti *x1); // [y0,y1] in [x0,x1]
int irin_pm(rmulti *y0, rmulti *y1, rmulti *x);           // [y0,y1] in [-abs(x),abs(x)]

#endif
