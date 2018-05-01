#ifndef ISYS_ICMULTI_H
#define ISYS_ICMULTI_H

#include<is_cmulti.h>
#include<is_func.h>

/**
 @file  is_icmulti.h
 @brief 多倍長精度実数型cmultiの機械区間演算に関する関数の宣言.
*/

/*
 * setting
 */
void icset_z(cmulti *y0, cmulti *y1, dcomplex x);                                       // [y0,y1]=x
void icset_zz(cmulti *y0, cmulti *y1, dcomplex x0, dcomplex x1);                        // [y0,y1]=[x0,x1]
void icset_dd(cmulti *y0, cmulti *y1, double xr, double xi);                            // [y0,y1]=x
void icset_d(cmulti *y0, cmulti *y1, double x);                                         // [y0,y1]=x
void icset_bigint(cmulti *z0, cmulti *z1, bigint *x);                                   // [z0,x1]=x.num/x.den
void iccopy(cmulti *y0, cmulti *y1, cmulti *x0, cmulti *x1);                            // [y0,y1]=[x0,x1]
void iccopy_r(cmulti *y0, cmulti *y1, rmulti *x0, rmulti *x1);                          // [y0,y1]=[x0,x1]
void iccopy_rr(cmulti *c0, cmulti *c1, rmulti *a0, rmulti *a1, rmulti *b0, rmulti *b1); // [c0,c1]=[a0,a1]+[b0,b1]i

/*
 * operatior of one argument
 */
void icconj(cmulti *y0, cmulti *y1, cmulti *x0, cmulti *x1);  // [y0,y1]=conj([x0,x1])
void icneg(cmulti *y0, cmulti *y1, cmulti *x0, cmulti *x1);   // [y0,y1]=-[x0,x1]
void icpm(cmulti *y0, cmulti *y1, cmulti *x);                 // [y0,y1]=[-x,x]
void icmid(cmulti *mid, cmulti *x0, cmulti *x1);              // [m-r,m+r]=[x0,x1]
void icrad(cmulti *rad, cmulti *x0, cmulti *x1);              // [m-r,m+r]=[x0,x1]
void icmr(cmulti *mid, cmulti *rad, cmulti *x0, cmulti *x1);  // [m-r,m+r]=[x0,x1]
void icabs2(rmulti *y0, rmulti *y1, cmulti *x0, cmulti *x1);  // [y0,y1]=abs([x0,x1])^2
void icsqrt(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1);  // [z0,z1]=sqrt([x0,x1])
void icexp(cmulti *y0, cmulti *y1, cmulti *x0, cmulti *x1);   // [y0,y1]=exp([x0,x1])
void iclog(cmulti *y0, cmulti *y1, cmulti *x0, cmulti *x1);   // [y0,y1]=log([x0,x1])
void icsin(cmulti *y0, cmulti *y1, cmulti *x0, cmulti *x1);   // [y0,y1]=sin([x0,x1])
void iccos(cmulti *y0, cmulti *y1, cmulti *x0, cmulti *x1);   // [y0,y1]=cos([x0,x1])
void ictan(cmulti *y0, cmulti *y1, cmulti *x0, cmulti *x1);   // [y0,y1]=tan([x0,x1])
void icasin(cmulti *y0, cmulti *y1, cmulti *x0, cmulti *x1);  // [y0,y1]=asin([x0,x1])
void icacos(cmulti *y0, cmulti *y1, cmulti *x0, cmulti *x1);  // [y0,y1]=acos([x0,x1])
void icatan(cmulti *y0, cmulti *y1, cmulti *x0, cmulti *x1);  // [y0,y1]=atan([x0,x1])
void icsinh(cmulti *y0, cmulti *y1, cmulti *x0, cmulti *x1);  // [y0,y1]=sinh([x0,x1])
void iccosh(cmulti *y0, cmulti *y1, cmulti *x0, cmulti *x1);  // [y0,y1]=cosh([x0,x1])
void ictanh(cmulti *y0, cmulti *y1, cmulti *x0, cmulti *x1);  // [y0,y1]=tanh([x0,x1])
void icasinh(cmulti *y0, cmulti *y1, cmulti *x0, cmulti *x1); // [y0,y1]=asinh([x0,x1])
void icacosh(cmulti *y0, cmulti *y1, cmulti *x0, cmulti *x1); // [y0,y1]=acosh([x0,x1])
void icatanh(cmulti *y0, cmulti *y1, cmulti *x0, cmulti *x1); // [y0,y1]=atanh([x0,x1])

/*
 * operatior of two arguments
 */
void icadd(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, cmulti *y0, cmulti *y1);     // [z0,z1]=[x0,x1]+[y0,y1]
void icadd_r(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, rmulti *y0, rmulti *y1);   // [z0,z1]=[x0,x1]+[y0,y1]
void icadd_pm(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, cmulti *y);               // [z0,z1]=[x0,x1]+[-y,y]
void icsub(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, cmulti *y0, cmulti *y1);     // [z0,z1]=[x0,x1]-[y0,y1]
void icsub_ws(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, cmulti *y0, cmulti *y1, int n_rws, rmulti **rws); // n_rws>=2
void icsub_r(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, rmulti *y0, rmulti *y1);   // [z0,z1]=[x0,x1]-[y0,y1]
void irsub_c(cmulti *z0, cmulti *z1, rmulti *x0, rmulti *x1, cmulti *y0, cmulti *y1);   // [z0,z1]=[x0,x1]-[y0,y1]
void icsub_d2(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, double y);                // [z0,z1]=[x0,x1]-[y,y]
void icmul(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, cmulti *y0, cmulti *y1);     // [z0,z1]=[x0,x1]*[y0,y1]
void icmul_ws(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, cmulti *y0, cmulti *y1, int n_rws, rmulti **rws, int n_cws, cmulti **cws); // n_rws>=8,n_cws>=2
void icmul_r1(cmulti *z0, cmulti *z1, rmulti *x0, rmulti *x1, cmulti *y0, cmulti *y1);  // [z0,z1]=[x0,x1]*[y0,y1]
void icmul_r2(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, rmulti *y0, rmulti *y1);  // [z0,z1]=[x0,x1]*[y0,y1]
void icmul_d(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, double y);                 // [z0,z1]=[x0,x1]*[y,y]
void icdot(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, cmulti *y0, cmulti *y1);     // [z0,z1]=conj([x0,x1])*[y0,y1]
void icdiv(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, cmulti *y0, cmulti *y1);     // [z0,x1]=[x0,x1]/[y0,y1]
void icdiv_r1(cmulti *z0, cmulti *z1, rmulti *x0, rmulti *x1, cmulti *y0, cmulti *y1);  // [z0,x1]=[x0,x1]/[y0,y1]
void icdiv_r2(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, rmulti *y0, rmulti *y1);  // [z0,x1]=[x0,x1]/[y0,y1]
void icinv(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1);                             // [z0,z1]=[1,1]/[x0,x1]
void icadd_mul(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, cmulti *y0, cmulti *y1); // [z0,z1]+=[x0,x1]*[y0,y1]
void icadd_mul_r1(cmulti *z0, cmulti *z1, rmulti *x0, rmulti *x1, cmulti *y0, cmulti *y1); // [z0,z1]+=[x0,x1]*[y0,y1]
void icadd_mul_r2(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, rmulti *y0, rmulti *y1); // [z0,z1]+=[x0,x1]*[y0,y1]
void icsub_mul(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, cmulti *y0, cmulti *y1); // [z0,z1]-=[x0,x1]*[y0,y1]
void icsub_mul_ws(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, cmulti *y0, cmulti *y1, int n_rws, rmulti **rws, int n_cws, cmulti **cws); // n_rws>=8,n_cws>=4
void icadd_dot(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, cmulti *y0, cmulti *y1); // [z0,z1]+=conj([x0,x1])*[y0,y1]
void icsub_dot(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, cmulti *y0, cmulti *y1); // [z0,z1]-=conj([x0,x1])*[y0,y1]
void icadd_abs2(rmulti *y0, rmulti *y1, cmulti *x0, cmulti *x1);                        // [y0,y1]=[y0,y1]+abs([x0,x1])^2
void icabs(rmulti *y0, rmulti *y1, cmulti *x0, cmulti *x1);                             // [y0,y1]=abs([x0,x1])
void icabs_ws(rmulti *y0, rmulti *y1, cmulti *x0, cmulti *x1, int n_cws, cmulti **cws); // n_cws>=2
void icabs_sub(rmulti *z0, rmulti *z1, cmulti *x0, cmulti *x1, cmulti *y0, cmulti *y1); // [z0,z1]=abs([x0,x1]-[y0,y1])
void icpow_si(cmulti *y0, cmulti *y1, cmulti *x0, cmulti *x1, long n);                  // [y0,y1]=[x0,x1]^n
void icget_polar(rmulti *r0, rmulti *r1, rmulti *theta0, rmulti *theta1, cmulti *z0, cmulti *z1); // [z0,z1]=[r0,r1]*exp(i*[theta0,theta1])
void icset_polar(cmulti *z0, cmulti *z1, rmulti *r0, rmulti *r1, rmulti *theta0, rmulti *theta1); // [z0,z1]=[r0,r1]*exp(i*[theta0,theta1])

/*
 * comparison
 */
int icin(cmulti *y0, cmulti *y1, cmulti *x0, cmulti *x1); // [y0,y1] in [x0,x1]
int icin_pm(cmulti *y0, cmulti *y1, cmulti *x);           // [y0,y1] in [-abs(x),abs(x)]

#endif
