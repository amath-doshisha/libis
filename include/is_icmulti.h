#ifndef ISYS_ICMULTI_H
#define ISYS_ICMULTI_H

#include<is_cmulti.h>
#include<is_func.h>

/**
 @file  is_icmulti.h
 @brief 多倍長精度実数型cmultiの機械区間演算に関する関数の宣言.
*/

int icset_z(cmulti *y0, cmulti *y1, dcomplex x);                                       // [y0,y1]=x
int icset_dd(cmulti *y0, cmulti *y1, double xr, double xi);                            // [y0,y1]=x
int icset_d(cmulti *y0, cmulti *y1, double x);                                         // [y0,y1]=x
int icset_bigint(cmulti *z0, cmulti *z1, bigint *x);                                   // [z0,x1]=x.num/x.den
int iccopy(cmulti *y0, cmulti *y1, cmulti *x0, cmulti *x1);                            // [y0,y1]=[x0,x1]
int icconj(cmulti *y0, cmulti *y1, cmulti *x0, cmulti *x1);                            // [y0,y1]=conj([x0,x1])
int icneg(cmulti *y0, cmulti *y1, cmulti *x0, cmulti *x1);                             // [y0,y1]=-[x0,x1]
int icpm(cmulti *y0, cmulti *y1, cmulti *x);                                           // [y0,y1]=[-x,x]
int icadd(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, cmulti *y0, cmulti *y1);     // [z0,z1]=[x0,x1]+[y0,y1]
int icadd_pm(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, cmulti *y);               // [z0,z1]=[x0,x1]+[-y,y]
int icmr(cmulti *mid, cmulti *rad, cmulti *x0, cmulti *x1);                            // [m-r,m+r]=[x0,x1]
int icsub(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, cmulti *y0, cmulti *y1);     // [z0,z1]=[x0,x1]-[y0,y1]
int icsub_d2(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, double y);                // [z0,z1]=[x0,x1]-[y,y]
int icmul(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, cmulti *y0, cmulti *y1);     // [z0,z1]=[x0,x1]*[y0,y1]
int icmul_d(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, double y);                 // [z0,z1]=[x0,x1]*[y,y]
int icdot(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, cmulti *y0, cmulti *y1);     // [z0,z1]=conj([x0,x1])*[y0,y1]
int icdiv(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, cmulti *y0, cmulti *y1);     // [z0,x1]=[x0,x1]/[y0,y1]
int icdiv_r2(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, rmulti *y0, rmulti *y1);  // [z0,x1]=[x0,x1]/[y0,y1]
int icinv(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1);                             // [z0,z1]=[1,1]/[x0,x1]
int icadd_mul(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, cmulti *y0, cmulti *y1); // [z0,z1]+=[x0,x1]*[y0,y1]
int icsub_mul(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, cmulti *y0, cmulti *y1); // [z0,z1]-=[x0,x1]*[y0,y1]
int icadd_dot(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, cmulti *y0, cmulti *y1); // [z0,z1]+=conj([x0,x1])*[y0,y1]
int icsub_dot(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, cmulti *y0, cmulti *y1); // [z0,z1]-=conj([x0,x1])*[y0,y1]
int icabs2(rmulti *y0, rmulti *y1, cmulti *x0, cmulti *x1);                            // [y0,y1]=abs([x0,x1])^2
int icadd_abs2(rmulti *y0, rmulti *y1, cmulti *x0, cmulti *x1);                        // [y0,y1]=[y0,y1]+abs([x0,x1])^2
int icabs(rmulti *y0, rmulti *y1, cmulti *x0, cmulti *x1);                             // [y0,y1]=abs([x0,x1])
int icabs_sub(rmulti *z0, rmulti *z1, cmulti *x0, cmulti *x1, cmulti *y0, cmulti *y1); // [z0,z1]=abs([x0,x1]-[y0,y1])
int icpow_si(cmulti *y0, cmulti *y1, cmulti *x0, cmulti *x1, long n);                  // [y0,y1]=[x0,x1]^n
int icget_polar(rmulti *r0, rmulti *r1, rmulti *theta0, rmulti *theta1, cmulti *z0, cmulti *z1); // [z0,z1]=[r0,r1]*exp(i*[theta0,theta1])
int icset_polar(cmulti *z0, cmulti *z1, rmulti *r0, rmulti *r1, rmulti *theta0, rmulti *theta1); // [z0,z1]=[r0,r1]*exp(i*[theta0,theta1])
int icsqrt(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1);  // [z0,z1]=sqrt([x0,x1])
int icexp(cmulti *y0, cmulti *y1, cmulti *x0, cmulti *x1);   // [y0,y1]=exp([x0,x1])
int iclog(cmulti *y0, cmulti *y1, cmulti *x0, cmulti *x1);   // [y0,y1]=log([x0,x1])
int icsin(cmulti *y0, cmulti *y1, cmulti *x0, cmulti *x1);   // [y0,y1]=sin([x0,x1])
int iccos(cmulti *y0, cmulti *y1, cmulti *x0, cmulti *x1);   // [y0,y1]=cos([x0,x1])
int ictan(cmulti *y0, cmulti *y1, cmulti *x0, cmulti *x1);   // [y0,y1]=tan([x0,x1])
int icasin(cmulti *y0, cmulti *y1, cmulti *x0, cmulti *x1);  // [y0,y1]=asin([x0,x1])
int icacos(cmulti *y0, cmulti *y1, cmulti *x0, cmulti *x1);  // [y0,y1]=acos([x0,x1])
int icatan(cmulti *y0, cmulti *y1, cmulti *x0, cmulti *x1);  // [y0,y1]=atan([x0,x1])
int icsinh(cmulti *y0, cmulti *y1, cmulti *x0, cmulti *x1);  // [y0,y1]=sinh([x0,x1])
int iccosh(cmulti *y0, cmulti *y1, cmulti *x0, cmulti *x1);  // [y0,y1]=cosh([x0,x1])
int ictanh(cmulti *y0, cmulti *y1, cmulti *x0, cmulti *x1);  // [y0,y1]=tanh([x0,x1])
int icasinh(cmulti *y0, cmulti *y1, cmulti *x0, cmulti *x1); // [y0,y1]=asinh([x0,x1])
int icacosh(cmulti *y0, cmulti *y1, cmulti *x0, cmulti *x1); // [y0,y1]=acosh([x0,x1])
int icatanh(cmulti *y0, cmulti *y1, cmulti *x0, cmulti *x1); // [y0,y1]=atanh([x0,x1])


int icin(cmulti *y0, cmulti *y1, cmulti *x0, cmulti *x1); // [y0,y1] in [x0,x1]
int icin_pm(cmulti *y0, cmulti *y1, cmulti *x);           // [y0,y1] in [-abs(x),abs(x)]

#endif
