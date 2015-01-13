#ifndef ISYS_IRMULTI_H
#define ISYS_IRMULTI_H

#include<is_rmulti.h>
#include<is_rvec.h>
#include<is_func.h>

/**
 @file  is_irmulti.h
 @brief 多倍長精度実数型rmultiの機械区間演算に関する関数の宣言.
*/

int irset_d(rmulti *y0, rmulti *y1, double x);              // [y0,y1]=[x,x]
int irset_bigint(rmulti *z0, rmulti *z1, bigint *x);        // [z0,z1]=x.num/x.den
int ircopy(rmulti *y0, rmulti *y1, rmulti *x0, rmulti *x1); // [y0,y1]=[x0,x1]

int irneg(rmulti *y0, rmulti *y1, rmulti *x0, rmulti *x1);                             // [y0,y1]=-[x0,x1]
int irpm(rmulti *y0, rmulti *y1, rmulti *x);                                           // [y0,y1]=[-abs(x),abs(x)]
int iradd(rmulti *z0, rmulti *z1, rmulti *x0, rmulti *x1, rmulti *y0, rmulti *y1);     // [z0,z1]=[x0,x1]+[y0,y1]
int iradd_pm(rmulti *z0, rmulti *z1, rmulti *x0, rmulti *x1, rmulti *y);               // [z0,z1]=[x0,x1]+[-y,y]
int irmr(rmulti *mid, rmulti *rad, rmulti *x0, rmulti *x1);                            // [m-r,m+r]=[x0,x1]
int irsub(rmulti *z0, rmulti *z1, rmulti *x0, rmulti *x1, rmulti *y0, rmulti *y1);     // [z0,z1]=[x0,x1]-[y0,y1]
int irsub_d2(rmulti *z0, rmulti *z1, rmulti *x0, rmulti *x1, double y);                // [z0,z1]=[x0,x1]-[y,y]
int irmul(rmulti *z0, rmulti *z1, rmulti *x0, rmulti *x1, rmulti *y0, rmulti *y1);     // [z0,z1]=[x0,x1]*[y0,y1]
int irmul_d(rmulti *z0, rmulti *z1, rmulti *x0, rmulti *x1, double y);                 // [z0,z1]=[x0,x1]*y
int irdiv(rmulti *z0, rmulti *z1, rmulti *x0, rmulti *x1, rmulti *y0, rmulti *y1);     // [z0,z1]=[x0,x1]/[y0,y1]
int iradd_mul(rmulti *z0, rmulti *z1, rmulti *x0, rmulti *x1, rmulti *y0, rmulti *y1); // [z0,z1]=[z0,z1]+[x0,x1]*[y0,y1]
int irsub_mul(rmulti *z0, rmulti *z1, rmulti *x0, rmulti *x1, rmulti *y0, rmulti *y1); // [z0,z1]=[z0,z1]-[x0,x1]*[y0,y1]

int irabs(rmulti *y0, rmulti *y1, rmulti *x0, rmulti *x1);                             // [y0,y1]=abs([x0,x1])
int irabs_sub(rmulti *z0, rmulti *z1, rmulti *x0, rmulti *x1, rmulti *y0, rmulti *y1); // [z0,z1]=abs([x0,x1]-[y0,y1])
int irdiv_abs(rmulti *z0, rmulti *z1, rmulti *x0, rmulti *x1, rmulti *y0, rmulti *y1); // [z0,z1]=[x0,x1]/abs([y0,y1])
int irdiv_abs_c(rmulti *z0, rmulti *z1, rmulti *x0, rmulti *x1, cmulti *y0, cmulti *y1); // [z0,z1]=[x0,x1]/abs([y0,y1])
int irpow_si(rmulti *y0, rmulti *y1, rmulti *x0, rmulti *x1, long n);                // [y0,y1]=[x0,x1]^n
int irsqrt(rmulti *y0, rmulti *y1, rmulti *x0, rmulti *x1);                          // [y0,y1]=sqrt([x0,x1])
int irexp(rmulti *y0, rmulti *y1, rmulti *x0, rmulti *x1);                           // [y0,y1]=exp([x0,x1])
int irlog(rmulti *y0, rmulti *y1, rmulti *x0, rmulti *x1);                           // [y0,y1]=log([x0,x1])
int irsin(rmulti *y0, rmulti *y1, rmulti *x0, rmulti *x1);                           // [y0,y1]=sin([x0,x1])
int ircos(rmulti *y0, rmulti *y1, rmulti *x0, rmulti *x1);                           // [y0,y1]=cos([x0,x1])
int irtan(rmulti *y0, rmulti *y1, rmulti *x0, rmulti *x1);                           // [y0,y1]=tan([x0,x1])
int irasin(rmulti *y0, rmulti *y1, rmulti *x0, rmulti *x1);                          // [y0,y1]=asin([x0,x1])
int iracos(rmulti *y0, rmulti *y1, rmulti *x0, rmulti *x1);                          // [y0,y1]=acos([x0,x1])
int iratan(rmulti *y0, rmulti *y1, rmulti *x0, rmulti *x1);                          // [y0,y1]=atan([x0,x1])
int iratan2(rmulti *z0, rmulti *z1, rmulti *x0, rmulti *x1, rmulti *y0, rmulti *y1); // [z0,z1]=atan2([x0,x1]/[y0,y1])
int irsinh(rmulti *y0, rmulti *y1, rmulti *x0, rmulti *x1);                          // [y0,y1]=sinh([x0,x1])
int ircosh(rmulti *y0, rmulti *y1, rmulti *x0, rmulti *x1);                          // [y0,y1]=cosh([x0,x1])
int irtanh(rmulti *y0, rmulti *y1, rmulti *x0, rmulti *x1);                          // [y0,y1]=tanh([x0,x1])
int irasinh(rmulti *y0, rmulti *y1, rmulti *x0, rmulti *x1);                         // [y0,y1]=asinh([x0,x1])
int iracosh(rmulti *y0, rmulti *y1, rmulti *x0, rmulti *x1);                         // [y0,y1]=acosh([x0,x1])
int iratanh(rmulti *y0, rmulti *y1, rmulti *x0, rmulti *x1);                         // [y0,y1]=atanh([x0,x1])

int irin(rmulti *y0, rmulti *y1, rmulti *x0, rmulti *x1); // [y0,y1] in [x0,x1]
int irin_pm(rmulti *y0, rmulti *y1, rmulti *x);           // [y0,y1] in [-abs(x),abs(x)]

#endif
