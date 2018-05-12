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
// [y0,y1]=[x0,x1]
void icset(cmulti *y0, cmulti *y1, cmulti *x0, cmulti *x1);
void icset_d(cmulti *y0, cmulti *y1, double x0, double x1);
void icset_z(cmulti *y0, cmulti *y1, dcomplex x0, dcomplex x1);
void icset_Z(cmulti *y0, cmulti *y1, double x0_r, double x0_i, double x1_r, double x1_i);
void icset_r(cmulti *y0, cmulti *y1, rmulti *x0, rmulti *x1);
void icset_zz(cmulti *y0, cmulti *y1, dcomplex x0, dcomplex x1);
void icset_zd(cmulti *y0, cmulti *y1, dcomplex x0, double x1);
void icset_dz(cmulti *y0, cmulti *y1, double x0, dcomplex x1);
void icset_rc(cmulti *y0, cmulti *y1, rmulti *x0, cmulti *x1);
void icset_cr(cmulti *y0, cmulti *y1, cmulti *x0, rmulti *x1);
void icset_s(cmulti *x0, cmulti *x1, char *s);
// [z0,x1]=x.num/x.den
void icset_bigint(cmulti *z0, cmulti *z1, bigint *x);

// [y0,y1]=[x0,x1]
void iccopy(cmulti *y0, cmulti *y1, cmulti *x0, cmulti *x1);
void iccopy_r(cmulti *y0, cmulti *y1, rmulti *x0, rmulti *x1);
void iccopy_rc(cmulti *y0, cmulti *y1, rmulti *x0, cmulti *x1);
void iccopy_cr(cmulti *y0, cmulti *y1, cmulti *x0, rmulti *x1);
// [c0,c1]=[a0,a1]+[b0,b1]i
void iccopy_C(cmulti *c0, cmulti *c1, rmulti *a0, rmulti *a1, rmulti *b0, rmulti *b1);


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
// [z0,z1]=[x0,x1]+[y0,y1]
void icadd_c(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, cmulti *y0, cmulti *y1);
void icadd_r(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, rmulti *y0, rmulti *y1);
void iradd_c(cmulti *z0, cmulti *z1, rmulti *x0, rmulti *x1, cmulti *y0, cmulti *y1);
void icadd_z(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, dcomplex y0, dcomplex y1);
void izadd_c(cmulti *z0, cmulti *z1, dcomplex x0, dcomplex x1, cmulti *y0, cmulti *y1);
void icadd_d(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, double y0, double y1);
void idadd_c(cmulti *z0, cmulti *z1, double x0, double x1, cmulti *y0, cmulti *y1);
void iradd_z(cmulti *z0, cmulti *z1, rmulti *x0, rmulti *x1, dcomplex y0, dcomplex y1);
void izadd_r(cmulti *z0, cmulti *z1, dcomplex x0, dcomplex x1, rmulti *y0, rmulti *y1);
// [z0,z1]=[x0,x1]-[y0,y1]
void icsub_ws(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, cmulti *y0, cmulti *y1, int n_rws, rmulti **rws); // n_rws>=2
void icsub_c(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, cmulti *y0, cmulti *y1);
void icsub_r(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, rmulti *y0, rmulti *y1);
void irsub_c(cmulti *z0, cmulti *z1, rmulti *x0, rmulti *x1, cmulti *y0, cmulti *y1);
void icsub_z(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, dcomplex y0, dcomplex y1);
void izsub_c(cmulti *z0, cmulti *z1, dcomplex x0, dcomplex x1, cmulti *y0, cmulti *y1);
void icsub_d(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, double y0, double y1);
void idsub_c(cmulti *z0, cmulti *z1, double x0, double x1, cmulti *y0, cmulti *y1);
void irsub_z(cmulti *z0, cmulti *z1, rmulti *x0, rmulti *x1, dcomplex y0, dcomplex y1);
void izsub_r(cmulti *z0, cmulti *z1, dcomplex x0, dcomplex x1, rmulti *y0, rmulti *y1);
// [z0,z1]=[x0,x1]*[y0,y1]
void icmul_c_ws(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, cmulti *y0, cmulti *y1, int n_rws, rmulti **rws, int n_cws, cmulti **cws); // n_rws>=8,n_cws>=2
void icmul_c(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, cmulti *y0, cmulti *y1);
void icmul_r(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, rmulti *y0, rmulti *y1);
void irmul_c(cmulti *z0, cmulti *z1, rmulti *x0, rmulti *x1, cmulti *y0, cmulti *y1);
void icmul_z(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, dcomplex y0, dcomplex y1);
void izmul_c(cmulti *z0, cmulti *z1, dcomplex x0, dcomplex x1, cmulti *y0, cmulti *y1);
void icmul_d(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, double y0, double y1);
void idmul_c(cmulti *z0, cmulti *z1, double x0, double x1, cmulti *y0, cmulti *y1);
void irmul_z(cmulti *z0, cmulti *z1, rmulti *x0, rmulti *x1, dcomplex y0, dcomplex y1);
void izmul_r(cmulti *z0, cmulti *z1, dcomplex x0, dcomplex x1, rmulti *y0, rmulti *y1);
// [z0,x1]=[x0,x1]/[y0,y1]
void icdiv_c(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, cmulti *y0, cmulti *y1);
void icdiv_r(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, rmulti *y0, rmulti *y1);
void irdiv_c(cmulti *z0, cmulti *z1, rmulti *x0, rmulti *x1, cmulti *y0, cmulti *y1);
void icdiv_z(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, dcomplex y0, dcomplex y1);
void izdiv_c(cmulti *z0, cmulti *z1, dcomplex x0, dcomplex x1, cmulti *y0, cmulti *y1);
void icdiv_d(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, double y0, double y1);
void iddiv_c(cmulti *z0, cmulti *z1, double x0, double x1, cmulti *y0, cmulti *y1);
void irdiv_z(cmulti *z0, cmulti *z1, rmulti *x0, rmulti *x1, dcomplex y0, dcomplex y1);
void izdiv_r(cmulti *z0, cmulti *z1, dcomplex x0, dcomplex x1, rmulti *y0, rmulti *y1);
// [z0,z1]=[x0,x1]+[-y,y]
void icadd_pm(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, cmulti *y);
void icdot(cmulti *z0, cmulti *z1, cmulti *x0, cmulti *x1, cmulti *y0, cmulti *y1);     // [z0,z1]=conj([x0,x1])*[y0,y1]
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
