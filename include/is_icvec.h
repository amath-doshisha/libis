#ifndef ISYS_ICVEC_H
#define ISYS_ICVEC_H

#include<is_dcomplex.h>
#include<is_cmulti.h>

/**
 @file  is_icvec.h
 @brief 多倍長精度複素数型cmultiの機械区間演算のベクトルに関する関数の宣言.
*/

/*
 * I/O
 */
void icvec_print(int n, cmulti **x0, cmulti **x1, char *name, char format, int digits);

/*
 * casting
 */
void icvec_set_cvec(int n, cmulti **y0, cmulti **y1, cmulti **x0, cmulti **x1);                                          // [cmulti,cmulti] <- [cmulti,cmulti]
void icvec_set_rvec_rvec(int n, cmulti **y0, cmulti **y1, rmulti **x0_r, rmulti **x0_i, rmulti **x1_r, rmulti **x1_i);   // [cmulti,cmulti] <- [(rmulti,rmulti),(rmulti,rmulti)]
void icvec_set_cvec_rvec(int n, cmulti **y0, cmulti **y1, cmulti **x0, rmulti **x1);                                     // [cmulti,cmulti] <- [cmulti,rmulti]
void icvec_set_rvec_cvec(int n, cmulti **y0, cmulti **y1, rmulti **x0, cmulti **x1);                                     // [cmulti,cmulti] <- [rmulti,cmulti]
void icvec_set_rvec(int n, cmulti **y0, cmulti **y1, rmulti **x0, rmulti **x1);                                          // [cmulti,cmulti] <- [rmulti,rmulti]
void irvec_set_cvec(int n, rmulti **y0, rmulti **y1, cmulti **x0, cmulti **x1);                                          // [rmulti,rmulti] <- [cmulti,cmulti]
void icvec_set_zvec(int n, cmulti **y0, cmulti **y1, dcomplex *x0, dcomplex *x1);                                        // [cmulti,cmulti] <- [dcomplex,dcomplex]
void icvec_set_dvec(int n, cmulti **y0, cmulti **y1, double *x0, double *x1);                                            // [cmulti,cmulti] <- [double,double]
void icvec_set_ivec(int n, cmulti **y0, cmulti **y1, int *x);                                                            // [cmulti,cmulti] <- [int,int]
void icvec_set_dvec_zvec(int n, cmulti **y0, cmulti **y1, double *x0, dcomplex *x1);                                     // [cmulti,cmulti] <- [double,dcomplex]
void icvec_set_zvec_dvec(int n, cmulti **y0, cmulti **y1, dcomplex *x0, double *x1);                                     // [cmulti,cmulti] <- [dcomplex,double]
void icvec_get_cvec(int n, cmulti **y, cmulti **x0, cmulti **x1);  // [cmulti,cmulti] -> cmulti
void irvec_get_cvec(int n, cmulti **y, rmulti **x0, rmulti **x1);  // [rmulti,rmulti] -> cmulti
void icvec_get_zvec(int n, dcomplex *y, cmulti **x0, cmulti **x1); // [cmulti,cmulti] -> dcomplex
void icvec_get_rvec(int n, rmulti **y, cmulti **x0, cmulti **x1);  // [cmulti,cmulti] -> rmulti
void icvec_get_dvec(int n, double *y, cmulti **x0, cmulti **x1);   // [cmulti,cmulti] -> double
void icvec_get_ivec(int n, int *y, cmulti **x0, cmulti **x1);      // [cmulti,cmulti] -> int
void cvec_set_icvec(int n, cmulti **y, cmulti **x0, cmulti **x1);  // [cmulti,cmulti] -> cmulti
void cvec_set_irvec(int n, cmulti **y, rmulti **x0, rmulti **x1);  // [rmulti,rmulti] -> cmulti
void zvec_set_icvec(int n, dcomplex *y, cmulti **x0, cmulti **x1); // [cmulti,cmulti] -> dcomplex
void rvec_set_icvec(int n, rmulti **y, cmulti **x0, cmulti **x1);  // [cmulti,cmulti] -> rmulti
void dvec_set_icvec(int n, double *y, cmulti **x0, cmulti **x1);   // [cmulti,cmulti] -> double
void ivec_set_icvec(int n, int *y, cmulti **x0, cmulti **x1);      // [cmulti,cmulti] -> int
void icvec_set_svec(int n, cmulti **y0, cmulti **y1, char **x);                           // [cmulti,cmulti] <- char
void icvec_get_svec(int n, char **y, cmulti **x0, cmulti **x1, char format, int digits);  // [cmulti,cmulti] -> char

/*
 * operations
 */
// y=x
void icvec_copy_cvec(int n, cmulti **y0, cmulti **y1, cmulti **x0, cmulti **x1);      // y=x
void icvec_copy_rvec(int n, cmulti **y0, cmulti **y1, rmulti **x0, rmulti **x1);      // y=x
void icvec_copy_cvec_rvec(int n, cmulti **y0, cmulti **y1, cmulti **x0, rmulti **x1); // y=x
void icvec_copy_rvec_cvec(int n, cmulti **y0, cmulti **y1, rmulti **x0, cmulti **x1); // y=x
void irvec_real_cvec(int n, rmulti **y0, rmulti **y1, cmulti **x0, cmulti **x1);      // y=real(x)
void irvec_imag_cvec(int n, rmulti **y0, rmulti **y1, cmulti **x0, cmulti **x1);      // y=imag(x)
void irvec_abs_cvec(int n, rmulti **y0, rmulti **y1, cmulti **x0, cmulti **x1);       // y=abs(x)
void icvec_conj_cvec(int n, cmulti **y0, cmulti **y1, cmulti **x0, cmulti **x1);      // y=conj(x)
void icvec_absc_cvec(int n, cmulti **y0, cmulti **y1, cmulti **x0, cmulti **x1);      // y=abs(real(x))+i*abs(imag(x))
void icvec_neg_cvec(int n, cmulti **y0, cmulti **y1, cmulti **x0, cmulti **x1);       // y=-x
void icvec_pow2_cvec(int n, cmulti **y0, cmulti **y1, cmulti **x0, cmulti **x1);      // y=x^2
void icsum_cvec(cmulti *y0, cmulti *y1, int n, cmulti **x0, cmulti **x1);             // y=sum(x)
void icaverage_cvec(cmulti *y0, cmulti *y1, int n, cmulti **x0, cmulti **x1);         // y=sum(x)/n
void irnorm1_cvec(rmulti *y0, rmulti *y1, int n, cmulti **x0, cmulti **x1);           // y=sum(abs(x))
void irsum_abs_cvec(rmulti *y0, rmulti *y1, int n, cmulti **x0, cmulti **x1);         // y=sum(abs(x))
void irsum_pow2_abs_cvec(rmulti *y0, rmulti *y1, int n, cmulti **x0, cmulti **x1);    // y=sum(abs(x)^2)
void irnorm1_cvec(rmulti *y0, rmulti *y1, int n, cmulti **x0, cmulti **x1);           // y=sum(abs(x))
void irnorm2_cvec(rmulti *y0, rmulti *y1, int n, cmulti **x0, cmulti **x1);           // y=sqrt(sum(abs(x)^2))
void irnorm_max_cvec(rmulti *y0, rmulti *y1, int n, cmulti **x0, cmulti **x1);        // y=max(abs(x))
void icmax_cvec(cmulti *y0, cmulti *y1, int n, cmulti **x0, cmulti **x1);             // [y0,y1]=[max(x0),max(x1)]
void icmin_cvec(cmulti *y0, cmulti *y1, int n, cmulti **x0, cmulti **x1);             // [y0,y1]=[min(x0),min(x1)]
void icmax_up_cvec(cmulti *y0, cmulti *y1, int n, cmulti **x0, cmulti **x1);          // [y0,y1]=[x0,max(x1)]
void icmin_up_cvec(cmulti *y0, cmulti *y1, int n, cmulti **x0, cmulti **x1);          // [y0,y1]=[x0,min(x1)]
void icmax_down_cvec(cmulti *y0, cmulti *y1, int n, cmulti **x0, cmulti **x1);        // [y0,y1]=[max(x0),x1]
void icmin_down_cvec(cmulti *y0, cmulti *y1, int n, cmulti **x0, cmulti **x1);        // [y0,y1]=[min(x0),x1]
void irmax_up_abs_cvec(rmulti *y0, rmulti *y1, int n, cmulti **x0, cmulti **x1);      // y=max(abs(x))
void irmin_down_abs_cvec(rmulti *y0, rmulti *y1, int n, cmulti **x0, cmulti **x1);    // y=min(abs(x))
void icvec_pm_cvec(int n, cmulti **y0, cmulti **y1, cmulti **x);                      // [y0,y1]=[-absc(x),absc(x)]
void icvec_mid_rad_cvec(int n, cmulti **mid, cmulti **rad, cmulti **x0, cmulti **x1); // mid=(x1+x0)/2, rad=x1-x0
void icvec_mr_cvec(int n, cmulti **mid, cmulti **rad, cmulti **x0, cmulti **x1);      // [m-r,m+r]=[x0,x1]
void icvec_mid_cvec(int n, cmulti **mid, cmulti **x0, cmulti **x1);                   // [m-r,m+r]=[x0,x1]
void icvec_rad_cvec(int n, cmulti **rad, cmulti **x0, cmulti **x1);                   // [m-r,m+r]=[x0,x1]

/*
 * z=f(x,y)
 */
// z=x+y
void icvec_add_cvec_cvec(int n, cmulti **z0, cmulti **z1, cmulti **x0, cmulti **x1, cmulti **y0, cmulti **y1);
void icvec_add_cvec_rvec(int n, cmulti **z0, cmulti **z1, cmulti **x0, cmulti **x1, rmulti **y0, rmulti **y1);
void icvec_add_rvec_cvec(int n, cmulti **z0, cmulti **z1, rmulti **x0, rmulti **x1, cmulti **y0, cmulti **y1);
void icvec_add_cvec_zvec(int n, cmulti **z0, cmulti **z1, cmulti **x0, cmulti **x1, dcomplex *y0, dcomplex *y1);
void icvec_add_zvec_cvec(int n, cmulti **z0, cmulti **z1, dcomplex *x0, dcomplex *x1, cmulti **y0, cmulti **y1);
void icvec_add_cvec_dvec(int n, cmulti **z0, cmulti **z1, cmulti **x0, cmulti **x1, double *y0, double *y1);
void icvec_add_dvec_cvec(int n, cmulti **z0, cmulti **z1, double *x0, double *x1, cmulti **y0, cmulti **y1);
void icvec_add_zvec_rvec(int n, cmulti **z0, cmulti **z1, dcomplex *x0, dcomplex *x1, rmulti **y0, rmulti **y1);
void icvec_add_rvec_zvec(int n, cmulti **z0, cmulti **z1, rmulti **x0, rmulti **x1, dcomplex *y0, dcomplex *y1);
void icvec_add_cvec_cscalar(int n, cmulti **z0, cmulti **z1, cmulti **x0, cmulti **x1, cmulti *y0, cmulti *y1);
void icvec_add_cvec_rscalar(int n, cmulti **z0, cmulti **z1, cmulti **x0, cmulti **x1, rmulti *y0, rmulti *y1);
void icvec_add_rvec_cscalar(int n, cmulti **z0, cmulti **z1, rmulti **x0, rmulti **x1, cmulti *y0, cmulti *y1);
void icvec_add_cvec_zscalar(int n, cmulti **z0, cmulti **z1, cmulti **x0, cmulti **x1, dcomplex y0, dcomplex y1);
void icvec_add_zvec_cscalar(int n, cmulti **z0, cmulti **z1, dcomplex *x0, dcomplex *x1, cmulti *y0, cmulti *y1);
void icvec_add_cvec_dscalar(int n, cmulti **z0, cmulti **z1, cmulti **x0, cmulti **x1, double y0, double y1);
void icvec_add_dvec_cscalar(int n, cmulti **z0, cmulti **z1, double *x0, double *x1, cmulti *y0, cmulti *y1);
void icvec_add_rvec_zscalar(int n, cmulti **z0, cmulti **z1, rmulti **x0, rmulti **x1, dcomplex y0, dcomplex y1);
void icvec_add_zvec_rscalar(int n, cmulti **z0, cmulti **z1, dcomplex *x0, dcomplex *x1, rmulti *y0, rmulti *y1);
void icvec_add_cscalar_cvec(int n, cmulti **z0, cmulti **z1, cmulti *x0, cmulti *x1, cmulti **y0, cmulti **y1);
void icvec_add_cscalar_rvec(int n, cmulti **z0, cmulti **z1, cmulti *x0, cmulti *x1, rmulti **y0, rmulti **y1);
void icvec_add_rscalar_cvec(int n, cmulti **z0, cmulti **z1, rmulti *x0, rmulti *x1, cmulti **y0, cmulti **y1);
void icvec_add_cscalar_zvec(int n, cmulti **z0, cmulti **z1, cmulti *x0, cmulti *x1, dcomplex *y0, dcomplex *y1);
void icvec_add_zscalar_cvec(int n, cmulti **z0, cmulti **z1, dcomplex x0, dcomplex x1, cmulti **y0, cmulti **y1);
void icvec_add_cscalar_dvec(int n, cmulti **z0, cmulti **z1, cmulti *x0, cmulti *x1, double *y0, double *y1);
void icvec_add_dscalar_cvec(int n, cmulti **z0, cmulti **z1, double x0, double x1, cmulti **y0, cmulti **y1);
void icvec_add_rscalar_zvec(int n, cmulti **z0, cmulti **z1, rmulti *x0, rmulti *x1, dcomplex *y0, dcomplex *y1);
void icvec_add_zscalar_rvec(int n, cmulti **z0, cmulti **z1, dcomplex x0, dcomplex x1, rmulti **y0, rmulti **y1);

// z=x-y
void icvec_sub_cvec_cvec(int n, cmulti **z0, cmulti **z1, cmulti **x0, cmulti **x1, cmulti **y0, cmulti **y1);
void icvec_sub_cvec_rvec(int n, cmulti **z0, cmulti **z1, cmulti **x0, cmulti **x1, rmulti **y0, rmulti **y1);
void icvec_sub_rvec_cvec(int n, cmulti **z0, cmulti **z1, rmulti **x0, rmulti **x1, cmulti **y0, cmulti **y1);
void icvec_sub_cvec_zvec(int n, cmulti **z0, cmulti **z1, cmulti **x0, cmulti **x1, dcomplex *y0, dcomplex *y1);
void icvec_sub_zvec_cvec(int n, cmulti **z0, cmulti **z1, dcomplex *x0, dcomplex *x1, cmulti **y0, cmulti **y1);
void icvec_sub_cvec_dvec(int n, cmulti **z0, cmulti **z1, cmulti **x0, cmulti **x1, double *y0, double *y1);
void icvec_sub_dvec_cvec(int n, cmulti **z0, cmulti **z1, double *x0, double *x1, cmulti **y0, cmulti **y1);
void icvec_sub_zvec_rvec(int n, cmulti **z0, cmulti **z1, dcomplex *x0, dcomplex *x1, rmulti **y0, rmulti **y1);
void icvec_sub_rvec_zvec(int n, cmulti **z0, cmulti **z1, rmulti **x0, rmulti **x1, dcomplex *y0, dcomplex *y1);
void icvec_sub_cvec_cscalar(int n, cmulti **z0, cmulti **z1, cmulti **x0, cmulti **x1, cmulti *y0, cmulti *y1);
void icvec_sub_cvec_rscalar(int n, cmulti **z0, cmulti **z1, cmulti **x0, cmulti **x1, rmulti *y0, rmulti *y1);
void icvec_sub_rvec_cscalar(int n, cmulti **z0, cmulti **z1, rmulti **x0, rmulti **x1, cmulti *y0, cmulti *y1);
void icvec_sub_cvec_zscalar(int n, cmulti **z0, cmulti **z1, cmulti **x0, cmulti **x1, dcomplex y0, dcomplex y1);
void icvec_sub_zvec_cscalar(int n, cmulti **z0, cmulti **z1, dcomplex *x0, dcomplex *x1, cmulti *y0, cmulti *y1);
void icvec_sub_cvec_dscalar(int n, cmulti **z0, cmulti **z1, cmulti **x0, cmulti **x1, double y0, double y1);
void icvec_sub_dvec_cscalar(int n, cmulti **z0, cmulti **z1, double *x0, double *x1, cmulti *y0, cmulti *y1);
void icvec_sub_rvec_zscalar(int n, cmulti **z0, cmulti **z1, rmulti **x0, rmulti **x1, dcomplex y0, dcomplex y1);
void icvec_sub_zvec_rscalar(int n, cmulti **z0, cmulti **z1, dcomplex *x0, dcomplex *x1, rmulti *y0, rmulti *y1);
void icvec_sub_cscalar_cvec(int n, cmulti **z0, cmulti **z1, cmulti *x0, cmulti *x1, cmulti **y0, cmulti **y1);
void icvec_sub_cscalar_rvec(int n, cmulti **z0, cmulti **z1, cmulti *x0, cmulti *x1, rmulti **y0, rmulti **y1);
void icvec_sub_rscalar_cvec(int n, cmulti **z0, cmulti **z1, rmulti *x0, rmulti *x1, cmulti **y0, cmulti **y1);
void icvec_sub_cscalar_zvec(int n, cmulti **z0, cmulti **z1, cmulti *x0, cmulti *x1, dcomplex *y0, dcomplex *y1);
void icvec_sub_zscalar_cvec(int n, cmulti **z0, cmulti **z1, dcomplex x0, dcomplex x1, cmulti **y0, cmulti **y1);
void icvec_sub_cscalar_dvec(int n, cmulti **z0, cmulti **z1, cmulti *x0, cmulti *x1, double *y0, double *y1);
void icvec_sub_dscalar_cvec(int n, cmulti **z0, cmulti **z1, double x0, double x1, cmulti **y0, cmulti **y1);
void icvec_sub_rscalar_zvec(int n, cmulti **z0, cmulti **z1, rmulti *x0, rmulti *x1, dcomplex *y0, dcomplex *y1);
void icvec_sub_zscalar_rvec(int n, cmulti **z0, cmulti **z1, dcomplex x0, dcomplex x1, rmulti **y0, rmulti **y1);
// z=x*y
void icvec_mul_cvec_cvec(int n, cmulti **z0, cmulti **z1, cmulti **x0, cmulti **x1, cmulti **y0, cmulti **y1);
void icvec_mul_cvec_rvec(int n, cmulti **z0, cmulti **z1, cmulti **x0, cmulti **x1, rmulti **y0, rmulti **y1);
void icvec_mul_rvec_cvec(int n, cmulti **z0, cmulti **z1, rmulti **x0, rmulti **x1, cmulti **y0, cmulti **y1);
void icvec_mul_cvec_zvec(int n, cmulti **z0, cmulti **z1, cmulti **x0, cmulti **x1, dcomplex *y0, dcomplex *y1);
void icvec_mul_zvec_cvec(int n, cmulti **z0, cmulti **z1, dcomplex *x0, dcomplex *x1, cmulti **y0, cmulti **y1);
void icvec_mul_cvec_dvec(int n, cmulti **z0, cmulti **z1, cmulti **x0, cmulti **x1, double *y0, double *y1);
void icvec_mul_dvec_cvec(int n, cmulti **z0, cmulti **z1, double *x0, double *x1, cmulti **y0, cmulti **y1);
void icvec_mul_zvec_rvec(int n, cmulti **z0, cmulti **z1, dcomplex *x0, dcomplex *x1, rmulti **y0, rmulti **y1);
void icvec_mul_rvec_zvec(int n, cmulti **z0, cmulti **z1, rmulti **x0, rmulti **x1, dcomplex *y0, dcomplex *y1);
void icvec_mul_cvec_cscalar(int n, cmulti **z0, cmulti **z1, cmulti **x0, cmulti **x1, cmulti *y0, cmulti *y1);
void icvec_mul_cvec_rscalar(int n, cmulti **z0, cmulti **z1, cmulti **x0, cmulti **x1, rmulti *y0, rmulti *y1);
void icvec_mul_rvec_cscalar(int n, cmulti **z0, cmulti **z1, rmulti **x0, rmulti **x1, cmulti *y0, cmulti *y1);
void icvec_mul_cvec_zscalar(int n, cmulti **z0, cmulti **z1, cmulti **x0, cmulti **x1, dcomplex y0, dcomplex y1);
void icvec_mul_zvec_cscalar(int n, cmulti **z0, cmulti **z1, dcomplex *x0, dcomplex *x1, cmulti *y0, cmulti *y1);
void icvec_mul_cvec_dscalar(int n, cmulti **z0, cmulti **z1, cmulti **x0, cmulti **x1, double y0, double y1);
void icvec_mul_dvec_cscalar(int n, cmulti **z0, cmulti **z1, double *x0, double *x1, cmulti *y0, cmulti *y1);
void icvec_mul_rvec_zscalar(int n, cmulti **z0, cmulti **z1, rmulti **x0, rmulti **x1, dcomplex y0, dcomplex y1);
void icvec_mul_zvec_rscalar(int n, cmulti **z0, cmulti **z1, dcomplex *x0, dcomplex *x1, rmulti *y0, rmulti *y1);
void icvec_mul_cscalar_cvec(int n, cmulti **z0, cmulti **z1, cmulti *x0, cmulti *x1, cmulti **y0, cmulti **y1);
void icvec_mul_cscalar_rvec(int n, cmulti **z0, cmulti **z1, cmulti *x0, cmulti *x1, rmulti **y0, rmulti **y1);
void icvec_mul_rscalar_cvec(int n, cmulti **z0, cmulti **z1, rmulti *x0, rmulti *x1, cmulti **y0, cmulti **y1);
void icvec_mul_cscalar_zvec(int n, cmulti **z0, cmulti **z1, cmulti *x0, cmulti *x1, dcomplex *y0, dcomplex *y1);
void icvec_mul_zscalar_cvec(int n, cmulti **z0, cmulti **z1, dcomplex x0, dcomplex x1, cmulti **y0, cmulti **y1);
void icvec_mul_cscalar_dvec(int n, cmulti **z0, cmulti **z1, cmulti *x0, cmulti *x1, double *y0, double *y1);
void icvec_mul_dscalar_cvec(int n, cmulti **z0, cmulti **z1, double x0, double x1, cmulti **y0, cmulti **y1);
void icvec_mul_rscalar_zvec(int n, cmulti **z0, cmulti **z1, rmulti *x0, rmulti *x1, dcomplex *y0, dcomplex *y1);
void icvec_mul_zscalar_rvec(int n, cmulti **z0, cmulti **z1, dcomplex x0, dcomplex x1, rmulti **y0, rmulti **y1);
// z=x/y
void icvec_div_cvec_cvec(int n, cmulti **z0, cmulti **z1, cmulti **x0, cmulti **x1, cmulti **y0, cmulti **y1);
void icvec_div_cvec_rvec(int n, cmulti **z0, cmulti **z1, cmulti **x0, cmulti **x1, rmulti **y0, rmulti **y1);
void icvec_div_rvec_cvec(int n, cmulti **z0, cmulti **z1, rmulti **x0, rmulti **x1, cmulti **y0, cmulti **y1);
void icvec_div_cvec_zvec(int n, cmulti **z0, cmulti **z1, cmulti **x0, cmulti **x1, dcomplex *y0, dcomplex *y1);
void icvec_div_zvec_cvec(int n, cmulti **z0, cmulti **z1, dcomplex *x0, dcomplex *x1, cmulti **y0, cmulti **y1);
void icvec_div_cvec_dvec(int n, cmulti **z0, cmulti **z1, cmulti **x0, cmulti **x1, double *y0, double *y1);
void icvec_div_dvec_cvec(int n, cmulti **z0, cmulti **z1, double *x0, double *x1, cmulti **y0, cmulti **y1);
void icvec_div_zvec_rvec(int n, cmulti **z0, cmulti **z1, dcomplex *x0, dcomplex *x1, rmulti **y0, rmulti **y1);
void icvec_div_rvec_zvec(int n, cmulti **z0, cmulti **z1, rmulti **x0, rmulti **x1, dcomplex *y0, dcomplex *y1);
void icvec_div_cvec_cscalar(int n, cmulti **z0, cmulti **z1, cmulti **x0, cmulti **x1, cmulti *y0, cmulti *y1);
void icvec_div_cvec_rscalar(int n, cmulti **z0, cmulti **z1, cmulti **x0, cmulti **x1, rmulti *y0, rmulti *y1);
void icvec_div_rvec_cscalar(int n, cmulti **z0, cmulti **z1, rmulti **x0, rmulti **x1, cmulti *y0, cmulti *y1);
void icvec_div_cvec_zscalar(int n, cmulti **z0, cmulti **z1, cmulti **x0, cmulti **x1, dcomplex y0, dcomplex y1);
void icvec_div_zvec_cscalar(int n, cmulti **z0, cmulti **z1, dcomplex *x0, dcomplex *x1, cmulti *y0, cmulti *y1);
void icvec_div_cvec_dscalar(int n, cmulti **z0, cmulti **z1, cmulti **x0, cmulti **x1, double y0, double y1);
void icvec_div_dvec_cscalar(int n, cmulti **z0, cmulti **z1, double *x0, double *x1, cmulti *y0, cmulti *y1);
void icvec_div_rvec_zscalar(int n, cmulti **z0, cmulti **z1, rmulti **x0, rmulti **x1, dcomplex y0, dcomplex y1);
void icvec_div_zvec_rscalar(int n, cmulti **z0, cmulti **z1, dcomplex *x0, dcomplex *x1, rmulti *y0, rmulti *y1);
void icvec_div_cscalar_cvec(int n, cmulti **z0, cmulti **z1, cmulti *x0, cmulti *x1, cmulti **y0, cmulti **y1);
void icvec_div_cscalar_rvec(int n, cmulti **z0, cmulti **z1, cmulti *x0, cmulti *x1, rmulti **y0, rmulti **y1);
void icvec_div_rscalar_cvec(int n, cmulti **z0, cmulti **z1, rmulti *x0, rmulti *x1, cmulti **y0, cmulti **y1);
void icvec_div_cscalar_zvec(int n, cmulti **z0, cmulti **z1, cmulti *x0, cmulti *x1, dcomplex *y0, dcomplex *y1);
void icvec_div_zscalar_cvec(int n, cmulti **z0, cmulti **z1, dcomplex x0, dcomplex x1, cmulti **y0, cmulti **y1);
void icvec_div_cscalar_dvec(int n, cmulti **z0, cmulti **z1, cmulti *x0, cmulti *x1, double *y0, double *y1);
void icvec_div_dscalar_cvec(int n, cmulti **z0, cmulti **z1, double x0, double x1, cmulti **y0, cmulti **y1);
void icvec_div_rscalar_zvec(int n, cmulti **z0, cmulti **z1, rmulti *x0, rmulti *x1, dcomplex *y0, dcomplex *y1);
void icvec_div_zscalar_rvec(int n, cmulti **z0, cmulti **z1, dcomplex x0, dcomplex x1, rmulti **y0, rmulti **y1);
// [z0,z1]=[x0,x1]+[-absc(y),absc(y)]
void icvec_add_pm_cvec(int n, cmulti **z0, cmulti **z1, cmulti **x0, cmulti **x1, cmulti **y);
// z=abs(x-y)
void irvec_abs_sub_cvec_cvec(int n, rmulti **z0, rmulti **z1, cmulti **x0, cmulti **x1, cmulti **y0, cmulti **y1);
void irvec_abs_sub_cvec_cscalar(int n, rmulti **z0, rmulti **z1, cmulti **x0, cmulti **x1, cmulti *y0, cmulti *y1);
void irvec_abs_sub_cvec_rscalar(int n, rmulti **z0, rmulti **z1, cmulti **x0, cmulti **x1, rmulti *y0, rmulti *y1);
void irvec_abs_sub_cvec_zscalar(int n, rmulti **z0, rmulti **z1, cmulti **x0, cmulti **x1, dcomplex y0, dcomplex y1);
void irvec_abs_sub_cvec_dscalar(int n, rmulti **z0, rmulti **z1, cmulti **x0, cmulti **x1, double y0, double y1);


////////////////////////////////
void icvec_mul_cmat_cvec(int m, int n, cmulti **y0, cmulti **y1, cmulti **A0, int LDA0, cmulti **A1, int LDA1, cmulti **x0, cmulti **x1);        // [y0,y1]=[A0,A1]*[x0,x1]
void icvec_add_lintr(int m, int n, cmulti **y0, cmulti **y1, cmulti **A0, int LDA0, cmulti **A1, int LDA1, cmulti **x0, cmulti **x1);    // [y0,y1]+=[A0,A1]*[x0,x1]
void icvec_sub_lintr(int m, int n, cmulti **y0, cmulti **y1, cmulti **A0, int LDA0, cmulti **A1, int LDA1, cmulti **x0, cmulti **x1);    // [y0,y1]-=[A0,A1]*[x0,x1]
void icvec_lintr_t(int m, int n, cmulti **y0, cmulti **y1, cmulti **A0, int LDA0, cmulti **A1, int LDA1, cmulti **x0, cmulti **x1);      // [y0,y1]=[A0,A1]'*[x0,x1]
void icvec_add_lintr_t(int m, int n, cmulti **y0, cmulti **y1, cmulti **A0, int LDA0, cmulti **A1, int LDA1, cmulti **x0, cmulti **x1);  // [y0,y1]+=[A0,A1]'*[x0,x1]
void icvec_sub_lintr_t(int m, int n, cmulti **y0, cmulti **y1, cmulti **A0, int LDA0, cmulti **A1, int LDA1, cmulti **x0, cmulti **x1);  // [y0,y1]-=[A0,A1]'*[x0,x1]
void icvec_lintr_ct(int m, int n, cmulti **y0, cmulti **y1, cmulti **A0, int LDA0, cmulti **A1, int LDA1, cmulti **x0, cmulti **x1);     // [y0,y1]=[A0,A1]'*[x0,x1]
void icvec_add_lintr_ct(int m, int n, cmulti **y0, cmulti **y1, cmulti **A0, int LDA0, cmulti **A1, int LDA1, cmulti **x0, cmulti **x1); // [y0,y1]+=[A0,A1]'*[x0,x1]
void icvec_sub_lintr_ct(int m, int n, cmulti **y0, cmulti **y1, cmulti **A0, int LDA0, cmulti **A1, int LDA1, cmulti **x0, cmulti **x1); // [y0,y1]-=[A0,A1]'*[x0,x1]

#endif
