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
void irmax_up_abc_cvec(rmulti *y0, rmulti *y1, int n, cmulti **x0, cmulti **x1);      // y=max(abs(x))
void irmin_down_abc_cvec(rmulti *y0, rmulti *y1, int n, cmulti **x0, cmulti **x1);    // y=min(abs(x))

/*
 * z=f(x,y)
 */
// [z0,z1]=[x0,z1]+[y0,y1]
void icvec_add_cvec(int n, cmulti **z0, cmulti **z1, cmulti **x0, cmulti **x1, cmulti **y0, cmulti **y1);
void icvec_add_rvec(int n, cmulti **z0, cmulti **z1, cmulti **x0, cmulti **x1, rmulti **y0, rmulti **y1);
void irvec_add_cvec(int n, cmulti **z0, cmulti **z1, rmulti **x0, rmulti **x1, cmulti **y0, cmulti **y1);
void icvec_add_zvec(int n, cmulti **z0, cmulti **z1, cmulti **x0, cmulti **x1, dcomplex *y0, dcomplex *y1);
void izvec_add_cvec(int n, cmulti **z0, cmulti **z1, dcomplex *x0, dcomplex *x1, cmulti **y0, cmulti **y1);
void icvec_add_dvec(int n, cmulti **z0, cmulti **z1, cmulti **x0, cmulti **x1, double *y0, double *y1);
void idvec_add_cvec(int n, cmulti **z0, cmulti **z1, double *x0, double *x1, cmulti **y0, cmulti **y1);
void izvec_add_rvec(int n, cmulti **z0, cmulti **z1, dcomplex *x0, dcomplex *x1, rmulti **y0, rmulti **y1);
void irvec_add_zvec(int n, cmulti **z0, cmulti **z1, rmulti **x0, rmulti **x1, dcomplex *y0, dcomplex *y1);
void icvec_add_cscalar(int n, cmulti **z0, cmulti **z1, cmulti **x0, cmulti **x1, cmulti *y0, cmulti *y1);
void icvec_add_rscalar(int n, cmulti **z0, cmulti **z1, cmulti **x0, cmulti **x1, rmulti *y0, rmulti *y1);
void irvec_add_cscalar(int n, cmulti **z0, cmulti **z1, rmulti **x0, rmulti **x1, cmulti *y0, cmulti *y1);
void icvec_add_zscalar(int n, cmulti **z0, cmulti **z1, cmulti **x0, cmulti **x1, dcomplex y0, dcomplex y1);
void izvec_add_cscalar(int n, cmulti **z0, cmulti **z1, dcomplex *x0, dcomplex *x1, cmulti *y0, cmulti *y1);
void icvec_add_dscalar(int n, cmulti **z0, cmulti **z1, cmulti **x0, cmulti **x1, double y0, double y1);
void idvec_add_cscalar(int n, cmulti **z0, cmulti **z1, double *x0, double *x1, cmulti *y0, cmulti *y1);
void irvec_add_zscalar(int n, cmulti **z0, cmulti **z1, rmulti **x0, rmulti **x1, dcomplex y0, dcomplex y1);
void izvec_add_rscalar(int n, cmulti **z0, cmulti **z1, dcomplex *x0, dcomplex *x1, rmulti *y0, rmulti *y1);
void icscalar_add_cvec(int n, cmulti **z0, cmulti **z1, cmulti *x0, cmulti *x1, cmulti **y0, cmulti **y1);
void icscalar_add_rvec(int n, cmulti **z0, cmulti **z1, cmulti *x0, cmulti *x1, rmulti **y0, rmulti **y1);
void irscalar_add_cvec(int n, cmulti **z0, cmulti **z1, rmulti *x0, rmulti *x1, cmulti **y0, cmulti **y1);
void icscalar_add_zvec(int n, cmulti **z0, cmulti **z1, cmulti *x0, cmulti *x1, dcomplex *y0, dcomplex *y1);
void izscalar_add_cvec(int n, cmulti **z0, cmulti **z1, dcomplex x0, dcomplex x1, cmulti **y0, cmulti **y1);
void icscalar_add_dvec(int n, cmulti **z0, cmulti **z1, cmulti *x0, cmulti *x1, double *y0, double *y1);
void idscalar_add_cvec(int n, cmulti **z0, cmulti **z1, double x0, double x1, cmulti **y0, cmulti **y1);
void irscalar_add_zvec(int n, cmulti **z0, cmulti **z1, rmulti *x0, rmulti *x1, dcomplex *y0, dcomplex *y1);
void izscalar_add_rvec(int n, cmulti **z0, cmulti **z1, dcomplex x0, dcomplex x1, rmulti **y0, rmulti **y1);
// [z0,z1]=[x0,z1]-[y0,y1]
void icvec_sub_cvec(int n, cmulti **z0, cmulti **z1, cmulti **x0, cmulti **x1, cmulti **y0, cmulti **y1);
void icvec_sub_rvec(int n, cmulti **z0, cmulti **z1, cmulti **x0, cmulti **x1, rmulti **y0, rmulti **y1);
void irvec_sub_cvec(int n, cmulti **z0, cmulti **z1, rmulti **x0, rmulti **x1, cmulti **y0, cmulti **y1);
void icvec_sub_zvec(int n, cmulti **z0, cmulti **z1, cmulti **x0, cmulti **x1, dcomplex *y0, dcomplex *y1);
void izvec_sub_cvec(int n, cmulti **z0, cmulti **z1, dcomplex *x0, dcomplex *x1, cmulti **y0, cmulti **y1);
void icvec_sub_dvec(int n, cmulti **z0, cmulti **z1, cmulti **x0, cmulti **x1, double *y0, double *y1);
void idvec_sub_cvec(int n, cmulti **z0, cmulti **z1, double *x0, double *x1, cmulti **y0, cmulti **y1);
void izvec_sub_rvec(int n, cmulti **z0, cmulti **z1, dcomplex *x0, dcomplex *x1, rmulti **y0, rmulti **y1);
void irvec_sub_zvec(int n, cmulti **z0, cmulti **z1, rmulti **x0, rmulti **x1, dcomplex *y0, dcomplex *y1);
void icvec_sub_cscalar(int n, cmulti **z0, cmulti **z1, cmulti **x0, cmulti **x1, cmulti *y0, cmulti *y1);
void icvec_sub_rscalar(int n, cmulti **z0, cmulti **z1, cmulti **x0, cmulti **x1, rmulti *y0, rmulti *y1);
void irvec_sub_cscalar(int n, cmulti **z0, cmulti **z1, rmulti **x0, rmulti **x1, cmulti *y0, cmulti *y1);
void icvec_sub_zscalar(int n, cmulti **z0, cmulti **z1, cmulti **x0, cmulti **x1, dcomplex y0, dcomplex y1);
void izvec_sub_cscalar(int n, cmulti **z0, cmulti **z1, dcomplex *x0, dcomplex *x1, cmulti *y0, cmulti *y1);
void icvec_sub_dscalar(int n, cmulti **z0, cmulti **z1, cmulti **x0, cmulti **x1, double y0, double y1);
void idvec_sub_cscalar(int n, cmulti **z0, cmulti **z1, double *x0, double *x1, cmulti *y0, cmulti *y1);
void irvec_sub_zscalar(int n, cmulti **z0, cmulti **z1, rmulti **x0, rmulti **x1, dcomplex y0, dcomplex y1);
void izvec_sub_rscalar(int n, cmulti **z0, cmulti **z1, dcomplex *x0, dcomplex *x1, rmulti *y0, rmulti *y1);
void icscalar_sub_cvec(int n, cmulti **z0, cmulti **z1, cmulti *x0, cmulti *x1, cmulti **y0, cmulti **y1);
void icscalar_sub_rvec(int n, cmulti **z0, cmulti **z1, cmulti *x0, cmulti *x1, rmulti **y0, rmulti **y1);
void irscalar_sub_cvec(int n, cmulti **z0, cmulti **z1, rmulti *x0, rmulti *x1, cmulti **y0, cmulti **y1);
void icscalar_sub_zvec(int n, cmulti **z0, cmulti **z1, cmulti *x0, cmulti *x1, dcomplex *y0, dcomplex *y1);
void izscalar_sub_cvec(int n, cmulti **z0, cmulti **z1, dcomplex x0, dcomplex x1, cmulti **y0, cmulti **y1);
void icscalar_sub_dvec(int n, cmulti **z0, cmulti **z1, cmulti *x0, cmulti *x1, double *y0, double *y1);
void idscalar_sub_cvec(int n, cmulti **z0, cmulti **z1, double x0, double x1, cmulti **y0, cmulti **y1);
void irscalar_sub_zvec(int n, cmulti **z0, cmulti **z1, rmulti *x0, rmulti *x1, dcomplex *y0, dcomplex *y1);
void izscalar_sub_rvec(int n, cmulti **z0, cmulti **z1, dcomplex x0, dcomplex x1, rmulti **y0, rmulti **y1);
// [z0,z1]=[x0,x1]*[y0,y1]
void icvec_mul_cvec(int n, cmulti **z0, cmulti **z1, cmulti **x0, cmulti **x1, cmulti **y0, cmulti **y1);
void icvec_mul_rvec(int n, cmulti **z0, cmulti **z1, cmulti **x0, cmulti **x1, rmulti **y0, rmulti **y1);
void irvec_mul_cvec(int n, cmulti **z0, cmulti **z1, rmulti **x0, rmulti **x1, cmulti **y0, cmulti **y1);
void icvec_mul_zvec(int n, cmulti **z0, cmulti **z1, cmulti **x0, cmulti **x1, dcomplex *y0, dcomplex *y1);
void izvec_mul_cvec(int n, cmulti **z0, cmulti **z1, dcomplex *x0, dcomplex *x1, cmulti **y0, cmulti **y1);
void icvec_mul_dvec(int n, cmulti **z0, cmulti **z1, cmulti **x0, cmulti **x1, double *y0, double *y1);
void idvec_mul_cvec(int n, cmulti **z0, cmulti **z1, double *x0, double *x1, cmulti **y0, cmulti **y1);
void izvec_mul_rvec(int n, cmulti **z0, cmulti **z1, dcomplex *x0, dcomplex *x1, rmulti **y0, rmulti **y1);
void irvec_mul_zvec(int n, cmulti **z0, cmulti **z1, rmulti **x0, rmulti **x1, dcomplex *y0, dcomplex *y1);
void icvec_mul_cscalar(int n, cmulti **z0, cmulti **z1, cmulti **x0, cmulti **x1, cmulti *y0, cmulti *y1);
void icvec_mul_rscalar(int n, cmulti **z0, cmulti **z1, cmulti **x0, cmulti **x1, rmulti *y0, rmulti *y1);
void irvec_mul_cscalar(int n, cmulti **z0, cmulti **z1, rmulti **x0, rmulti **x1, cmulti *y0, cmulti *y1);
void icvec_mul_zscalar(int n, cmulti **z0, cmulti **z1, cmulti **x0, cmulti **x1, dcomplex y0, dcomplex y1);
void izvec_mul_cscalar(int n, cmulti **z0, cmulti **z1, dcomplex *x0, dcomplex *x1, cmulti *y0, cmulti *y1);
void icvec_mul_dscalar(int n, cmulti **z0, cmulti **z1, cmulti **x0, cmulti **x1, double y0, double y1);
void idvec_mul_cscalar(int n, cmulti **z0, cmulti **z1, double *x0, double *x1, cmulti *y0, cmulti *y1);
void irvec_mul_zscalar(int n, cmulti **z0, cmulti **z1, rmulti **x0, rmulti **x1, dcomplex y0, dcomplex y1);
void izvec_mul_rscalar(int n, cmulti **z0, cmulti **z1, dcomplex *x0, dcomplex *x1, rmulti *y0, rmulti *y1);
void icscalar_mul_cvec(int n, cmulti **z0, cmulti **z1, cmulti *x0, cmulti *x1, cmulti **y0, cmulti **y1);
void icscalar_mul_rvec(int n, cmulti **z0, cmulti **z1, cmulti *x0, cmulti *x1, rmulti **y0, rmulti **y1);
void irscalar_mul_cvec(int n, cmulti **z0, cmulti **z1, rmulti *x0, rmulti *x1, cmulti **y0, cmulti **y1);
void icscalar_mul_zvec(int n, cmulti **z0, cmulti **z1, cmulti *x0, cmulti *x1, dcomplex *y0, dcomplex *y1);
void izscalar_mul_cvec(int n, cmulti **z0, cmulti **z1, dcomplex x0, dcomplex x1, cmulti **y0, cmulti **y1);
void icscalar_mul_dvec(int n, cmulti **z0, cmulti **z1, cmulti *x0, cmulti *x1, double *y0, double *y1);
void idscalar_mul_cvec(int n, cmulti **z0, cmulti **z1, double x0, double x1, cmulti **y0, cmulti **y1);
void irscalar_mul_zvec(int n, cmulti **z0, cmulti **z1, rmulti *x0, rmulti *x1, dcomplex *y0, dcomplex *y1);
void izscalar_mul_rvec(int n, cmulti **z0, cmulti **z1, dcomplex x0, dcomplex x1, rmulti **y0, rmulti **y1);
// [z0,z1]=[x0,x1]/[y0,y1]
void icvec_div_cvec(int n, cmulti **z0, cmulti **z1, cmulti **x0, cmulti **x1, cmulti **y0, cmulti **y1);
void icvec_div_rvec(int n, cmulti **z0, cmulti **z1, cmulti **x0, cmulti **x1, rmulti **y0, rmulti **y1);
void irvec_div_cvec(int n, cmulti **z0, cmulti **z1, rmulti **x0, rmulti **x1, cmulti **y0, cmulti **y1);
void icvec_div_zvec(int n, cmulti **z0, cmulti **z1, cmulti **x0, cmulti **x1, dcomplex *y0, dcomplex *y1);
void izvec_div_cvec(int n, cmulti **z0, cmulti **z1, dcomplex *x0, dcomplex *x1, cmulti **y0, cmulti **y1);
void icvec_div_dvec(int n, cmulti **z0, cmulti **z1, cmulti **x0, cmulti **x1, double *y0, double *y1);
void idvec_div_cvec(int n, cmulti **z0, cmulti **z1, double *x0, double *x1, cmulti **y0, cmulti **y1);
void izvec_div_rvec(int n, cmulti **z0, cmulti **z1, dcomplex *x0, dcomplex *x1, rmulti **y0, rmulti **y1);
void irvec_div_zvec(int n, cmulti **z0, cmulti **z1, rmulti **x0, rmulti **x1, dcomplex *y0, dcomplex *y1);
void icvec_div_cscalar(int n, cmulti **z0, cmulti **z1, cmulti **x0, cmulti **x1, cmulti *y0, cmulti *y1);
void icvec_div_rscalar(int n, cmulti **z0, cmulti **z1, cmulti **x0, cmulti **x1, rmulti *y0, rmulti *y1);
void irvec_div_cscalar(int n, cmulti **z0, cmulti **z1, rmulti **x0, rmulti **x1, cmulti *y0, cmulti *y1);
void icvec_div_zscalar(int n, cmulti **z0, cmulti **z1, cmulti **x0, cmulti **x1, dcomplex y0, dcomplex y1);
void izvec_div_cscalar(int n, cmulti **z0, cmulti **z1, dcomplex *x0, dcomplex *x1, cmulti *y0, cmulti *y1);
void icvec_div_dscalar(int n, cmulti **z0, cmulti **z1, cmulti **x0, cmulti **x1, double y0, double y1);
void idvec_div_cscalar(int n, cmulti **z0, cmulti **z1, double *x0, double *x1, cmulti *y0, cmulti *y1);
void irvec_div_zscalar(int n, cmulti **z0, cmulti **z1, rmulti **x0, rmulti **x1, dcomplex y0, dcomplex y1);
void izvec_div_rscalar(int n, cmulti **z0, cmulti **z1, dcomplex *x0, dcomplex *x1, rmulti *y0, rmulti *y1);
void icscalar_div_cvec(int n, cmulti **z0, cmulti **z1, cmulti *x0, cmulti *x1, cmulti **y0, cmulti **y1);
void icscalar_div_rvec(int n, cmulti **z0, cmulti **z1, cmulti *x0, cmulti *x1, rmulti **y0, rmulti **y1);
void irscalar_div_cvec(int n, cmulti **z0, cmulti **z1, rmulti *x0, rmulti *x1, cmulti **y0, cmulti **y1);
void icscalar_div_zvec(int n, cmulti **z0, cmulti **z1, cmulti *x0, cmulti *x1, dcomplex *y0, dcomplex *y1);
void izscalar_div_cvec(int n, cmulti **z0, cmulti **z1, dcomplex x0, dcomplex x1, cmulti **y0, cmulti **y1);
void icscalar_div_dvec(int n, cmulti **z0, cmulti **z1, cmulti *x0, cmulti *x1, double *y0, double *y1);
void idscalar_div_cvec(int n, cmulti **z0, cmulti **z1, double x0, double x1, cmulti **y0, cmulti **y1);
void irscalar_div_zvec(int n, cmulti **z0, cmulti **z1, rmulti *x0, rmulti *x1, dcomplex *y0, dcomplex *y1);
void izscalar_div_rvec(int n, cmulti **z0, cmulti **z1, dcomplex x0, dcomplex x1, rmulti **y0, rmulti **y1);

// xc=(x1+x0)/2, xr=x1-x0
void icvec_center_radius(int n, cmulti **xc, cmulti **xr, cmulti **x0, cmulti **x1);
void icvec_pm(int n, cmulti **y0, cmulti **y1, cmulti **x);                                          // [y0,y1]=[-x,x]





// [z0,z1]=[x0,x1]+[-y,y]
void icvec_add_pm(int n, cmulti **z0, cmulti **z1, cmulti **x0, cmulti **x1, cmulti **y);
void icvec_mid(int n, cmulti **mid, cmulti **x0, cmulti **x1);                                       // [m-r,m+r]=[x0,x1]
void icvec_rad(int n, cmulti **rad, cmulti **x0, cmulti **x1);                                       // [m-r,m+r]=[x0,x1]
void icvec_mr(int n, cmulti **mid, cmulti **rad, cmulti **x0, cmulti **x1);                          // [m-r,m+r]=[x0,x1]
void icvec_add_cscalar(int n, cmulti **z0, cmulti **z1, cmulti **x0, cmulti **x1, cmulti *y0, cmulti *y1); // [z0,z1]=[x0,z1]+[y0,y1]
void icvec_abs_sub(int n, rmulti **z0, rmulti **z1, cmulti **x0, cmulti **x1, cmulti **y0, cmulti **y1); // [z0,z1]=abs([x0,x1]-[y0,y1])
void icvec_lintr(int m, int n, cmulti **y0, cmulti **y1, cmulti **A0, int LDA0, cmulti **A1, int LDA1, cmulti **x0, cmulti **x1);        // [y0,y1]=[A0,A1]*[x0,x1]
void icvec_add_lintr(int m, int n, cmulti **y0, cmulti **y1, cmulti **A0, int LDA0, cmulti **A1, int LDA1, cmulti **x0, cmulti **x1);    // [y0,y1]+=[A0,A1]*[x0,x1]
void icvec_sub_lintr(int m, int n, cmulti **y0, cmulti **y1, cmulti **A0, int LDA0, cmulti **A1, int LDA1, cmulti **x0, cmulti **x1);    // [y0,y1]-=[A0,A1]*[x0,x1]
void icvec_lintr_t(int m, int n, cmulti **y0, cmulti **y1, cmulti **A0, int LDA0, cmulti **A1, int LDA1, cmulti **x0, cmulti **x1);      // [y0,y1]=[A0,A1]'*[x0,x1]
void icvec_add_lintr_t(int m, int n, cmulti **y0, cmulti **y1, cmulti **A0, int LDA0, cmulti **A1, int LDA1, cmulti **x0, cmulti **x1);  // [y0,y1]+=[A0,A1]'*[x0,x1]
void icvec_sub_lintr_t(int m, int n, cmulti **y0, cmulti **y1, cmulti **A0, int LDA0, cmulti **A1, int LDA1, cmulti **x0, cmulti **x1);  // [y0,y1]-=[A0,A1]'*[x0,x1]
void icvec_lintr_ct(int m, int n, cmulti **y0, cmulti **y1, cmulti **A0, int LDA0, cmulti **A1, int LDA1, cmulti **x0, cmulti **x1);     // [y0,y1]=[A0,A1]'*[x0,x1]
void icvec_add_lintr_ct(int m, int n, cmulti **y0, cmulti **y1, cmulti **A0, int LDA0, cmulti **A1, int LDA1, cmulti **x0, cmulti **x1); // [y0,y1]+=[A0,A1]'*[x0,x1]
void icvec_sub_lintr_ct(int m, int n, cmulti **y0, cmulti **y1, cmulti **A0, int LDA0, cmulti **A1, int LDA1, cmulti **x0, cmulti **x1); // [y0,y1]-=[A0,A1]'*[x0,x1]

#endif
