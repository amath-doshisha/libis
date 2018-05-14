#ifndef ISYS_IRVEC_H
#define ISYS_IRVEC_H

#include<is_dcomplex.h>
#include<is_rmulti.h>

/**
 @file  is_irvec.h
 @brief 多倍長精度実数型rmultiの機械区間演算のベクトルに関する関数の宣言.
*/

/*
 * I/O
 */
void irvec_print(int n, rmulti **x0, rmulti **x1, char *name, char format, int digits);

/*
 * casting
 */
void irvec_set_rvec(int n, rmulti **y0, rmulti **y1, rmulti **x0, rmulti **x1);          // [rmulti,rmulti] <- [rmulti,rmulti]
void irvec_set_dvec(int n, rmulti **y0, rmulti **y1, double *x0, double *x1);            // [rmulti,rmulti] <- [double,double]
void irvec_set_ivec(int n, rmulti **y0, rmulti **y1, int *x0, int *x1);                  // [rmulti,rmulti] <- [int,int]
void irvec_set_zvec(int n, rmulti **y0, rmulti **y1, dcomplex *x0, dcomplex *x1);        // [rmulti,rmulti] <- [dcomplex,dcomplex]
void rvec_set_irvec(int n, rmulti **y, rmulti **x0, rmulti **x1);                        // [rmulti,rmulti] -> rmulti
void dvec_set_irvec(int n, double *y, rmulti **x0, rmulti **x1);                         // [rmulti,rmulti] -> double
void ivec_set_irvec(int n, int *y, rmulti **x0, rmulti **x);                             // [rmulti,rmulti] -> int
void zvec_set_irvec(int n, dcomplex *y, rmulti **x0, rmulti **x1);                       // [rmulti,rmulti] -> dcomplex
void irvec_get_rvec(int n, rmulti **y, rmulti **x0, rmulti **x1);                        // [rmulti,rmulti] -> rmulti
void irvec_get_dvec(int n, double *y, rmulti **x0, rmulti **x1);                         // [rmulti,rmulti] -> double
void irvec_get_ivec(int n, int *y, rmulti **x0, rmulti **x);                             // [rmulti,rmulti] -> int
void irvec_set_svec(int n, rmulti **y0, rmulti **y1, char **x);                          // [rmulti,rmulti] <- char
void irvec_get_svec(int n, char **y, rmulti **x0, rmulti **x1, char format, int digits); // [rmulti,rmulti] -> char

/*
 * z=f(x,y)
 */
// [z0,z1]=[x0,z1]+[y0,y1]
void irvec_add_rvec(int n, rmulti **z0, rmulti **z1, rmulti **x0, rmulti **x1, rmulti **y0, rmulti **y1);
void irvec_add_dvec(int n, rmulti **z0, rmulti **z1, rmulti **x0, rmulti **x1, double *y0, double *y1);
void idvec_add_rvec(int n, rmulti **z0, rmulti **z1, double *x0, double *x1, rmulti **y0, rmulti **y1);
void irvec_add_rscalar(int n, rmulti **z0, rmulti **z1, rmulti **x0, rmulti **x1, rmulti *y0, rmulti *y1);
void irvec_add_dscalar(int n, rmulti **z0, rmulti **z1, rmulti **x0, rmulti **x1, double y0, double y1);
void idvec_add_rscalar(int n, rmulti **z0, rmulti **z1, double *x0, double *x1, rmulti *y0, rmulti *y1);
void irscalar_add_rvec(int n, rmulti **z0, rmulti **z1, rmulti *x0, rmulti *x1, rmulti **y0, rmulti **y1);
void irscalar_add_dvec(int n, rmulti **z0, rmulti **z1, rmulti *x0, rmulti *x1, double *y0, double *y1);
void idscalar_add_rvec(int n, rmulti **z0, rmulti **z1, double x0, double x1, rmulti **y0, rmulti **y1);
// [z0,z1]=[x0,z1]-[y0,y1]
void irvec_sub_rvec(int n, rmulti **z0, rmulti **z1, rmulti **x0, rmulti **x1, rmulti **y0, rmulti **y1);
void irvec_sub_dvec(int n, rmulti **z0, rmulti **z1, rmulti **x0, rmulti **x1, double *y0, double *y1);
void idvec_sub_rvec(int n, rmulti **z0, rmulti **z1, double *x0, double *x1, rmulti **y0, rmulti **y1);
void irvec_sub_rscalar(int n, rmulti **z0, rmulti **z1, rmulti **x0, rmulti **x1, rmulti *y0, rmulti *y1);
void irvec_sub_dscalar(int n, rmulti **z0, rmulti **z1, rmulti **x0, rmulti **x1, double y0, double y1);
void idvec_sub_rscalar(int n, rmulti **z0, rmulti **z1, double *x0, double *x1, rmulti *y0, rmulti *y1);
void irscalar_sub_rvec(int n, rmulti **z0, rmulti **z1, rmulti *x0, rmulti *x1, rmulti **y0, rmulti **y1);
void irscalar_sub_dvec(int n, rmulti **z0, rmulti **z1, rmulti *x0, rmulti *x1, double *y0, double *y1);
void idscalar_sub_rvec(int n, rmulti **z0, rmulti **z1, double x0, double x1, rmulti **y0, rmulti **y1);
// [z0,z1]=[x0,z1]*[y0,y1]
void irvec_mul_rvec(int n, rmulti **z0, rmulti **z1, rmulti **x0, rmulti **x1, rmulti **y0, rmulti **y1);
void irvec_mul_dvec(int n, rmulti **z0, rmulti **z1, rmulti **x0, rmulti **x1, double *y0, double *y1);
void idvec_mul_rvec(int n, rmulti **z0, rmulti **z1, double *x0, double *x1, rmulti **y0, rmulti **y1);
void irvec_mul_rscalar(int n, rmulti **z0, rmulti **z1, rmulti **x0, rmulti **x1, rmulti *y0, rmulti *y1);
void irvec_mul_dscalar(int n, rmulti **z0, rmulti **z1, rmulti **x0, rmulti **x1, double y0, double y1);
void idvec_mul_rscalar(int n, rmulti **z0, rmulti **z1, double *x0, double *x1, rmulti *y0, rmulti *y1);
void irscalar_mul_rvec(int n, rmulti **z0, rmulti **z1, rmulti *x0, rmulti *x1, rmulti **y0, rmulti **y1);
void irscalar_mul_dvec(int n, rmulti **z0, rmulti **z1, rmulti *x0, rmulti *x1, double *y0, double *y1);
void idscalar_mul_rvec(int n, rmulti **z0, rmulti **z1, double x0, double x1, rmulti **y0, rmulti **y1);
// [z0,z1]=[x0,x1]/[y0,y1]
void irvec_div_rvec(int n, rmulti **z0, rmulti **z1, rmulti **x0, rmulti **x1, rmulti **y0, rmulti **y1);
void irvec_div_dvec(int n, rmulti **z0, rmulti **z1, rmulti **x0, rmulti **x1, double *y0, double *y1);
void idvec_div_rvec(int n, rmulti **z0, rmulti **z1, double *x0, double *x1, rmulti **y0, rmulti **y1);
void irvec_div_rscalar(int n, rmulti **z0, rmulti **z1, rmulti **x0, rmulti **x1, rmulti *y0, rmulti *y1);
void irvec_div_dscalar(int n, rmulti **z0, rmulti **z1, rmulti **x0, rmulti **x1, double y0, double y1);
void idvec_div_rscalar(int n, rmulti **z0, rmulti **z1, double *x0, double *x1, rmulti *y0, rmulti *y1);
void irscalar_div_rvec(int n, rmulti **z0, rmulti **z1, rmulti *x0, rmulti *x1, rmulti **y0, rmulti **y1);
void irscalar_div_dvec(int n, rmulti **z0, rmulti **z1, rmulti *x0, rmulti *x1, double *y0, double *y1);
void idscalar_div_rvec(int n, rmulti **z0, rmulti **z1, double x0, double x1, rmulti **y0, rmulti **y1);

/*
 * operations
 */
// [y0,y1]=[x0,x1]
void irvec_copy(int n, rmulti **y0, rmulti **y1, rmulti **x0, rmulti **x1);
// y=real(x)
void irvec_center_radius(int n, rmulti **xc, rmulti **xr, rmulti **x0, rmulti **x1);                     // xc=(x1+x0)/2, xr=x1-x0
void irvec_neg(int n, rmulti **y0, rmulti **y1, rmulti **x0, rmulti **x1);                               // [y0,y1]=-[x0,x1]
void irvec_pm(int n, rmulti **y0, rmulti **y1, rmulti **x);                                              // [y0,y1]=[-x,x]
void irvec_add_pm(int n, rmulti **z0, rmulti **z1, rmulti **x0, rmulti **x1, rmulti **y);                // [z0,z1]=[x0,z1]+[-y,y]
void irvec_mid(int n, rmulti **mid, rmulti **x0, rmulti **x1);                                           // [m-r,m+r]=[x0,x1]
void irvec_rad(int n, rmulti **rad, rmulti **x0, rmulti **x1);                                           // [m-r,m+r]=[x0,x1]
void irvec_mr(int n, rmulti **mid, rmulti **rad, rmulti **x0, rmulti **x1);                              // [m-r,m+r]=[x0,x1]
void irvec_sum_pow2(rmulti *y0, rmulti *y1, int n, rmulti **x0, rmulti **x1);                            // y=sum(x.^2)
void irvec_abs(int n, rmulti **y0, rmulti **y1, rmulti **x0, rmulti **x1);                               // [y0,y1]=abs([x0,x1]
void irvec_abs_sub(int n, rmulti **z0, rmulti **z1, rmulti **x0, rmulti **x1, rmulti **y0, rmulti **y1);        // [z0,z1]=abs([x0,x1]-[y0,y1])
void irvec_max(rmulti *y0, rmulti *y1, int n, rmulti **x0, rmulti **x1);      // [y0,y1]=[max(x0),max(x1)]
void irvec_umax(rmulti *y0, rmulti *y1, int n, rmulti **x0, rmulti **x1);     // [y0,y1]=[x0,max(x1)]
void irvec_dmax(rmulti *y0, rmulti *y1, int n, rmulti **x0, rmulti **x1);     // [y0,y1]=[max(x0),x1]
void irvec_min(rmulti *y0, rmulti *y1, int n, rmulti **x0, rmulti **x1);      // [y0,y1]=[min(x0),min(x1)]
void irvec_umin(rmulti *y0, rmulti *y1, int n, rmulti **x0, rmulti **x1);     // [y0,y1]=[x0,min(x1)]
void irvec_dmin(rmulti *y0, rmulti *y1, int n, rmulti **x0, rmulti **x1);     // [y0,y1]=[min(x0),x1]
void irvec_umax_abs(rmulti *y0, rmulti *y1, int n, rmulti **x0, rmulti **x1); // value=max(abs(x))
void irvec_dmin_abs(rmulti *y0, rmulti *y1, int n, rmulti **x0, rmulti **x1); // value=min(abs(x))
void irvec_sum(rmulti *y0, rmulti *y1, int n, rmulti **x0, rmulti **x1);      // value=sum(x)
void irvec_lintr(int m, int n, rmulti **y0, rmulti **y1, rmulti **A0, int LDA0, rmulti **A1, int LDA1, rmulti **x0, rmulti **x1);        // [y0,y1]=[A0,A1]*[x0,x1]
void irvec_add_lintr(int m, int n, rmulti **y0, rmulti **y1, rmulti **A0, int LDA0, rmulti **A1, int LDA1, rmulti **x0, rmulti **x1);    // [y0,y1]+=[A0,A1]*[x0,x1]
void irvec_sub_lintr(int m, int n, rmulti **y0, rmulti **y1, rmulti **A0, int LDA0, rmulti **A1, int LDA1, rmulti **x0, rmulti **x1);    // [y0,y1]-=[A0,A1]*[x0,x1]
void irvec_lintr_t(int m, int n, rmulti **y0, rmulti **y1, rmulti **A0, int LDA0, rmulti **A1, int LDA1, rmulti **x0, rmulti **x1);      // [y0,y1]=[A0,A1]'*[x0,x1]
void irvec_add_lintr_t(int m, int n, rmulti **y0, rmulti **y1, rmulti **A0, int LDA0, rmulti **A1, int LDA1, rmulti **x0, rmulti **x1);  // [y0,y1]+=[A0,A1]'*[x0,x1]
void irvec_sub_lintr_t(int m, int n, rmulti **y0, rmulti **y1, rmulti **A0, int LDA0, rmulti **A1, int LDA1, rmulti **x0, rmulti **x1);  // [y0,y1]-=[A0,A1]'*[x0,x1]

#endif

