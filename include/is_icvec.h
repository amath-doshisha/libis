#ifndef ISYS_ICVEC_H
#define ISYS_ICVEC_H

#include<is_cmulti.h>


/**
 @file  is_icvec.h
 @brief 多倍長精度実数型cmultiの機械区間演算のベクトルに関する関数の宣言.
*/

/*
 * setting
 */
// [y0,y1]=x
void icvec_set_s(int n, cmulti **y0, cmulti **y1, char **x);
void icvec_set_si(int n, cmulti **y0, cmulti **y1, int *x);
void icvec_set(int n, cmulti **y0, cmulti **y1, cmulti **x0, cmulti **x1);
void icvec_set_d(int n, cmulti **y0, cmulti **y1, double *x0, double *x1);
void icvec_set_z(int n, cmulti **y0, cmulti **y1, dcomplex *x0, dcomplex *x1);
void icvec_set_r(int n, cmulti **y0, cmulti **y1, rmulti **x0, rmulti **x1);
void icvec_set_dz(int n, cmulti **y0, cmulti **y1, double *x0, dcomplex *x1);
void icvec_set_zd(int n, cmulti **y0, cmulti **y1, dcomplex *x0, double *x1);
void icvec_set_cr(int n, cmulti **y0, cmulti **y1, cmulti **x0, rmulti **x1);
void icvec_set_rc(int n, cmulti **y0, cmulti **y1, rmulti **x0, cmulti **x1);
void dvec_set_ic(int n, double *y, cmulti **x0, cmulti **x1);
void zvec_set_ic(int n, dcomplex *y, cmulti **x0, cmulti **x1);

/*
 * I/O
 */
void icvec_print(int n, cmulti **x0, cmulti **x1, char *name, char format, int digits);

/*
 * casting
 */
void icvec_get_si(int n, int *y, cmulti **x0, cmulti **x1);
void icvec_get_z(int n, dcomplex *y, cmulti **x0, cmulti **x1);
void icvec_get_s(int n, char **y, cmulti **x0, cmulti **x1, char format, int digits);
void icvec_get_c(int n, cmulti **y, cmulti **x0, cmulti **x1);

/*
 * operations
 */
// [y0,y1]=[x0,x1]
void icvec_copy(int n, cmulti **y0, cmulti **y1, cmulti **x0, cmulti **x1);
void icvec_copy_r(int n, cmulti **y0, cmulti **y1, rmulti **x0, rmulti **x1);
void icvec_copy_cr(int n, cmulti **y0, cmulti **y1, cmulti **x0, rmulti **x1);
void icvec_copy_rc(int n, cmulti **y0, cmulti **y1, rmulti **x0, cmulti **x1);
void icvec_copy_C(int n, cmulti **y0, cmulti **y1, rmulti **x0r, rmulti **x0i, rmulti **x1r, rmulti **x1i);
// [y0,y1]=real([x0,x1])
void icvec_real(int n, rmulti **y0, rmulti **y1, cmulti **x0, cmulti **x1);
// [y0,y1]=conj([x0,x1])
void icvec_conj(int n, cmulti **y0, cmulti **y1, cmulti **x0, cmulti **x1);
// xc=(x1+x0)/2, xr=x1-x0
void icvec_center_radius(int n, cmulti **xc, cmulti **xr, cmulti **x0, cmulti **x1);
void icvec_neg(int n, cmulti **y0, cmulti **y1, cmulti **x0, cmulti **x1);                           // [y0,y1]=-[x0,x1]
void icvec_pm(int n, cmulti **y0, cmulti **y1, cmulti **x);                                          // [y0,y1]=[-x,x]
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
// [z0,z1]=[x0,x1]+[-y,y]
void icvec_add_pm(int n, cmulti **z0, cmulti **z1, cmulti **x0, cmulti **x1, cmulti **y);
void icvec_mid(int n, cmulti **mid, cmulti **x0, cmulti **x1);                                       // [m-r,m+r]=[x0,x1]
void icvec_rad(int n, cmulti **rad, cmulti **x0, cmulti **x1);                                       // [m-r,m+r]=[x0,x1]
void icvec_mr(int n, cmulti **mid, cmulti **rad, cmulti **x0, cmulti **x1);                          // [m-r,m+r]=[x0,x1]
void icvec_add_cscalar(int n, cmulti **z0, cmulti **z1, cmulti **x0, cmulti **x1, cmulti *y0, cmulti *y1); // [z0,z1]=[x0,z1]+[y0,y1]
void icvec_mul_c(int n, cmulti **z0, cmulti **z1, cmulti **x0, cmulti **x1, cmulti *y0, cmulti *y1); // [z0,z1]=[x0,x1]*[y0,y1]
void icvec_sum_abs2(rmulti *y0, rmulti *y1, int n, cmulti **x0, cmulti **x1);                        // [y0,y1]=sum(abs([x0,x1]).^2)
void icvec_abs(int n, rmulti **y0, rmulti **y1, cmulti **x0, cmulti **x1);                           // [y0,y1]=abs([x0,x1])
void icvec_abs_sub(int n, rmulti **z0, rmulti **z1, cmulti **x0, cmulti **x1, cmulti **y0, cmulti **y1); // [z0,z1]=abs([x0,x1]-[y0,y1])
void icvec_max(cmulti *y0, cmulti *y1, int n, cmulti **x0, cmulti **x1);                             // [y0,y1]=[max(x0),max(x1)]
void icvec_umax(cmulti *y0, cmulti *y1, int n, cmulti **x0, cmulti **x1);                            // [y0,y1]=[x0,max(x1)]
void icvec_dmax(cmulti *y0, cmulti *y1, int n, cmulti **x0, cmulti **x1);                            // [y0,y1]=[max(x0),x1]
void icvec_min(cmulti *y0, cmulti *y1, int n, cmulti **x0, cmulti **x1);                             // [y0,y1]=[min(x0),min(x1)]
void icvec_umin(cmulti *y0, cmulti *y1, int n, cmulti **x0, cmulti **x1);                            // [y0,y1]=[x0,min(x1)]
void icvec_dmin(cmulti *y0, cmulti *y1, int n, cmulti **x0, cmulti **x1);                            // [y0,y1]=[min(x0),x1]
void icvec_umax_abs(rmulti *y0, rmulti *y1, int n, cmulti **x0, cmulti **x1);                        // value=max(abs(x))
void icvec_dmin_abs(rmulti *y0, rmulti *y1, int n, cmulti **x0, cmulti **x1);                        // value=min(abs(x))
void icvec_sum(cmulti *y0, cmulti *y1, int n, cmulti **x0, cmulti **x1);                             // value=sum(x)
void icvec_lintr(int m, int n, cmulti **y0, cmulti **y1, cmulti **A0, int LDA0, cmulti **A1, int LDA1, cmulti **x0, cmulti **x1);        // [y0,y1]=[A0,A1]*[x0,x1]
void icvec_add_lintr(int m, int n, cmulti **y0, cmulti **y1, cmulti **A0, int LDA0, cmulti **A1, int LDA1, cmulti **x0, cmulti **x1);    // [y0,y1]+=[A0,A1]*[x0,x1]
void icvec_sub_lintr(int m, int n, cmulti **y0, cmulti **y1, cmulti **A0, int LDA0, cmulti **A1, int LDA1, cmulti **x0, cmulti **x1);    // [y0,y1]-=[A0,A1]*[x0,x1]
void icvec_lintr_t(int m, int n, cmulti **y0, cmulti **y1, cmulti **A0, int LDA0, cmulti **A1, int LDA1, cmulti **x0, cmulti **x1);      // [y0,y1]=[A0,A1]'*[x0,x1]
void icvec_add_lintr_t(int m, int n, cmulti **y0, cmulti **y1, cmulti **A0, int LDA0, cmulti **A1, int LDA1, cmulti **x0, cmulti **x1);  // [y0,y1]+=[A0,A1]'*[x0,x1]
void icvec_sub_lintr_t(int m, int n, cmulti **y0, cmulti **y1, cmulti **A0, int LDA0, cmulti **A1, int LDA1, cmulti **x0, cmulti **x1);  // [y0,y1]-=[A0,A1]'*[x0,x1]
void icvec_lintr_ct(int m, int n, cmulti **y0, cmulti **y1, cmulti **A0, int LDA0, cmulti **A1, int LDA1, cmulti **x0, cmulti **x1);     // [y0,y1]=[A0,A1]'*[x0,x1]
void icvec_add_lintr_ct(int m, int n, cmulti **y0, cmulti **y1, cmulti **A0, int LDA0, cmulti **A1, int LDA1, cmulti **x0, cmulti **x1); // [y0,y1]+=[A0,A1]'*[x0,x1]
void icvec_sub_lintr_ct(int m, int n, cmulti **y0, cmulti **y1, cmulti **A0, int LDA0, cmulti **A1, int LDA1, cmulti **x0, cmulti **x1); // [y0,y1]-=[A0,A1]'*[x0,x1]
void icvec_func(cmulti *y0, cmulti *y1, func_t *f, int n, cmulti **x0, cmulti **x1);                // [y0,y1]=f(x)
void icvec_func_list(int m, cmulti **y0, cmulti **y1, func_t *f, int n, cmulti **x0, cmulti **x1);  // [y0,y1]=f(x)

#endif
