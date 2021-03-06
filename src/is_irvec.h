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
 * y=f(x)
 */
void irvec_copy_rvec(int n, rmulti **y0, rmulti **y1, rmulti **x0, rmulti **x1);      // y=x
void irvec_abs_rvec(int n, rmulti **y0, rmulti **y1, rmulti **x0, rmulti **x1);        // y=abs(x)
void irvec_neg_rvec(int n, rmulti **y0, rmulti **y1, rmulti **x0, rmulti **x1);       // y=-x
void irvec_pow2_rvec(int n, rmulti **y0, rmulti **y1, rmulti **x0, rmulti **x1);      // y=x^2
void irsum_rvec(rmulti *y0, rmulti *y1, int n, rmulti **x0, rmulti **x1);             // y=sum(x)
void iraverage_rvec(rmulti *y0, rmulti *y1, int n, rmulti **x0, rmulti **x1);         // y=sum(x)/n
void irsum_abs_rvec(rmulti *y0, rmulti *y1, int n, rmulti **x0, rmulti **x1);         // y=sum(abs(x))
void irsum_pow2_abs_rvec(rmulti *y0, rmulti *y1, int n, rmulti **x0, rmulti **x1);    // y=sum(abs(x)^2)
void irnorm1_rvec(rmulti *y0, rmulti *y1, int n, rmulti **x0, rmulti **x1);           // y=sum(abs(x))
void irnorm2_rvec(rmulti *y0, rmulti *y1, int n, rmulti **x0, rmulti **x1);           // y=sqrt(sum(abs(x)^2))
void irnorm_max_rvec(rmulti *y0, rmulti *y1, int n, rmulti **x0, rmulti **x1);        // y=max(abs(x))
void irmax_rvec(rmulti *y0, rmulti *y1, int n, rmulti **x0, rmulti **x1);             // [y0,y1]=[max(x0),max(x1)]
void irmin_rvec(rmulti *y0, rmulti *y1, int n, rmulti **x0, rmulti **x1);             // [y0,y1]=[min(x0),min(x1)]
void irmax_up_rvec(rmulti *y0, rmulti *y1, int n, rmulti **x0, rmulti **x1);          // [y0,y1]=[x0,max(x1)]
void irmin_up_rvec(rmulti *y0, rmulti *y1, int n, rmulti **x0, rmulti **x1);          // [y0,y1]=[x0,min(x1)]
void irmax_down_rvec(rmulti *y0, rmulti *y1, int n, rmulti **x0, rmulti **x1);        // [y0,y1]=[max(x0),x1]
void irmin_down_rvec(rmulti *y0, rmulti *y1, int n, rmulti **x0, rmulti **x1);        // [y0,y1]=[min(x0),x1]
void irmax_up_abs_rvec(rmulti *y0, rmulti *y1, int n, rmulti **x0, rmulti **x1);      // y=max(abs(x))
void irmin_down_abs_rvec(rmulti *y0, rmulti *y1, int n, rmulti **x0, rmulti **x1);    // y=min(abs(x))
void irvec_pm_rvec(int n, rmulti **y0, rmulti **y1, rmulti **x);                      // [y0,y1]=[-abs(x),abs(x)]
void irvec_mid_rad_rvec(int n, rmulti **mid, rmulti **rad, rmulti **x0, rmulti **x1); // mid=(x1+x0)/2, rad=x1-x0
void irvec_mr_rvec(int n, rmulti **mid, rmulti **rad, rmulti **x0, rmulti **x1);      // [m-r,m+r]=[x0,x1]
void irvec_mid_rvec(int n, rmulti **mid, rmulti **x0, rmulti **x1);                   // [m-r,m+r]=[x0,x1]
void irvec_rad_rvec(int n, rmulti **rad, rmulti **x0, rmulti **x1);                   // [m-r,m+r]=[x0,x1]
void irvec_sqrt_rvec(int n, rmulti **y0, rmulti **y1, rmulti **x0, rmulti **x1);      // [y0,y1]=sqrt([x0,x1])

/*
 * z=f(x,y)
 */
// z=x+y
void irvec_add_rvec_rvec(int n, rmulti **z0, rmulti **z1, rmulti **x0, rmulti **x1, rmulti **y0, rmulti **y1);
void irvec_add_rvec_dvec(int n, rmulti **z0, rmulti **z1, rmulti **x0, rmulti **x1, double *y0, double *y1);
void irvec_add_dvec_rvec(int n, rmulti **z0, rmulti **z1, double *x0, double *x1, rmulti **y0, rmulti **y1);
void irvec_add_rvec_rscalar(int n, rmulti **z0, rmulti **z1, rmulti **x0, rmulti **x1, rmulti *y0, rmulti *y1);
void irvec_add_rvec_dscalar(int n, rmulti **z0, rmulti **z1, rmulti **x0, rmulti **x1, double y0, double y1);
void irvec_add_dvec_rscalar(int n, rmulti **z0, rmulti **z1, double *x0, double *x1, rmulti *y0, rmulti *y1);
void irvec_add_rscalar_rvec(int n, rmulti **z0, rmulti **z1, rmulti *x0, rmulti *x1, rmulti **y0, rmulti **y1);
void irvec_add_rscalar_dvec(int n, rmulti **z0, rmulti **z1, rmulti *x0, rmulti *x1, double *y0, double *y1);
void irvec_add_dscalar_rvec(int n, rmulti **z0, rmulti **z1, double x0, double x1, rmulti **y0, rmulti **y1);
// z=x-y
void irvec_sub_rvec_rvec(int n, rmulti **z0, rmulti **z1, rmulti **x0, rmulti **x1, rmulti **y0, rmulti **y1);
void irvec_sub_rvec_dvec(int n, rmulti **z0, rmulti **z1, rmulti **x0, rmulti **x1, double *y0, double *y1);
void irvec_sub_dvec_rvec(int n, rmulti **z0, rmulti **z1, double *x0, double *x1, rmulti **y0, rmulti **y1);
void irvec_sub_rvec_rscalar(int n, rmulti **z0, rmulti **z1, rmulti **x0, rmulti **x1, rmulti *y0, rmulti *y1);
void irvec_sub_rvec_dscalar(int n, rmulti **z0, rmulti **z1, rmulti **x0, rmulti **x1, double y0, double y1);
void irvec_sub_dvec_rscalar(int n, rmulti **z0, rmulti **z1, double *x0, double *x1, rmulti *y0, rmulti *y1);
void irvec_sub_rscalar_rvec(int n, rmulti **z0, rmulti **z1, rmulti *x0, rmulti *x1, rmulti **y0, rmulti **y1);
void irvec_sub_rscalar_dvec(int n, rmulti **z0, rmulti **z1, rmulti *x0, rmulti *x1, double *y0, double *y1);
void irvec_sub_dscalar_rvec(int n, rmulti **z0, rmulti **z1, double x0, double x1, rmulti **y0, rmulti **y1);
// z=x*y
void irvec_mul_rvec_rvec(int n, rmulti **z0, rmulti **z1, rmulti **x0, rmulti **x1, rmulti **y0, rmulti **y1);
void irvec_mul_rvec_dvec(int n, rmulti **z0, rmulti **z1, rmulti **x0, rmulti **x1, double *y0, double *y1);
void irvec_mul_dvec_rvec(int n, rmulti **z0, rmulti **z1, double *x0, double *x1, rmulti **y0, rmulti **y1);
void irvec_mul_rvec_rscalar(int n, rmulti **z0, rmulti **z1, rmulti **x0, rmulti **x1, rmulti *y0, rmulti *y1);
void irvec_mul_rvec_dscalar(int n, rmulti **z0, rmulti **z1, rmulti **x0, rmulti **x1, double y0, double y1);
void irvec_mul_dvec_rscalar(int n, rmulti **z0, rmulti **z1, double *x0, double *x1, rmulti *y0, rmulti *y1);
void irvec_mul_rscalar_rvec(int n, rmulti **z0, rmulti **z1, rmulti *x0, rmulti *x1, rmulti **y0, rmulti **y1);
void irvec_mul_rscalar_dvec(int n, rmulti **z0, rmulti **z1, rmulti *x0, rmulti *x1, double *y0, double *y1);
void irvec_mul_dscalar_rvec(int n, rmulti **z0, rmulti **z1, double x0, double x1, rmulti **y0, rmulti **y1);
// z=x/y
void irvec_div_rvec_rvec(int n, rmulti **z0, rmulti **z1, rmulti **x0, rmulti **x1, rmulti **y0, rmulti **y1);
void irvec_div_rvec_dvec(int n, rmulti **z0, rmulti **z1, rmulti **x0, rmulti **x1, double *y0, double *y1);
void irvec_div_dvec_rvec(int n, rmulti **z0, rmulti **z1, double *x0, double *x1, rmulti **y0, rmulti **y1);
void irvec_div_rvec_rscalar(int n, rmulti **z0, rmulti **z1, rmulti **x0, rmulti **x1, rmulti *y0, rmulti *y1);
void irvec_div_rvec_dscalar(int n, rmulti **z0, rmulti **z1, rmulti **x0, rmulti **x1, double y0, double y1);
void irvec_div_dvec_rscalar(int n, rmulti **z0, rmulti **z1, double *x0, double *x1, rmulti *y0, rmulti *y1);
void irvec_div_rscalar_rvec(int n, rmulti **z0, rmulti **z1, rmulti *x0, rmulti *x1, rmulti **y0, rmulti **y1);
void irvec_div_rscalar_dvec(int n, rmulti **z0, rmulti **z1, rmulti *x0, rmulti *x1, double *y0, double *y1);
void irvec_div_dscalar_rvec(int n, rmulti **z0, rmulti **z1, double x0, double x1, rmulti **y0, rmulti **y1);
// [z0,z1]=[x0,z1]+[-abs(y),abs(y)]
void irvec_add_pm_rvec(int n, rmulti **z0, rmulti **z1, rmulti **x0, rmulti **x1, rmulti **y);
// z=abs(x-y)
void irvec_abs_sub_rvec_rvec(int n, rmulti **z0, rmulti **z1, rmulti **x0, rmulti **x1, rmulti **y0, rmulti **y1);
void irvec_abs_sub_rvec_rscalar(int n, rmulti **z0, rmulti **z1, rmulti **x0, rmulti **x1, rmulti *y0, rmulti *y1);
void irvec_abs_sub_rvec_dscalar(int n, rmulti **z0, rmulti **z1, rmulti **x0, rmulti **x1, double y0, double y1);

/*
 * operations
 */
// y=A*x
void irvec_mul_rmat_rvec(int m, int n, rmulti **y0, rmulti **y1, rmulti **A0, int LDA0, rmulti **A1, int LDA1, rmulti **x0, rmulti **x1);
// y=y+A*x
void irvec_add_lintr(int m, int n, rmulti **y0, rmulti **y1, rmulti **A0, int LDA0, rmulti **A1, int LDA1, rmulti **x0, rmulti **x1);
// y=y-A*x
void irvec_sub_lintr(int m, int n, rmulti **y0, rmulti **y1, rmulti **A0, int LDA0, rmulti **A1, int LDA1, rmulti **x0, rmulti **x1);
// y=A'*x
void irvec_lintr_t(int m, int n, rmulti **y0, rmulti **y1, rmulti **A0, int LDA0, rmulti **A1, int LDA1, rmulti **x0, rmulti **x1);
// y=y+A'*x
void irvec_add_lintr_t(int m, int n, rmulti **y0, rmulti **y1, rmulti **A0, int LDA0, rmulti **A1, int LDA1, rmulti **x0, rmulti **x1);
// y=y-A'*x
void irvec_sub_lintr_t(int m, int n, rmulti **y0, rmulti **y1, rmulti **A0, int LDA0, rmulti **A1, int LDA1, rmulti **x0, rmulti **x1);

#endif

