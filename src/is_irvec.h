#ifndef ISYS_IRVEC_H
#define ISYS_IRVEC_H

#include<is_rmulti.h>
#include<is_rvec.h>
#include<is_func.h>

/**
 @file  is_irvec.h
 @brief 多倍長精度実数型rmultiの機械区間演算のベクトルに関する関数の宣言.
*/

int irvec_set_d(int n, rmulti **y0, rmulti **y1, double *x);                                            // [y0,y1]=x
int irvec_copy(int n, rmulti **y0, rmulti **y1, rmulti **x0, rmulti **x1);                              // [y0,y1]=[x0,x1]
void irvec_print(int n, rmulti **x0, rmulti **x1, const char *name, const char *f, int digits);         // output
int irvec_center_radius(int n, rmulti **xc, rmulti **xr, rmulti **x0, rmulti **x1);                     // xc=(x1+x0)/2, xr=x1-x0
int irvec_neg(int n, rmulti **y0, rmulti **y1, rmulti **x0, rmulti **x1);                               // [y0,y1]=-[x0,x1]
int irvec_pm(int n, rmulti **y0, rmulti **y1, rmulti **x);                                              // [y0,y1]=[-x,x]
int irvec_add_pm(int n, rmulti **z0, rmulti **z1, rmulti **x0, rmulti **x1, rmulti **y);                // [z0,z1]=[x0,z1]+[-y,y]
int irvec_mid(int n, rmulti **mid, rmulti **x0, rmulti **x1);                                           // [m-r,m+r]=[x0,x1]
int irvec_rad(int n, rmulti **rad, rmulti **x0, rmulti **x1);                                           // [m-r,m+r]=[x0,x1]
int irvec_mr(int n, rmulti **mid, rmulti **rad, rmulti **x0, rmulti **x1);                              // [m-r,m+r]=[x0,x1]
int irvec_add(int n, rmulti **z0, rmulti **z1, rmulti **x0, rmulti **x1, rmulti **y0, rmulti **y1);     // [z0,z1]=[x0,z1]+[y0,y1]

int irvec_add_r(int n, rmulti **z0, rmulti **z1, rmulti **x0, rmulti **x1, rmulti *y0, rmulti *y1);     // [z0,z1]=[x0,z1]+[y0,y1]
int irvec_sub(int n, rmulti **z0, rmulti **z1, rmulti **x0, rmulti **x1, rmulti **y0, rmulti **y1);     // [z0,z1]=[x0,z1]-[y0,y1]
int irvec_mul_r(int n, rmulti **z0, rmulti **z1, rmulti **x0, rmulti **x1, rmulti *y0, rmulti *y1);     // [z0,z1]=[x0,z1]*[y0,y1]
int irvec_div(int n, rmulti **z0, rmulti **z1, rmulti **x0, rmulti **x1, rmulti **y0, rmulti **y1);       // [z0,z1]=[x0,x1]/[y0,y1]
int irvec_sum_pow2(rmulti *y0, rmulti *y1, int n, rmulti **x0, rmulti **x1);                            // y=sum(x.^2)

int irvec_abs(int n, rmulti **y0, rmulti **y1, rmulti **x0, rmulti **x1);                               // [y0,y1]=abs([x0,x1]
int irvec_abs_sub(int n, rmulti **z0, rmulti **z1, rmulti **x0, rmulti **x1, rmulti **y0, rmulti **y1);        // [z0,z1]=abs([x0,x1]-[y0,y1])
int irvec_umax(rmulti *y0, rmulti *y1, int n, rmulti **x0, rmulti **x1);     // [y0,y1]=[x0,max(x1)]
int irvec_dmax(rmulti *y0, rmulti *y1, int n, rmulti **x0, rmulti **x1);     // [y0,y1]=[max(x0),x1]
int irvec_umin(rmulti *y0, rmulti *y1, int n, rmulti **x0, rmulti **x1);     // [y0,y1]=[x0,min(x1)]
int irvec_dmin(rmulti *y0, rmulti *y1, int n, rmulti **x0, rmulti **x1);     // [y0,y1]=[min(x0),x1]
int irvec_umax_abs(rmulti *y0, rmulti *y1, int n, rmulti **x0, rmulti **x1); // value=max(abs(x))
int irvec_dmin_abs(rmulti *y0, rmulti *y1, int n, rmulti **x0, rmulti **x1); // value=min(abs(x))
//追加
int irvec_sum(rmulti *y0, rmulti *y1, int n, rmulti **x0, rmulti **x1);      // value=sum(x)
//ここまで

int irvec_lintr(int m, int n, rmulti **y0, rmulti **y1, rmulti **A0, int LDA0, rmulti **A1, int LDA1, rmulti **x0, rmulti **x1);        // [y0,y1]=[A0,A1]*[x0,x1]
int irvec_add_lintr(int m, int n, rmulti **y0, rmulti **y1, rmulti **A0, int LDA0, rmulti **A1, int LDA1, rmulti **x0, rmulti **x1);    // [y0,y1]+=[A0,A1]*[x0,x1]
int irvec_sub_lintr(int m, int n, rmulti **y0, rmulti **y1, rmulti **A0, int LDA0, rmulti **A1, int LDA1, rmulti **x0, rmulti **x1);    // [y0,y1]-=[A0,A1]*[x0,x1]
int irvec_lintr_t(int m, int n, rmulti **y0, rmulti **y1, rmulti **A0, int LDA0, rmulti **A1, int LDA1, rmulti **x0, rmulti **x1);      // [y0,y1]=[A0,A1]'*[x0,x1]
int irvec_add_lintr_t(int m, int n, rmulti **y0, rmulti **y1, rmulti **A0, int LDA0, rmulti **A1, int LDA1, rmulti **x0, rmulti **x1);  // [y0,y1]+=[A0,A1]'*[x0,x1]
int irvec_sub_lintr_t(int m, int n, rmulti **y0, rmulti **y1, rmulti **A0, int LDA0, rmulti **A1, int LDA1, rmulti **x0, rmulti **x1);  // [y0,y1]-=[A0,A1]'*[x0,x1]
int irvec_func(rmulti *y0, rmulti *y1, func_t *f, int n, rmulti **x0, rmulti **x1);               // [y0,y1]=f([x0,x1])
int irvec_func_list(int m, rmulti **y0, rmulti **y1, func_t *f, int n, rmulti **x0, rmulti **x1); // [y0,y1]=f([x0,x1])

#endif

