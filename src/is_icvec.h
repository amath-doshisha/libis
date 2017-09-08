#ifndef ISYS_ICVEC_H
#define ISYS_ICVEC_H

#include<is_cmulti.h>

int icvec_copy(int n, cmulti **y0, cmulti **y1, cmulti **x0, cmulti **x1);                          // [y0,y1]=[x0,x1]
void icvec_print(int n, cmulti **x0, cmulti **x1, const char *name, const char *f, int digits);     // output
int icvec_center_radius(int n, cmulti **xc, cmulti **xr, cmulti **x0, cmulti **x1);                 // xc=(x1+x0)/2, xr=x1-x0
int icvec_neg(int n, cmulti **y0, cmulti **y1, cmulti **x0, cmulti **x1);                           // [y0,y1]=-[x0,x1]
int icvec_pm(int n, cmulti **y0, cmulti **y1, cmulti **x);                                          // [y0,y1]=[-x,x]
int icvec_add(int n, cmulti **z0, cmulti **z1, cmulti **x0, cmulti **x1, cmulti **y0, cmulti **y1); // [z0,z1]=[x0,z1]+[y0,y1]
int icvec_add_pm(int n, cmulti **z0, cmulti **z1, cmulti **x0, cmulti **x1, cmulti **y);            // [z0,z1]=[x0,x1]+[-y,y]
int icvec_mid(int n, cmulti **mid, cmulti **x0, cmulti **x1);                                       // [m-r,m+r]=[x0,x1]
int icvec_rad(int n, cmulti **rad, cmulti **x0, cmulti **x1);                                       // [m-r,m+r]=[x0,x1]
int icvec_mr(int n, cmulti **mid, cmulti **rad, cmulti **x0, cmulti **x1);                          // [m-r,m+r]=[x0,x1]
int icvec_add_c(int n, cmulti **z0, cmulti **z1, cmulti **x0, cmulti **x1, cmulti *y0, cmulti *y1); // [z0,z1]=[x0,z1]+[y0,y1]
int icvec_sub(int n, cmulti **z0, cmulti **z1, cmulti **x0, cmulti **x1, cmulti **y0, cmulti **y1); // [z0,z1]=[x0,z1]-[y0,y1]
int icvec_mul_c(int n, cmulti **z0, cmulti **z1, cmulti **x0, cmulti **x1, cmulti *y0, cmulti *y1); // [z0,z1]=[x0,x1]*[y0,y1]
int icvec_sum_abs2(rmulti *y0, rmulti *y1, int n, cmulti **x0, cmulti **x1);                        // [y0,y1]=sum(abs([x0,x1]).^2)

int icvec_abs(int n, rmulti **y0, rmulti **y1, cmulti **x0, cmulti **x1);                           // [y0,y1]=abs([x0,x1])
int icvec_abs_sub(int n, rmulti **z0, rmulti **z1, cmulti **x0, cmulti **x1, cmulti **y0, cmulti **y1); // [z0,z1]=abs([x0,x1]-[y0,y1])
int icvec_umax(cmulti *y0, cmulti *y1, int n, cmulti **x0, cmulti **x1);                            // [y0,y1]=[x0,max(x1)]
int icvec_dmax(cmulti *y0, cmulti *y1, int n, cmulti **x0, cmulti **x1);                            // [y0,y1]=[max(x0),x1]
int icvec_umin(rmulti *y0, rmulti *y1, int n, rmulti **x0, rmulti **x1);                            // [y0,y1]=[x0,min(x1)]
int icvec_dmin(cmulti *y0, cmulti *y1, int n, cmulti **x0, cmulti **x1);                            // [y0,y1]=[min(x0),x1]
int icvec_umax_abs(rmulti *y0, rmulti *y1, int n, cmulti **x0, cmulti **x1);                        // value=max(abs(x))
int icvec_dmin_abs(rmulti *y0, rmulti *y1, int n, cmulti **x0, cmulti **x1);                        // value=min(abs(x))

int icvec_lintr(int m, int n, cmulti **y0, cmulti **y1, cmulti **A0, int LDA0, cmulti **A1, int LDA1, cmulti **x0, cmulti **x1);        // [y0,y1]=[A0,A1]*[x0,x1]
int icvec_add_lintr(int m, int n, cmulti **y0, cmulti **y1, cmulti **A0, int LDA0, cmulti **A1, int LDA1, cmulti **x0, cmulti **x1);    // [y0,y1]+=[A0,A1]*[x0,x1]
int icvec_sub_lintr(int m, int n, cmulti **y0, cmulti **y1, cmulti **A0, int LDA0, cmulti **A1, int LDA1, cmulti **x0, cmulti **x1);    // [y0,y1]-=[A0,A1]*[x0,x1]
int icvec_lintr_t(int m, int n, cmulti **y0, cmulti **y1, cmulti **A0, int LDA0, cmulti **A1, int LDA1, cmulti **x0, cmulti **x1);      // [y0,y1]=[A0,A1]'*[x0,x1]
int icvec_add_lintr_t(int m, int n, cmulti **y0, cmulti **y1, cmulti **A0, int LDA0, cmulti **A1, int LDA1, cmulti **x0, cmulti **x1);  // [y0,y1]+=[A0,A1]'*[x0,x1]
int icvec_sub_lintr_t(int m, int n, cmulti **y0, cmulti **y1, cmulti **A0, int LDA0, cmulti **A1, int LDA1, cmulti **x0, cmulti **x1);  // [y0,y1]-=[A0,A1]'*[x0,x1]
int icvec_lintr_ct(int m, int n, cmulti **y0, cmulti **y1, cmulti **A0, int LDA0, cmulti **A1, int LDA1, cmulti **x0, cmulti **x1);     // [y0,y1]=[A0,A1]'*[x0,x1]
int icvec_add_lintr_ct(int m, int n, cmulti **y0, cmulti **y1, cmulti **A0, int LDA0, cmulti **A1, int LDA1, cmulti **x0, cmulti **x1); // [y0,y1]+=[A0,A1]'*[x0,x1]
int icvec_sub_lintr_ct(int m, int n, cmulti **y0, cmulti **y1, cmulti **A0, int LDA0, cmulti **A1, int LDA1, cmulti **x0, cmulti **x1); // [y0,y1]-=[A0,A1]'*[x0,x1]
int icvec_func(cmulti *y0, cmulti *y1, func_t *f, int n, cmulti **x0, cmulti **x1);                // [y0,y1]=f(x)
int icvec_func_list(int m, cmulti **y0, cmulti **y1, func_t *f, int n, cmulti **x0, cmulti **x1);  // [y0,y1]=f(x)

#endif
