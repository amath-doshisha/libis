#ifndef IS_ZVEC_H
#define IS_ZVEC_H

#include<is_dcomplex.h>

/**
 @file  is_zvec.h
 @brief 倍精度複素数型dcomplexのベクトルに関する関数の宣言
 */

/*
 * allocation
 */
dcomplex *zvec_allocate(int n);
dcomplex *zvec_allocate_clone_zvec(int n, dcomplex *y);
dcomplex *zvec_allocate_clone_dvec(int n, double *y);
dcomplex* zvec_free(dcomplex *x);

/*
 * I/O
 */
void zvec_print(int n, dcomplex *x, char *name, char format, int digits);
void zvec_save(int n, dcomplex *x, int offset, char* fmt, ...);
void zvec_save_cplane(int n, dcomplex *x, char* fmt, ...);
void zvec_save_cplane_label(int n, dcomplex *x, int offset, int *label, char* fmt, ...);
void zvec_bin_save(int n, dcomplex *x, char* fmt, ...);
dcomplex *zvec_bin_load(int *n, char* fmt, ...);

/*
 * initialization
 */
void zvec_set_nan(int n, dcomplex *x);                            // x=nan(n,1)
void zvec_set_inf(int n, dcomplex *x, int rsgn, int isgn);        // x=sgn*inf(n,1)
void zvec_set_zeros(int n, dcomplex *x);                          // x=zeros(n,1)
void zvec_set_ones(int n, dcomplex *x);                           // x=ones(n,1)
void zvec_set_all_z(int n, dcomplex *x, dcomplex a);              // x=ones(n,1)*a
void zvec_set_all_dd(int n, dcomplex *x, double a_r, double a_i); // x=ones(n,1)*a
void zvec_set_all_d(int n, dcomplex *x, double a);                // x=ones(n,1)*a
void zvec_set_unit(int n, dcomplex *x, int k);                    // x=zeros(n,1); x[k]=1
void zvec_set_grid(int n, dcomplex *x);                           // x=0:(n-1)
void zvec_set_rand(int n, dcomplex *x, double a, double b);       // x=rand(n,1)*a+b

/*
 * casting
 */
void zvec_set_zvec(int n, dcomplex *y, dcomplex *x);                       // dcomplex <- dcomplex
void zvec_set_dvec_dvec(int n, dcomplex *y, double *x_r, double *x_i);     // dcomplex <- (double,double)
void zvec_set_dvec(int n, dcomplex *y, double *x);                         // dcomplex <- double
void zvec_get_dvec(int n, double *y, dcomplex *x);                         // dcomplex -> double
void dvec_set_zvec(int n, double *y, dcomplex *x);                         // dcomplex -> double
void zvec_set_ivec(int n, dcomplex *y, int *x);                            // dcomplex <- int
void zvec_get_ivec(int n, int *y, dcomplex *x);                            // dcomplex -> int
void ivec_set_zvec(int n, int *y, dcomplex *x);                            // dcomplex -> int
void zvec_get_svec(int n, char **y, dcomplex *x, char format, int digits); // dcomplex -> char
void zvec_set_svec(int n, dcomplex *y, char **x);                          // dcomplex <- char

/*
 * convert its elements
 */
void zvec_swap(int n, dcomplex *x, dcomplex *y);  // x <=> y
void zvec_reverse(int n, dcomplex *x);            // x[0]=x[n-1]; x[1]=x[n-2]; x[2]=x[n-3]; ...
void zvec_swap_at(dcomplex *x, int i, int j);     // x[i] <=> x[j]
void zvec_swap_index(int n, dcomplex *x, int *I); // x[i] <=> x[I[i]]
// sort x as x[0] <= x[1] <= x[2] <= ... <= x[n-1]
// if I==NULL, then I is not ussed.
// if I!=NULL, then I is stored with sorted indexes
void zvec_sort(int n, dcomplex *x, int *I);
// Don't call this function directly!
void zvec_quick_sort(int n, dcomplex *x, int *I, int left, int right);
// store the list of indexes, x is not destroyed
void zvec_sort_index(int *I, int n, dcomplex *X);

/*
 * y=f(x)
 */
void zvec_copy(int n, dcomplex *y, dcomplex *x);                // y=x
void zvec_copy_index(int n, dcomplex *Y, dcomplex *X, int *I);  // Y[i]=X[I[i]], 0<=i<n
void zvec_neg(int n, dcomplex *y, dcomplex *x);                 // y=-x
void zvec_real(int n, double *y, dcomplex *x);                  // y=real(x)
void zvec_conj(int n, dcomplex *y, dcomplex *x);                // y=conj(x)
double zvec_norm1(int n, dcomplex *x);                          // y=sum(abs(x))
double zvec_norm2(int n, dcomplex *x);                          // y=sqrt(sum(abs(x).^2))
dcomplex zvec_max(int n, dcomplex *x);                          // y=max(x)
dcomplex zvec_min(int n, dcomplex *x);                          // y=min(x)
dcomplex zvec_sum(int n, dcomplex *x);                          // y=sum(x)
dcomplex zvec_average(int n, dcomplex *x);                      // y=sum(x)/n
double zvec_norm_max(int n, dcomplex *x);                       // y=max(abs(x))
double zvec_sum_pow2_abs(int n, dcomplex *x);                   // y=sum(abs(x).^2)
double zvec_max_abs(int n, dcomplex *x);                        // y=max(abs(x))
double zvec_min_abs(int n, dcomplex *x);                        // y=min(abs(x))
double zvec_max_abs_index(int n, dcomplex *x, int *I);          // y=max(abs(x)), (*I)=k where k=max{ k | abs(x[k]) } unless I==NULL
double zvec_min_abs_index(int n, dcomplex *x, int *I);          // y=min(abs(x)), (*I)=k where k=min{ k | abs(x[k]) } unless I==NULL

/*
 * z=f(x,y)
 */
// z=x+y
void zvec_add_zvec(int n, dcomplex *z, dcomplex *x, dcomplex *y);
void zvec_add_dvec(int n, dcomplex *z, dcomplex *x, double *y);
void dvec_add_zvec(int n, dcomplex *z, double *x, dcomplex *y);
void zvec_add_zscalar(int n, dcomplex *z, dcomplex *x, dcomplex y);
void zvec_add_dscalar(int n, dcomplex *z, dcomplex *x, double y);
void dvec_add_zscalar(int n, dcomplex *z, double *x, dcomplex y);
void zscalar_add_zvec(int n, dcomplex *z, dcomplex x, dcomplex *y);
void zscalar_add_dvec(int n, dcomplex *z, dcomplex x, double *y);
void dscalar_add_zvec(int n, dcomplex *z, double x, dcomplex *y);
// z=x-y
void zvec_sub_zvec(int n, dcomplex *z, dcomplex *x, dcomplex *y);
void zvec_sub_dvec(int n, dcomplex *z, dcomplex *x, double *y);
void dvec_sub_zvec(int n, dcomplex *z, double *x, dcomplex *y);
void zvec_sub_zscalar(int n, dcomplex *z, dcomplex *x, dcomplex y);
void zvec_sub_dscalar(int n, dcomplex *z, dcomplex *x, double y);
void dvec_sub_zscalar(int n, dcomplex *z, double *x, dcomplex y);
void zscalar_sub_zvec(int n, dcomplex *z, dcomplex x, dcomplex *y);
void zscalar_sub_dvec(int n, dcomplex *z, dcomplex x, double *y);
void dscalar_sub_zvec(int n, dcomplex *z, double x, dcomplex *y);
// z=x*y
void zvec_mul_zvec(int n, dcomplex *z, dcomplex *x, dcomplex *y);
void zvec_mul_dvec(int n, dcomplex *z, dcomplex *x, double *y);
void dvec_mul_zvec(int n, dcomplex *z, double *x, dcomplex *y);
void zvec_mul_zscalar(int n, dcomplex *z, dcomplex *x, dcomplex y);
void zvec_mul_dscalar(int n, dcomplex *z, dcomplex *x, double y);
void dvec_mul_zscalar(int n, dcomplex *z, double *x, dcomplex y);
void zscalar_mul_zvec(int n, dcomplex *z, dcomplex x, dcomplex *y);
void zscalar_mul_dvec(int n, dcomplex *z, dcomplex x, double *y);
void dscalar_mul_zvec(int n, dcomplex *z, double x, dcomplex *y);
// z=x/y
void zvec_div_zvec(int n, dcomplex *z, dcomplex *x, dcomplex *y);
void zvec_div_dvec(int n, dcomplex *z, dcomplex *x, double *y);
void dvec_div_zvec(int n, dcomplex *z, double *x, dcomplex *y);
void zvec_div_zscalar(int n, dcomplex *z, dcomplex *x, dcomplex y);
void zvec_div_dscalar(int n, dcomplex *z, dcomplex *x, double y);
void dvec_div_zscalar(int n, dcomplex *z, double *x, dcomplex y);
void zscalar_div_zvec(int n, dcomplex *z, dcomplex x, dcomplex *y);
void zscalar_div_dvec(int n, dcomplex *z, dcomplex x, double *y);
void dscalar_div_zvec(int n, dcomplex *z, double x, dcomplex *y);
// others
dcomplex zvec_dot(int n, dcomplex *x, dcomplex *y);                  // z=x'*y
double zvec_dcos_abs(int n, dcomplex *x, dcomplex *y);               // z=abs(x'*y)/\sqrt(x'*x)/sqrt(y'*y)
double zvec_angle(int n, dcomplex *x, dcomplex *y);                  // z=acos(abs(x'*y)/\sqrt(x'*x)/sqrt(y'*y))
double zvec_angle_deg(int n, dcomplex *x, dcomplex *y);              // z=(180/PI)*acos(abs(x'*y)/\sqrt(x'*x)/sqrt(y'*y))
double zvec_dist_norm1(int n, dcomplex *x, dcomplex *y);             // z=sum(abs(x-y))
double zvec_dist_norm_max(int n, dcomplex *x, dcomplex *y);          // z=max(|x-y|)
void zvec_abs_sub(int n, double *z, dcomplex *x, dcomplex *y);       // z=abs(x-y)
void zvec_abs_sub_scalar(int n, double *z, dcomplex *x, dcomplex y); // z=abs(x-y)
void zvec_normalize(int n, dcomplex *y, dcomplex *x);                // y=x/sqrt(x'*x)
void zvec_normalize_sgn(int n, dcomplex *y, dcomplex *x);            // y=x/sqrt(x'*x) where x[i]/|x[i]|=1, abs(x[i]) is max

/*
 * y=y+f(x)
 */
void zvec_add_scaled(int n, dcomplex *y, dcomplex a, dcomplex *x);                    // y=y+a*x
void zvec_sub_scaled(int n, dcomplex *y, dcomplex a, dcomplex *x);                    // y=y-a*x
void zvec_orthogonalize(int n, dcomplex *y, dcomplex *x);                             // y=y-(x'*y)*x where x'*x=1



/*
 * 作業中：dmat.cに移動する
 */
void zvec_lintr(int m, int n, dcomplex *y, dcomplex *A, int LDA, dcomplex *x);        // y=A*x
void zvec_add_lintr(int m, int n, dcomplex *y, dcomplex *A, int LDA, dcomplex *x);    // y=y+A*x
void zvec_sub_lintr(int m, int n, dcomplex *y, dcomplex *A, int LDA, dcomplex *x);    // y=y-A*x
void zvec_lintr_t(int m, int n, dcomplex *y, dcomplex *A, int LDA, dcomplex *x);      // y=A^T*x
void zvec_add_lintr_t(int m, int n, dcomplex *y, dcomplex *A, int LDA, dcomplex *x);  // y=y+A^T*x
void zvec_sub_lintr_t(int m, int n, dcomplex *y, dcomplex *A, int LDA, dcomplex *x);  // y=y-A^T*x
void zvec_lintr_ct(int m, int n, dcomplex *y, dcomplex *A, int LDA, dcomplex *x);     // y=A'*x
void zvec_add_lintr_ct(int m, int n, dcomplex *y, dcomplex *A, int LDA, dcomplex *x); // y=y+A'*x
void zvec_sub_lintr_ct(int m, int n, dcomplex *y, dcomplex *A, int LDA, dcomplex *x); // y=y-A'*x


#endif
