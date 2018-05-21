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
 * setting
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
void zvec_set_dvec_dvec(int n, dcomplex *y, double *xr, double *xi);       // dcomplex <- (double,double)
void zvec_set_dvec(int n, dcomplex *y, double *x);                         // dcomplex <- double
void zvec_get_dvec(int n, double *y, dcomplex *x);                         // dcomplex -> double
void dvec_set_zvec(int n, double *y, dcomplex *x);                         // dcomplex -> double
void zvec_set_ivec(int n, dcomplex *y, int *x);                            // dcomplex <- int
void zvec_get_ivec(int n, int *y, dcomplex *x);                            // dcomplex -> int
void ivec_set_zvec(int n, int *y, dcomplex *x);                            // dcomplex -> int
void zvec_get_svec(int n, char **y, dcomplex *x, char format, int digits); // dcomplex -> char
void zvec_set_svec(int n, dcomplex *y, char **x);                          // dcomplex <- char

/*
 * convert itself
 */
void zvec_swap(int n, dcomplex *x, dcomplex *y);  // x <=> y
void zvec_reverse(int n, dcomplex *x);            // x[0]=x[n-1]; x[1]=x[n-2]; x[2]=x[n-3]; ...
void zvec_swap_at(dcomplex *x, int i, int j);     // x[i] <=> x[j]
void zvec_swap_index(int n, dcomplex *x, int *I); // x[i] <=> x[I[i]]
void zvec_sort(int n, dcomplex *x, int *I);       // sort x as x[0] <= x[1] <= x[2] <= ... <= x[n-1]. If I==NULL, then I is not ussed. If I!=NULL, then I is stored with sorted indexes.
void zvec_quick_sort(int n, dcomplex *x, int *I, int left, int right); // Don't call this function directly!
void zvec_sort_index(int *I, int n, dcomplex *X); // store the list of indexes, x is not destroyed

/*
 * y=f(x)
 */
void zvec_copy_zvec(int n, dcomplex *y, dcomplex *x);                // y=x
void zvec_copy_zvec_index(int n, dcomplex *y, dcomplex *x, int *I);  // y=x(I)
void zvec_index_copy_zvec(int n, dcomplex *y, int *I, dcomplex *x);  // y(I)=x
void dvec_real_zvec(int n, double *y, dcomplex *x);                  // y=real(x)
void dvec_imag_zvec(int n, double *y, dcomplex *x);                  // y=imag(x)
void dvec_abs_zvec(int n, double *y, dcomplex *x);                   // y=abs(x)
void dvec_arg_zvec(int n, double *y, dcomplex *x);                   // y=arg(x)
void zvec_conj_zvec(int n, dcomplex *y, dcomplex *x);                // y=conj(x)
void zvec_absc_zvec(int n, dcomplex *y, dcomplex *x);                // y=abs(real(x))+i*abs(imag(x))
void zvec_neg_zvec(int n, dcomplex *y, dcomplex *x);                 // y=-x
void zvec_pow2_zvec(int n, dcomplex *y, dcomplex *x);                // y=x^2
void dvec_log2_abs_zvec(int n, double *y, dcomplex *x);              // y=log2(abs(x))
void zvec_normalize_zvec(int n, dcomplex *y, dcomplex *x);           // y=x/sqrt(x'*x)
void zvec_normalize_sgn_zvec(int n, dcomplex *y, dcomplex *x);       // y=x/sqrt(x'*x) where x[i]/|x[i]|=1, abs(x[i]) is max
dcomplex zsum_zvec(int n, dcomplex *x);                              // y=sum(x)
dcomplex zaverage_zvec(int n, dcomplex *x);                          // y=sum(x)/n
double dsum_abs_zvec(int n, dcomplex *x);                            // y=sum(abs(x))
double dsum_pow2_abs_zvec(int n, dcomplex *x);                       // y=sum(abs(x)^2)
double dnorm1_zvec(int n, dcomplex *x);                              // y=sum(abs(x))
double dnorm2_zvec(int n, dcomplex *x);                              // y=sqrt(sum(abs(x)^2))
double dnorm_max_zvec(int n, dcomplex *x);                           // y=max(abs(x))
dcomplex zmax_zvec(int n, dcomplex *x);                              // y=max(x)
dcomplex zmin_zvec(int n, dcomplex *x);                              // y=min(x)
double dmax_abs_zvec(int n, dcomplex *x);                            // y=max(abs(x))
double dmin_abs_zvec(int n, dcomplex *x);                            // y=min(abs(x))
double dmax_abs_zvec_index(int n, dcomplex *x, int *I);              // y=max(abs(x)), (*I)=k where k=max{ k | abs(x[k]) } unless I==NULL
double dmin_abs_zvec_index(int n, dcomplex *x, int *I);              // y=min(abs(x)), (*I)=k where k=min{ k | abs(x[k]) } unless I==NULL
double dmax_max_absc_zvec(int n, dcomplex *x);                       // y=max(real(absc(x)),imag(absc(x)))

/*
 * z=f(x,y)
 */
// z=x+y
void zvec_add_zvec_zvec(int n, dcomplex *z, dcomplex *x, dcomplex *y);
void zvec_add_zvec_dvec(int n, dcomplex *z, dcomplex *x, double *y);
void zvec_add_dvec_zvec(int n, dcomplex *z, double *x, dcomplex *y);
void zvec_add_zvec_zscalar(int n, dcomplex *z, dcomplex *x, dcomplex y);
void zvec_add_zvec_dscalar(int n, dcomplex *z, dcomplex *x, double y);
void zvec_add_dvec_zscalar(int n, dcomplex *z, double *x, dcomplex y);
void zvec_add_zscalar_zvec(int n, dcomplex *z, dcomplex x, dcomplex *y);
void zvec_add_zscalar_dvec(int n, dcomplex *z, dcomplex x, double *y);
void zvec_add_dscalar_zvec(int n, dcomplex *z, double x, dcomplex *y);
// z=x-y
void zvec_sub_zvec_zvec(int n, dcomplex *z, dcomplex *x, dcomplex *y);
void zvec_sub_zvec_dvec(int n, dcomplex *z, dcomplex *x, double *y);
void zvec_sub_dvec_zvec(int n, dcomplex *z, double *x, dcomplex *y);
void zvec_sub_zvec_zscalar(int n, dcomplex *z, dcomplex *x, dcomplex y);
void zvec_sub_zvec_dscalar(int n, dcomplex *z, dcomplex *x, double y);
void zvec_sub_dvec_zscalar(int n, dcomplex *z, double *x, dcomplex y);
void zvec_sub_zscalar_zvec(int n, dcomplex *z, dcomplex x, dcomplex *y);
void zvec_sub_zscalar_dvec(int n, dcomplex *z, dcomplex x, double *y);
void zvec_sub_dscalar_zvec(int n, dcomplex *z, double x, dcomplex *y);
// z=x*y
void zvec_mul_zvec_zvec(int n, dcomplex *z, dcomplex *x, dcomplex *y);
void zvec_mul_zvec_dvec(int n, dcomplex *z, dcomplex *x, double *y);
void zvec_mul_dvec_zvec(int n, dcomplex *z, double *x, dcomplex *y);
void zvec_mul_zvec_zscalar(int n, dcomplex *z, dcomplex *x, dcomplex y);
void zvec_mul_zvec_dscalar(int n, dcomplex *z, dcomplex *x, double y);
void zvec_mul_dvec_zscalar(int n, dcomplex *z, double *x, dcomplex y);
void zvec_mul_zscalar_zvec(int n, dcomplex *z, dcomplex x, dcomplex *y);
void zvec_mul_zscalar_dvec(int n, dcomplex *z, dcomplex x, double *y);
void zvec_mul_dscalar_zvec(int n, dcomplex *z, double x, dcomplex *y);
// z=x/y
void zvec_div_zvec_zvec(int n, dcomplex *z, dcomplex *x, dcomplex *y);
void zvec_div_zvec_dvec(int n, dcomplex *z, dcomplex *x, double *y);
void zvec_div_dvec_zvec(int n, dcomplex *z, double *x, dcomplex *y);
void zvec_div_zvec_zscalar(int n, dcomplex *z, dcomplex *x, dcomplex y);
void zvec_div_zvec_dscalar(int n, dcomplex *z, dcomplex *x, double y);
void zvec_div_dvec_zscalar(int n, dcomplex *z, double *x, dcomplex y);
void zvec_div_zscalar_zvec(int n, dcomplex *z, dcomplex x, dcomplex *y);
void zvec_div_zscalar_dvec(int n, dcomplex *z, dcomplex x, double *y);
void zvec_div_dscalar_zvec(int n, dcomplex *z, double x, dcomplex *y);
// z=x^y
void zvec_pow_zvec_zvec(int n, dcomplex *z, dcomplex *x, dcomplex *y);
void zvec_pow_zvec_dvec(int n, dcomplex *z, dcomplex *x, double *y);
void zvec_pow_zvec_zscalar(int n, dcomplex *z, dcomplex *x, dcomplex y);
void zvec_pow_zvec_dscalar(int n, dcomplex *z, dcomplex *x, double y);
void zvec_pow_zvec_iscalar(int n, dcomplex *z, dcomplex *x, int y);
// z=abs(x-y)
void dvec_abs_sub_zvec_zvec(int n, double *z, dcomplex *x, dcomplex *y);
void dvec_abs_sub_zvec_zscalar(int n, double *z, dcomplex *x, dcomplex y);
// z=sum(abs(x-y))
double dsum_abs_sub_zvec_zvec(int n, dcomplex *x, dcomplex *y);
// z=max(|x-y|)
double dmax_abs_sub_zvec_zvec(int n, dcomplex *x, dcomplex *y);
// z=max(|x-y|^2)
double dmax_pow2_abs_sub_zvec_zvec(int n, dcomplex *x, dcomplex *y);
// z=x'*y
dcomplex zinnprod_zvec_zvec(int n, dcomplex *x, dcomplex *y);
// z=(x'*y)/\sqrt(x'*x)/sqrt(y'*y)
dcomplex zdcos_zvec_zvec(int n, dcomplex *x, dcomplex *y);
// z=abs(x'*y)/\sqrt(x'*x)/sqrt(y'*y)
double ddcos_abs_zvec_zvec(int n, dcomplex *x, dcomplex *y);
// z=acos(abs(x'*y)/\sqrt(x'*x)/sqrt(y'*y))
double dangle_zvec_zvec(int n, dcomplex *x, dcomplex *y);
// z=(180/PI)*acos(abs(x'*y)/\sqrt(x'*x)/sqrt(y'*y))
double ddeg_angle_zvec_zvec(int n, dcomplex *x, dcomplex *y);
// z=max(abs(x)/abs(y))
double dmax_div_abs_zvec_zvec(int n, dcomplex *x, dcomplex *y);

/*
 * y=y+f(x)
 */
void zvec_add_scaled(int n, dcomplex *y, dcomplex a, dcomplex *x);                    // y=y+a*x
void zvec_sub_scaled(int n, dcomplex *y, dcomplex a, dcomplex *x);                    // y=y-a*x
void zvec_orthogonalize(int n, dcomplex *y, dcomplex *x);                             // y=y-(x'*y)*x where x'*x=1



/*
 * 作業中：dmat.cに移動する
 */
void zvec_mul_zmat_zvec(int m, int n, dcomplex *y, dcomplex *A, int LDA, dcomplex *x);        // y=A*x
void zvec_mul_zmat_t_zvec(int m, int n, dcomplex *y, dcomplex *A, int LDA, dcomplex *x);      // y=A^T*x
void zvec_mul_zmat_ct_zvec(int m, int n, dcomplex *y, dcomplex *A, int LDA, dcomplex *x);     // y=A'*x


void zvec_add_mul_zmat_zvec(int m, int n, dcomplex *y, dcomplex *A, int LDA, dcomplex *x);    // y=y+A*x
void zvec_add_mul_zmat_t_zvec(int m, int n, dcomplex *y, dcomplex *A, int LDA, dcomplex *x);  // y=y+A^T*x
void zvec_add_mul_zmat_ct_zvec(int m, int n, dcomplex *y, dcomplex *A, int LDA, dcomplex *x); // y=y+A'*x

void zvec_sub_zmat_zvec(int m, int n, dcomplex *y, dcomplex *A, int LDA, dcomplex *x);    // y=y-A*x
void zvec_sub_zmat_t_zvec(int m, int n, dcomplex *y, dcomplex *A, int LDA, dcomplex *x);  // y=y-A^T*x
void zvec_sub_zmat_ct_zvec(int m, int n, dcomplex *y, dcomplex *A, int LDA, dcomplex *x); // y=y-A'*x


#endif
