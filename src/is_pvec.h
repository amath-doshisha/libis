#ifndef IS_PVEC_H
#define IS_PVEC_H
#include<is_poly.h>

/*
 * allocation
 */
poly* pvec_allocate(int n, int p_m, int p_n);
poly* pvec_allocate_from_str(int *N_ret, char **str);
poly* pvec_allocate_from_file(int *N_ret, char* fmt, ...);
poly* pvec_free(int n, poly *x);

/*
 * input and output
 */
void pvec_print(int n, poly *x, char *name);
void pvec_put(int n, poly *x, char *name);

#endif
