#ifndef IS_PMAT_H
#define IS_PMAT_H

/*
 * allocation
 */
poly* pmat_allocate(int m, int n, int p_m, int p_n);
poly* pmat_allocate_as_jacobian(int n, poly *F);
poly* pmat_free(int LDA, int n, poly *A);

/*
 * input and output
 */
void pmat_print(int m, int n, poly *A, int LDA, char *name);
void pmat_put(int m, int n, poly *A, int LDA, char *name);

#endif
