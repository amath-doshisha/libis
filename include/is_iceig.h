#ifndef ISYS_ICEIG_H
#define ISYS_ICEIG_H

// [F]=[A*X-lambda*X]
void iceig_residual(int n, cmulti **F0, cmulti **F1, cmulti **A0, int LDA0, cmulti **A1, int LDA1, cmulti **X0, cmulti **X1, cmulti *lambda0, cmulti *lambda1);

// Input: x,lambda,A,n
// Output: e
// Return: 0 if success, 1 if fail
int iceig_1pair_krawczyk(int n, cmulti **e, cmulti **A, int LDA, cmulti **x, cmulti *lambda, int debug);
int iceig_krawczyk(int n, int k, cmulti **E_vec, int LDE, cmulti **E_val, cmulti **A, int LDA, cmulti **X, int LDX, cmulti **lambda, int debug);
#endif
