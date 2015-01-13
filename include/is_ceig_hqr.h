#ifndef IS_CEIG_HQR_H
#define IS_CEIG_HQR_H

#include<is_rmulti.h>
#include<is_cmulti.h>

// QR method for computing eigenvalues of the Hessenberg-type matrix
void ceig_hqr(int n, cmulti **lambda, cmulti **A, int LDA, int debug);

#endif
