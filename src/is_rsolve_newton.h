#ifndef ISYS_RSOLVE_NEWTON_H
#define ISYS_RSOLVE_NEWTON_H

#include<is_rmulti.h>

enum { RSOLVE_NEWTON_SUCCESS=0,
       RSOLVE_NEWTON_FAILED=1,
       RSOLVE_NEWTON_FAILED_MORESTEP,
       RSOLVE_NEWTON_FAILED_STEPMAX
};

void rsolve_func_residual(int m, rmulti *norm, rmulti **F, rmulti **x, func_t *fF);
int rsolve_newton(int m, rmulti **x, func_t *fF, int step_max, int debug);
int rsolve_krawczyk(int m, rmulti **e, rmulti **x, func_t *fF, int debug);

#endif
