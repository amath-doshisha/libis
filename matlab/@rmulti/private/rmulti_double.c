#include "mex.h"
#include <isys.h>
#include "mx2rmat.c"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  int n=1,m=1,i,j;
  long prec=53;
  rmulti **A=NULL;
    
  (void) prhs;
    
  /* Check for proper number of input and  output arguments */
  if(nrhs>=1){
    if(mxIsStruct(prhs[0])){}else{ mexErrMsgIdAndTxt("MATLAB:rmulti_double:maxrhs","The arg1 should be Struct."); }
  }
  if(nrhs>1){ mexErrMsgIdAndTxt("MATLAB:rmulti_double:maxrhs","Too many input arguments."); }
  if(nlhs>1){ mexErrMsgIdAndTxt("MATLAB:rmulti_double:maxlhs","Too many output arguments."); }
    
  // allocate
  m=mxGetM(prhs[0]);
  n=mxGetN(prhs[0]);
  A=rmat_allocate_prec(m,n,prec);
  rmat_set_zeros(m,n,A,m);

  // get rmat from struct
  mx2rmat(m,n,A,m,prhs[0]);
  rmat_round(m,n,A,m,53);

  // convert double
  plhs[0]=mxCreateDoubleMatrix(m,n,mxREAL);
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      mxGetPr(plhs[0])[j*m+i]=rget_d(MAT(A,i,j,m));
    }
  }
       
  // free
  A=rmat_free(m,n,A);
}
