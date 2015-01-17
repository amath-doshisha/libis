#include "mex.h"
#include <isys.h>
#include "mx2cmat.c"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  int n=1,m=1,i,j;
  long prec=53;
  cmulti **A=NULL;
  dcomplex z;
    
  (void) prhs;
    
  /* Check for proper number of input and  output arguments */
  if(nrhs>=1){
    if(mxIsStruct(prhs[0])){}else{ mexErrMsgIdAndTxt("MATLAB:cmulti_double:maxrhs","The arg1 should be Struct."); }
  }
  if(nrhs>1){ mexErrMsgIdAndTxt("MATLAB:cmulti_double:maxrhs","Too many input arguments."); }
  if(nlhs>1){ mexErrMsgIdAndTxt("MATLAB:cmulti_double:maxlhs","Too many output arguments."); }
    
  // allocate
  m=mxGetM(prhs[0]);
  n=mxGetN(prhs[0]);
  A=cmat_allocate_prec(m,n,prec);
  cmat_set_zeros(m,n,A,m);

  // get cmat from struct
  mx2cmat(m,n,A,m,prhs[0]);
  cmat_round(m,n,A,m,53);

  // convert double
  plhs[0]=mxCreateDoubleMatrix(m,n,mxCOMPLEX);
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      z=cget_z(MAT(A,i,j,m));
      mxGetPr(plhs[0])[j*m+i]=Z_R(z);
      mxGetPi(plhs[0])[j*m+i]=Z_I(z);
    }
  }
       
  // free
  A=cmat_free(m,n,A);
}
