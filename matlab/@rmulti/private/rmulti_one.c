#include "mex.h"
#include <string.h>
#include <isys.h>
#include"rmat2mx.c"
#include"mx2rmat.c"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  int mX=-1,nX=-1,mY=-1,nY=-1,j;
  long prec=53;
  rmulti **X=NULL,**Y=NULL;
  char *cmd=NULL;
    
  (void) prhs;
    
  /* Check for proper number of input and  output arguments */
  if(nrhs>=1){
    if(mxIsChar(prhs[0])){ cmd=mxArrayToString(prhs[0]); }
    else{ mexErrMsgIdAndTxt("MATLAB:rmulti_one","The arg1 should be Char."); }
  }
  if(nrhs>=2){
    if(mxIsStruct(prhs[1])){}
    else{ mexErrMsgIdAndTxt("MATLAB:rmulti_one","The arg2 should be Struct."); }
  }
  if(nrhs>2){ mexErrMsgIdAndTxt("MATLAB:rmulti_one:maxrhs","Too many input arguments."); }
  if(nlhs>1){ mexErrMsgIdAndTxt("MATLAB:rmulti_one:maxlhs","Too many output arguments."); }

  // allocate
  mX=mxGetM(prhs[1]);
  nX=mxGetN(prhs[1]);
  X=rmat_allocate_prec(mX,nX,prec);
  mx2rmat(mX,nX,X,mX,prhs[1]);
  prec=rmat_get_prec_max(mX,nX,X,mX);

  // command
  if(char_eq(cmd,"abs")){
    // y=abs(x)
    mY=mX;
    nY=nX;
    Y=rmat_allocate_prec(mY,nY,prec);
    rmat_abs(mY,nY,Y,mY,X,mX);
  }else if(char_eq(cmd,"max")){
    // y=max(x)
    if(mX==1){
      mY=1; nY=1;
      Y=rmat_allocate_prec(mY,nY,prec);
      rvec_max(Y[0],nX,X);
    }else{
      mY=1;
      nY=nX;
      Y=rmat_allocate_prec(mY,nY,prec);
      for(j=0; j<nX; j++){ rvec_max(Y[j],mX,&COL(X,j,mX)); }
    }
  }else{ mexErrMsgIdAndTxt("MATLAB:rmulti_one","The command is illeagal."); }

  // done
  plhs[0]=rmat2mx(mY,nY,Y,mY);
       
  // free
  X=rmat_free(mX,nX,X);
  Y=rmat_free(mY,nY,Y);
}
