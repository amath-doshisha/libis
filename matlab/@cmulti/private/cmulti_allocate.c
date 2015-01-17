#include "mex.h"
#include <string.h>
#include <isys.h>
#include"../../@rmulti/private/mxGetDataLong.c"
#include"cmat2mx.c"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int n=1,m=1;
    long prec=53;
    cmulti **A=NULL;
    char *cmd=NULL;
    
    (void) prhs;
    
    /* Check for proper number of input and  output arguments */
    if(nrhs<1){ mexErrMsgIdAndTxt("MATLAB:cmulti_allocate:maxrhs","Too few input arguments."); }
    // get command
    if(nrhs>=1){
      if(mxIsChar(prhs[0])){ cmd=mxArrayToString(prhs[0]); }
      else{ mexErrMsgIdAndTxt("MATLAB:cmulti_allocate","The arg1 should be Char."); }
    }
    // get precision
    if(nrhs>=2){
      if(mxIsNumeric(prhs[1]) && mxGetM(prhs[1])==1 && mxGetN(prhs[1])==1){ prec=mxGetDataLong(prhs[1]); }
      else{ mexErrMsgIdAndTxt("MATLAB:cmulti_allocate","The arg2 should be Integer(1x1)."); }
    }
    // get options
    if(char_eq(cmd,"set")){
      if(nrhs>3){ mexErrMsgIdAndTxt("MATLAB:cmulti_allocate:maxrhs","Too many input arguments."); }
      if(mxIsDouble(prhs[2])){ m=mxGetM(prhs[2]); n=mxGetN(prhs[2]); }
      else{ mexErrMsgIdAndTxt("MATLAB:cmulti_allocate","The arg3 should be Double(m,n)."); }
    }else{
      if(nrhs>4){ mexErrMsgIdAndTxt("MATLAB:cmulti_allocate:maxrhs","Too many input arguments."); }
      if(nrhs>=3){
	if(mxIsNumeric(prhs[2]) && mxGetM(prhs[2])==1 && mxGetN(prhs[2])==1){ m=mxGetDataLong(prhs[2]); }
	else{ mexErrMsgIdAndTxt("MATLAB:cmulti_allocate","The arg3 should be Integer(1x1)."); }
      }
      if(nrhs>=4){
	if(mxIsNumeric(prhs[3]) && mxGetM(prhs[3])==1 && mxGetN(prhs[3])==1){ n=mxGetDataLong(prhs[3]); }
	else{ mexErrMsgIdAndTxt("MATLAB:cmulti_allocate","The arg4 should be Integer(1x1)."); }
      }
    }
    if(nlhs>1){ mexErrMsgIdAndTxt("MATLAB:cmulti_allocate:maxlhs","Too many output arguments."); }
    
    // allocate
    A=cmat_allocate_prec(m,n,prec);
    cmat_set_zeros(m,n,A,m);
    if(char_eq(cmd,"zeros")){
      cmat_set_zeros(m,n,A,m);
    }else if(char_eq(cmd,"ones")){
      cmat_set_ones(m,n,A,m);
    }else if(char_eq(cmd,"eye")){
      cmat_set_eye(m,n,A,m);
    }else if(char_eq(cmd,"rand")){
      cmat_set_rand(m,n,A,m,1,0);
    }else if(char_eq(cmd,"set")) {
      if(mxIsComplex(prhs[2])){
	cmat_set_dd(m,n,A,m,mxGetPr(prhs[2]),m,mxGetPi(prhs[2]),m);
      }else{
	cmat_set_d(m,n,A,m,mxGetPr(prhs[2]),m);
      }
    }else{
      cmat_set_zeros(m,n,A,m);
    }
    plhs[0]=cmat2mx(m,n,A,m);
       
    // free
    A=cmat_free(m,n,A);
}
