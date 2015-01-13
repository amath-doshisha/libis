#include "mex.h"
#include <string.h>
#include <isys.h>
#include"mxGetDataLong.c"
#include"rmat2mx.c"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int n=1,m=1;
    long prec=53;
    rmulti **A=NULL;
    char *cmd=NULL;
    
    (void) prhs;
    
    /* Check for proper number of input and  output arguments */
    if(nrhs<1){ mexErrMsgIdAndTxt("MATLAB:rmulti_allocate:maxrhs","Too few input arguments."); }
    if(nrhs>=1){
      if(mxIsChar(prhs[0])){ cmd=mxArrayToString(prhs[0]); }
      else{ mexErrMsgIdAndTxt("MATLAB:rmulti_allocate","The arg1 should be Char."); }
    }
    if(nrhs>=2){
      if(mxIsNumeric(prhs[1]) && mxGetM(prhs[1])==1 && mxGetN(prhs[1])==1){ prec=mxGetDataLong(prhs[1]); }
      else{ mexErrMsgIdAndTxt("MATLAB:rmulti_allocate","The arg2 should be Integer(1x1)."); }
    }
    if(char_eq(cmd,"set")){
      if(nrhs>3){ mexErrMsgIdAndTxt("MATLAB:rmulti_allocate:maxrhs","Too many input arguments."); }
      if(mxIsDouble(prhs[2])){ m=mxGetM(prhs[2]); n=mxGetN(prhs[2]); }
      else{ mexErrMsgIdAndTxt("MATLAB:rmulti_allocate","The arg3 should be Double(m,n)."); }
    }else{
      if(nrhs>4){ mexErrMsgIdAndTxt("MATLAB:rmulti_allocate:maxrhs","Too many input arguments."); }
      if(nrhs>=3){
	if(mxIsNumeric(prhs[2]) && mxGetM(prhs[2])==1 && mxGetN(prhs[2])==1){ m=mxGetDataLong(prhs[2]); }
	else{ mexErrMsgIdAndTxt("MATLAB:rmulti_allocate","The arg3 should be Integer(1x1)."); }
      }
      if(nrhs>=4){
	if(mxIsNumeric(prhs[3]) && mxGetM(prhs[3])==1 && mxGetN(prhs[3])==1){ n=mxGetDataLong(prhs[3]); }
	else{ mexErrMsgIdAndTxt("MATLAB:rmulti_allocate","The arg4 should be Integer(1x1)."); }
      }
    }
    if(nlhs>1){ mexErrMsgIdAndTxt("MATLAB:rmulti_allocate:maxlhs","Too many output arguments."); }
    
    // allocate
    A=rmat_allocate_prec(m,n,prec);
    rmat_set_zeros(m,n,A,m);
    if     (char_eq(cmd,"zeros")){ rmat_set_zeros(m,n,A,m); }
    else if(char_eq(cmd,"ones")) { rmat_set_ones(m,n,A,m); }
    else if(char_eq(cmd,"eye"))  { rmat_set_eye(m,n,A,m); }
    else if(char_eq(cmd,"rand")) { rmat_set_rand(m,n,A,m,1,0); }
    else if(char_eq(cmd,"set"))  { rmat_set_d(m,n,A,m,mxGetPr(prhs[2]),m); }
    else                         { rmat_set_zeros(m,n,A,m); }
    plhs[0]=rmat2mx(m,n,A,m);
       
    // free
    A=rmat_free(m,n,A);
}
