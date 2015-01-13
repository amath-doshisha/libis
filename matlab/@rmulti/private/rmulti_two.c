#include "mex.h"
#include <string.h>
#include <isys.h>
#include"rmat2mx.c"
#include"mx2rmat.c"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  int mX=-1,nX=-1,mY=-1,nY=-1,mZ=-1,nZ=-1,info;
  long prec=53,prec0,prec1;
  rmulti **X=NULL,**Y=NULL,**Z=NULL;
  char *cmd=NULL;
    
  (void) prhs;
    
  /* Check for proper number of input and  output arguments */
  if(nrhs>=1){
    if(mxIsChar(prhs[0])){ cmd=mxArrayToString(prhs[0]); }
    else{ mexErrMsgIdAndTxt("MATLAB:rmulti_two","The arg1 should be Char."); }
  }
  if(nrhs>=2){
    if(mxIsStruct(prhs[1])){}
    else{ mexErrMsgIdAndTxt("MATLAB:rmulti_two","The arg2 should be Struct."); }
  }
  if(nrhs>=3){
    if(mxIsStruct(prhs[2])){}
    else{ mexErrMsgIdAndTxt("MATLAB:rmulti_two","The arg3 should be Struct."); }
  }
  if(nrhs>3){ mexErrMsgIdAndTxt("MATLAB:rmulti_two:maxrhs","Too many input arguments."); }
  if(nlhs>1){ mexErrMsgIdAndTxt("MATLAB:rmulti_two:maxlhs","Too many output arguments."); }

  // allocate
  mX=mxGetM(prhs[1]); nX=mxGetN(prhs[1]); X=rmat_allocate_prec(mX,nX,prec); mx2rmat(mX,nX,X,mX,prhs[1]); prec0=rmat_get_prec_max(mX,nX,X,mX);
  mY=mxGetM(prhs[2]); nY=mxGetN(prhs[2]); Y=rmat_allocate_prec(mY,nY,prec); mx2rmat(mY,nY,Y,mY,prhs[2]); prec1=rmat_get_prec_max(mY,nY,Y,mY);
  prec=MAX2(prec0,prec1);

  // command
  if(char_eq(cmd,"plus")){
    // x+y
    if(mY==1 && nY==1){
      mZ=mX;
      nZ=nX;
      Z=rmat_allocate_prec(mZ,nZ,prec);
      rmat_add_r(mZ,nZ,Z,mZ,X,mX,Y[0]);
    }else if(mX==1 && nX==1){
      mZ=mY;
      nZ=nY;
      Z=rmat_allocate_prec(mZ,nZ,prec);
      rmat_add_r(mZ,nZ,Z,mZ,Y,mY,X[0]);
    }else if(mX==mY && nX==nY){
      mZ=mX;
      nZ=nX;
      Z=rmat_allocate_prec(mZ,nZ,prec);
      rmat_add(mZ,nZ,Z,mZ,X,mX,Y,mY);
    }else{ mexErrMsgTxt("rmulti_two, Z=X+Y. The sizes of X,Y should be same."); }
  }else if(char_eq(cmd,"minus")){
    // x-y
    if(mY==1 && nY==1){
      mZ=mX;
      nZ=nX;
      Z=rmat_allocate_prec(mZ,nZ,prec);
      rmat_sub_r2(mZ,nZ,Z,mZ,X,mX,Y[0]);
    }else if(mX==1 && nX==1){
      mZ=mY;
      nZ=nY;
      Z=rmat_allocate_prec(mZ,nZ,prec);
      rmat_sub_r1(mZ,nZ,Z,mZ,X[0],Y,mY);
    }else if(mX==mY && nX==nY){
      mZ=mX;
      nZ=nX;
      Z=rmat_allocate_prec(mZ,nZ,prec);
      rmat_sub(mZ,nZ,Z,mZ,X,mX,Y,mY);
    }else{ mexErrMsgTxt("rmulti_two, Z=X-Y. The sizes of X,Y should be same."); }
  }else if(char_eq(cmd,"times")){
    // x.*y
    if(mY==1 && nY==1){
      mZ=mX;
      nZ=nX;
      Z=rmat_allocate_prec(mZ,nZ,prec);
      rmat_mul_r(mZ,nZ,Z,mZ,X,mX,Y[0]);
    }else if(mX==1 && nX==1){
      mZ=mY;
      nZ=nY;
      Z=rmat_allocate_prec(mZ,nZ,prec);
      rmat_mul_r(mZ,nZ,Z,mZ,Y,mY,X[0]);
    }else if(mX==mY && nX==nY){
      mZ=mX;
      nZ=nX;
      Z=rmat_allocate_prec(mZ,nZ,prec);
      rmat_mul(mZ,nZ,Z,mZ,X,mX,Y,mY);
    }else{ mexErrMsgTxt("rmulti_two, Z=X.*Y. The sizes of X,Y should be same."); }
  }else if(char_eq(cmd,"mtimes")){
    // x*y
    if(mY==1 && nY==1){
      mZ=mX;
      nZ=nX;
      Z=rmat_allocate_prec(mZ,nZ,prec);
      rmat_mul_r(mZ,nZ,Z,mZ,X,mX,Y[0]);
    }else if(mX==1 && nX==1){
      mZ=mY;
      nZ=nY;
      Z=rmat_allocate_prec(mZ,nZ,prec);
      rmat_mul_r(mZ,nZ,Z,mZ,Y,mY,X[0]);
    }else if(nX==mY){
      mZ=mX;
      nZ=nY;
      Z=rmat_allocate_prec(mZ,nZ,prec);
      rmat_prod(mX,nX,nY,Z,mZ,X,mX,Y,mY);
    }else{ mexErrMsgTxt("rmulti_two, Z=X*Y. The sizes of X,Y should be (m,n) and (n,k)."); }
  }else if(char_eq(cmd,"rdivide")){
    // x./y
    if(mX!=mY || nX!=nY){ mexErrMsgTxt("rmulti_two, Z=X./Y. The sizes of X,Y should be same."); }
    mZ=mX;
    nZ=nX;
    Z=rmat_allocate_prec(mZ,nZ,prec);
    rmat_div(mZ,nZ,Z,mZ,X,mX,Y,mY);
  }else if(char_eq(cmd,"mldivide")){
    // z=x\y
    if(mX!=nX || nX!=mY){ mexErrMsgTxt("rmulti_two, Z=X\\Y. The sizes of X,Y should be (n,n) and (n,k)."); }
    mZ=mY;
    nZ=nY;
    Z=rmat_allocate_prec(mZ,nZ,prec); 
    rmat_copy(mZ,nZ,Z,mZ,Y,mY);
    rsolve(mZ,nZ,Z,mZ,X,mX,&info);
  }else{ mexErrMsgTxt("MATLAB:rmulti_two. The command is illeagal."); }

  // done
  plhs[0]=rmat2mx(mZ,nZ,Z,mZ);
       
  // free
  X=rmat_free(mX,nX,X);
  Y=rmat_free(mY,nY,Y);
  Z=rmat_free(mZ,nZ,Z);
}
