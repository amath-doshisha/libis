#include "mex.h"
#include <string.h>
#include <isys.h>
#include"./mxGetDataLong.c"
#include"./rmat2mx.c"
#include"./mx2rmat.c"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
#define CHECK_NUMERIC_SCALAR(X) (mxIsNumeric(X) && mxGetM(X)==1 && mxGetN(X)==1)
#define ERROR(X,Y)              (mexErrMsgIdAndTxt(X,Y))
  int n=0,m=0,LDA=0,i,j;
  int mX=0,nX=0,LDX=0,mY=0,nY=0,LDY=0,mZ=0,nZ=0,LDZ=0,info;
  long prec=53;
  rmulti **rA=NULL;
  rmulti **rX=NULL,**rY=NULL,**rZ=NULL;
  char *mode=NULL,*cmd=NULL;

  (void) prhs;
    
  /*
   * check number of argments
   */
  if(nrhs<3){
    ERROR("MATLAB:rmulti_data:maxrhs","Too few input arguments.");
  }

  /*
   * get precision
   */
  if(CHECK_NUMERIC_SCALAR(prhs[0])){
    prec=mxGetDataLong(prhs[0]);
  }else{
    ERROR("MATLAB:rmulti_data","The arg1 should be Integer(1x1).");
  }

  /*
   * get mode
   */
  if(mxIsChar(prhs[1])){
    mode=mxArrayToString(prhs[1]);
  }else{
    ERROR("MATLAB:rmulti_data","The arg2 should be Char.");
  }

  /*
   * get command
   */
  if(mxIsChar(prhs[2])){
    cmd=mxArrayToString(prhs[2]);
  }else{
    ERROR("MATLAB:rmulti_data","The arg3 should be Char.");
  }


  /*
   * command
   */ 
  if(char_eq(mode,"init") && char_eq(cmd,"set")){
    /*
     * allocate rmulti by double
     */ 
    if(nrhs>4){ ERROR("MATLAB:rmulti_data:maxrhs","Too many input arguments."); }
    if(nlhs>1){ ERROR("MATLAB:rmulti_data:maxlhs","Too many output arguments."); }
    if(!mxIsDouble(prhs[3])){ ERROR("MATLAB:rmulti_data","The arg4 should be Double(m,n)."); }
    m=mxGetM(prhs[3]); n=mxGetN(prhs[3]); LDA=m; rA=rmat_allocate_prec(LDA,n,prec);
    rmat_set_d(m,n,rA,LDA,mxGetPr(prhs[3]),m);
    plhs[0]=rmat2mx(m,n,rA,LDA);
    rA=rmat_free(LDA,n,rA);
  }else if(char_eq(mode,"init") && char_eq(cmd,"zeros")){
    /*
     * A=zeros(m,n)
     */ 
    if(nrhs>5){ ERROR("MATLAB:rmulti_data:maxrhs","Too many input arguments."); }
    if(nlhs>1){ ERROR("MATLAB:rmulti_data:maxlhs","Too many output arguments."); }
    if(CHECK_NUMERIC_SCALAR(prhs[3])){ m=mxGetDataLong(prhs[3]); }else{ ERROR("MATLAB:rmulti_data","The arg4 should be Integer(1,1)."); }
    if(CHECK_NUMERIC_SCALAR(prhs[4])){ n=mxGetDataLong(prhs[4]); }else{ ERROR("MATLAB:rmulti_data","The arg5 should be Integer(1,1)."); }
    LDA=m; rA=rmat_allocate_prec(LDA,n,prec);
    rmat_set_zeros(m,n,rA,LDA);
    plhs[0]=rmat2mx(m,n,rA,LDA);
    rA=rmat_free(LDA,n,rA);
  }else if(char_eq(mode,"init") && char_eq(cmd,"ones")){
    /*
     * A=ones(m,n)
     */ 
    if(nrhs>5){ ERROR("MATLAB:rmulti_data:maxrhs","Too many input arguments."); }
    if(nlhs>1){ ERROR("MATLAB:rmulti_data:maxlhs","Too many output arguments."); }
    if(CHECK_NUMERIC_SCALAR(prhs[3])){ m=mxGetDataLong(prhs[3]); }else{ ERROR("MATLAB:rmulti_data","The arg4 should be Integer(1,1)."); }
    if(CHECK_NUMERIC_SCALAR(prhs[4])){ n=mxGetDataLong(prhs[4]); }else{ ERROR("MATLAB:rmulti_data","The arg5 should be Integer(1,1)."); }
    LDA=m; rA=rmat_allocate_prec(LDA,n,prec);
    rmat_set_ones(m,n,rA,LDA);
    plhs[0]=rmat2mx(m,n,rA,LDA);
    rA=rmat_free(LDA,n,rA);
  }else if(char_eq(mode,"init") && char_eq(cmd,"eye")){
    /*
     * A=eye(m,n)
     */ 
    if(nrhs>5){ ERROR("MATLAB:rmulti_data:maxrhs","Too many input arguments."); }
    if(nlhs>1){ ERROR("MATLAB:rmulti_data:maxlhs","Too many output arguments."); }
    if(CHECK_NUMERIC_SCALAR(prhs[3])){ m=mxGetDataLong(prhs[3]); }else{ ERROR("MATLAB:rmulti_data","The arg4 should be Integer(1,1)."); }
    if(CHECK_NUMERIC_SCALAR(prhs[4])){ n=mxGetDataLong(prhs[4]); }else{ ERROR("MATLAB:rmulti_data","The arg5 should be Integer(1,1)."); }
    LDA=m; rA=rmat_allocate_prec(LDA,n,prec);
    rmat_set_eye(m,n,rA,LDA); plhs[0]=rmat2mx(m,n,rA,LDA);
    rA=rmat_free(LDA,n,rA);
  }else if(char_eq(mode,"init") && char_eq(cmd,"rand")){
    /*
     * A=rand(m,n)
     */ 
    if(nrhs>5){ ERROR("MATLAB:rmulti_data:maxrhs","Too many input arguments."); }
    if(nlhs>1){ ERROR("MATLAB:rmulti_data:maxlhs","Too many output arguments."); }
    if(CHECK_NUMERIC_SCALAR(prhs[3])){ m=mxGetDataLong(prhs[3]); }else{ ERROR("MATLAB:rmulti_data","The arg4 should be Integer(1,1)."); }
    if(CHECK_NUMERIC_SCALAR(prhs[4])){ n=mxGetDataLong(prhs[4]); }else{ ERROR("MATLAB:rmulti_data","The arg5 should be Integer(1,1)."); }
    LDA=m; rA=rmat_allocate_prec(LDA,n,prec);
    rmat_set_rand(m,n,rA,LDA,1,0);
    plhs[0]=rmat2mx(m,n,rA,LDA);
    rA=rmat_free(LDA,n,rA);
  }else if(char_eq(mode,"one") && char_eq(cmd,"double")){
    /*
     * Y=double(X)
     */ 
    if(nrhs>4){ ERROR("MATLAB:rmulti_data:maxrhs","Too many input arguments."); }
    if(nlhs>1){ ERROR("MATLAB:rmulti_data:maxlhs","Too many output arguments."); }
    if(!mxIsStruct(prhs[3])){ ERROR("MATLAB:rmulti_data","The arg4 should be Struct."); }
    m=mxGetM(prhs[3]); n=mxGetN(prhs[3]); LDA=m; rA=rmat_allocate_prec(LDA,n,prec);
    mx2rmat(m,n,rA,LDA,prhs[3]);
    rmat_round(m,n,rA,LDA,53);
    plhs[0]=mxCreateDoubleMatrix(m,n,mxREAL);
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	mxGetPr(plhs[0])[j*m+i]=rget_d(MAT(rA,i,j,LDA));
      }
    }
    rA=rmat_free(LDA,n,rA);
  }else if(char_eq(mode,"one") && char_eq(cmd,"prec")){
    /*
     * get to prec
     */ 
    if(nrhs>4){ ERROR("MATLAB:rmulti_data:maxrhs","Too many input arguments."); }
    if(nlhs>1){ ERROR("MATLAB:rmulti_data:maxlhs","Too many output arguments."); }
    if(!mxIsStruct(prhs[3])){ ERROR("MATLAB:rmulti_data","The arg4 should be Struct."); }
    m=mxGetM(prhs[3]); n=mxGetN(prhs[3]); LDA=m; rA=rmat_allocate_prec(LDA,n,prec);
    mx2rmat(m,n,rA,LDA,prhs[3]);
    plhs[0]=mxCreateDoubleMatrix(m,n,mxREAL);
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	mxGetPr(plhs[0])[j*m+i]=rget_prec(MAT(rA,i,j,LDA));
      }
    }
    rA=rmat_free(LDA,n,rA);
  }else if(char_eq(mode,"one") && char_eq(cmd,"exp")){
    /*
     * get to exp
     */ 
    if(nrhs>4){ ERROR("MATLAB:rmulti_data:maxrhs","Too many input arguments."); }
    if(nlhs>1){ ERROR("MATLAB:rmulti_data:maxlhs","Too many output arguments."); }
    if(!mxIsStruct(prhs[3])){ ERROR("MATLAB:rmulti_data","The arg4 should be Struct."); }
    m=mxGetM(prhs[3]); n=mxGetN(prhs[3]); LDA=m; rA=rmat_allocate_prec(LDA,n,prec);
    mx2rmat(m,n,rA,LDA,prhs[3]);
    plhs[0]=mxCreateDoubleMatrix(m,n,mxREAL);
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	mxGetPr(plhs[0])[j*m+i]=rget_exp(MAT(rA,i,j,LDA));
      }
    }
    rA=rmat_free(LDA,n,rA);
  }else if(char_eq(mode,"one") && char_eq(cmd,"sign")){
    /*
     * get to sign
     */ 
    if(nrhs>4){ ERROR("MATLAB:rmulti_data:maxrhs","Too many input arguments."); }
    if(nlhs>1){ ERROR("MATLAB:rmulti_data:maxlhs","Too many output arguments."); }
    if(!mxIsStruct(prhs[3])){ ERROR("MATLAB:rmulti_data","The arg4 should be Struct."); }
    m=mxGetM(prhs[3]); n=mxGetN(prhs[3]); LDA=m; rA=rmat_allocate_prec(LDA,n,prec);
    mx2rmat(m,n,rA,LDA,prhs[3]);
    plhs[0]=mxCreateDoubleMatrix(m,n,mxREAL);
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	mxGetPr(plhs[0])[j*m+i]=rget_sgn(MAT(rA,i,j,LDA));
      }
    }
    rA=rmat_free(LDA,n,rA);
  }else if(char_eq(mode,"one") && char_eq(cmd,"round")){
    /*
     * y=round(x,prec)
     */ 
    if(nrhs>4){ ERROR("MATLAB:rmulti_data:maxrhs","Too many input arguments."); }
    if(nlhs>1){ ERROR("MATLAB:rmulti_data:maxlhs","Too many output arguments."); }
    if(!mxIsStruct(prhs[3])){ ERROR("MATLAB:rmulti_data","The arg4 should be Struct."); }
    m=mxGetM(prhs[3]); n=mxGetN(prhs[3]); LDA=m; rA=rmat_allocate_prec(LDA,n,prec);
    mx2rmat(m,n,rA,LDA,prhs[3]);
    rmat_round(m,n,rA,LDA,prec);
    plhs[0]=rmat2mx(m,n,rA,LDA);
    rA=rmat_free(LDA,n,rA);
  }else if(char_eq(mode,"one") && char_eq(cmd,"abs")){
    /*
     * y=abs(x)
     */ 
    if(nrhs>4){ ERROR("MATLAB:rmulti_data:maxrhs","Too many input arguments."); }
    if(nlhs>1){ ERROR("MATLAB:rmulti_data:maxlhs","Too many output arguments."); }
    if(!mxIsStruct(prhs[3])){ ERROR("MATLAB:rmulti_data","The arg4 should be Struct."); }
    m=mxGetM(prhs[3]);
    n=mxGetN(prhs[3]);
    LDA=m; rA=rmat_allocate_prec(LDA,n,prec);
    mx2rmat(m,n,rA,LDA,prhs[3]);
    rmat_abs(m,n,rA,LDA,rA,LDA);
    rmat_round(m,n,rA,LDA,prec);
    plhs[0]=rmat2mx(m,n,rA,LDA);
    rA=rmat_free(LDA,n,rA);
  }else if(char_eq(mode,"one") && char_eq(cmd,"max")){
    /*
     * y=max(x)
     */ 
    if(nrhs>4){ ERROR("MATLAB:rmulti_data:maxrhs","Too many input arguments."); }
    if(nlhs>1){ ERROR("MATLAB:rmulti_data:maxlhs","Too many output arguments."); }
    if(!mxIsStruct(prhs[3])){ ERROR("MATLAB:rmulti_data","The arg4 should be Struct."); }
    // allocate
    mX=mxGetM(prhs[3]);
    nX=mxGetN(prhs[3]);
    LDX=mX; rX=rmat_allocate_prec(LDX,nX,prec);
    mx2rmat(mX,nX,rX,LDX,prhs[3]);
    if(mX==1){
      mY=1; nY=1;
      LDY=mY; rY=rmat_allocate_prec(LDY,nY,prec);
      rvec_max(rY[0],nX,rX);
    }else{
      mY=1; nY=nX;
      LDY=mY; rY=rmat_allocate_prec(LDY,nY,prec);
      for(j=0; j<nX; j++){ rvec_max(rY[j],mX,&COL(rX,j,LDX)); }
    }
    // done
    plhs[0]=rmat2mx(mY,nY,rY,LDY);
    // free
    rX=rmat_free(LDX,nX,rX);
    rY=rmat_free(LDY,nY,rY);
  }else if(char_eq(mode,"two") && char_eq(cmd,"plus")){
    /*
     * z=x+y
     */
    if(nrhs>5){ ERROR("MATLAB:rmulti_data:maxrhs","Too many input arguments."); }
    if(nlhs>1){ ERROR("MATLAB:rmulti_data:maxlhs","Too many output arguments."); }
    if(!mxIsStruct(prhs[3])){ ERROR("MATLAB:rmulti_data","The arg4 should be Struct."); }
    if(!mxIsStruct(prhs[4])){ ERROR("MATLAB:rmulti_data","The arg5 should be Struct."); }
    // allocate and compute
    mX=mxGetM(prhs[3]); nX=mxGetN(prhs[3]); LDX=mX; rX=rmat_allocate_prec(LDX,nX,prec); mx2rmat(mX,nX,rX,LDX,prhs[3]);
    mY=mxGetM(prhs[4]); nY=mxGetN(prhs[4]); LDY=mY; rY=rmat_allocate_prec(LDY,nY,prec); mx2rmat(mY,nY,rY,LDY,prhs[4]);
    if     (mY==1  && nY==1) { mZ=mX; nZ=nX; LDZ=mZ; rZ=rmat_allocate_prec(LDZ,nZ,prec); rmat_add_r(mZ,nZ,rZ,LDZ,rX,LDX,rY[0]); }
    else if(mX==1  && nX==1) { mZ=mY; nZ=nY; LDZ=mZ; rZ=rmat_allocate_prec(LDZ,nZ,prec); rmat_add_r(mZ,nZ,rZ,LDZ,rY,LDY,rX[0]); }
    else if(mX==mY && nX==nY){ mZ=mX; nZ=nX; LDZ=mZ; rZ=rmat_allocate_prec(LDZ,nZ,prec); rmat_add  (mZ,nZ,rZ,LDZ,rX,LDX,rY,LDY); }
    else                     { ERROR("MATLAB:rmulti_data","Z=X+Y. The sizes of X,Y should be same."); }
    // done
    plhs[0]=rmat2mx(mZ,nZ,rZ,LDZ);
    // free
    rX=rmat_free(LDX,nX,rX);
    rY=rmat_free(LDY,nY,rY);
    rZ=rmat_free(LDZ,nZ,rZ);
  }else if(char_eq(mode,"two") && char_eq(cmd,"minus")){
    /*
     * z=x-y
     */
    if(nrhs>5){ ERROR("MATLAB:rmulti_data:maxrhs","Too many input arguments."); }
    if(nlhs>1){ ERROR("MATLAB:rmulti_data:maxlhs","Too many output arguments."); }
    if(!mxIsStruct(prhs[3])){ ERROR("MATLAB:rmulti_data","The arg4 should be Struct."); }
    if(!mxIsStruct(prhs[4])){ ERROR("MATLAB:rmulti_data","The arg5 should be Struct."); }
    // allocate and compute
    mX=mxGetM(prhs[3]); nX=mxGetN(prhs[3]); LDX=mX; rX=rmat_allocate_prec(LDX,nX,prec); mx2rmat(mX,nX,rX,LDX,prhs[3]);
    mY=mxGetM(prhs[4]); nY=mxGetN(prhs[4]); LDY=mY; rY=rmat_allocate_prec(LDY,nY,prec); mx2rmat(mY,nY,rY,LDY,prhs[4]);
    if     (mY==1  && nY==1) { mZ=mX; nZ=nX; LDZ=mZ; rZ=rmat_allocate_prec(LDZ,nZ,prec); rmat_sub_r2(mZ,nZ,rZ,LDZ,rX,LDX,rY[0]);  }
    else if(mX==1  && nX==1) { mZ=mY; nZ=nY; LDZ=mZ; rZ=rmat_allocate_prec(LDZ,nZ,prec); rmat_sub_r1(mZ,nZ,rZ,LDZ,rX[0], rY,LDY); }
    else if(mX==mY && nX==nY){ mZ=mX; nZ=nX; LDZ=mZ; rZ=rmat_allocate_prec(LDZ,nZ,prec); rmat_sub   (mZ,nZ,rZ,LDZ,rX,LDX,rY,LDY); }
    else                     { ERROR("MATLAB:rmulti_data","Z=X-Y. The sizes of X,Y should be same."); }
    // done
    plhs[0]=rmat2mx(mZ,nZ,rZ,LDZ);
    // free
    rX=rmat_free(LDX,nX,rX);
    rY=rmat_free(LDY,nY,rY);
    rZ=rmat_free(LDZ,nZ,rZ);
  }else if(char_eq(mode,"two") && char_eq(cmd,"times")){
    /*
     * z=x.*y
     */
    if(nrhs>5){ ERROR("MATLAB:rmulti_data:maxrhs","Too many input arguments."); }
    if(nlhs>1){ ERROR("MATLAB:rmulti_data:maxlhs","Too many output arguments."); }
    if(!mxIsStruct(prhs[3])){ ERROR("MATLAB:rmulti_data","The arg4 should be Struct."); }
    if(!mxIsStruct(prhs[4])){ ERROR("MATLAB:rmulti_data","The arg5 should be Struct."); }
    // allocate and compute
    mX=mxGetM(prhs[3]); nX=mxGetN(prhs[3]); LDX=mX; rX=rmat_allocate_prec(LDX,nX,prec); mx2rmat(mX,nX,rX,LDX,prhs[3]);
    mY=mxGetM(prhs[4]); nY=mxGetN(prhs[4]); LDY=mY; rY=rmat_allocate_prec(LDY,nY,prec); mx2rmat(mY,nY,rY,LDY,prhs[4]);
    if     (mY==1  && nY==1) { mZ=mX; nZ=nX; LDZ=mZ; rZ=rmat_allocate_prec(LDZ,nZ,prec); rmat_mul_r(mZ,nZ,rZ,LDZ,rX,LDX,rY[0]);  }
    else if(mX==1  && nX==1) { mZ=mY; nZ=nY; LDZ=mZ; rZ=rmat_allocate_prec(LDZ,nZ,prec); rmat_mul_r(mZ,nZ,rZ,LDZ,rY,LDY,rX[0]);  }
    else if(mX==mY && nX==nY){ mZ=mX; nZ=nX; LDZ=mZ; rZ=rmat_allocate_prec(LDZ,nZ,prec); rmat_mul  (mZ,nZ,rZ,LDZ,rX,LDX,rY,LDY); }
    else                     { ERROR("MATLAB:rmulti_data","Z=X.*Y. The sizes of X,Y should be same."); }
    // done
    plhs[0]=rmat2mx(mZ,nZ,rZ,LDZ);
    // free
    rX=rmat_free(LDX,nX,rX);
    rY=rmat_free(LDY,nY,rY);
    rZ=rmat_free(LDZ,nZ,rZ);
  }else if(char_eq(mode,"two") && char_eq(cmd,"mtimes")){
    /*
     * z=x*y
     */
    if(nrhs>5){ ERROR("MATLAB:rmulti_data:maxrhs","Too many input arguments."); }
    if(nlhs>1){ ERROR("MATLAB:rmulti_data:maxlhs","Too many output arguments."); }
    if(!mxIsStruct(prhs[3])){ ERROR("MATLAB:rmulti_data","The arg4 should be Struct."); }
    if(!mxIsStruct(prhs[4])){ ERROR("MATLAB:rmulti_data","The arg5 should be Struct."); }
    // allocate and compute
    mX=mxGetM(prhs[3]); nX=mxGetN(prhs[3]); LDX=mX; rX=rmat_allocate_prec(LDX,nX,prec); mx2rmat(mX,nX,rX,LDX,prhs[3]);
    mY=mxGetM(prhs[4]); nY=mxGetN(prhs[4]); LDY=mY; rY=rmat_allocate_prec(LDY,nY,prec); mx2rmat(mY,nY,rY,LDY,prhs[4]);
    if     (mY==1 && nY==1){ mZ=mX; nZ=nX; LDZ=mZ; rZ=rmat_allocate_prec(LDZ,nZ,prec); rmat_mul_r(mZ,nZ,rZ,LDZ,rX,LDX,rY[0]); }
    else if(mX==1 && nX==1){ mZ=mY; nZ=nY; LDZ=mZ; rZ=rmat_allocate_prec(LDZ,nZ,prec); rmat_mul_r(mZ,nZ,rZ,LDZ,rY,LDY,rX[0]); }
    else if(nX==mY)        { mZ=mX; nZ=nY; LDZ=mZ; rZ=rmat_allocate_prec(LDZ,nZ,prec); rmat_prod(mX,nX,nY,rZ,LDZ,rX,LDX,rY,LDY); }
    else                   { ERROR("MATLAB:rmulti_data","Z=X*Y. The sizes of X,Y should be (m,n) and (n,k)."); }
    // done
    plhs[0]=rmat2mx(mZ,nZ,rZ,LDZ);
    // free
    rX=rmat_free(LDX,nX,rX);
    rY=rmat_free(LDY,nY,rY);
    rZ=rmat_free(LDZ,nZ,rZ);
  }else if(char_eq(mode,"two") && char_eq(cmd,"rdivide")){
    /*
     * z=x./y
     */
    if(nrhs>5){ ERROR("MATLAB:rmulti_data:maxrhs","Too many input arguments."); }
    if(nlhs>1){ ERROR("MATLAB:rmulti_data:maxlhs","Too many output arguments."); }
    if(!mxIsStruct(prhs[3])){ ERROR("MATLAB:rmulti_data","The arg4 should be Struct."); }
    if(!mxIsStruct(prhs[4])){ ERROR("MATLAB:rmulti_data","The arg5 should be Struct."); }
    // allocate and compute
    mX=mxGetM(prhs[3]); nX=mxGetN(prhs[3]); LDX=mX; rX=rmat_allocate_prec(LDX,nX,prec); mx2rmat(mX,nX,rX,LDX,prhs[3]);
    mY=mxGetM(prhs[4]); nY=mxGetN(prhs[4]); LDY=mY; rY=rmat_allocate_prec(LDY,nY,prec); mx2rmat(mY,nY,rY,LDY,prhs[4]);
    if     (mY==1  && nY==1) { mZ=mX; nZ=nX; LDZ=mZ; rZ=rmat_allocate_prec(LDZ,nZ,prec); rmat_div_r2(mZ,nZ,rZ,LDZ,rX,LDX,rY[0]);  }
    else if(mX==1  && nX==1) { mZ=mY; nZ=nY; LDZ=mZ; rZ=rmat_allocate_prec(LDZ,nZ,prec); rmat_div_r1(mZ,nZ,rZ,LDZ,rX[0],rY,LDY);  }
    else if(mX==mY && nX==nY){ mZ=mX; nZ=nX; LDZ=mZ; rZ=rmat_allocate_prec(LDZ,nZ,prec); rmat_div   (mZ,nZ,rZ,LDZ,rX,LDX,rY,LDY); }
    else                     { ERROR("MATLAB:rmulti_data","Z=X./Y. The sizes of X,Y should be same."); }
    // done
    plhs[0]=rmat2mx(mZ,nZ,rZ,LDZ);
    // free
    rX=rmat_free(LDX,nX,rX);
    rY=rmat_free(LDY,nY,rY);
    rZ=rmat_free(LDZ,nZ,rZ);
  }else if(char_eq(mode,"two") && char_eq(cmd,"mldivide")){
    /*
     * z=x\y
     */
    if(nrhs>5){ ERROR("MATLAB:rmulti_data:maxrhs","Too many input arguments."); }
    if(nlhs>1){ ERROR("MATLAB:rmulti_data:maxlhs","Too many output arguments."); }
    if(!mxIsStruct(prhs[3])){ ERROR("MATLAB:rmulti_data","The arg4 should be Struct."); }
    if(!mxIsStruct(prhs[4])){ ERROR("MATLAB:rmulti_data","The arg5 should be Struct."); }
    // allocate and compute
    mX=mxGetM(prhs[3]); nX=mxGetN(prhs[3]); LDX=mX; rX=rmat_allocate_prec(LDX,nX,prec); mx2rmat(mX,nX,rX,LDX,prhs[3]);
    mY=mxGetM(prhs[4]); nY=mxGetN(prhs[4]); LDY=mY; rY=rmat_allocate_prec(LDY,nY,prec); mx2rmat(mY,nY,rY,LDY,prhs[4]);
    // 
    if(mX!=nX || nX!=mY){ ERROR("MATLAB:rmulti_data","Z=X\\Y. The sizes of X,Y should be (n,n) and (n,k)."); }
    mZ=mY; nZ=nY; LDZ=mZ; rZ=rmat_allocate_prec(LDZ,nZ,prec); 
    rmat_copy(mZ,nZ,rZ,LDZ,rY,LDY); rsolve(mZ,nZ,rZ,LDZ,rX,LDX,&info);
    // done
    plhs[0]=rmat2mx(mZ,nZ,rZ,LDZ);
    // free
    rX=rmat_free(LDX,nX,rX);
    rY=rmat_free(LDY,nY,rY);
    rZ=rmat_free(LDZ,nZ,rZ);
  }else{
    ERROR("MATLAB:rmulti_data","Illeagal command.");
  }
        
}
