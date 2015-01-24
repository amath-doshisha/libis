#include "mex.h"
#include <string.h>
#include <isys.h>
#include"../../@rmulti/private/mxGetDataLong.c"
#include"../../@rmulti/private/mx2rmat.c"
#include"../../@rmulti/private/rmat2mx.c"
#include"./cmat2mx.c"
#include"./mx2cmat.c"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
#define CHECK_NUMERIC_SCALAR(X) (mxIsNumeric(X) && mxGetM(X)==1 && mxGetN(X)==1)
#define ERROR(X,Y)              (mexErrMsgIdAndTxt(X,Y))
  int n=0,m=0,LDA=0,i,j,mX=0,nX=0,LDX=0,mY=0,nY=0,LDY=0,mZ=0,nZ=0,LDZ=0,info,LDB=0;
  long prec=53;
  cmulti **cA=NULL,**cX=NULL,**cY=NULL,**cZ=NULL;
  rmulti **rB=NULL;
  char *mode=NULL,*cmd=NULL;
  dcomplex z;
    
  (void) prhs;
    
  /*
   * check number of argments
   */
  if(nrhs<3){
    ERROR("MATLAB:cmulti_data:maxrhs","Too few input arguments.");
  }

  /*
   * get precision
   */
  if(CHECK_NUMERIC_SCALAR(prhs[0])){
    prec=mxGetDataLong(prhs[0]);
  }else{
    ERROR("MATLAB:cmulti_data","The arg1 should be Integer(1x1).");
  }

  /*
   * get mode
   */
  if(mxIsChar(prhs[1])){
    mode=mxArrayToString(prhs[1]);
  }else{
    ERROR("MATLAB:cmulti_data","The arg2 should be Char.");
  }

  /*
   * get command
   */
  if(mxIsChar(prhs[2])){
    cmd=mxArrayToString(prhs[2]);
  }else{
    ERROR("MATLAB:cmulti_data","The arg3 should be Char.");
  }

  /*
   * command
   */ 
  if(char_eq(mode,"init") && char_eq(cmd,"set")){
    /*
     * allocate cmulti by double
     */ 
    if(nrhs>4){ ERROR("MATLAB:cmulti_data:maxrhs","Too many input arguments."); }
    if(nlhs>1){ ERROR("MATLAB:cmulti_data:maxlhs","Too many output arguments."); }
    if(!mxIsDouble(prhs[3])){ ERROR("MATLAB:cmulti_data","The arg4 should be Double(m,n)."); }
    m=mxGetM(prhs[3]); n=mxGetN(prhs[3]); LDA=m; cA=cmat_allocate_prec(LDA,n,prec);
    if(mxIsComplex(prhs[3])){ cmat_set_dd(m,n,cA,LDA,mxGetPr(prhs[3]),m,mxGetPi(prhs[3]),m); }
    else                    { cmat_set_d (m,n,cA,LDA,mxGetPr(prhs[3]),m); }
    plhs[0]=cmat2mx(m,n,cA,LDA);
    cA=cmat_free(LDA,n,cA);
  }else if(char_eq(mode,"init") && char_eq(cmd,"set_r")){
    /*
     * allocate cmulti by rmulti
     */ 
    if(nrhs>4){ ERROR("MATLAB:cmulti_data:maxrhs","Too many input arguments."); }
    if(nlhs>1){ ERROR("MATLAB:cmulti_data:maxlhs","Too many output arguments."); }
    if(!mxIsStruct(prhs[3])){ ERROR("MATLAB:cmulti_data","The arg4 should be Struct(m,n)."); }
    m=mxGetM(prhs[3]); n=mxGetN(prhs[3]);
    LDB=m; rB=rmat_allocate_prec(LDB,n,prec); mx2rmat(m,n,rB,LDB,prhs[3]);
    LDA=m; cA=cmat_allocate_prec(LDA,n,prec); cmat_copy_r(m,n,cA,LDA,rB,LDB);
    plhs[0]=cmat2mx(m,n,cA,LDA);
    cA=cmat_free(LDA,n,cA);
    rB=rmat_free(LDB,n,rB);
  }else if(char_eq(mode,"init") && char_eq(cmd,"zeros")){
    /*
     * A=zeros(m,n)
     */ 
    if(nrhs>5){ ERROR("MATLAB:cmulti_data:maxrhs","Too many input arguments."); }
    if(nlhs>1){ ERROR("MATLAB:cmulti_data:maxlhs","Too many output arguments."); }
    if(CHECK_NUMERIC_SCALAR(prhs[3])){ m=mxGetDataLong(prhs[3]); }else{ ERROR("MATLAB:cmulti_data","The arg4 should be Integer(1,1)."); }
    if(CHECK_NUMERIC_SCALAR(prhs[4])){ n=mxGetDataLong(prhs[4]); }else{ ERROR("MATLAB:cmulti_data","The arg5 should be Integer(1,1)."); }
    LDA=m; cA=cmat_allocate_prec(LDA,n,prec);
    cmat_set_zeros(m,n,cA,LDA);
    plhs[0]=cmat2mx(m,n,cA,LDA);
    cA=cmat_free(LDA,n,cA);
  }else if(char_eq(mode,"init") && char_eq(cmd,"ones")){
    /*
     * A=ones(m,n)
     */ 
    if(nrhs>5){ ERROR("MATLAB:cmulti_data:maxrhs","Too many input arguments."); }
    if(nlhs>1){ ERROR("MATLAB:cmulti_data:maxlhs","Too many output arguments."); }
    if(CHECK_NUMERIC_SCALAR(prhs[3])){ m=mxGetDataLong(prhs[3]); }else{ ERROR("MATLAB:cmulti_data","The arg4 should be Integer(1,1)."); }
    if(CHECK_NUMERIC_SCALAR(prhs[4])){ n=mxGetDataLong(prhs[4]); }else{ ERROR("MATLAB:cmulti_data","The arg5 should be Integer(1,1)."); }
    LDA=m; cA=cmat_allocate_prec(LDA,n,prec);
    cmat_set_ones(m,n,cA,LDA);
    plhs[0]=cmat2mx(m,n,cA,LDA);
    cA=cmat_free(LDA,n,cA);
  }else if(char_eq(mode,"init") && char_eq(cmd,"eye")){
    /*
     * A=eye(m,n)
     */ 
    if(nrhs>5){ ERROR("MATLAB:cmulti_data:maxrhs","Too many input arguments."); }
    if(nlhs>1){ ERROR("MATLAB:cmulti_data:maxlhs","Too many output arguments."); }
    if(CHECK_NUMERIC_SCALAR(prhs[3])){ m=mxGetDataLong(prhs[3]); }else{ ERROR("MATLAB:cmulti_data","The arg4 should be Integer(1,1)."); }
    if(CHECK_NUMERIC_SCALAR(prhs[4])){ n=mxGetDataLong(prhs[4]); }else{ ERROR("MATLAB:cmulti_data","The arg5 should be Integer(1,1)."); }
    LDA=m; cA=cmat_allocate_prec(LDA,n,prec);
    cmat_set_eye(m,n,cA,LDA); plhs[0]=cmat2mx(m,n,cA,LDA);
    cA=cmat_free(LDA,n,cA);
  }else if(char_eq(mode,"init") && char_eq(cmd,"rand")){
    /*
     * A=rand(m,n)
     */ 
    if(nrhs>5){ ERROR("MATLAB:cmulti_data:maxrhs","Too many input arguments."); }
    if(nlhs>1){ ERROR("MATLAB:cmulti_data:maxlhs","Too many output arguments."); }
    if(CHECK_NUMERIC_SCALAR(prhs[3])){ m=mxGetDataLong(prhs[3]); }else{ ERROR("MATLAB:cmulti_data","The arg4 should be Integer(1,1)."); }
    if(CHECK_NUMERIC_SCALAR(prhs[4])){ n=mxGetDataLong(prhs[4]); }else{ ERROR("MATLAB:cmulti_data","The arg5 should be Integer(1,1)."); }
    LDA=m; cA=cmat_allocate_prec(LDA,n,prec);
    cmat_set_rand(m,n,cA,LDA,1,0);
    plhs[0]=cmat2mx(m,n,cA,LDA);
    cA=cmat_free(LDA,n,cA);
  }else if(char_eq(mode,"one") && char_eq(cmd,"double")){
    /*
     * Y=double(X)
     */ 
    if(nrhs>4){ ERROR("MATLAB:cmulti_data:maxrhs","Too many input arguments."); }
    if(nlhs>1){ ERROR("MATLAB:cmulti_data:maxlhs","Too many output arguments."); }
    if(!mxIsStruct(prhs[3])){ ERROR("MATLAB:cmulti_data","The arg4 should be Struct."); }
    m=mxGetM(prhs[3]); n=mxGetN(prhs[3]); LDA=m; cA=cmat_allocate_prec(LDA,n,prec);
    cmat_set_zeros(m,n,cA,LDA);
    mx2cmat(m,n,cA,LDA,prhs[3]);
    cmat_round(m,n,cA,LDA,53);
    plhs[0]=mxCreateDoubleMatrix(m,n,mxCOMPLEX);
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	z=cget_z(MAT(cA,i,j,m));
	mxGetPr(plhs[0])[j*m+i]=Z_R(z);
	mxGetPi(plhs[0])[j*m+i]=Z_I(z);
      }
    }
    cA=cmat_free(LDA,n,cA);
  }else if(char_eq(mode,"one") && char_eq(cmd,"prec")){
    /*
     * get to prec
     */ 
    if(nrhs>4){ ERROR("MATLAB:cmulti_data:maxrhs","Too many input arguments."); }
    if(nlhs>1){ ERROR("MATLAB:cmulti_data:maxlhs","Too many output arguments."); }
    if(!mxIsStruct(prhs[3])){ ERROR("MATLAB:cmulti_data","The arg4 should be Struct."); }
    m=mxGetM(prhs[3]); n=mxGetN(prhs[3]); LDA=m; cA=cmat_allocate_prec(LDA,n,prec);
    mx2cmat(m,n,cA,LDA,prhs[3]);
    plhs[0]=mxCreateDoubleMatrix(m,2*n,mxREAL);
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	mxGetPr(plhs[0])[(2*j  )*m+i]=rget_prec(C_R(MAT(cA,i,j,LDA)));
	mxGetPr(plhs[0])[(2*j+1)*m+i]=rget_prec(C_I(MAT(cA,i,j,LDA)));
      }
    }
    cA=cmat_free(LDA,n,cA);
  }else if(char_eq(mode,"one") && char_eq(cmd,"exp")){
    /*
     * get to exp
     */ 
    if(nrhs>4){ ERROR("MATLAB:cmulti_data:maxrhs","Too many input arguments."); }
    if(nlhs>1){ ERROR("MATLAB:cmulti_data:maxlhs","Too many output arguments."); }
    if(!mxIsStruct(prhs[3])){ ERROR("MATLAB:cmulti_data","The arg4 should be Struct."); }
    m=mxGetM(prhs[3]); n=mxGetN(prhs[3]); LDA=m; cA=cmat_allocate_prec(LDA,n,prec);
    mx2cmat(m,n,cA,LDA,prhs[3]);
    plhs[0]=mxCreateDoubleMatrix(m,2*n,mxREAL);
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	mxGetPr(plhs[0])[(2*j  )*m+i]=rget_exp(C_R(MAT(cA,i,j,LDA)));
	mxGetPr(plhs[0])[(2*j+1)*m+i]=rget_exp(C_I(MAT(cA,i,j,LDA)));
      }
    }
    cA=cmat_free(LDA,n,cA);
  }else if(char_eq(mode,"one") && char_eq(cmd,"sign")){
    /*
     * get to sign
     */ 
    if(nrhs>4){ ERROR("MATLAB:cmulti_data:maxrhs","Too many input arguments."); }
    if(nlhs>1){ ERROR("MATLAB:cmulti_data:maxlhs","Too many output arguments."); }
    if(!mxIsStruct(prhs[3])){ ERROR("MATLAB:cmulti_data","The arg4 should be Struct."); }
    m=mxGetM(prhs[3]); n=mxGetN(prhs[3]); LDA=m; cA=cmat_allocate_prec(LDA,n,prec);
    mx2cmat(m,n,cA,LDA,prhs[3]);
    plhs[0]=mxCreateDoubleMatrix(m,2*n,mxREAL);
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	mxGetPr(plhs[0])[(2*j  )*m+i]=rget_sgn(C_R(MAT(cA,i,j,LDA)));
	mxGetPr(plhs[0])[(2*j+1)*m+i]=rget_sgn(C_I(MAT(cA,i,j,LDA)));
      }
    }
    cA=cmat_free(LDA,n,cA);
  }else if(char_eq(mode,"one") && char_eq(cmd,"round")){
    /*
     * y=round(x,prec)
     */ 
    if(nrhs>4){ ERROR("MATLAB:cmulti_data:maxrhs","Too many input arguments."); }
    if(nlhs>1){ ERROR("MATLAB:cmulti_data:maxlhs","Too many output arguments."); }
    if(!mxIsStruct(prhs[3])){ ERROR("MATLAB:cmulti_data","The arg4 should be Struct."); }
    m=mxGetM(prhs[3]); n=mxGetN(prhs[3]); LDA=m; cA=cmat_allocate_prec(LDA,n,prec);
    mx2cmat(m,n,cA,LDA,prhs[3]);
    cmat_round(m,n,cA,LDA,prec);
    plhs[0]=cmat2mx(m,n,cA,LDA);
    cA=cmat_free(LDA,n,cA);
  }else if(char_eq(mode,"one") && char_eq(cmd,"real")){
    /*
     * y=real(x)
     */ 
    if(nrhs>4){ ERROR("MATLAB:cmulti_data:maxrhs","Too many input arguments."); }
    if(nlhs>1){ ERROR("MATLAB:cmulti_data:maxlhs","Too many output arguments."); }
    if(!mxIsStruct(prhs[3])){ ERROR("MATLAB:cmulti_data","The arg4 should be Struct."); }
    m=mxGetM(prhs[3]);
    n=mxGetN(prhs[3]);
    LDA=m; cA=cmat_allocate_prec(LDA,n,prec);
    LDB=m; rB=rmat_allocate_prec(LDB,n,prec);
    mx2cmat(m,n,cA,LDA,prhs[3]);
    cmat_real(m,n,rB,LDB,cA,LDA);
    plhs[0]=rmat2mx(m,n,rB,LDB);
    cA=cmat_free(LDA,n,cA);
    rB=rmat_free(LDB,n,rB);
  }else if(char_eq(mode,"one") && char_eq(cmd,"imag")){
    /*
     * y=imag(x)
     */ 
    if(nrhs>4){ ERROR("MATLAB:cmulti_data:maxrhs","Too many input arguments."); }
    if(nlhs>1){ ERROR("MATLAB:cmulti_data:maxlhs","Too many output arguments."); }
    if(!mxIsStruct(prhs[3])){ ERROR("MATLAB:cmulti_data","The arg4 should be Struct."); }
    m=mxGetM(prhs[3]);
    n=mxGetN(prhs[3]);
    LDA=m; cA=cmat_allocate_prec(LDA,n,prec);
    LDB=m; rB=rmat_allocate_prec(LDB,n,prec);
    mx2cmat(m,n,cA,LDA,prhs[3]);
    cmat_imag(m,n,rB,LDB,cA,LDA);
    plhs[0]=rmat2mx(m,n,rB,LDB);
    cA=cmat_free(LDA,n,cA);
    rB=rmat_free(LDB,n,rB);
  }else if(char_eq(mode,"one") && char_eq(cmd,"abs")){
    /*
     * y=abs(x)
     */ 
    if(nrhs>4){ ERROR("MATLAB:cmulti_data:maxrhs","Too many input arguments."); }
    if(nlhs>1){ ERROR("MATLAB:cmulti_data:maxlhs","Too many output arguments."); }
    if(!mxIsStruct(prhs[3])){ ERROR("MATLAB:cmulti_data","The arg4 should be Struct."); }
    m=mxGetM(prhs[3]);
    n=mxGetN(prhs[3]);
    LDA=m; cA=cmat_allocate_prec(LDA,n,prec);
    LDB=m; rB=rmat_allocate_prec(LDB,n,prec);
    mx2cmat(m,n,cA,LDA,prhs[3]);
    cmat_abs(m,n,rB,LDB,cA,LDA);
    plhs[0]=rmat2mx(m,n,rB,LDB);
    cA=cmat_free(LDA,n,cA);
    rB=rmat_free(LDB,n,rB);
  }else if(char_eq(mode,"one") && char_eq(cmd,"max")){
    /*
     * y=max(x)
     */ 
    if(nrhs>4){ ERROR("MATLAB:cmulti_data:maxrhs","Too many input arguments."); }
    if(nlhs>1){ ERROR("MATLAB:cmulti_data:maxlhs","Too many output arguments."); }
    if(!mxIsStruct(prhs[3])){ ERROR("MATLAB:cmulti_data","The arg4 should be Struct."); }
    // allocate
    mX=mxGetM(prhs[3]);
    nX=mxGetN(prhs[3]);
    LDX=mX; cX=cmat_allocate_prec(LDX,nX,prec);
    mx2cmat(mX,nX,cX,LDX,prhs[3]);
    if(mX==1){
      mY=1; nY=1;
      LDY=mY; cY=cmat_allocate_prec(LDY,nY,prec);
      cvec_max(cY[0],nX,cX);
    }else{
      mY=1; nY=nX;
      LDY=mY; cY=cmat_allocate_prec(LDY,nY,prec);
      for(j=0; j<nX; j++){ cvec_max(cY[j],mX,&COL(cX,j,LDX)); }
    }
    // done
    plhs[0]=cmat2mx(mY,nY,cY,LDY);
    // free
    cX=cmat_free(LDX,nX,cX);
    cY=cmat_free(LDY,nY,cY);
  }else if(char_eq(mode,"two") && char_eq(cmd,"plus")){
    /*
     * z=x+y
     */
    if(nrhs>5){ ERROR("MATLAB:cmulti_data:maxrhs","Too many input arguments."); }
    if(nlhs>1){ ERROR("MATLAB:cmulti_data:maxlhs","Too many output arguments."); }
    if(!mxIsStruct(prhs[3])){ ERROR("MATLAB:cmulti_data","The arg4 should be Struct."); }
    if(!mxIsStruct(prhs[4])){ ERROR("MATLAB:cmulti_data","The arg5 should be Struct."); }
    // allocate and compute
    mX=mxGetM(prhs[3]); nX=mxGetN(prhs[3]); LDX=mX; cX=cmat_allocate_prec(LDX,nX,prec); mx2cmat(mX,nX,cX,LDX,prhs[3]);
    mY=mxGetM(prhs[4]); nY=mxGetN(prhs[4]); LDY=mY; cY=cmat_allocate_prec(LDY,nY,prec); mx2cmat(mY,nY,cY,LDY,prhs[4]);
    if     (mY==1  && nY==1) { mZ=mX; nZ=nX; LDZ=mZ; cZ=cmat_allocate_prec(LDZ,nZ,prec); cmat_add_c(mZ,nZ,cZ,LDZ,cX,LDX,cY[0]); }
    else if(mX==1  && nX==1) { mZ=mY; nZ=nY; LDZ=mZ; cZ=cmat_allocate_prec(LDZ,nZ,prec); cmat_add_c(mZ,nZ,cZ,LDZ,cY,LDY,cX[0]); }
    else if(mX==mY && nX==nY){ mZ=mX; nZ=nX; LDZ=mZ; cZ=cmat_allocate_prec(LDZ,nZ,prec); cmat_add  (mZ,nZ,cZ,LDZ,cX,LDX,cY,LDY); }
    else                     { ERROR("MATLAB:cmulti_data","Z=X+Y. The sizes of X,Y should be same."); }
    // done
    plhs[0]=cmat2mx(mZ,nZ,cZ,LDZ);
    // free
    cX=cmat_free(LDX,nX,cX);
    cY=cmat_free(LDY,nY,cY);
    cZ=cmat_free(LDZ,nZ,cZ);
  }else if(char_eq(mode,"two") && char_eq(cmd,"minus")){
    /*
     * z=x-y
     */
    if(nrhs>5){ ERROR("MATLAB:cmulti_data:maxrhs","Too many input arguments."); }
    if(nlhs>1){ ERROR("MATLAB:cmulti_data:maxlhs","Too many output arguments."); }
    if(!mxIsStruct(prhs[3])){ ERROR("MATLAB:cmulti_data","The arg4 should be Struct."); }
    if(!mxIsStruct(prhs[4])){ ERROR("MATLAB:cmulti_data","The arg5 should be Struct."); }
    // allocate and compute
    mX=mxGetM(prhs[3]); nX=mxGetN(prhs[3]); LDX=mX; cX=cmat_allocate_prec(LDX,nX,prec); mx2cmat(mX,nX,cX,LDX,prhs[3]);
    mY=mxGetM(prhs[4]); nY=mxGetN(prhs[4]); LDY=mY; cY=cmat_allocate_prec(LDY,nY,prec); mx2cmat(mY,nY,cY,LDY,prhs[4]);
    if     (mY==1  && nY==1) { mZ=mX; nZ=nX; LDZ=mZ; cZ=cmat_allocate_prec(LDZ,nZ,prec); cmat_sub_c2(mZ,nZ,cZ,LDZ,cX,LDX,cY[0]);  }
    else if(mX==1  && nX==1) { mZ=mY; nZ=nY; LDZ=mZ; cZ=cmat_allocate_prec(LDZ,nZ,prec); cmat_sub_c1(mZ,nZ,cZ,LDZ,cX[0], cY,LDY); }
    else if(mX==mY && nX==nY){ mZ=mX; nZ=nX; LDZ=mZ; cZ=cmat_allocate_prec(LDZ,nZ,prec); cmat_sub   (mZ,nZ,cZ,LDZ,cX,LDX,cY,LDY); }
    else                     { ERROR("MATLAB:cmulti_data","Z=X-Y. The sizes of X,Y should be same."); }
    // done
    plhs[0]=cmat2mx(mZ,nZ,cZ,LDZ);
    // free
    cX=cmat_free(LDX,nX,cX);
    cY=cmat_free(LDY,nY,cY);
    cZ=cmat_free(LDZ,nZ,cZ);
  }else if(char_eq(mode,"two") && char_eq(cmd,"times")){
    /*
     * z=x.*y
     */
    if(nrhs>5){ ERROR("MATLAB:cmulti_data:maxrhs","Too many input arguments."); }
    if(nlhs>1){ ERROR("MATLAB:cmulti_data:maxlhs","Too many output arguments."); }
    if(!mxIsStruct(prhs[3])){ ERROR("MATLAB:cmulti_data","The arg4 should be Struct."); }
    if(!mxIsStruct(prhs[4])){ ERROR("MATLAB:cmulti_data","The arg5 should be Struct."); }
    // allocate and compute
    mX=mxGetM(prhs[3]); nX=mxGetN(prhs[3]); LDX=mX; cX=cmat_allocate_prec(LDX,nX,prec); mx2cmat(mX,nX,cX,LDX,prhs[3]);
    mY=mxGetM(prhs[4]); nY=mxGetN(prhs[4]); LDY=mY; cY=cmat_allocate_prec(LDY,nY,prec); mx2cmat(mY,nY,cY,LDY,prhs[4]);
    if     (mY==1  && nY==1) { mZ=mX; nZ=nX; LDZ=mZ; cZ=cmat_allocate_prec(LDZ,nZ,prec); cmat_mul_c(mZ,nZ,cZ,LDZ,cX,LDX,cY[0]);  }
    else if(mX==1  && nX==1) { mZ=mY; nZ=nY; LDZ=mZ; cZ=cmat_allocate_prec(LDZ,nZ,prec); cmat_mul_c(mZ,nZ,cZ,LDZ,cY,LDY,cX[0]);  }
    else if(mX==mY && nX==nY){ mZ=mX; nZ=nX; LDZ=mZ; cZ=cmat_allocate_prec(LDZ,nZ,prec); cmat_mul  (mZ,nZ,cZ,LDZ,cX,LDX,cY,LDY); }
    else                     { ERROR("MATLAB:cmulti_data","Z=X.*Y. The sizes of X,Y should be same."); }
    // done
    plhs[0]=cmat2mx(mZ,nZ,cZ,LDZ);
    // free
    cX=cmat_free(LDX,nX,cX);
    cY=cmat_free(LDY,nY,cY);
    cZ=cmat_free(LDZ,nZ,cZ);
  }else if(char_eq(mode,"two") && char_eq(cmd,"mtimes")){
    /*
     * z=x*y
     */
    if(nrhs>5){ ERROR("MATLAB:cmulti_data:maxrhs","Too many input arguments."); }
    if(nlhs>1){ ERROR("MATLAB:cmulti_data:maxlhs","Too many output arguments."); }
    if(!mxIsStruct(prhs[3])){ ERROR("MATLAB:cmulti_data","The arg4 should be Struct."); }
    if(!mxIsStruct(prhs[4])){ ERROR("MATLAB:cmulti_data","The arg5 should be Struct."); }
    // allocate and compute
    mX=mxGetM(prhs[3]); nX=mxGetN(prhs[3]); LDX=mX; cX=cmat_allocate_prec(LDX,nX,prec); mx2cmat(mX,nX,cX,LDX,prhs[3]);
    mY=mxGetM(prhs[4]); nY=mxGetN(prhs[4]); LDY=mY; cY=cmat_allocate_prec(LDY,nY,prec); mx2cmat(mY,nY,cY,LDY,prhs[4]);
    if     (mY==1 && nY==1){ mZ=mX; nZ=nX; LDZ=mZ; cZ=cmat_allocate_prec(LDZ,nZ,prec); cmat_mul_c(mZ,nZ,cZ,LDZ,cX,LDX,cY[0]); }
    else if(mX==1 && nX==1){ mZ=mY; nZ=nY; LDZ=mZ; cZ=cmat_allocate_prec(LDZ,nZ,prec); cmat_mul_c(mZ,nZ,cZ,LDZ,cY,LDY,cX[0]); }
    else if(nX==mY)        { mZ=mX; nZ=nY; LDZ=mZ; cZ=cmat_allocate_prec(LDZ,nZ,prec); cmat_prod(mX,nX,nY,cZ,LDZ,cX,LDX,cY,LDY); }
    else                   { ERROR("MATLAB:cmulti_data","Z=X*Y. The sizes of X,Y should be (m,n) and (n,k)."); }
    // done
    plhs[0]=cmat2mx(mZ,nZ,cZ,LDZ);
    // free
    cX=cmat_free(LDX,nX,cX);
    cY=cmat_free(LDY,nY,cY);
    cZ=cmat_free(LDZ,nZ,cZ);
  }else if(char_eq(mode,"two") && char_eq(cmd,"rdivide")){
    /*
     * z=x./y
     */
    if(nrhs>5){ ERROR("MATLAB:cmulti_data:maxrhs","Too many input arguments."); }
    if(nlhs>1){ ERROR("MATLAB:cmulti_data:maxlhs","Too many output arguments."); }
    if(!mxIsStruct(prhs[3])){ ERROR("MATLAB:cmulti_data","The arg4 should be Struct."); }
    if(!mxIsStruct(prhs[4])){ ERROR("MATLAB:cmulti_data","The arg5 should be Struct."); }
    // allocate and compute
    mX=mxGetM(prhs[3]); nX=mxGetN(prhs[3]); LDX=mX; cX=cmat_allocate_prec(LDX,nX,prec); mx2cmat(mX,nX,cX,LDX,prhs[3]);
    mY=mxGetM(prhs[4]); nY=mxGetN(prhs[4]); LDY=mY; cY=cmat_allocate_prec(LDY,nY,prec); mx2cmat(mY,nY,cY,LDY,prhs[4]);
    if     (mY==1  && nY==1) { mZ=mX; nZ=nX; LDZ=mZ; cZ=cmat_allocate_prec(LDZ,nZ,prec); cmat_div_c2(mZ,nZ,cZ,LDZ,cX,LDX,cY[0]);  }
    else if(mX==1  && nX==1) { mZ=mY; nZ=nY; LDZ=mZ; cZ=cmat_allocate_prec(LDZ,nZ,prec); cmat_div_c1(mZ,nZ,cZ,LDZ,cX[0],cY,LDY);  }
    else if(mX==mY && nX==nY){ mZ=mX; nZ=nX; LDZ=mZ; cZ=cmat_allocate_prec(LDZ,nZ,prec); cmat_div   (mZ,nZ,cZ,LDZ,cX,LDX,cY,LDY); }
    else                     { ERROR("MATLAB:cmulti_data","Z=X./Y. The sizes of X,Y should be same."); }
    // done
    plhs[0]=cmat2mx(mZ,nZ,cZ,LDZ);
    // free
    cX=cmat_free(LDX,nX,cX);
    cY=cmat_free(LDY,nY,cY);
    cZ=cmat_free(LDZ,nZ,cZ);
  }else if(char_eq(mode,"two") && char_eq(cmd,"mldivide")){
    /*
     * z=x\y
     */
    if(nrhs>5){ ERROR("MATLAB:cmulti_data:maxrhs","Too many input arguments."); }
    if(nlhs>1){ ERROR("MATLAB:cmulti_data:maxlhs","Too many output arguments."); }
    if(!mxIsStruct(prhs[3])){ ERROR("MATLAB:cmulti_data","The arg4 should be Struct."); }
    if(!mxIsStruct(prhs[4])){ ERROR("MATLAB:cmulti_data","The arg5 should be Struct."); }
    // allocate and compute
    mX=mxGetM(prhs[3]); nX=mxGetN(prhs[3]); LDX=mX; cX=cmat_allocate_prec(LDX,nX,prec); mx2cmat(mX,nX,cX,LDX,prhs[3]);
    mY=mxGetM(prhs[4]); nY=mxGetN(prhs[4]); LDY=mY; cY=cmat_allocate_prec(LDY,nY,prec); mx2cmat(mY,nY,cY,LDY,prhs[4]);
    // 
    if(mX!=nX || nX!=mY){ ERROR("MATLAB:cmulti_data","Z=X\\Y. The sizes of X,Y should be (n,n) and (n,k)."); }
    mZ=mY; nZ=nY; LDZ=mZ; cZ=cmat_allocate_prec(LDZ,nZ,prec); 
    cmat_copy(mZ,nZ,cZ,LDZ,cY,LDY); csolve(mZ,nZ,cZ,LDZ,cX,LDX,&info);
    // done
    plhs[0]=cmat2mx(mZ,nZ,cZ,LDZ);
    // free
    cX=cmat_free(LDX,nX,cX);
    cY=cmat_free(LDY,nY,cY);
    cZ=cmat_free(LDZ,nZ,cZ);
  }else{
    ERROR("MATLAB:cmulti_data","Illeagal command.");
  }

}
