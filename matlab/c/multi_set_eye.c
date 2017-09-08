/**
 * @breif 1行列の生成 y=ones(m,n,l)
 */
void multi_set_eye(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  char type='r',*buf=NULL;
  int m=1,n=1,l=1;
  multi *A=NULL,*s=NULL;
  if(nlhs>1){ mexErrMsgIdAndTxt("MATLAB:multi_mex:maxlhs","Too many output arguments."); }
  if(!(IS_CHAR(nrhs,prhs,N0))){ MATLAB_ERROR("multi_set_eye: The 1st-arg should be Char."); }
  if(IS_CHAR(nrhs,prhs,N0)){ buf=mxArrayToString(prhs[N0]); type=buf[0]; }
  if(IS_DUBL(nrhs,prhs,N0+1) && IS_ROW(nrhs,prhs,N0+1)){
    s=multi_allocate_mxArray(prhs[N0+1]);
    if(_N(s)>0){ m=_D(s)[0]; }
    if(_N(s)>1){ n=_D(s)[1]; }
    if(_N(s)>2){ l=_D(s)[2]; }
  }else if(IS_NUMR(nrhs,prhs,N0+1)){
    m=GET_DOUBLE(prhs[N0+1])[0];
  }else{ MATLAB_ERROR("multi_set_eye: The 2nd-arg should be scalar Double or row double vector."); }  
  if(IS_NUMR(nrhs,prhs,N0+2)){ n=GET_DOUBLE(prhs[N0+2])[0]; }
  if(IS_NUMR(nrhs,prhs,N0+3)){ l=GET_DOUBLE(prhs[N0+3])[0]; }
  if(nrhs-N0==2 && IS_NUMR(nrhs,prhs,N0+1)){ n=m; l=1; }
  if(l!=1){ MATLAB_ERROR("multi_set_eye: This function support only a matrix."); }
  A=multi_allocate(type,m,n,l);
       if(_T(A)=='r'){ rmat_set_eye(_M(A),_N(A),_R(A),_LD1(A)); }
  else if(_T(A)=='c'){ cmat_set_eye(_M(A),_N(A),_C(A),_LD1(A)); }
  else if(_T(A)=='R'){ rmat_set_eye(_M(A),_N(A),_R0(A),_LD1(A));
                       rmat_set_eye(_M(A),_N(A),_R1(A),_LD1(A)); }
  else if(_T(A)=='C'){ cmat_set_eye(_M(A),_N(A),_C0(A),_LD1(A));
                       cmat_set_eye(_M(A),_N(A),_C1(A),_LD1(A)); }
  else{ MATLAB_ERROR("multi_set_eye: Unkown type"); }
  plhs[0]=mxCreateStructMulti(A);
  A=multi_free(A);
  s=multi_free(s);
  if(buf!=NULL){ mxFree(buf); }
  return;
}

//EOF
