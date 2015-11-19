/**
 * @breif 1行列の生成
 */
void multi_set_ones(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  char type='r',*buf=NULL;
  int m=1,n=1,l=1;
  multi *A=NULL;
  if(!(IS_CHAR(nrhs,prhs,2))){ MATLAB_ERROR("multi_set_zeros: The arg2 should be Char."); }
  if(IS_CHAR(nrhs,prhs,2)){ buf=mxArrayToString(prhs[2]); type=buf[0]; }
  if(IS_NUMR(nrhs,prhs,3)){ m=GET_DOUBLE(prhs[3])[0]; }
  if(IS_NUMR(nrhs,prhs,4)){ n=GET_DOUBLE(prhs[4])[0]; }
  if(IS_NUMR(nrhs,prhs,5)){ l=GET_DOUBLE(prhs[5])[0]; }
  A=multi_allocate(type,m,n,l);
       if(_T(A)=='r'){ rmat3_set_ones(_M(A),_N(A),_L(A),_R(A),_LD1(A),_LD2(A)); }
  else if(_T(A)=='c'){ cmat3_set_ones(_M(A),_N(A),_L(A),_C(A),_LD1(A),_LD2(A)); }
  else if(_T(A)=='R'){ rmat3_set_ones(_M(A),_N(A),_L(A),_R0(A),_LD1(A),_LD2(A)); rmat3_set_ones(_M(A),_N(A),_L(A),_R1(A),_LD1(A),_LD2(A)); }
  else if(_T(A)=='C'){ cmat3_set_ones(_M(A),_N(A),_L(A),_C0(A),_LD1(A),_LD2(A)); cmat3_set_ones(_M(A),_N(A),_L(A),_C1(A),_LD1(A),_LD2(A)); }
  else{ MATLAB_ERROR("multi_set_zeros: Unkown type"); }
  plhs[0]=mxCreateStructMulti(A);
  if(nlhs>1){ mexErrMsgIdAndTxt("MATLAB:multi_mex:maxlhs","Too many output arguments."); }
  A=multi_free(A);
  if(buf!=NULL){ mxFree(buf); }
  return;
}

//EOF
