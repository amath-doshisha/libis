//追加
/**
 * @brief [C0,C1]=imulti(A,B)
 */
void multi_icopy2(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  multi *x0=NULL,*x1=NULL,*y=NULL;
  if(nlhs>1){ mexErrMsgIdAndTxt("MATLAB:multi_mex:maxlhs","Too many output arguments."); }
  if(nrhs>N0+2){ mexErrMsgIdAndTxt("MATLAB:multi_mex:maxrhs","Too many input arguments."); }
  if(!(IS_STRT(nrhs,prhs,N0))){ MATLAB_ERROR("multi_icopy2: The 1st-arg should be Struct."); }
  if(!(IS_STRT(nrhs,prhs,N0+1))){ MATLAB_ERROR("multi_icopy2: The 2nd-arg should be Struct."); }
  // allocate by clone
  x0=multi_allocate_mxArray(prhs[N0]);
  x1=multi_allocate_mxArray(prhs[N0+1]);
  if(!(_M(x0)==_M(x1) && _N(x0)==_N(x1) && _L(x0)==_L(x1))){ MATLAB_ERROR("multi_icopy2(x0,x1): size(x0)!=size(x1) error!"); }
  // allocate by default precision
       if(_T(x0)=='r' && _T(x1)=='r'){ y=multi_allocate('R',_M(x0),_N(x0),_L(x0)); }
  else if(_T(x0)=='c' && _T(x1)=='c'){ y=multi_allocate('C',_M(x0),_N(x0),_L(x0)); }
  else{ MATLAB_ERROR("multi_icopy2: Unkown type"); }
  // copy
       if(_T(y)=='R'){ irmat3_copy(_M(y),_N(y),_L(y),_R0(y),_R1(y),_LD1(y),_LD2(y),_R(x0),_R(x1),_LD1(x0),_LD2(x0)); }
  else if(_T(y)=='C'){ icmat3_copy(_M(y),_N(y),_L(y),_C0(y),_C1(y),_LD1(y),_LD2(y),_C(x0),_C(x1),_LD1(x0),_LD2(x0)); }
  else{ MATLAB_ERROR("multi_icopy2: Unkown type"); }
  // done
  plhs[0]=mxCreateStructMulti(y);
  x0=multi_free(x0);
  x1=multi_free(x1);
  y=multi_free(y);
  return;
}

//EOF
