/**
 * @brief y=down(x)
 */
void multi_down(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  multi *x=NULL,*y=NULL;
  if(nlhs>1){ mexErrMsgIdAndTxt("MATLAB:multi_mex:maxlhs","Too many output arguments."); }
  if(!(IS_STRT(nrhs,prhs,N0))){ MATLAB_ERROR("multi_down: The 1st-arg should be Struct."); }
  // allocate by clone
  x=multi_allocate_mxArray(prhs[N0]);
  // allocate by default precision
       if(_T(x)=='r' || _T(x)=='R'){ y=multi_allocate('r',_M(x),_N(x),_L(x)); }
  else if(_T(x)=='c' || _T(x)=='C'){ y=multi_allocate('c',_M(x),_N(x),_L(x)); }
  // copy
       if(_T(x)=='r' && _T(y)=='r'){ rmat3_copy(_M(y),_N(y),_L(y),_R(y),_LD1(y),_LD2(y),_R(x), _LD1(x),_LD2(x)); }
  else if(_T(x)=='R' && _T(y)=='r'){ rmat3_copy(_M(y),_N(y),_L(y),_R(y),_LD1(y),_LD2(y),_R0(x),_LD1(x),_LD2(x)); }
  else if(_T(x)=='c' && _T(y)=='c'){ cmat3_copy(_M(y),_N(y),_L(y),_C(y),_LD1(y),_LD2(y),_C(x), _LD1(x),_LD2(x)); }
  else if(_T(x)=='C' && _T(y)=='c'){ cmat3_copy(_M(y),_N(y),_L(y),_C(y),_LD1(y),_LD2(y),_C0(x),_LD1(x),_LD2(x)); }
  else{ MATLAB_ERROR("multi_down: Unkown type"); }
  // done
  plhs[0]=mxCreateStructMulti(y);
  x=multi_free(x);
  y=multi_free(y);
  return;
}

//EOF
