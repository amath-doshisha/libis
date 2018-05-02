/**
 * @brief y=get_exp10(x) or y=get_exp10(x,offset)
 */
void multi_get_exp10(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  multi *x=NULL,*y=NULL;
  int offset=0;
  if(nlhs>1){ mexErrMsgIdAndTxt("MATLAB:multi_mex:maxlhs","Too many output arguments."); }
  if(nrhs>N0+2){ mexErrMsgIdAndTxt("MATLAB:multi_mex:maxrhs","Too many input arguments."); }
  if(!(IS_STRT(nrhs,prhs,N0))){ MATLAB_ERROR("multi_get_exp10: The 1st-arg should be Struct."); }
  if(!(IS_NUMR(nrhs,prhs,N0+1))){ MATLAB_ERROR("multi_get_exp10: The 2nd-arg should be scalar Double."); }
  // allocate by clone
  x=multi_allocate_mxArray(prhs[N0]);
  // operation
  if(IS_NUMR(nrhs,prhs,N0+1)){ offset=GET_DOUBLE(prhs[N0+1])[0];}
  // allocate by clone
       if(_T(x)=='r'){ y=multi_allocate('r',_M(x),_N(x),_L(x)); rmat3_exp10_floor_log10_abs_sub(_M(x),_N(x),_L(x),_R(y),_LD1(y),_LD2(y),_R(x),_LD1(x),_LD2(x),offset); }
  else if(_T(x)=='c'){ y=multi_allocate('c',_M(x),_N(x),_L(x)); cmat3_exp10_floor_log10_abs_sub(_M(x),_N(x),_L(x),_C(y),_LD1(y),_LD2(y),_C(x),_LD1(x),_LD2(x),offset); }
  else{ MATLAB_ERROR("multi_get_exp10: error"); }
  // done
  plhs[0]=mxCreateStructMulti(y);
  x=multi_free(x);
  y=multi_free(y);
  return;
}
  //EOF

