/**
 * @brief y=get_exp(x)
 */
void multi_get_exp(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  array *x=NULL,*y=NULL;
  if(nlhs>1){ mexErrMsgIdAndTxt("MATLAB:multi_mex:maxlhs","Too many output arguments."); }
  if(!(IS_STRT(nrhs,prhs,N0))){ MATLAB_ERROR("multi_get_exp: The 1st-arg should be Struct."); }
  // allocate by clone
  x=mxArray_to_array(prhs[N0]);
  // allocate by default precision
  y=multi_allocate('i',_M(x),_N(x),_L(x));
  // exp
  if(_T(x)=='r'){ rmat3_get_exp(_M(y),_N(y),_L(y),_I(y),_LD1(y),_LD2(y),_R(x),_LD1(x),_LD2(x)); }
  else{ MATLAB_ERROR("multi_get_exp: Unkown type"); }
  // done
  plhs[0]=array_to_mxArray(y);
  x=array_free(x);
  y=array_free(y);
  return;
}

//EOF
