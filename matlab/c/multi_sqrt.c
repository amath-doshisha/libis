/**
 * @brief y=sqrt(x)
 */
void multi_sqrt(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  array *x=NULL,*y=NULL;
  // intput
  if(nlhs>1){ mexErrMsgIdAndTxt("MATLAB:multi_mex:maxlhs","Too many output arguments."); }
  x=mxArray_to_array(prhs[N0]);
  // y=sqrt(x)
  y=array_sqrt(x);
  // output
  plhs[0]=array_to_mxArray(y);
  x=array_free(x);
  y=array_free(y);
  return;
}

/*
void multi_sqrt(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  array *x=NULL,*y=NULL;
  int info=-1;
  if(nlhs>1){ mexErrMsgIdAndTxt("MATLAB:multi_mex:maxlhs","Too many output arguments."); }
  if(!(IS_STRT(nrhs,prhs,N0))){ MATLAB_ERROR("multi_sqrt: The 1st-arg should be Struct."); }
  // allocate by clone
  x=mxArray_to_array(prhs[N0]);
  // operation
  // 正のみ
       if(_T(x)=='r'){ y=multi_allocate('r',_M(x),_N(x),_L(x));  rmat3_sqrt(_M(y),_N(y),_L(y),_R(y),_LD1(y),_LD2(y),_R(x),_LD1(x),_LD2(x)); }
  else if(_T(x)=='R'){ y=multi_allocate('R',_M(x),_N(x),_L(x)); irmat3_sqrt(_M(y),_N(y),_L(y),_R0(y),_R1(y),_LD1(y),_LD2(y),_R0(x),_R1(x),_LD1(x),_LD2(x)); }
  else{ MATLAB_ERROR("multi_sqrt: Unkown type"); }
  // done
  plhs[0]=array_to_mxArray(y);
  x=array_free(x);
  y=array_free(y);
  return;
}
*/

//EOF
