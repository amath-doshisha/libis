/**
 * @brief double型へのキャスト y=double(x)
 */
void multi_get_d(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  array *x=NULL,*y=NULL;
  // input
  if(nlhs>1){ mexErrMsgIdAndTxt("MATLAB:multi_mex:maxlhs","Too many output arguments."); }
  if(!(IS_STRT(nrhs,prhs,N0))){ MATLAB_ERROR("multi_get_d: The 1st-arg should be Struct."); }
  x=mxArray_to_array(prhs[N0]);
  // y=get_double(x)
  y=array_get_double(x);
  // output
  plhs[0]=array_to_mxArray(y);  
  x=array_free(x);
  y=array_free(y);
}

//EOF
