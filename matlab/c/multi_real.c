/**
 * @brief y=real(x)
 */
void multi_real(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  array *x=NULL,*y=NULL;
  // input
  if(nlhs>1){ mexErrMsgIdAndTxt("MATLAB:multi_mex:maxlhs","Too many output arguments."); }
  x=mxArray_to_array(prhs[N0]);
  // y=get_real(x)
  y=array_get_real(x);
  // output
  plhs[0]=array_to_mxArray(y);
  x=array_free(x);
  y=array_free(y);
  return;
}

//EOF
