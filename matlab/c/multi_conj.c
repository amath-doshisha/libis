/**
 * @brief y=conj(x)
 */
void multi_conj(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  array *x=NULL,*y=NULL;
  // intput
  if(nlhs>1){ mexErrMsgIdAndTxt("MATLAB:multi_mex:maxlhs","Too many output arguments."); }
  x=mxArray_to_array(prhs[N0]);
  // y=conj(x)
  y=array_conj(x);
  // output
  plhs[0]=array_to_mxArray(y);
  x=array_free(x);
  y=array_free(y);
  return;
}

//EOF
