/**
 * @brief z=plus(x,y)
 */
void multi_plus(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  array *x=NULL,*y=NULL,*z=NULL;
  // input
  if(nlhs>1){ mexErrMsgIdAndTxt("MATLAB:multi_mex:maxlhs","Too many output arguments."); }
  if(nrhs-N0!=2){ MATLAB_ERROR("The function z=multi_plus(x,y) needs two input arguments."); }
  x=mxArray_to_array(prhs[N0]);
  y=mxArray_to_array(prhs[N0+1]);
  // z=plus(x,y)
  z=array_add(x,y);
  // error
  if(z==NULL){
    MATLAB_ERROR("z=multi_plus(x,y): not supported size.");
  }
  // output
  plhs[0]=array_to_mxArray(z);  
  x=array_free(x);
  y=array_free(y);
  z=array_free(z);
  return;
}

//EOF
