/**
 * @brief 行列の生成 y=multi(x)
 */
void multi_get(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  char type='r',*buf=NULL;
  array *x=NULL,*y=NULL;
  // input-1
  if(nlhs>1){ mexErrMsgIdAndTxt("MATLAB:multi_mex:maxlhs","Too many output arguments."); }
  x=mxArray_to_array(prhs[N0]);
  // input-2
  if(!(IS_CHAR(nrhs,prhs,N0+1))){ MATLAB_ERROR("multi_get: The 2nd-arg should be Char."); }
  if(IS_CHAR(nrhs,prhs,N0+1)){
    buf=mxArrayToString(prhs[N0+1]);
    type=buf[0];
  }
  // y=get(type,x)
  y=array_get(type,x);
  if(y==NULL){ MATLAB_ERROR("Not supported type."); }
  // output
  plhs[0]=array_to_mxArray(y);
  x=array_free(x);
  y=array_free(y);
  if(buf!=NULL){ mxFree(buf); }
  return;
}

//EOF
