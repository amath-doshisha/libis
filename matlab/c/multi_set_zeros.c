/**
 * @brief 零行列の生成 x=zeros(type,size)
 */
void multi_set_zeros(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  char type='r',*buf=NULL;
  int i,ndim,*dim=NULL;
  array *x=NULL;
  // input-1
  if(nlhs>1){ mexErrMsgIdAndTxt("MATLAB:multi_mex:maxlhs","Too many output arguments."); }
  // input-2
  if(!(IS_CHAR(nrhs,prhs,N0))){ MATLAB_ERROR("multi_set_zeros: The 1st-arg should be Char."); }
  if(IS_CHAR(nrhs,prhs,N0)){ buf=mxArrayToString(prhs[N0]); type=buf[0]; }
  // input-3
  if(nrhs-N0==2 && IS_NUMR(nrhs,prhs,N0+1)){
    // only a number
    ndim=2;
    dim=ivec_allocate(ndim);    
    dim[0]=GET_DOUBLE(prhs[N0+1])[0];
    dim[1]=dim[0];
  }else if(nrhs-N0==2 && IS_DUBL(nrhs,prhs,N0+1) && IS_ROW(nrhs,prhs,N0+1)){
    // only a row vector
    ndim=mxGetNumberOfDimensions(prhs[N0+1]);
    dim=ivec_allocate(ndim);
    for(i=0; i<ndim; i++){ dim[i]=GET_DOUBLE(prhs[N0+1])[i]; }
  }else if(nrhs-N0>=3 && IS_NUMR(nrhs,prhs,N0+1)){
    ndim=nrhs-N0-1;
    dim=ivec_allocate(ndim);    
    for(i=0; i<ndim; i++){
      if(!IS_NUMR(nrhs,prhs,N0+1+i)){ MATLAB_ERROR("multi_set_zeros: The 2nd, 3rd, ... args should be scalar."); }
      dim[i]=GET_DOUBLE(prhs[N0+1+i])[0];
    }
  }else{ MATLAB_ERROR("multi_set_zeros: The 2nd-arg should be scalar Double or row double vector."); }  
  // set_zeros(x,type,size)
  x=array_allocate(type,ndim,dim);
  array_set_zeros(x);
  // output
  plhs[0]=array_to_mxArray(x);
  x=array_free(x);
  dim=ivec_free(dim);
  if(buf!=NULL){ mxFree(buf); }
  return;
}

//EOF
