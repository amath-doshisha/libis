/**
 * @brief chare型のセルへのキャスト y=get_s(x,format,digits)
 */
void multi_get_s(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  int i,digits;
  char *format=NULL;
  mwSize ndim=3,*dim=NULL;
  mxArray *s=NULL;
  array *x=NULL,*y=NULL;
  // [input-1] array x
  if(nlhs>1){ mexErrMsgIdAndTxt("MATLAB:multi_mex:maxlhs","Too many output arguments."); }
  if(!(IS_STRT(nrhs,prhs,N0))){ MATLAB_ERROR("multi_get_s: The 1st-arg should be Struct."); }
  x=mxArray_to_array(prhs[N0]);
  // [input-2] char format
  if(!(IS_CHAR(nrhs,prhs,N0+1))){ MATLAB_ERROR("multi_get_s: The 2nd-arg should be Char."); }
  format=mxArrayToString(prhs[N0+1]);  
  // [input-3] int digits
  if(!IS_NUMR(nrhs,prhs,N0+2)){ MATLAB_ERROR("multi_get_s: The 3rd-arg should be scalar."); }
  digits=GET_DOUBLE(prhs[N0+2])[0];
  // y=get_char(x,format,digits)
  y=array_get_char(x,format[0],digits);  
  // output
  ndim=ARRAY_NDIM(y);
  dim=(mwSize*)mxMalloc(sizeof(mwSize)*ndim);
  for(i=0; i<ndim; i++){ dim[i]=ARRAY_DIM(y,i); } 
  plhs[0]=mxCreateCellArray(ndim,dim);
  for(i=0; i<ARRAY_SIZE(y); i++){
    s=mxCreateString(ARRAY_SVEC(y,i));
    mxSetCell(plhs[0],i,s);
  }
  // done
  if(format!=NULL){ mxFree(format); }
  if(dim!=NULL){ mxFree(dim); }
  x=array_free(x);
  y=array_free(y);
}

//EOF
