/**
 * @brief 行列の生成 y=imulti(x) or y=imulti(x0,x1)
 */
void multi_get_imulti(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  array *x=NULL,*x0=NULL,*x1=NULL,*y=NULL;
  if(nlhs>1){ mexErrMsgIdAndTxt("MATLAB:multi_mex:maxlhs","Too many output arguments."); }
  if(nrhs-N0>=3){ mexErrMsgIdAndTxt("MATLAB:multi_mex:maxrhs","Too many input arguments."); }
  // Case: y=imulti(x)
  if(nrhs-N0==1){
    // input
    x=mxArray_to_array(prhs[N0]);
    // y=get_imulti(x)
    y=array_get_imulti(x);
    // output
    plhs[0]=array_to_mxArray(y);
    x=array_free(x);
    y=array_free(y);
    return;
  }
  // Case: y=get_imulti(x,y)
  if(nrhs-N0==2){
    // input
    x0=mxArray_to_array(prhs[N0]);
    x1=mxArray_to_array(prhs[N0+1]);
    // y=get_imulti2(x)
    y=array_get_imulti2(x0,x1);
    // output
    plhs[0]=array_to_mxArray(y);
    x0=array_free(x0);
    x1=array_free(x1);
    y=array_free(y);    
    return;
  }
}

//EOF
