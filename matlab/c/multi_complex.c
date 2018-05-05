/**
 * @brief multi型行列から複素数行列の生成 y=complex(x) or z=complex(x,y)
 */
void multi_complex(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  array *x=NULL,*y=NULL,*z=NULL;
  if(nlhs>1){ mexErrMsgIdAndTxt("MATLAB:multi_mex:maxlhs","Too many output arguments."); }
  if(nrhs-N0>=3){ mexErrMsgIdAndTxt("MATLAB:multi_mex:maxrhs","Too many input arguments."); }
  // Case: y=complex(x)
  if(nrhs-N0==1){
    // input
    x=mxArray_to_array(prhs[N0]);
    // y=get_complex(x)
    y=array_get_complex(x);
    // output
    plhs[0]=array_to_mxArray(y);
    x=array_free(x);
    y=array_free(y);
    return;
  }
  // Case: z=complex(x,y)
  if(nrhs-N0==2){
    // input
    x=mxArray_to_array(prhs[N0]);
    y=mxArray_to_array(prhs[N0+1]);
    // y=get_complex(x)
    z=array_get_complex2(x,y);
    // output
    plhs[0]=array_to_mxArray(z);
    x=array_free(x);
    y=array_free(y);    
    z=array_free(z);
  }
}

//EOF
