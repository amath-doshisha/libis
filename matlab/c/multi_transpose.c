/**
 * @brief y=x.'
 */
void multi_transpose(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  array *x=NULL,*y=NULL;
  if(nlhs>1){ mexErrMsgIdAndTxt("MATLAB:multi_mex:maxlhs","Too many output arguments."); }
  if(!(IS_STRT(nrhs,prhs,N0))){ MATLAB_ERROR("multi_transpose: The 1st-arg should be Struct."); }
  // allocate by clone
  x=mxArray_to_array(prhs[N0]);
  // check size
  if(_L(x)!=1){ MATLAB_ERROR("multi_transpose: 3D-Array cannot be transposed."); }
  // transpose
  y=multi_allocate(_T(x),_N(x),_M(x),_L(x));
       if(_T(y)=='r'){ rmat_copy_t(_M(x),_N(x),_R(y),_LD1(y),_R(x),_LD1(x)); }
  else if(_T(y)=='c'){ cmat_copy_t(_M(x),_N(x),_C(y),_LD1(y),_C(x),_LD1(x)); }
       //追加
  else if(_T(y)=='R'){ irmat_copy_t(_M(x),_N(x),_R0(y),_LD1(y),_R1(y),_LD1(y),_R0(x),_LD1(x),_R1(x),_LD1(x)); }
       //else if(_T(y)=='c'){ cmat_copy_t(_M(x),_N(x),_C(y),_LD1(y),_C(x),_LD1(x)); }
       //ここまで
  else{ MATLAB_ERROR("multi_transpose: Unkown type"); }
  // done
  plhs[0]=array_to_mxArray(y);
  x=array_free(x);
  y=array_free(y);
  return;
}

//EOF
