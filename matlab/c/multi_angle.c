
/**
 * @brief y=angle(x)
 */
void multi_angle(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  array *x=NULL,*y=NULL;
  int info=-1;
  if(nlhs>1){ mexErrMsgIdAndTxt("MATLAB:multi_mex:maxlhs","Too many output arguments."); }
  if(!(IS_STRT(nrhs,prhs,N0))){ MATLAB_ERROR("multi_angle: The 1st-arg should be Struct."); }
  // allocate by clone
  x=mxArray_to_array(prhs[N0]);
  y=multi_allocate('r',_M(x),_N(x),_L(x)); 
  // operation
       if(_T(x)=='r'){ rmat3_set_zeros(_M(y),_N(y),_L(y),_R(y),_LD1(y),_LD2(y)); }
  else if(_T(x)=='c'){ cmat3_arg(_M(y),_N(y),_L(y),_R(y),_LD1(y),_LD2(y),_C(x),_LD1(x),_LD2(x)); }
  else{ MATLAB_ERROR("multi_angle: Unkown type"); }
  // done
  plhs[0]=array_to_mxArray(y);
  x=array_free(x);
  y=array_free(y);
  return;
}

//EOF
