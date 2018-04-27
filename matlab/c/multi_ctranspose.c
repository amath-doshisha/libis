/**
 * @brief y=x'
 */
void multi_ctranspose(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  multi *x=NULL,*y=NULL;
  if(nlhs>1){ mexErrMsgIdAndTxt("MATLAB:multi_mex:maxlhs","Too many output arguments."); }
  if(!(IS_STRT(nrhs,prhs,N0))){ MATLAB_ERROR("multi_ctranspose: The 1st-arg should be Struct."); }
  // allocate by clone
  x=multi_allocate_mxArray(prhs[N0]);
  // check size
  if(_L(x)!=1){ MATLAB_ERROR("multi_ctranspose: 3D-Array cannot be ctransposed."); }
  // ctranspose
  y=multi_allocate(_T(x),_N(x),_M(x),_L(x));
       if(_T(y)=='r'){ rmat_copy_t (_M(x),_N(x),_R(y),_LD1(y),_R(x),_LD1(x)); }
  else if(_T(y)=='c'){ cmat_copy_ct(_M(x),_N(x),_C(y),_LD1(y),_C(x),_LD1(x)); }
  else if(_T(y)=='R'){ irmat_copy_t(_M(x),_N(x),_R0(y),_LD1(y),_R1(y),_LD1(y),_R0(x),_LD1(x),_R1(x),_LD1(x)); }
  else if(_T(y)=='C'){ icmat_copy_ct(_M(x),_N(x),_C0(y),_LD1(y),_C1(y),_LD1(y),_C0(x),_LD1(x),_C1(x),_LD1(x)); }
  else{ MATLAB_ERROR("multi_ctranspose: Unkown type"); }
  // done
  plhs[0]=mxCreateStructMulti(y);
  x=multi_free(x);
  y=multi_free(y);
  return;
}

//EOF
