/**
 * @brief y=mag(x)
 */
// y=magnitude([x])=max(abs(x0),abs(x1))
// ir,ic型をrmulti型で返す
void multi_mag(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  multi *x=NULL,*y=NULL,*u=NULL;
  int info=-1;
  if(nlhs>1){ mexErrMsgIdAndTxt("MATLAB:multi_mex:maxlhs","Too many output arguments."); }
  if(!(IS_STRT(nrhs,prhs,N0))){ MATLAB_ERROR("multi_mag: The 1st-arg should be Struct."); }
  // allocate by clone
  x=multi_allocate_mxArray(prhs[N0]);
  // operation
       if(_T(x)=='R'){ y=multi_allocate('r',_M(x),_N(x),_L(x)); u=multi_allocate('R',_M(x),_N(x),_L(x)); irmat3_abs(_M(u),_N(u),_L(u),_R0(u),_R1(u),_LD1(u),_LD2(u),_R0(x),_R1(x),_LD1(x),_LD2(x)); rmat3_copy(_M(y),_N(y),_L(y),_R(y),_LD1(y),_LD2(y),_R1(u),_LD1(u),_LD2(u)); }
  else if(_T(x)=='C'){ y=multi_allocate('r',_M(x),_N(x),_L(x)); u=multi_allocate('R',_M(x),_N(x),_L(x)); icmat3_abs(_M(u),_N(u),_L(u),_R0(u),_R1(u),_LD1(u),_LD2(u),_C0(x),_C1(x),_LD1(x),_LD2(x)); rmat3_copy(_M(y),_N(y),_L(y),_R(y),_LD1(y),_LD2(y),_R1(u),_LD1(u),_LD2(u)); }
  else{ MATLAB_ERROR("multi_mag: Unkown type"); }
  // done
  plhs[0]=mxCreateStructMulti(y);
  x=multi_free(x);
  y=multi_free(y);
  u=multi_free(u);
  return;
}

//EOF
