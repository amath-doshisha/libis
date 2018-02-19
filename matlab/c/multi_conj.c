/**
 * @breif y=conj(x)
 */
void multi_conj(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  multi *x=NULL,*y=NULL;
  int info=-1;
  if(nlhs>1){ mexErrMsgIdAndTxt("MATLAB:multi_mex:maxlhs","Too many output arguments."); }
  if(!(IS_STRT(nrhs,prhs,N0))){ MATLAB_ERROR("multi_conj: The 1st-arg should be Struct."); }
  // allocate by clone
  x=multi_allocate_mxArray(prhs[N0]);
  // operation
       if(_T(x)=='r'){ y=multi_allocate('r',_M(x),_N(x),_L(x)); rmat3_copy(_M(y),_N(y),_L(y),_R(y),_LD1(y),_LD2(y),_R(x),_LD1(x),_LD2(x)); }
  else if(_T(x)=='c'){ y=multi_allocate('c',_M(x),_N(x),_L(x)); cmat3_conj(_M(y),_N(y),_L(y),_C(y),_LD1(y),_LD2(y),_C(x),_LD1(x),_LD2(x)); }
  else if(_T(x)=='R'){ y=multi_allocate('R',_M(x),_N(x),_L(x)); irmat3_copy(_M(y),_N(y),_L(y),_R0(y),_R1(y),_LD1(y),_LD2(y),_R0(x),_R1(x),_LD1(x),_LD2(x)); }
  else if(_T(x)=='C'){ y=multi_allocate('C',_M(x),_N(x),_L(x)); icmat3_conj(_M(y),_N(y),_L(y),_C0(y),_C1(y),_LD1(y),_LD2(y),_C0(x),_C1(x),_LD1(x),_LD2(x)); }
  else{ MATLAB_ERROR("multi_conj: Unkown type"); }
  // done
  plhs[0]=mxCreateStructMulti(y);
  x=multi_free(x);
  y=multi_free(y);
  return;
}

//EOF
