/**
 * @breif uminus
 */
void multi_uminus(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  multi *x=NULL,*y=NULL;
  if(nlhs>1){ mexErrMsgIdAndTxt("MATLAB:multi_mex:maxlhs","Too many output arguments."); }
  if(!(IS_STRT(nrhs,prhs,N0))){ MATLAB_ERROR("multi_uminus: The 1st-arg should be Struct."); }
  // allocate by clone
  x=multi_allocate_mxArray(prhs[N0]);
  // allocate by default precision
  y=multi_allocate(_T(x),_M(x),_N(x),_L(x));
  // uminus
       if(_T(y)=='r'){ rmat3_neg(_M(y),_N(y),_L(y),_R(y),_LD1(y),_LD2(y),_R(x),_LD1(x),_LD2(x)); }
  else if(_T(y)=='c'){ cmat3_neg(_M(y),_N(y),_L(y),_C(y),_LD1(y),_LD2(y),_C(x),_LD1(x),_LD2(x)); }
       //編集済み
  else if(_T(y)=='R'){ irmat3_neg(_M(y),_N(y),_L(y),_R0(y),_R1(y),_LD1(y),_LD2(y),_R0(x),_R1(x),_LD1(x),_LD2(x)); }
  else if(_T(y)=='C'){ icmat3_neg(_M(y),_N(y),_L(y),_C0(y),_C1(y),_LD1(y),_LD2(y),_C0(x),_C1(x),_LD1(x),_LD2(x)); }
  else{ MATLAB_ERROR("multi_uminus: Unkown type"); }
  // done
  plhs[0]=mxCreateStructMulti(y);
  x=multi_free(x);
  y=multi_free(y);
  return;
}

//EOF
