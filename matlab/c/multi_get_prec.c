/**
 * @brief y=get_prec(x)
 */
void multi_get_prec(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  multi *x=NULL,*y=NULL;
  if(nlhs>1){ mexErrMsgIdAndTxt("MATLAB:multi_mex:maxlhs","Too many output arguments."); }
  if(!(IS_STRT(nrhs,prhs,N0))){ MATLAB_ERROR("multi_get_prec: The 1st-arg should be Struct."); }
  // allocate by clone
  x=multi_allocate_mxArray(prhs[N0]);
  // allocate by default precision
  y=multi_allocate('i',_M(x),_N(x),_L(x));
  // prec
       if(_T(x)=='r'){ rmat3_get_prec(_M(y),_N(y),_L(y),_I(y),_LD1(y),_LD2(y),_R(x), _LD1(x),_LD2(x)); }
  else if(_T(x)=='c'){ cmat3_get_prec(_M(y),_N(y),_L(y),_I(y),_LD1(y),_LD2(y),_C(x), _LD1(x),_LD2(x)); }
  else if(_T(x)=='R'){ rmat3_get_prec(_M(y),_N(y),_L(y),_I(y),_LD1(y),_LD2(y),_R0(x),_LD1(x),_LD2(x)); }
  else if(_T(x)=='C'){ cmat3_get_prec(_M(y),_N(y),_L(y),_I(y),_LD1(y),_LD2(y),_C0(x),_LD1(x),_LD2(x)); }
  else{ MATLAB_ERROR("multi_get_prec: Unkown type"); }
  // done
  plhs[0]=mxCreateStructMulti(y);
  x=multi_free(x);
  y=multi_free(y);
  return;
}

//EOF
