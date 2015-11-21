/**
 * @breif double型行列から行列の生成
 */
void multi_set_d(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  multi *x=NULL,*y=NULL;
  if(nlhs>1){ mexErrMsgIdAndTxt("MATLAB:multi_mex:maxlhs","Too many output arguments."); }
  if(!(IS_DUBL(nrhs,prhs,2))){ MATLAB_ERROR("multi_set_d: The arg2 should be Double."); }
  // allocate
  x=multi_allocate_mxArray(prhs[2]);
       if(_T(x)=='d'){ y=multi_allocate('r',_M(x),_N(x),_L(x)); }
  else if(_T(x)=='z'){ y=multi_allocate('c',_M(x),_N(x),_L(x)); }
  else               { MATLAB_ERROR("multi_set_d: Unkown type"); }
  // set_d
       if(_T(y)=='r'){ rmat3_set_d(_M(y),_N(y),_L(y),_R(y),_LD1(y),_LD2(y),_D(x),_LD1(x),_LD2(x)); }
  else if(_T(y)=='c'){ cmat3_set_z(_M(y),_N(y),_L(y),_C(y),_LD1(y),_LD2(y),_Z(x),_LD1(x),_LD2(x)); }
  else               { MATLAB_ERROR("multi_set_d: Unkown type"); }
  // done
  plhs[0]=mxCreateStructMulti(y);
  x=multi_free(x);
  y=multi_free(y);
  return;
}

//EOF
