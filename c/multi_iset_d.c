/**
 * @breif double型行列から行列の生成 [y0,y1]=multi(x)
 */
void multi_iset_d(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  multi *x=NULL,*y=NULL;
  if(nlhs>1){ mexErrMsgIdAndTxt("MATLAB:multi_mex:maxlhs","Too many output arguments."); }
  if(!(IS_DUBL(nrhs,prhs,N0))){ MATLAB_ERROR("multi_iset_d: The 1st-arg should be Double."); }
  // allocate
  x=multi_allocate_mxArray(prhs[N0]);
       if(_T(x)=='d'){ y=multi_allocate('R',_M(x),_N(x),_L(x)); }
  else if(_T(x)=='z'){ y=multi_allocate('C',_M(x),_N(x),_L(x)); }
  else               { MATLAB_ERROR("multi_iset_d: Unkown type of x"); }
  // iset_d
       if(_T(y)=='R'){ irmat3_set_d(_M(y),_N(y),_L(y),_R0(y),_R1(y),_LD1(y),_LD2(y),_D(x),_LD1(x),_LD2(x)); }
  else if(_T(y)=='C'){ icmat3_set_z(_M(y),_N(y),_L(y),_C0(y),_C1(y),_LD1(y),_LD2(y),_Z(x),_LD1(x),_LD2(x)); }
  else               { MATLAB_ERROR("multi_iset_d: Unkown type of y"); }
  // done
  plhs[0]=mxCreateStructMulti(y);
  x=multi_free(x);
  y=multi_free(y);
  return;
}

//EOF
