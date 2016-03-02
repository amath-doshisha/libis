/**
 * @breif 文字列のセル型行列から行列の生成 y=multi(x)
 */
void multi_set_s(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  multi *x=NULL,*y=NULL;
  if(nlhs>1){ mexErrMsgIdAndTxt("MATLAB:multi_mex:maxlhs","Too many output arguments."); }
  if(!(IS_CELL(nrhs,prhs,N0))){ MATLAB_ERROR("multi_set_s: The 1st-arg should be Cell."); }
  // allocate
  x=multi_allocate_mxArray(prhs[N0]);
  if(_T(x)=='s'){ y=multi_allocate('r',_M(x),_N(x),_L(x)); }
  else          { MATLAB_ERROR("multi_set_s: Unkown type"); }
  // set_s
  rmat3_set_s(_M(y),_N(y),_L(y),_R(y),_LD1(y),_LD2(y),_S(x),_LD1(x),_LD2(x));
  // done
  plhs[0]=mxCreateStructMulti(y);
  x=multi_free(x);
  y=multi_free(y);
  return;
}

//EOF
