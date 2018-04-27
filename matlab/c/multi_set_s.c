
int multi_s_isreal(multi *x)
{
  int i,j,k,ret=1;
  if(_T(x)!='s'){ MATLAB_ERROR("multi_s_isreal: Unkown type"); }
  for(k=0; k<_L(x); k++){
    for(j=0; j<_N(x); j++){
      for(i=0; i<_M(x); i++){
	if(str_has_any_char(MAT3(_S(x),i,j,k,_LD1(x),_LD2(x)),"iIjJ")){ ret=0; };
      }
    }
  }
  return ret;
}

/**
 * @brief 文字列のセル型行列から行列の生成 y=multi(x)
 */
void multi_set_s(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  multi *x=NULL,*y=NULL;
  if(nlhs>1){ mexErrMsgIdAndTxt("MATLAB:multi_mex:maxlhs","Too many output arguments."); }
  if(!(IS_CELL(nrhs,prhs,N0))){ MATLAB_ERROR("multi_set_s: The 1st-arg should be Cell."); }
  // allocate
  x=multi_allocate_mxArray(prhs[N0]);
  if(_T(x)!='s'){ MATLAB_ERROR("multi_set_s: Unkown type"); }
  if(multi_s_isreal(x)){ y=multi_allocate('r',_M(x),_N(x),_L(x)); }
  else                 { y=multi_allocate('c',_M(x),_N(x),_L(x)); }
  // set_s
       if(_T(y)=='r'){ rmat3_set_s(_M(y),_N(y),_L(y),_R(y),_LD1(y),_LD2(y),_S(x),_LD1(x),_LD2(x)); }
  else if(_T(y)=='c'){ cmat3_set_s(_M(y),_N(y),_L(y),_C(y),_LD1(y),_LD2(y),_S(x),_LD1(x),_LD2(x)); }
  else               { MATLAB_ERROR("multi_set_s: Unkown type"); }
  // done
  plhs[0]=mxCreateStructMulti(y);
  x=multi_free(x);
  y=multi_free(y);
  return;
}

//EOF
