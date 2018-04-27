/**
 * @brief double型行列から行列の生成 [y0,y1]=imulti(x0,x1)
 */
void multi_iset_dd(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  multi *x0=NULL,*x1=NULL,*y=NULL;
  if(nlhs>1){ mexErrMsgIdAndTxt("MATLAB:multi_mex:maxlhs","Too many output arguments."); }
  if(nrhs>N0+2){ mexErrMsgIdAndTxt("MATLAB:multi_mex:maxrhs","Too many input arguments."); }
  if(!(IS_DUBL(nrhs,prhs,N0))){ MATLAB_ERROR("multi_iset_dd: The 1st-arg should be Double."); }
  if(!(IS_DUBL(nrhs,prhs,N0+1))){ MATLAB_ERROR("multi_iset_dd: The 2nd-arg should be Double."); }
  // allocate
  x0=multi_allocate_mxArray(prhs[N0]);
  x1=multi_allocate_mxArray(prhs[N0+1]);
  if(!(_M(x0)==_M(x1) && _N(x0)==_N(x1) && _L(x0)==_L(x1))){ MATLAB_ERROR("multi_iset_dd(x0,x1): size(x0)!=size(x1) error!"); }
  
  // iset_dd
  if(_T(x0)=='d' && _T(x1)=='d'){ y=multi_allocate('R',_M(x0),_N(x0),_L(x0)); irmat3_set_dd(_M(y),_N(y),_L(y),_R0(y),_R1(y),_LD1(y),_LD2(y),_D(x0),_LD1(x0),_LD2(x0),_D(x1),_LD1(x1),_LD2(x1)); }
       //  else if(_T(x0)=='z' && _T(x1)=='d'){ y=multi_allocate('C',_M(x0),_N(x0),_L(x0));    }
       //  else if(_T(x0)=='d' && _T(x1)=='z'){ y=multi_allocate('C',_M(x0),_N(x0),_L(x0));    }
  else if(_T(x0)=='z' && _T(x1)=='z'){ y=multi_allocate('C',_M(x0),_N(x0),_L(x0)); icmat3_set_zz(_M(y),_N(y),_L(y),_C0(y),_C1(y),_LD1(y),_LD2(y),_Z(x0),_LD1(x0),_LD2(x0),_Z(x1),_LD1(x1),_LD2(x1));   }
  else{ MATLAB_ERROR("multi_iset_dd: error"); }
  // done
  plhs[0]=mxCreateStructMulti(y);
  x0=multi_free(x0);
  x1=multi_free(x1);
  y=multi_free(y);
  return;
}

//EOF
