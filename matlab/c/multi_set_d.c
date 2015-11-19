/**
 * @breif double型行列から行列の生成
 */
void multi_set_d(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  char type='r';
  int m=1,n=1,l=1,i,j,k,t;
  multi *A=NULL;
  const mxArray *x=NULL;
  dcomplex z;
  double d;
  if(!(IS_DUBL(nrhs,prhs,2))){ MATLAB_ERROR("multi_set_d: The arg2 should be Double."); }
  x=prhs[2];
  // get type
  if(mxIsComplex(x)){ type='c'; }
  // get size
  m=mxGetM(x);
  n=mxGetN(x);
  if(mxGetNumberOfDimensions(x)>=3){ l=mxGetDimensions(x)[2]; }
  if(mxGetNumberOfDimensions(x)>=4){ MATLAB_ERROR("multi_set_d: Number of dimensions is too big."); }  
  // allocate
  A=multi_allocate(type,m,n,l);
  // set
  for(k=0; k<_L(A); k++){
    for(j=0; j<_N(A); j++){
      for(i=0; i<_M(A); i++){
	t=i+j*_M(A)+k*_M(A)*_N(A);
             if(_T(A)=='r'){ d=mxGetPr(x)[t];                      rset_d(MAT3(_R(A),i,j,k,_LD1(A),_LD2(A)),d); }
	else if(_T(A)=='c'){ Z_SET(z,mxGetPr(x)[t],mxGetPi(x)[t]); cset_z(MAT3(_C(A),i,j,k,_LD1(A),_LD2(A)),z); }
	else{ MATLAB_ERROR("multi_set_d: Unknow type."); }
      }
    }
  }
  // done
  plhs[0]=mxCreateStructMulti(A);
  if(nlhs>1){ mexErrMsgIdAndTxt("MATLAB:multi_mex:maxlhs","Too many output arguments."); }
  A=multi_free(A);
  return;
}

//EOF
