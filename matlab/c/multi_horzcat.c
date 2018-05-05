/**
 * @brief y=[x1 x2 ...]
 */
void multi_horzcat(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  char type='r';
  int num,i=0,m,n,l;
  array **x=NULL,*y=NULL;
  if(nlhs>1){ mexErrMsgIdAndTxt("MATLAB:multi_mex:maxlhs","Too many output arguments."); }
  // get number of args
  num=nrhs-N0;
  // allocate
  x=malloc(num*sizeof(array*));
  for(i=0; i<num; i++){
    if(!(IS_STRT(nrhs,prhs,N0+i))){ MATLAB_ERROR("multi_plus: The arg should be Struct."); }
    x[i]=mxArray_to_array(prhs[N0+i]);
  }
  // check size
  type=_T(x[0]);
  m=_M(x[0]);
  n=_N(x[0]);
  l=_L(x[0]);  
  for(i=1; i<num; i++){
    if(_T(x[i])=='c'){ type='c'; }    
    if(m!=_M(x[i])){ MATLAB_ERROR("multi_horzcat: The first dimension of all arguments should be same."); }
    if(l!=_L(x[i])){ MATLAB_ERROR("multi_horzcat: The third dimension of all arguments should be same."); }
    n=n+_N(x[i]);
  }
  // allocate
  y=multi_allocate(type,m,n,l);
  // concat
  for(i=0, n=0; i<num; i++){
         if(_T(y)=='r' && _T(x[i])=='r'){ rmat3_copy      (_M(x[i]),_N(x[i]),_L(x[i]),&MAT3(_R(y),0,n,0,_LD1(x[i]),_LD2(x[i])),_LD1(y),_LD2(y),_R(x[i]),_LD1(x[i]),_LD2(x[i])); }
    else if(_T(y)=='c' && _T(x[i])=='r'){ cmat3_copy_rmat3(_M(x[i]),_N(x[i]),_L(x[i]),&MAT3(_C(y),0,n,0,_LD1(x[i]),_LD2(x[i])),_LD1(y),_LD2(y),_R(x[i]),_LD1(x[i]),_LD2(x[i])); }
    else if(_T(y)=='c' && _T(x[i])=='c'){ cmat3_copy      (_M(x[i]),_N(x[i]),_L(x[i]),&MAT3(_C(y),0,n,0,_LD1(x[i]),_LD2(x[i])),_LD1(y),_LD2(y),_C(x[i]),_LD1(x[i]),_LD2(x[i])); }
    else{ MATLAB_ERROR("multi_horzcat: Unkown type"); }
    n=n+_N(x[i]);    
  }
  // done
  plhs[0]=array_to_mxArray(y);
  // free
  for(i=0; i<num; i++){ x[i]=array_free(x[i]); }
  free(x); x=NULL;
  y=array_free(y);
  return;
}

//EOF
