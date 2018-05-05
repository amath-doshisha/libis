/**
 * @brief z=mrdivide(x,y)
 */
void multi_mrdivide(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  array *x=NULL,*y=NULL,*z=NULL,*u=NULL,*v=NULL;
  int info=-1;
  if(nlhs>1){ mexErrMsgIdAndTxt("MATLAB:multi_mex:maxlhs","Too many output arguments."); }
  if(!(IS_STRT(nrhs,prhs,N0))){ MATLAB_ERROR("multi_mrdivide: The 1st-arg should be Struct."); }
  if(!(IS_STRT(nrhs,prhs,N0+1))){ MATLAB_ERROR("multi_mrdivide: The 2nd-arg should be Struct."); }
  // allocate by clone
  x=mxArray_to_array(prhs[N0]);
  y=mxArray_to_array(prhs[N0+1]);
  // operation
  if(_M(y)==1 && _N(y)==1 && _L(y)==1){
    // In the case of y is scalar
         if(_T(x)=='r' && _T(y)=='r'){ z=multi_allocate('r',_M(x),_N(x),_L(x)); rmat3_div_r2(_M(z),_N(z),_L(z),_R(z),_LD1(z),_LD2(z),_R(x),_LD1(x),_LD2(x),MAT3(_R(y),0,0,0,_LD1(y),_LD2(y))); }
    else if(_T(x)=='r' && _T(y)=='c'){ z=multi_allocate('c',_M(x),_N(x),_L(x)); rmat3_div_c2(_M(z),_N(z),_L(z),_C(z),_LD1(z),_LD2(z),_R(x),_LD1(x),_LD2(x),MAT3(_C(y),0,0,0,_LD1(y),_LD2(y))); }
    else if(_T(x)=='c' && _T(y)=='r'){ z=multi_allocate('c',_M(x),_N(x),_L(x)); cmat3_div_r2(_M(z),_N(z),_L(z),_C(z),_LD1(z),_LD2(z),_C(x),_LD1(x),_LD2(x),MAT3(_R(y),0,0,0,_LD1(y),_LD2(y))); }
    else if(_T(x)=='c' && _T(y)=='c'){ z=multi_allocate('c',_M(x),_N(x),_L(x)); cmat3_div_c2(_M(z),_N(z),_L(z),_C(z),_LD1(z),_LD2(z),_C(x),_LD1(x),_LD2(x),MAT3(_C(y),0,0,0,_LD1(y),_LD2(y))); }
	 //追加
    else if(_T(x)=='R' && _T(y)=='R'){ z=multi_allocate('R',_M(x),_N(x),_L(x));irmat3_div_r2(_M(z),_N(z),_L(z),_R0(z),_R1(z),_LD1(z),_LD2(z),_R0(x),_R1(x),_LD1(x),_LD2(x),MAT3(_R0(y),0,0,0,_LD1(y),_LD2(y)),MAT3(_R1(y),0,0,0,_LD1(y),_LD2(y))); }
    else if(_T(x)=='r' && _T(y)=='R'){ z=multi_allocate('R',_M(x),_N(x),_L(x));irmat3_div_r2(_M(z),_N(z),_L(z),_R0(z),_R1(z),_LD1(z),_LD2(z),_R(x), _R(x), _LD1(x),_LD2(x),MAT3(_R0(y),0,0,0,_LD1(y),_LD2(y)),MAT3(_R1(y),0,0,0,_LD1(y),_LD2(y))); }
    else if(_T(x)=='R' && _T(y)=='r'){ z=multi_allocate('R',_M(x),_N(x),_L(x));irmat3_div_r2(_M(z),_N(z),_L(z),_R0(z),_R1(z),_LD1(z),_LD2(z),_R0(x),_R1(x),_LD1(x),_LD2(x),MAT3(_R(y),0,0,0, _LD1(y),_LD2(y)),MAT3(_R(y),0,0,0, _LD1(y),_LD2(y))); }	 
    else if(_T(x)=='C' && _T(y)=='C'){ z=multi_allocate('C',_M(x),_N(x),_L(x));icmat3_div_c2(_M(z),_N(z),_L(z),_C0(z),_C1(z),_LD1(z),_LD2(z),_C0(x),_C1(x),_LD1(x),_LD2(x),MAT3(_C0(y),0,0,0,_LD1(y),_LD2(y)),MAT3(_C1(y),0,0,0,_LD1(y),_LD2(y))); }
    else if(_T(x)=='c' && _T(y)=='C'){ z=multi_allocate('C',_M(x),_N(x),_L(x));icmat3_div_c2(_M(z),_N(z),_L(z),_C0(z),_C1(z),_LD1(z),_LD2(z),_C(x), _C(x), _LD1(x),_LD2(x),MAT3(_C0(y),0,0,0,_LD1(y),_LD2(y)),MAT3(_C1(y),0,0,0,_LD1(y),_LD2(y))); }
    else if(_T(x)=='C' && _T(y)=='c'){ z=multi_allocate('C',_M(x),_N(x),_L(x));icmat3_div_c2(_M(z),_N(z),_L(z),_C0(z),_C1(z),_LD1(z),_LD2(z),_C0(x),_C1(x),_LD1(x),_LD2(x),MAT3(_C(y),0,0,0, _LD1(y),_LD2(y)),MAT3(_C(y),0,0,0, _LD1(y),_LD2(y))); }	 
    else if(_T(x)=='R' && _T(y)=='C'){ z=multi_allocate('C',_M(x),_N(x),_L(x));irmat3_div_c2(_M(z),_N(z),_L(z),_C0(z),_C1(z),_LD1(z),_LD2(z),_R0(x),_R1(x),_LD1(x),_LD2(x),MAT3(_C0(y),0,0,0,_LD1(y),_LD2(y)),MAT3(_C1(y),0,0,0,_LD1(y),_LD2(y))); }
    else if(_T(x)=='r' && _T(y)=='C'){ z=multi_allocate('C',_M(x),_N(x),_L(x));irmat3_div_c2(_M(z),_N(z),_L(z),_C0(z),_C1(z),_LD1(z),_LD2(z),_R(x), _R(x), _LD1(x),_LD2(x),MAT3(_C0(y),0,0,0,_LD1(y),_LD2(y)),MAT3(_C1(y),0,0,0,_LD1(y),_LD2(y))); }
    else if(_T(x)=='R' && _T(y)=='c'){ z=multi_allocate('C',_M(x),_N(x),_L(x));irmat3_div_c2(_M(z),_N(z),_L(z),_C0(z),_C1(z),_LD1(z),_LD2(z),_R0(x),_R1(x),_LD1(x),_LD2(x),MAT3(_C(y),0,0,0, _LD1(y),_LD2(y)),MAT3(_C(y),0,0,0, _LD1(y),_LD2(y))); }	 
    else if(_T(x)=='C' && _T(y)=='R'){ z=multi_allocate('C',_M(x),_N(x),_L(x));icmat3_div_r2(_M(z),_N(z),_L(z),_C0(z),_C1(z),_LD1(z),_LD2(z),_C0(x),_C1(x),_LD1(x),_LD2(x),MAT3(_R0(y),0,0,0,_LD1(y),_LD2(y)),MAT3(_R1(y),0,0,0,_LD1(y),_LD2(y))); }
    else if(_T(x)=='c' && _T(y)=='R'){ z=multi_allocate('C',_M(x),_N(x),_L(x));icmat3_div_r2(_M(z),_N(z),_L(z),_C0(z),_C1(z),_LD1(z),_LD2(z),_C(x), _C(x), _LD1(x),_LD2(x),MAT3(_R0(y),0,0,0,_LD1(y),_LD2(y)),MAT3(_R1(y),0,0,0,_LD1(y),_LD2(y))); }
    else if(_T(x)=='C' && _T(y)=='r'){ z=multi_allocate('C',_M(x),_N(x),_L(x));icmat3_div_r2(_M(z),_N(z),_L(z),_C0(z),_C1(z),_LD1(z),_LD2(z),_C0(x),_C1(x),_LD1(x),_LD2(x),MAT3(_R(y),0,0,0, _LD1(y),_LD2(y)),MAT3(_R(y),0,0,0, _LD1(y),_LD2(y))); }
	 //ここまで
    else{ MATLAB_ERROR("multi_rdivide: Unkown type"); }
  }else if(_M(x)==1 && _N(x)==1 && _L(x)==1){
    // In the case of x is scalar
    MATLAB_ERROR("multi_mrdivide: z=x/y: Dimensions of x and y are NOT same.");
  }else if(_M(y)==_N(y) && _L(y)==1 && _M(y)==_N(x) && _M(x)==1 && _L(x)==1){
    // In the case where y is square matrix and x's colums is same as y's row
         if(_T(x)=='r' && _T(y)=='r'){ u=multi_allocate('r',_N(y),_M(y),1); v=multi_allocate('r',_N(x),_M(x),1); z=multi_allocate('r',_M(x),_N(x),1); rmat_copy_t     (_M(y),_N(y),_R(u),_LD1(u),_R(y),_LD1(y)); rmat_copy_t     (_M(x),_N(x),_R(v),_LD1(v),_R(x),_LD1(x)); rsolve(_M(v),_N(v),_R(v),_LD1(v),_R(u),_LD1(u),&info); rmat_copy_t(_M(v),_N(v),_R(z),_LD1(z),_R(v),_LD1(v)); }
    else if(_T(x)=='r' && _T(y)=='c'){ u=multi_allocate('c',_N(y),_M(y),1); v=multi_allocate('c',_N(x),_M(x),1); z=multi_allocate('c',_M(x),_N(x),1); cmat_copy_t     (_M(y),_N(y),_C(u),_LD1(u),_C(y),_LD1(y)); cmat_copy_rmat_t(_M(x),_N(x),_C(v),_LD1(v),_R(x),_LD1(x)); csolve(_M(v),_N(v),_C(v),_LD1(v),_C(u),_LD1(u),&info); cmat_copy_t(_M(v),_N(v),_C(z),_LD1(z),_C(v),_LD1(v)); }
    else if(_T(x)=='c' && _T(y)=='r'){ u=multi_allocate('c',_N(y),_M(y),1); v=multi_allocate('c',_N(x),_M(x),1); z=multi_allocate('c',_M(x),_N(x),1); cmat_copy_rmat_t(_M(y),_N(y),_C(u),_LD1(u),_R(y),_LD1(y)); cmat_copy_t     (_M(x),_N(x),_C(v),_LD1(v),_C(x),_LD1(x)); csolve(_M(v),_N(v),_C(v),_LD1(v),_C(u),_LD1(u),&info); cmat_copy_t(_M(v),_N(v),_C(z),_LD1(z),_C(v),_LD1(v)); }
    else if(_T(x)=='c' && _T(y)=='c'){ u=multi_allocate('c',_N(y),_M(y),1); v=multi_allocate('c',_N(x),_M(x),1); z=multi_allocate('c',_M(x),_N(x),1); cmat_copy_t     (_M(y),_N(y),_C(u),_LD1(u),_C(y),_LD1(y)); cmat_copy_t     (_M(x),_N(x),_C(v),_LD1(v),_C(x),_LD1(x)); csolve(_M(v),_N(v),_C(v),_LD1(v),_C(u),_LD1(u),&info); cmat_copy_t(_M(v),_N(v),_C(z),_LD1(z),_C(v),_LD1(v)); }
    else{ MATLAB_ERROR("multi_mtimes: Unkown type"); }
    if(info!=0){ MATLAB_ERROR("multi_mtimes: csolve error."); }
  }else{
    MATLAB_ERROR("multi_mrdivide: z=x/y: This function is not implemented");
  }
  // done
  plhs[0]=array_to_mxArray(z);
  x=array_free(x);
  y=array_free(y);
  z=array_free(z);
  u=array_free(u);
  v=array_free(v);
  return;
}

//EOF


