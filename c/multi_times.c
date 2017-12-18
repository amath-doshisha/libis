/**
 * @breif z=times(x,y)
 */
void multi_times(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  multi *x=NULL,*y=NULL,*z=NULL;
  if(nlhs>1){ mexErrMsgIdAndTxt("MATLAB:multi_mex:maxlhs","Too many output arguments."); }
  if(!(IS_STRT(nrhs,prhs,N0))){ MATLAB_ERROR("multi_times: The 1st-arg should be Struct."); }
  if(!(IS_STRT(nrhs,prhs,N0+1))){ MATLAB_ERROR("multi_times: The 2nd-arg should be Struct."); }
  // allocate by clone
  x=multi_allocate_mxArray(prhs[N0]);
  y=multi_allocate_mxArray(prhs[N0+1]);
  // operation
  if(_M(x)==_M(y) && _N(x)==_N(y) && _L(x)==_L(y)){  // In the case of same sizes
         if(_T(x)=='r' && _T(y)=='r'){ z=multi_allocate('r',_M(x),_N(x),_L(x)); rmat3_mul      (_M(z),_N(z),_L(z),_R(z),_LD1(z),_LD2(z),_R(x),_LD1(x),_LD2(x),_R(y),_LD1(y),_LD2(y)); }
    else if(_T(x)=='r' && _T(y)=='c'){ z=multi_allocate('c',_M(x),_N(x),_L(x)); cmat3_mul_rmat3(_M(z),_N(z),_L(z),_C(z),_LD1(z),_LD2(z),_C(y),_LD1(y),_LD2(y),_R(x),_LD1(x),_LD2(x)); }
    else if(_T(x)=='c' && _T(y)=='r'){ z=multi_allocate('c',_M(x),_N(x),_L(x)); cmat3_mul_rmat3(_M(z),_N(z),_L(z),_C(z),_LD1(z),_LD2(z),_C(x),_LD1(x),_LD2(x),_R(y),_LD1(y),_LD2(y)); }
    else if(_T(x)=='c' && _T(y)=='c'){ z=multi_allocate('c',_M(x),_N(x),_L(x)); cmat3_mul      (_M(z),_N(z),_L(z),_C(z),_LD1(z),_LD2(z),_C(x),_LD1(x),_LD2(x),_C(y),_LD1(y),_LD2(y)); }
	 //編集済み
    else if(_T(x)=='R' && _T(y)=='R'){ z=multi_allocate('R',_M(x),_N(x),_L(x)); irmat3_mul      (_M(z),_N(z),_L(z),_R0(z),_R1(z),_LD1(z),_LD2(z),_R0(x),_R1(x),_LD1(x),_LD2(x),_R0(y),_R1(y),_LD1(y),_LD2(y)); }
    else if(_T(x)=='r' && _T(y)=='R'){ z=multi_allocate('R',_M(x),_N(x),_L(x)); irmat3_mul      (_M(z),_N(z),_L(z),_R0(z),_R1(z),_LD1(z),_LD2(z),_R(x), _R(x), _LD1(x),_LD2(x),_R0(y),_R1(y),_LD1(y),_LD2(y)); }
    else if(_T(x)=='R' && _T(y)=='r'){ z=multi_allocate('R',_M(x),_N(x),_L(x)); irmat3_mul      (_M(z),_N(z),_L(z),_R0(z),_R1(z),_LD1(z),_LD2(z),_R0(x),_R1(x),_LD1(x),_LD2(x),_R(y), _R(y), _LD1(y),_LD2(y)); }	 
    else if(_T(x)=='C' && _T(y)=='C'){ z=multi_allocate('C',_M(x),_N(x),_L(x)); icmat3_mul      (_M(z),_N(z),_L(z),_C0(z),_C1(z),_LD1(z),_LD2(z),_C0(x),_C1(x),_LD1(x),_LD2(x),_C0(y),_C1(y),_LD1(y),_LD2(y)); }
    else if(_T(x)=='c' && _T(y)=='C'){ z=multi_allocate('C',_M(x),_N(x),_L(x)); icmat3_mul      (_M(z),_N(z),_L(z),_C0(z),_C1(z),_LD1(z),_LD2(z),_C(x), _C(x), _LD1(x),_LD2(x),_C0(y),_C1(y),_LD1(y),_LD2(y)); }
    else if(_T(x)=='C' && _T(y)=='c'){ z=multi_allocate('C',_M(x),_N(x),_L(x)); icmat3_mul      (_M(z),_N(z),_L(z),_C0(z),_C1(z),_LD1(z),_LD2(z),_C0(x),_C1(x),_LD1(x),_LD2(x),_C(y), _C(y), _LD1(y),_LD2(y)); }	 
    else if(_T(x)=='R' && _T(y)=='C'){ z=multi_allocate('C',_M(x),_N(x),_L(x)); icmat3_mul_rmat3(_M(z),_N(z),_L(z),_C0(z),_C1(z),_LD1(z),_LD2(z),_C0(y),_C1(y),_LD1(y),_LD2(y),_R0(x),_R1(x),_LD1(x),_LD2(x)); }
    else if(_T(x)=='r' && _T(y)=='C'){ z=multi_allocate('C',_M(x),_N(x),_L(x)); icmat3_mul_rmat3(_M(z),_N(z),_L(z),_C0(z),_C1(z),_LD1(z),_LD2(z),_C0(y),_C1(y),_LD1(y),_LD2(y),_R(x), _R(x), _LD1(x),_LD2(x)); }
    else if(_T(x)=='R' && _T(y)=='c'){ z=multi_allocate('C',_M(x),_N(x),_L(x)); icmat3_mul_rmat3(_M(z),_N(z),_L(z),_C0(z),_C1(z),_LD1(z),_LD2(z),_C(y), _C(y), _LD1(y),_LD2(y),_R0(x),_R1(x),_LD1(x),_LD2(x)); }	 
    else if(_T(x)=='C' && _T(y)=='R'){ z=multi_allocate('C',_M(x),_N(x),_L(x)); icmat3_mul_rmat3(_M(z),_N(z),_L(z),_C0(z),_C1(z),_LD1(z),_LD2(z),_C0(x),_C1(x),_LD1(x),_LD2(x),_R0(y),_R1(y),_LD1(y),_LD2(y)); }
    else if(_T(x)=='c' && _T(y)=='R'){ z=multi_allocate('C',_M(x),_N(x),_L(x)); icmat3_mul_rmat3(_M(z),_N(z),_L(z),_C0(z),_C1(z),_LD1(z),_LD2(z),_C(x), _C(x), _LD1(x),_LD2(x),_R0(y),_R1(y),_LD1(y),_LD2(y)); }
    else if(_T(x)=='C' && _T(y)=='r'){ z=multi_allocate('C',_M(x),_N(x),_L(x)); icmat3_mul_rmat3(_M(z),_N(z),_L(z),_C0(z),_C1(z),_LD1(z),_LD2(z),_C0(x),_C1(x),_LD1(x),_LD2(x),_R(y), _R(y), _LD1(y),_LD2(y)); }
    else{ MATLAB_ERROR("multi_times: Unkown type"); }
  }else if(_M(x)==1 && _N(x)==1 && _L(x)==1){    // In the case of x is scalar    
         if(_T(x)=='r' && _T(y)=='r'){ z=multi_allocate('r',_M(y),_N(y),_L(y)); rmat3_mul_r(_M(z),_N(z),_L(z),_R(z),_LD1(z),_LD2(z),_R(y),_LD1(y),_LD2(y),MAT3(_R(x),0,0,0,_LD1(x),_LD2(x))); }
    else if(_T(x)=='r' && _T(y)=='c'){ z=multi_allocate('c',_M(y),_N(y),_L(y)); cmat3_mul_r(_M(z),_N(z),_L(z),_C(z),_LD1(z),_LD2(z),_C(y),_LD1(y),_LD2(y),MAT3(_R(x),0,0,0,_LD1(x),_LD2(x))); }
    else if(_T(x)=='c' && _T(y)=='r'){ z=multi_allocate('c',_M(y),_N(y),_L(y)); rmat3_mul_c(_M(z),_N(z),_L(z),_C(z),_LD1(z),_LD2(z),_R(y),_LD1(y),_LD2(y),MAT3(_C(x),0,0,0,_LD1(x),_LD2(x))); }
    else if(_T(x)=='c' && _T(y)=='c'){ z=multi_allocate('c',_M(y),_N(y),_L(y)); cmat3_mul_c(_M(z),_N(z),_L(z),_C(z),_LD1(z),_LD2(z),_C(y),_LD1(y),_LD2(y),MAT3(_C(x),0,0,0,_LD1(x),_LD2(x))); }
    else if(_T(x)=='R' && _T(y)=='R'){ z=multi_allocate('R',_M(y),_N(y),_L(y)); irmat3_mul_r(_M(z),_N(z),_L(z),_R0(z),_R1(z),_LD1(z),_LD2(z),_R0(y),_R1(y),_LD1(y),_LD2(y),MAT3(_R0(x),0,0,0,_LD1(x),_LD2(x)),MAT3(_R1(x),0,0,0,_LD1(x),_LD2(x))); }
    else if(_T(x)=='r' && _T(y)=='R'){ z=multi_allocate('R',_M(y),_N(y),_L(y)); irmat3_mul_r(_M(z),_N(z),_L(z),_R0(z),_R1(z),_LD1(z),_LD2(z),_R0(y),_R1(y),_LD1(y),_LD2(y),MAT3(_R(x),0,0,0, _LD1(x),_LD2(x)),MAT3(_R(x),0,0,0, _LD1(x),_LD2(x))); }
    else if(_T(x)=='R' && _T(y)=='r'){ z=multi_allocate('R',_M(y),_N(y),_L(y)); irmat3_mul_r(_M(z),_N(z),_L(z),_R0(z),_R1(z),_LD1(z),_LD2(z),_R(y), _R(y), _LD1(y),_LD2(y),MAT3(_R0(x),0,0,0,_LD1(x),_LD2(x)),MAT3(_R1(x),0,0,0,_LD1(x),_LD2(x))); }	 
    else if(_T(x)=='C' && _T(y)=='C'){ z=multi_allocate('C',_M(y),_N(y),_L(y)); icmat3_mul_c(_M(z),_N(z),_L(z),_C0(z),_C1(z),_LD1(z),_LD2(z),_C0(y),_C1(y),_LD1(y),_LD2(y),MAT3(_C0(x),0,0,0,_LD1(x),_LD2(x)),MAT3(_C1(x),0,0,0,_LD1(x),_LD2(x))); }
    else if(_T(x)=='c' && _T(y)=='C'){ z=multi_allocate('C',_M(y),_N(y),_L(y)); icmat3_mul_c(_M(z),_N(z),_L(z),_C0(z),_C1(z),_LD1(z),_LD2(z),_C0(y),_C1(y),_LD1(y),_LD2(y),MAT3(_C(x),0,0,0, _LD1(x),_LD2(x)),MAT3(_C(x),0,0,0, _LD1(x),_LD2(x))); }
    else if(_T(x)=='C' && _T(y)=='c'){ z=multi_allocate('C',_M(y),_N(y),_L(y)); icmat3_mul_c(_M(z),_N(z),_L(z),_C0(z),_C1(z),_LD1(z),_LD2(z),_C(y), _C(y), _LD1(y),_LD2(y),MAT3(_C0(x),0,0,0,_LD1(x),_LD2(x)),MAT3(_C1(x),0,0,0,_LD1(x),_LD2(x))); }	 
    else if(_T(x)=='R' && _T(y)=='C'){ z=multi_allocate('C',_M(y),_N(y),_L(y)); icmat3_mul_r(_M(z),_N(z),_L(z),_C0(z),_C1(z),_LD1(z),_LD2(z),_C0(y),_C1(y),_LD1(y),_LD2(y),MAT3(_R0(x),0,0,0,_LD1(x),_LD2(x)),MAT3(_R1(x),0,0,0,_LD1(x),_LD2(x))); }
    else if(_T(x)=='r' && _T(y)=='C'){ z=multi_allocate('C',_M(y),_N(y),_L(y)); icmat3_mul_r(_M(z),_N(z),_L(z),_C0(z),_C1(z),_LD1(z),_LD2(z),_C0(y),_C1(y),_LD1(y),_LD2(y),MAT3(_R(x),0,0,0, _LD1(x),_LD2(x)),MAT3(_R(x),0,0,0, _LD1(x),_LD2(x))); }
    else if(_T(x)=='R' && _T(y)=='c'){ z=multi_allocate('C',_M(y),_N(y),_L(y)); icmat3_mul_r(_M(z),_N(z),_L(z),_C0(z),_C1(z),_LD1(z),_LD2(z),_C(y), _C(y), _LD1(y),_LD2(y),MAT3(_R0(x),0,0,0,_LD1(x),_LD2(x)),MAT3(_R1(x),0,0,0,_LD1(x),_LD2(x))); }	 
    else if(_T(x)=='C' && _T(y)=='R'){ z=multi_allocate('C',_M(y),_N(y),_L(y)); irmat3_mul_c(_M(z),_N(z),_L(z),_C0(z),_C1(z),_LD1(z),_LD2(z),_R0(y),_R1(y),_LD1(y),_LD2(y),MAT3(_C0(x),0,0,0,_LD1(x),_LD2(x)),MAT3(_C1(x),0,0,0,_LD1(x),_LD2(x))); }
    else if(_T(x)=='c' && _T(y)=='R'){ z=multi_allocate('C',_M(y),_N(y),_L(y)); irmat3_mul_c(_M(z),_N(z),_L(z),_C0(z),_C1(z),_LD1(z),_LD2(z),_R0(y),_R1(y),_LD1(y),_LD2(y),MAT3(_C(x),0,0,0, _LD1(x),_LD2(x)),MAT3(_C(x),0,0,0, _LD1(x),_LD2(x))); }
    else if(_T(x)=='C' && _T(y)=='r'){ z=multi_allocate('C',_M(y),_N(y),_L(y)); irmat3_mul_c(_M(z),_N(z),_L(z),_C0(z),_C1(z),_LD1(z),_LD2(z),_R(y), _R(y), _LD1(y),_LD2(y),MAT3(_C0(x),0,0,0,_LD1(x),_LD2(x)),MAT3(_C1(x),0,0,0,_LD1(x),_LD2(x))); }
    else{ MATLAB_ERROR("multi_times: Unkown type"); }
  }else if(_M(y)==1 && _N(y)==1 && _L(y)==1){    // In the case of y is scalar
         if(_T(x)=='r' && _T(y)=='r'){ z=multi_allocate('r',_M(x),_N(x),_L(x)); rmat3_mul_r(_M(z),_N(z),_L(z),_R(z),_LD1(z),_LD2(z),_R(x),_LD1(x),_LD2(x),MAT3(_R(y),0,0,0,_LD1(y),_LD2(y))); }
    else if(_T(x)=='r' && _T(y)=='c'){ z=multi_allocate('c',_M(x),_N(x),_L(x)); rmat3_mul_c(_M(z),_N(z),_L(z),_C(z),_LD1(z),_LD2(z),_R(x),_LD1(x),_LD2(x),MAT3(_C(y),0,0,0,_LD1(y),_LD2(y))); }
    else if(_T(x)=='c' && _T(y)=='r'){ z=multi_allocate('c',_M(x),_N(x),_L(x)); cmat3_mul_r(_M(z),_N(z),_L(z),_C(z),_LD1(z),_LD2(z),_C(x),_LD1(x),_LD2(x),MAT3(_R(y),0,0,0,_LD1(y),_LD2(y))); }
    else if(_T(x)=='c' && _T(y)=='c'){ z=multi_allocate('c',_M(x),_N(x),_L(x)); cmat3_mul_c(_M(z),_N(z),_L(z),_C(z),_LD1(z),_LD2(z),_C(x),_LD1(x),_LD2(x),MAT3(_C(y),0,0,0,_LD1(y),_LD2(y))); }
    else if(_T(x)=='R' && _T(y)=='R'){ z=multi_allocate('R',_M(x),_N(x),_L(x)); irmat3_mul_r(_M(z),_N(z),_L(z),_R0(z),_R1(z),_LD1(z),_LD2(z),_R0(x),_R1(x),_LD1(x),_LD2(x),MAT3(_R0(y),0,0,0,_LD1(y),_LD2(y)),MAT3(_R1(y),0,0,0,_LD1(y),_LD2(y))); }
    else if(_T(x)=='r' && _T(y)=='R'){ z=multi_allocate('R',_M(x),_N(x),_L(x)); irmat3_mul_r(_M(z),_N(z),_L(z),_R0(z),_R1(z),_LD1(z),_LD2(z),_R(x), _R(x), _LD1(x),_LD2(x),MAT3(_R0(y),0,0,0,_LD1(y),_LD2(y)),MAT3(_R1(y),0,0,0,_LD1(y),_LD2(y))); }
    else if(_T(x)=='R' && _T(y)=='r'){ z=multi_allocate('R',_M(x),_N(x),_L(x)); irmat3_mul_r(_M(z),_N(z),_L(z),_R0(z),_R1(z),_LD1(z),_LD2(z),_R0(x),_R1(x),_LD1(x),_LD2(x),MAT3(_R(y),0,0,0, _LD1(y),_LD2(y)),MAT3(_R(y),0,0,0, _LD1(y),_LD2(y))); }	 
    else if(_T(x)=='C' && _T(y)=='C'){ z=multi_allocate('C',_M(x),_N(x),_L(x)); icmat3_mul_c(_M(z),_N(z),_L(z),_C0(z),_C1(z),_LD1(z),_LD2(z),_C0(x),_C1(x),_LD1(x),_LD2(x),MAT3(_C0(y),0,0,0,_LD1(y),_LD2(y)),MAT3(_C1(y),0,0,0,_LD1(y),_LD2(y))); }
    else if(_T(x)=='c' && _T(y)=='C'){ z=multi_allocate('C',_M(x),_N(x),_L(x)); icmat3_mul_c(_M(z),_N(z),_L(z),_C0(z),_C1(z),_LD1(z),_LD2(z),_C(x), _C(x), _LD1(x),_LD2(x),MAT3(_C0(y),0,0,0,_LD1(y),_LD2(y)),MAT3(_C1(y),0,0,0,_LD1(y),_LD2(y))); }
    else if(_T(x)=='C' && _T(y)=='c'){ z=multi_allocate('C',_M(x),_N(x),_L(x)); icmat3_mul_c(_M(z),_N(z),_L(z),_C0(z),_C1(z),_LD1(z),_LD2(z),_C0(x),_C1(x),_LD1(x),_LD2(x),MAT3(_C(y),0,0,0, _LD1(y),_LD2(y)),MAT3(_C(y),0,0,0, _LD1(y),_LD2(y))); }	 
    else if(_T(x)=='R' && _T(y)=='C'){ z=multi_allocate('C',_M(x),_N(x),_L(x)); irmat3_mul_c(_M(z),_N(z),_L(z),_C0(z),_C1(z),_LD1(z),_LD2(z),_R0(x),_R1(x),_LD1(x),_LD2(x),MAT3(_C0(y),0,0,0,_LD1(y),_LD2(y)),MAT3(_C1(y),0,0,0,_LD1(y),_LD2(y))); }
    else if(_T(x)=='r' && _T(y)=='C'){ z=multi_allocate('C',_M(x),_N(x),_L(x)); irmat3_mul_c(_M(z),_N(z),_L(z),_C0(z),_C1(z),_LD1(z),_LD2(z),_R(x), _R(x), _LD1(x),_LD2(x),MAT3(_C0(y),0,0,0,_LD1(y),_LD2(y)),MAT3(_C1(y),0,0,0,_LD1(y),_LD2(y))); }
    else if(_T(x)=='R' && _T(y)=='c'){ z=multi_allocate('C',_M(x),_N(x),_L(x)); irmat3_mul_c(_M(z),_N(z),_L(z),_C0(z),_C1(z),_LD1(z),_LD2(z),_R0(x),_R1(x),_LD1(x),_LD2(x),MAT3(_C(y),0,0,0, _LD1(y),_LD2(y)),MAT3(_C(y),0,0,0, _LD1(y),_LD2(y))); }	 
    else if(_T(x)=='C' && _T(y)=='R'){ z=multi_allocate('C',_M(x),_N(x),_L(x)); icmat3_mul_r(_M(z),_N(z),_L(z),_C0(z),_C1(z),_LD1(z),_LD2(z),_C0(x),_C1(x),_LD1(x),_LD2(x),MAT3(_R0(y),0,0,0,_LD1(y),_LD2(y)),MAT3(_R1(y),0,0,0,_LD1(y),_LD2(y))); }
    else if(_T(x)=='c' && _T(y)=='R'){ z=multi_allocate('C',_M(x),_N(x),_L(x)); icmat3_mul_r(_M(z),_N(z),_L(z),_C0(z),_C1(z),_LD1(z),_LD2(z),_C(x), _C(x), _LD1(x),_LD2(x),MAT3(_R0(y),0,0,0,_LD1(y),_LD2(y)),MAT3(_R1(y),0,0,0,_LD1(y),_LD2(y))); }
    else if(_T(x)=='C' && _T(y)=='r'){ z=multi_allocate('C',_M(x),_N(x),_L(x)); icmat3_mul_r(_M(z),_N(z),_L(z),_C0(z),_C1(z),_LD1(z),_LD2(z),_C0(x),_C1(x),_LD1(x),_LD2(x),MAT3(_R(y),0,0,0, _LD1(y),_LD2(y)),MAT3(_R(y),0,0,0, _LD1(y),_LD2(y))); }	 
    else{ MATLAB_ERROR("multi_times: Unkown type"); }
  }else{ MATLAB_ERROR("multi_times: z=x.*y: Dimensions of x and y are NOT same."); }
  // done
  plhs[0]=mxCreateStructMulti(z);
  
  x=multi_free(x);
  y=multi_free(y);
  z=multi_free(z);
  return;
}

//EOF
