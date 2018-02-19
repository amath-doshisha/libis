/**
 * @breif multi型行列から複素数行列の生成 z=complex(x,y) or complex(x)
 */
void multi_complex(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  multi *x=NULL,*y=NULL,*z=NULL;
  if(nlhs>1){ mexErrMsgIdAndTxt("MATLAB:multi_mex:maxlhs","Too many output arguments."); }
  if(nrhs>N0+2){ mexErrMsgIdAndTxt("MATLAB:multi_mex:maxrhs","Too many input arguments."); }
  if(!(IS_STRT(nrhs,prhs,N0))){ MATLAB_ERROR("multi_complex: The 1st-arg should be Struct."); }
  if(nrhs==N0+1){
  // allocate
  x=multi_allocate_mxArray(prhs[N0]);
  // operation
       if(_T(x)=='r'){ z=multi_allocate('c',_M(x),_N(x),_L(x)); cmat3_copy_rmat3  (_M(z),_N(z),_L(z),_C(z),_LD1(z),_LD2(z),_R(x),_LD1(x),_LD2(x));  }
  else if(_T(x)=='c'){ z=multi_allocate('c',_M(x),_N(x),_L(x)); cmat3_copy        (_M(z),_N(z),_L(z),_C(z),_LD1(z),_LD2(z),_C(x),_LD1(x),_LD2(x));  }
  else if(_T(x)=='R'){ z=multi_allocate('C',_M(x),_N(x),_L(x)); icmat3_copy_irmat3(_M(z),_N(z),_L(z),_C0(z),_C1(z),_LD1(z),_LD2(z),_R0(x),_R1(x),_LD1(x),_LD2(x));  }
  else if(_T(x)=='C'){ z=multi_allocate('C',_M(x),_N(x),_L(x)); icmat3_copy       (_M(z),_N(z),_L(z),_C0(z),_C1(z),_LD1(z),_LD2(z),_C0(x),_C1(x),_LD1(x),_LD2(x));  }
  else{ MATLAB_ERROR("multi_complex: Unkown type"); }
  }else if(nrhs==N0+2){
  if(!(IS_STRT(nrhs,prhs,N0+1))){ MATLAB_ERROR("multi_complex: The 2nd-arg should be Struct."); }
  // allocate
  x=multi_allocate_mxArray(prhs[N0]);
  y=multi_allocate_mxArray(prhs[N0+1]);
  // operation
  if(_M(x)==_M(y) && _N(x)==_N(y) && _L(x)==_L(y)){  // In the case of same sizes
        if(_T(x)=='r' && _T(y)=='r'){ z=multi_allocate('c',_M(x),_N(x),_L(x)); cmat3_copy_rmat3_rmat3(_M(z),_N(z),_L(z),_C(z),_LD1(z),_LD2(z),_R(x),_LD1(x),_LD2(x),_R(y),_LD1(y),_LD2(y)); }
   else if(_T(x)=='r' && _T(y)=='R'){ z=multi_allocate('C',_M(x),_N(x),_L(x)); icmat3_copy_irmat3_irmat3(_M(z),_N(z),_L(z),_C0(z),_C1(z),_LD1(z),_LD2(z),_R(x), _R(x), _LD1(x),_LD2(x),_R0(y),_R1(y),_LD1(y),_LD2(y)); }
   else if(_T(x)=='R' && _T(y)=='r'){ z=multi_allocate('C',_M(x),_N(x),_L(x)); icmat3_copy_irmat3_irmat3(_M(z),_N(z),_L(z),_C0(z),_C1(z),_LD1(z),_LD2(z),_R0(x),_R1(x),_LD1(x),_LD2(x),_R(y), _R(y), _LD1(y),_LD2(y)); }
   else if(_T(x)=='R' && _T(y)=='R'){ z=multi_allocate('C',_M(x),_N(x),_L(x)); icmat3_copy_irmat3_irmat3(_M(z),_N(z),_L(z),_C0(z),_C1(z),_LD1(z),_LD2(z),_R0(x),_R1(x),_LD1(x),_LD2(x),_R0(y),_R1(y),_LD1(y),_LD2(y)); }
   else{ MATLAB_ERROR("multi_complex: Unkown type"); }
  }else if(_M(x)==1 && _N(x)==1 && _L(x)==1){    // In the case of x is scalar 
             if(_T(x)=='r' && _T(y)=='r'){ z=multi_allocate('c',_M(y),_N(y),_L(y)); cmat3_copy_r_rmat3(_M(z),_N(z),_L(z),_C(z),_LD1(z),_LD2(z),MAT3(_R(x),0,0,0,_LD1(x),_LD2(x)),_R(y),_LD1(y),_LD2(y)); }
	else if(_T(x)=='r' && _T(y)=='R'){ z=multi_allocate('C',_M(y),_N(y),_L(y)); icmat3_copy_ir_irmat3(_M(z),_N(z),_L(z),_C0(z),_C1(z),_LD1(z),_LD2(z),MAT3(_R(x),0,0,0,_LD1(x),_LD2(x)),MAT3( _R(x),0,0,0, _LD1(x),_LD2(x)),_R0(y),_R1(y),_LD1(y),_LD2(y));}
	else if(_T(x)=='R' && _T(y)=='r'){ z=multi_allocate('C',_M(y),_N(y),_L(y)); icmat3_copy_ir_irmat3(_M(z),_N(z),_L(z),_C0(z),_C1(z),_LD1(z),_LD2(z),MAT3(_R0(x),0,0,0,_LD1(x),_LD2(x)),MAT3(_R1(x),0,0,0,_LD1(x),_LD2(x)),_R(y), _R(y), _LD1(y),_LD2(y));}
	else if(_T(x)=='R' && _T(y)=='R'){ z=multi_allocate('C',_M(y),_N(y),_L(y)); icmat3_copy_ir_irmat3(_M(z),_N(z),_L(z),_C0(z),_C1(z),_LD1(z),_LD2(z),MAT3(_R0(x),0,0,0,_LD1(x),_LD2(x)),MAT3(_R1(x),0,0,0,_LD1(x),_LD2(x)),_R0(y),_R1(y),_LD1(y),_LD2(y));}
  
   else{ MATLAB_ERROR("multi_complex: Unkown type"); }
  }else if(_M(y)==1 && _N(y)==1 && _L(y)==1){    // In the case of y is scalar
             if(_T(x)=='r' && _T(y)=='r'){ z=multi_allocate('c',_M(x),_N(x),_L(x)); cmat3_copy_rmat3_r(_M(z),_N(z),_L(z),_C(z),_LD1(z),_LD2(z),_R(x),_LD1(x),_LD2(x),MAT3(_R(y),0,0,0,_LD1(y),_LD2(y))); }
	else if(_T(x)=='r' && _T(y)=='R'){ z=multi_allocate('C',_M(x),_N(x),_L(x)); icmat3_copy_irmat3_ir(_M(z),_N(z),_L(z),_C0(z),_C1(z),_LD1(z),_LD2(z),_R(x), _R(x), _LD1(x),_LD2(x),MAT3(_R0(y),0,0,0,_LD1(x),_LD2(x)),MAT3(_R1(y),0,0,0,_LD1(x),_LD2(x)));}
	else if(_T(x)=='R' && _T(y)=='r'){ z=multi_allocate('C',_M(x),_N(x),_L(x)); icmat3_copy_irmat3_ir(_M(z),_N(z),_L(z),_C0(z),_C1(z),_LD1(z),_LD2(z),_R0(x),_R1(x),_LD1(x),_LD2(x),MAT3(_R(y),0,0,0, _LD1(x),_LD2(x)),MAT3(_R(y),0,0,0, _LD1(x),_LD2(x)));}
	else if(_T(x)=='R' && _T(y)=='R'){ z=multi_allocate('C',_M(x),_N(x),_L(x)); icmat3_copy_irmat3_ir(_M(z),_N(z),_L(z),_C0(z),_C1(z),_LD1(z),_LD2(z),_R0(x),_R1(x),_LD1(x),_LD2(x),MAT3(_R0(y),0,0,0,_LD1(x),_LD2(x)),MAT3(_R1(y),0,0,0,_LD1(x),_LD2(x)));}
  
  }else{ MATLAB_ERROR("multi_complex: Unknown type"); }
 
  }else{ MATLAB_ERROR("multi_complex: error"); }
  // done
  plhs[0]=mxCreateStructMulti(z);
  x=multi_free(x);
  y=multi_free(y);
  z=multi_free(z);
  return;
}

//EOF
