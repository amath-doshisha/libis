/**
 * @brief y=max(x) or y=max(x,z)
 */
//区間の場合 [y0,y1]=[max(x0),max(x1)]
void multi_max(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  array *x=NULL,*y=NULL,*z=NULL;
  if(nlhs>1){ mexErrMsgIdAndTxt("MATLAB:multi_mex:maxlhs","Too many output arguments."); }
  if(nrhs>N0+2){ mexErrMsgIdAndTxt("MATLAB:multi_mex:maxrhs","Too many input arguments."); }
  if(!(IS_STRT(nrhs,prhs,N0))){ MATLAB_ERROR("multi_max: The 1st-arg should be Struct."); }
  // allocate by clone
  x=mxArray_to_array(prhs[N0]);
  // operation
  ccmp_set_abs_arg();
  if(nrhs==N0+1){
    if((_M(x)==1 && _L(x)==1) || (_N(x)==1 && _L(x)==1)){
           if(_T(x)=='r'){ y=multi_allocate('r',1,1,1); rvec_max(_R(y)[0],_M(x)*_N(x)*_L(x),_R(x)); }
      else if(_T(x)=='c'){ y=multi_allocate('c',1,1,1); cvec_max(_C(y)[0],_M(x)*_N(x)*_L(x),_C(x)); }
      else if(_T(x)=='R'){ y=multi_allocate('R',1,1,1); irvec_max(_R0(y)[0],_R1(y)[0],_M(x)*_N(x)*_L(x),_R0(x),_R1(x)); }
      else if(_T(x)=='C'){ y=multi_allocate('C',1,1,1); icvec_max(_C0(y)[0],_C1(y)[0],_M(x)*_N(x)*_L(x),_C0(x),_C1(x)); }
    }else if(_L(x)==1){
           if(_T(x)=='r'){ y=multi_allocate('r',1,_N(x),1); rvec_max_rmat(_M(x),_N(x),_R(y),_R(x),_LD1(x)); }
      else if(_T(x)=='c'){ y=multi_allocate('c',1,_N(x),1); cvec_max_cmat(_M(x),_N(x),_C(y),_C(x),_LD1(x)); }
      else if(_T(x)=='R'){ y=multi_allocate('R',1,_N(x),1); irvec_max_irmat(_M(x),_N(x),_R0(y),_R1(y),_R0(x),_R1(x),_LD1(x)); }
      else if(_T(x)=='C'){ y=multi_allocate('C',1,_N(x),1); icvec_max_icmat(_M(x),_N(x),_C0(y),_C1(y),_C0(x),_C1(x),_LD1(x)); }
    }else if(_L(x)>1 && _M(x)==1 && _N(x)==1){
           if(_T(x)=='r'){ y=multi_allocate('r',1,1,1); rvec_max(_R(y)[0],_M(x)*_N(x)*_L(x),_R(x)); }
      else if(_T(x)=='c'){ y=multi_allocate('c',1,1,1); cvec_max(_C(y)[0],_M(x)*_N(x)*_L(x),_C(x)); }
      else if(_T(x)=='R'){ y=multi_allocate('R',1,1,1); irvec_max(_R0(y)[0],_R1(y)[0],_M(x)*_N(x)*_L(x),_R0(x),_R1(x)); }
      else if(_T(x)=='C'){ y=multi_allocate('C',1,1,1); icvec_max(_C0(y)[0],_C1(y)[0],_M(x)*_N(x)*_L(x),_C0(x),_C1(x)); }
    }else if(_L(x)>1 && _M(x)>1 && _N(x)==1){
           if(_T(x)=='r'){ y=multi_allocate('r',1,1,_L(x)); rvec_max_rmat(_M(x),_L(x),_R(y),_R(x),_M(x)); }
      else if(_T(x)=='c'){ y=multi_allocate('c',1,1,_L(x)); cvec_max_cmat(_M(x),_L(x),_C(y),_C(x),_M(x)); }
      else if(_T(x)=='R'){ y=multi_allocate('R',1,1,_L(x)); irvec_max_irmat(_M(x),_L(x),_R0(y),_R1(y),_R0(x),_R1(x),_M(x)); }
      else if(_T(x)=='C'){ y=multi_allocate('C',1,1,_L(x)); icvec_max_icmat(_M(x),_L(x),_C0(y),_C1(y),_C0(x),_C1(x),_M(x)); }
    }else if(_L(x)>1 && _M(x)==1 && _N(x)>1){
           if(_T(x)=='r'){ y=multi_allocate('r',1,1,_L(x)); rvec_max_rmat(_N(x),_L(x),_R(y),_R(x),_N(x)); }
      else if(_T(x)=='c'){ y=multi_allocate('c',1,1,_L(x)); cvec_max_cmat(_N(x),_L(x),_C(y),_C(x),_N(x)); }
      else if(_T(x)=='R'){ y=multi_allocate('R',1,1,_L(x)); irvec_max_irmat(_N(x),_L(x),_R0(y),_R1(y),_R0(x),_R1(x),_N(x)); }
      else if(_T(x)=='C'){ y=multi_allocate('C',1,1,_L(x)); icvec_max_icmat(_N(x),_L(x),_C0(y),_C1(y),_C0(x),_C1(x),_N(x)); }
    }else if(_L(x)>1 && _M(x)>1 && _N(x)>1){
           if(_T(x)=='r'){ y=multi_allocate('r',1,_N(x),_L(x)); rmat3_max(_M(x),_N(x),_L(x),_R(y),_LD1(y),_LD2(y),_R(x),_LD1(x),_LD2(x)); }
      else if(_T(x)=='c'){ y=multi_allocate('c',1,_N(x),_L(x)); cmat3_max(_M(x),_N(x),_L(x),_C(y),_LD1(y),_LD2(y),_C(x),_LD1(x),_LD2(x)); }
      else if(_T(x)=='R'){ y=multi_allocate('R',1,_N(x),_L(x)); irmat3_max(_M(x),_N(x),_L(x),_R0(y),_R1(y),_LD1(y),_LD2(y),_R0(x),_R1(x),_LD1(x),_LD2(x)); }
      else if(_T(x)=='C'){ y=multi_allocate('C',1,_N(x),_L(x)); icmat3_max(_M(x),_N(x),_L(x),_C0(y),_C1(y),_LD1(y),_LD2(y),_C0(x),_C1(x),_LD1(x),_LD2(x)); }
    }else{ MATLAB_ERROR("multi_max: error"); }
  }else if(nrhs==N0+2){
    // allocate by clone
    z=mxArray_to_array(prhs[N0+1]);
    if(_M(x)==_M(z) && _N(x)==_N(z) && _L(x)==_L(z)){
      if(_T(x)=='r' && _T(z)=='r'){ y=multi_allocate('r',_M(x),_N(x),_L(x)); rmat3_max2(_M(x),_N(x),_L(x),_R(y),_LD1(y),_LD2(y),_R(x),_LD1(x),_LD2(x),_R(z),_LD1(z),_LD2(z)); }
    }else{ MATLAB_ERROR("multi_max: error not supported input combination");}
  }else{ MATLAB_ERROR("multi_max: error"); }
  // done
  plhs[0]=array_to_mxArray(y);
  x=array_free(x);
  y=array_free(y);
  z=array_free(z); 
  return;
}
  //EOF


