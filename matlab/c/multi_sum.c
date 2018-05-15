//作成
/**
 * @brief y=sum(x)
 */
void multi_sum(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  array *x=NULL,*y=NULL;
  if(nlhs>1){ mexErrMsgIdAndTxt("MATLAB:multi_mex:maxlhs","Too many output arguments."); }
  if(!(IS_STRT(nrhs,prhs,N0))){ MATLAB_ERROR("multi_max: The 1st-arg should be Struct."); }
  // allocate by clone
  x=mxArray_to_array(prhs[N0]);
  // operation
  if((_M(x)==1 && _L(x)==1) || (_N(x)==1 && _L(x)==1)){
        if(_T(x)=='r'){ y=multi_allocate('r',1,1,1); rsum_rvec(_R(y)[0],_M(x)*_N(x)*_L(x),_R(x)); }
   else if(_T(x)=='c'){ y=multi_allocate('c',1,1,1); csum_cvec(_C(y)[0],_M(x)*_N(x)*_L(x),_C(x)); }
   else if(_T(x)=='R'){ y=multi_allocate('R',1,1,1); irsum_rvec(_R0(y)[0],_R1(y)[0],_M(x)*_N(x)*_L(x),_R0(x),_R1(x)); }
   else if(_T(x)=='C'){ y=multi_allocate('C',1,1,1); icsum_cvec(_C0(y)[0],_C1(y)[0],_M(x)*_N(x)*_L(x),_C0(x),_C1(x)); }
  }else if(_L(x)==1){
        if(_T(x)=='r'){ y=multi_allocate('r',1,_N(x),1); rvec_sum_rmat(_M(x),_N(x),_R(y),_R(x),_LD1(x)); }
   else if(_T(x)=='c'){ y=multi_allocate('c',1,_N(x),1); cvec_sum_cmat(_M(x),_N(x),_C(y),_C(x),_LD1(x)); }
   else if(_T(x)=='R'){ y=multi_allocate('R',1,_N(x),1); irvec_sum_irmat(_M(x),_N(x),_R0(y),_R1(y),_R0(x),_R1(x),_LD1(x)); }
   else if(_T(x)=='C'){ y=multi_allocate('C',1,_N(x),1); icvec_sum_icmat(_M(x),_N(x),_C0(y),_C1(y),_C0(x),_C1(x),_LD1(x)); }
  }else if(_L(x)>1 && _M(x)==1 && _N(x)==1){
        if(_T(x)=='r'){ y=multi_allocate('r',1,1,1); rsum_rvec(_R(y)[0],_M(x)*_N(x)*_L(x),_R(x)); }
   else if(_T(x)=='c'){ y=multi_allocate('c',1,1,1); csum_cvec(_C(y)[0],_M(x)*_N(x)*_L(x),_C(x)); }
   else if(_T(x)=='R'){ y=multi_allocate('R',1,1,1); irsum_rvec(_R0(y)[0],_R1(y)[0],_M(x)*_N(x)*_L(x),_R0(x),_R1(x)); }
   else if(_T(x)=='C'){ y=multi_allocate('C',1,1,1); icsum_cvec(_C0(y)[0],_C1(y)[0],_M(x)*_N(x)*_L(x),_C0(x),_C1(x)); }
  }else if(_L(x)>1 && _M(x)>1 && _N(x)==1){
        if(_T(x)=='r'){ y=multi_allocate('r',1,1,_L(x)); rvec_sum_rmat(_M(x),_L(x),_R(y),_R(x),_M(x)); }
   else if(_T(x)=='c'){ y=multi_allocate('c',1,1,_L(x)); cvec_sum_cmat(_M(x),_L(x),_C(y),_C(x),_M(x)); }
   else if(_T(x)=='R'){ y=multi_allocate('R',1,1,_L(x)); irvec_sum_irmat(_M(x),_L(x),_R0(y),_R1(y),_R0(x),_R1(x),_M(x)); }
   else if(_T(x)=='C'){ y=multi_allocate('C',1,1,_L(x)); icvec_sum_icmat(_M(x),_L(x),_C0(y),_C1(y),_C0(x),_C1(x),_M(x)); }
  }else if(_L(x)>1 && _M(x)==1 && _N(x)>1){
        if(_T(x)=='r'){ y=multi_allocate('r',1,1,_L(x)); rvec_sum_rmat(_N(x),_L(x),_R(y),_R(x),_N(x)); }
   else if(_T(x)=='c'){ y=multi_allocate('c',1,1,_L(x)); cvec_sum_cmat(_N(x),_L(x),_C(y),_C(x),_N(x)); }
   else if(_T(x)=='R'){ y=multi_allocate('R',1,1,_L(x)); irvec_sum_irmat(_N(x),_L(x),_R0(y),_R1(y),_R0(x),_R1(x),_N(x)); }
   else if(_T(x)=='C'){ y=multi_allocate('C',1,1,_L(x)); icvec_sum_icmat(_N(x),_L(x),_C0(y),_C1(y),_C0(x),_C1(x),_N(x)); }
  }else if(_L(x)>1 && _M(x)>1 && _N(x)>1){
        if(_T(x)=='r'){ y=multi_allocate('r',1,_N(x),_L(x)); rmat3_sum(_M(x),_N(x),_L(x),_R(y),_LD1(y),_LD2(y),_R(x),_LD1(x),_LD2(x)); }
   else if(_T(x)=='c'){ y=multi_allocate('c',1,_N(x),_L(x)); cmat3_sum(_M(x),_N(x),_L(x),_C(y),_LD1(y),_LD2(y),_C(x),_LD1(x),_LD2(x)); }
   else if(_T(x)=='R'){ y=multi_allocate('R',1,_N(x),_L(x)); irmat3_sum(_M(x),_N(x),_L(x),_R0(y),_R1(y),_LD1(y),_LD2(y),_R0(x),_R1(x),_LD1(x),_LD2(x)); }
   else if(_T(x)=='C'){ y=multi_allocate('C',1,_N(x),_L(x)); icmat3_sum(_M(x),_N(x),_L(x),_C0(y),_C1(y),_LD1(y),_LD2(y),_C0(x),_C1(x),_LD1(x),_LD2(x)); }
  }else{ MATLAB_ERROR("multi_sum: error"); }
  // done
  plhs[0]=array_to_mxArray(y);
  x=array_free(x);
  y=array_free(y);
  return;
}

//EOF


