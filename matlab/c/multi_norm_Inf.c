/**
 * @brief y=norm_Inf(x)
 */
//区間行列の場合    [y0,y1]=norm_Inf([x])=max(sum(abs([x]')))
//区間ベクトルの場合 [y0,y1]=norm_Inf([x])=max(abs([x]))
void multi_norm_Inf(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  multi *x=NULL,*y=NULL,*u=NULL,*v=NULL,*w=NULL;
  if(nlhs>1){ mexErrMsgIdAndTxt("MATLAB:multi_mex:maxlhs","Too many output arguments."); }
  if(!(IS_STRT(nrhs,prhs,N0))){ MATLAB_ERROR("multi_norm_Inf: The 1st-arg should be Struct."); }
  // allocate by clone
  x=multi_allocate_mxArray(prhs[N0]);
  // operation
  // abs([x])
       if(_T(x)=='r'){ u=multi_allocate('r',_M(x),_N(x),_L(x)); rmat3_abs(_M(u),_N(u),_L(u),_R(u),_LD1(u),_LD2(u),_R(x),_LD1(x),_LD2(x)); }
  else if(_T(x)=='c'){ u=multi_allocate('r',_M(x),_N(x),_L(x)); cmat3_abs(_M(u),_N(u),_L(u),_R(u),_LD1(u),_LD2(u),_C(x),_LD1(x),_LD2(x)); }
  else if(_T(x)=='R'){ u=multi_allocate('R',_M(x),_N(x),_L(x)); irmat3_abs(_M(u),_N(u),_L(u),_R0(u),_R1(u),_LD1(u),_LD2(u),_R0(x),_R1(x),_LD1(x),_LD2(x)); }
  else if(_T(x)=='C'){ u=multi_allocate('R',_M(x),_N(x),_L(x)); icmat3_abs(_M(u),_N(u),_L(u),_R0(u),_R1(u),_LD1(u),_LD2(u),_C0(x),_C1(x),_LD1(x),_LD2(x)); }
  else{ MATLAB_ERROR("multi_abs: Unkown type"); }
  if((_M(u)==1 && _L(u)==1) || (_N(u)==1 && _L(u)==1)){
         if(_T(u)=='r'){ y=multi_allocate('r',1,1,1); rvec_max(_R(y)[0],_M(u)*_N(u)*_L(u),_R(u)); }
    else if(_T(u)=='R'){ y=multi_allocate('R',1,1,1); irvec_max(_R0(y)[0],_R1(y)[0],_M(u)*_N(u)*_L(u),_R0(u),_R1(u)); }
  }else if(_L(u)==1){
         if(_T(u)=='r'){ y=multi_allocate('r',1,1,1); v=multi_allocate(_T(u),_N(u),_M(u),_L(u));  rmat_copy_t(_M(u),_N(u),_R(v),_LD1(v),_R(u),_LD1(u));
	   w=multi_allocate('r',1,_N(v),1); rvec_sum_rmat(_M(v),_N(v),_R(w),_R(v),_LD1(v)); rvec_max(_R(y)[0],_N(w),_R(w)); }
    else if(_T(u)=='R'){ y=multi_allocate('R',1,1,1); v=multi_allocate(_T(u),_N(u),_M(u),_L(u)); irmat_copy_t(_M(u),_N(u),_R0(v),_LD1(v),_R1(v),_LD1(v),_R0(u),_LD1(u),_R1(u),_LD1(u));
      w=multi_allocate('R',1,_N(v),1); irvec_sum_irmat(_M(v),_N(v),_R0(w),_R1(w),_R0(v),_R1(v),_LD1(v)); irvec_max(_R0(y)[0],_R1(y)[0],_N(w),_R0(w),_R1(w)); }
   }else{ MATLAB_ERROR("multi_norm_Inf: error"); }
  // done
  plhs[0]=mxCreateStructMulti(y);
  x=multi_free(x);
  y=multi_free(y);
  u=multi_free(u);
  v=multi_free(v);
  w=multi_free(w);
  return;
}

//EOF

