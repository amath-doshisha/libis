/**
 * @brief z=mldivide(x,y)
 */
void multi_mldivide(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  multi *x=NULL,*y=NULL,*z=NULL,*u=NULL,*v=NULL;
  int i,info=-1;
  if(nlhs>1){ mexErrMsgIdAndTxt("MATLAB:multi_mex:maxlhs","Too many output arguments."); }
  if(!(IS_STRT(nrhs,prhs,N0))){ MATLAB_ERROR("multi_mldivide: The 1st-arg should be Struct."); }
  if(!(IS_STRT(nrhs,prhs,N0+1))){ MATLAB_ERROR("multi_mldivide: The 2nd-arg should be Struct."); }
  // allocate by clone
  x=multi_allocate_mxArray(prhs[N0]);
  y=multi_allocate_mxArray(prhs[N0+1]);
  // operation
  if(_M(x)==1 && _N(x)==1 && _L(x)==1){
    // In the case of x is scalar
         if(_T(x)=='r' && _T(y)=='r'){ z=multi_allocate('r',_M(x),_N(x),_L(x)); rmat3_div_r2(_M(z),_N(z),_L(z),_R(z),_LD1(z),_LD2(z),_R(y),_LD1(y),_LD2(y),MAT3(_R(x),0,0,0,_LD1(x),_LD2(x))); }
    else if(_T(x)=='r' && _T(y)=='c'){ z=multi_allocate('c',_M(x),_N(x),_L(x)); cmat3_div_r2(_M(z),_N(z),_L(z),_C(z),_LD1(z),_LD2(z),_C(y),_LD1(y),_LD2(y),MAT3(_R(x),0,0,0,_LD1(x),_LD2(x))); }
    else if(_T(x)=='c' && _T(y)=='r'){ z=multi_allocate('c',_M(x),_N(x),_L(x)); rmat3_div_c2(_M(z),_N(z),_L(z),_C(z),_LD1(z),_LD2(z),_R(y),_LD1(y),_LD2(y),MAT3(_C(x),0,0,0,_LD1(x),_LD2(x))); }
    else if(_T(x)=='c' && _T(y)=='c'){ z=multi_allocate('c',_M(x),_N(x),_L(x)); cmat3_div_c2(_M(z),_N(z),_L(z),_C(z),_LD1(z),_LD2(z),_C(y),_LD1(y),_LD2(y),MAT3(_C(x),0,0,0,_LD1(x),_LD2(x))); }
    else{ MATLAB_ERROR("multi_rdivide: Unkown type"); }
  }else if(_M(y)==1 && _N(y)==1 && _L(y)==1){
    // In the case of y is scalar
    MATLAB_ERROR("multi_mldivide: z=x\\y: Dimensions of x and y are NOT same.");
    //[時間の都合で実装を後に回す。]
    //もし、xが、N*1の行ベクトルであった場合、[例: [1 1 3]\1 ]
    //実数の場合は、行ベクトルのうちで、一番大きい値について、除算を行い、それ以外の要素については、零埋めする。
    //複素数の場合は、おそらく絶対値を比較して、大きい要素で除算を行う。絶対値が等しい場合は、おそらく位相角の一番大きい方or小さい方の除算を行っていると予想されるが、未確認。
    //全く同じ値の場合は、先に出現した値について除算を行い、それ以外の要素については、零埋めする。
  }else if(_M(x)==_N(x) && _M(x)==_M(y) && _L(x)==1 && _L(y)==1){
    // In the case where x is square matrix and y's colums is same as x's row
    //ランク落ちするとNaNを返す
         if(_T(x)=='r' && _T(y)=='r'){ z=multi_allocate('r',_N(x),_N(y),1); rsolve(_M(y),_N(y),_R(y),_LD1(y),_R(x),_LD1(x),&info); rmat_copy(_M(y),_N(y),_R(z),_LD1(z),_R(y),_LD1(y)); if(info!=0){rmat_set_nan(_M(y),_N(y),_R(z),_LD1(z));}}
    else if(_T(x)=='r' && _T(y)=='c'){ v=multi_allocate('c',_M(x),_N(x),1); cmat_copy_rmat(_M(x),_N(x),_C(v),_LD1(v),_R(x),_LD1(x)); z=multi_allocate('c',_N(x),_N(y),1); csolve(_M(y),_N(y),_C(y),_LD1(y),_C(v),_LD1(v),&info); cmat_copy(_M(y),_N(y),_C(z),_LD1(z),_C(y),_LD1(y)); if(info!=0){cmat_set_nan(_M(y),_N(y),_C(z),_LD1(z));}}
    else if(_T(x)=='c' && _T(y)=='r'){ u=multi_allocate('c',_M(y),_N(y),1); cmat_copy_rmat(_M(y),_N(y),_C(u),_LD1(u),_R(y),_LD1(y)); z=multi_allocate('c',_N(x),_N(y),1); csolve(_M(y),_N(y),_C(u),_LD1(u),_C(x),_LD1(x),&info); cmat_copy(_M(y),_N(y),_C(z),_LD1(z),_C(u),_LD1(u)); if(info!=0){cmat_set_nan(_M(y),_N(y),_C(z),_LD1(z));}}
    else if(_T(x)=='c' && _T(y)=='c'){ z=multi_allocate('c',_N(x),_N(y),1); csolve(_M(y),_N(y),_C(y),_LD1(y),_C(x),_LD1(x),&info); cmat_copy(_M(y),_N(y),_C(z),_LD1(z),_C(y),_LD1(y)); if(info!=0){cmat_set_nan(_M(y),_N(y),_C(z),_LD1(z));}}
    //追加
    else if(_T(x)=='R' && _T(y)=='R'){ z=multi_allocate('R',_N(x),_N(y),1); irsolve(_M(y),_N(y),_R0(y),_R1(y),_LD1(y),_R0(x),_R1(x),_LD1(x),&info); irmat_copy(_M(z),_N(z),_R0(z),_LD1(z),_R1(z),_LD1(z),_R0(y),_LD1(y),_R1(y),_LD1(y)); if(info!=0){irmat_set_nan(_M(z),_N(z),_R0(z),_LD1(z),_R1(z),_LD1(z)); }} 
    else if(_T(x)=='r' && _T(y)=='R'){ z=multi_allocate('R',_N(x),_N(y),1); v=multi_allocate('R',_M(x),_N(x),1); irmat_copy(_M(v),_N(v),_R0(v),_LD1(v),_R1(v),_LD1(v),_R(x),_LD1(x),_R(x),_LD1(x)); irsolve(_M(y),_N(y),_R0(y),_R1(y),_LD1(y),_R0(v),_R1(v),_LD1(v),&info); irmat_copy(_M(z),_N(z),_R0(z),_LD1(z),_R1(z),_LD1(z),_R0(y),_LD1(y),_R1(y),_LD1(y)); if(info!=0){irmat_set_nan(_M(z),_N(z),_R0(z),_LD1(z),_R1(z),_LD1(z)); }}   
    else if(_T(x)=='R' && _T(y)=='r'){ z=multi_allocate('R',_N(x),_N(y),1); u=multi_allocate('R',_M(y),_N(y),1); irmat_copy(_M(u),_N(u),_R0(u),_LD1(u),_R1(u),_LD1(u),_R(y),_LD1(y),_R(y),_LD1(y)); irsolve(_M(u),_N(u),_R0(u),_R1(u),_LD1(u),_R0(x),_R1(x),_LD1(x),&info); irmat_copy(_M(z),_N(z),_R0(z),_LD1(z),_R1(z),_LD1(z),_R0(u),_LD1(u),_R1(u),_LD1(u)); if(info!=0){irmat_set_nan(_M(z),_N(z),_R0(z),_LD1(z),_R1(z),_LD1(z)); }}

	
    else if(_T(x)=='C' && _T(y)=='C'){ z=multi_allocate('C',_N(x),_N(y),1); icsolve(_M(y),_N(y),_C0(y),_C1(y),_LD1(y),_C0(x),_C1(x),_LD1(x),&info); icmat_copy(_M(z),_N(z),_C0(z),_LD1(z),_C1(z),_LD1(z),_C0(y),_LD1(y),_C1(y),_LD1(y));
      if(info!=0){icmat_set_nan(_M(z),_N(z),_C0(z),_LD1(z),_C1(z),_LD1(z)); }} 	 
    else if(_T(x)=='c' && _T(y)=='C'){ z=multi_allocate('C',_N(x),_N(y),1); v=multi_allocate('C',_M(x),_N(x),1); icmat_copy(_M(v),_N(v),_C0(v),_LD1(v),_C1(v),_LD1(v),_C(x),_LD1(x),_C(x),_LD1(x));
      icsolve(_M(y),_N(y),_C0(y),_C1(y),_LD1(y),_C0(v),_C1(v),_LD1(v),&info); icmat_copy(_M(z),_N(z),_C0(z),_LD1(z),_C1(z),_LD1(z),_C0(y),_LD1(y),_C1(y),_LD1(y));
      if(info!=0){icmat_set_nan(_M(z),_N(z),_C0(z),_LD1(z),_C1(z),_LD1(z)); }}     
    else if(_T(x)=='C' && _T(y)=='c'){ z=multi_allocate('C',_N(x),_N(y),1); u=multi_allocate('C',_M(y),_N(y),1); icmat_copy(_M(u),_N(u),_C0(u),_LD1(u),_C1(u),_LD1(u),_C(y),_LD1(y),_C(y),_LD1(y));
      icsolve(_M(u),_N(u),_C0(u),_C1(u),_LD1(u),_C0(x),_C1(x),_LD1(x),&info); icmat_copy(_M(z),_N(z),_C0(z),_LD1(z),_C1(z),_LD1(z),_C0(u),_LD1(u),_C1(u),_LD1(u));
      if(info!=0){icmat_set_nan(_M(z),_N(z),_C0(z),_LD1(z),_C1(z),_LD1(z)); }} 


    else if(_T(x)=='R' && _T(y)=='C'){ z=multi_allocate('C',_N(x),_N(y),1); v=multi_allocate('C',_M(x),_N(x),1); icmat_copy_irmat(_M(v),_N(v),_C0(v),_LD1(v),_C1(v),_LD1(v),_R0(x),_LD1(x),_R1(x),_LD1(x));
      icsolve(_M(y),_N(y),_C0(y),_C1(y),_LD1(y),_C0(v),_C1(v),_LD1(v),&info); icmat_copy(_M(z),_N(z),_C0(z),_LD1(z),_C1(z),_LD1(z),_C0(y),_LD1(y),_C1(y),_LD1(y)); if(info!=0){icmat_set_nan(_M(z),_N(z),_C0(z),_LD1(z),_C1(z),_LD1(z)); }}
	 
    else if(_T(x)=='r' && _T(y)=='C'){ z=multi_allocate('C',_N(x),_N(y),1); v=multi_allocate('C',_M(x),_N(x),1); icmat_copy_rmat(_M(v),_N(v),_C0(v),_LD1(v),_C1(v),_LD1(v),_R(x),_LD1(x));
      icsolve(_M(y),_N(y),_C0(y),_C1(y),_LD1(y),_C0(v),_C1(v),_LD1(v),&info); icmat_copy(_M(z),_N(z),_C0(z),_LD1(z),_C1(z),_LD1(z),_C0(y),_LD1(y),_C1(y),_LD1(y)); if(info!=0){icmat_set_nan(_M(z),_N(z),_C0(z),_LD1(z),_C1(z),_LD1(z)); }}
	 
    else if(_T(x)=='R' && _T(y)=='c'){ z=multi_allocate('C',_N(x),_N(y),1); v=multi_allocate('C',_M(x),_N(x),1); u=multi_allocate('C',_M(y),_N(y),1);
      icmat_copy_irmat(_M(v),_N(v),_C0(v),_LD1(v),_C1(v),_LD1(v),_R0(x),_LD1(x),_R1(x),_LD1(x)); icmat_copy(_M(u),_N(u),_C0(u),_LD1(u),_C1(u),_LD1(u),_C(y),_LD1(y),_C(y),_LD1(y));
      icsolve(_M(u),_N(u),_C0(u),_C1(u),_LD1(u),_C0(v),_C1(v),_LD1(v),&info); icmat_copy(_M(z),_N(z),_C0(z),_LD1(z),_C1(z),_LD1(z),_C0(u),_LD1(u),_C1(u),_LD1(u)); if(info!=0){icmat_set_nan(_M(z),_N(z),_C0(z),_LD1(z),_C1(z),_LD1(z)); }}

	 

    else if(_T(x)=='C' && _T(y)=='R'){ z=multi_allocate('C',_N(x),_N(y),1); u=multi_allocate('C',_M(y),_N(y),1); icmat_copy_irmat(_M(u),_N(u),_C0(u),_LD1(u),_C1(u),_LD1(u),_R0(y),_LD1(y),_R1(y),_LD1(y));
      icsolve(_M(u),_N(u),_C0(u),_C1(u),_LD1(u),_C0(x),_C1(x),_LD1(x),&info); icmat_copy(_M(z),_N(z),_C0(z),_LD1(z),_C1(z),_LD1(z),_C0(u),_LD1(u),_C1(u),_LD1(u));
      if(info!=0){icmat_set_nan(_M(z),_N(z),_C0(z),_LD1(z),_C1(z),_LD1(z)); }}
    
    else if(_T(x)=='c' && _T(y)=='R'){ z=multi_allocate('C',_N(x),_N(y),1); v=multi_allocate('C',_M(x),_N(x),1); icmat_copy(_M(v),_N(v),_C0(v),_LD1(v),_C1(v),_LD1(v),_C(x),_LD1(x),_C(x),_LD1(x));
      u=multi_allocate('C',_M(y),_N(y),1); icmat_copy_irmat(_M(u),_N(u),_C0(u),_LD1(u),_C1(u),_LD1(u),_R0(y),_LD1(y),_R1(y),_LD1(y));
      icsolve(_M(u),_N(u),_C0(u),_C1(u),_LD1(u),_C0(v),_C1(v),_LD1(v),&info); icmat_copy(_M(z),_N(z),_C0(z),_LD1(z),_C1(z),_LD1(z),_C0(u),_LD1(u),_C1(u),_LD1(u)); if(info!=0){icmat_set_nan(_M(z),_N(z),_C0(z),_LD1(z),_C1(z),_LD1(z)); }}
      
    else if(_T(x)=='C' && _T(y)=='r'){ z=multi_allocate('C',_N(x),_N(y),1); u=multi_allocate('C',_M(y),_N(y),1); icmat_copy_rmat(_M(u),_N(u),_C0(u),_LD1(u),_C1(u),_LD1(u),_R(y),_LD1(y));
      icsolve(_M(u),_N(u),_C0(u),_C1(u),_LD1(u),_C0(x),_C1(x),_LD1(x),&info); icmat_copy(_M(z),_N(z),_C0(z),_LD1(z),_C1(z),_LD1(z),_C0(u),_LD1(u),_C1(u),_LD1(u)); if(info!=0){icmat_set_nan(_M(z),_N(z),_C0(z),_LD1(z),_C1(z),_LD1(z)); }} 
    
    //ここまで
    else{ MATLAB_ERROR("multi_mldivide: Unkown type"); }
    if(info!=0){ MATLAB_WARN("multi_mldivide: not full rank.");}
  }else{
    MATLAB_ERROR("multi_mldivide: z=x\\y: This function is not implemented");
  }
  // done
  plhs[0]=mxCreateStructMulti(z);
  x=multi_free(x);
  y=multi_free(y);
  z=multi_free(z);
  u=multi_free(u);
  v=multi_free(v);
  return;
}

//EOF


