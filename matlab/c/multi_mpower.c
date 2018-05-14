/**
 * @brief z=mpower(x,y)
 */
void multi_mpower(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  array *x=NULL,*y=NULL,*z=NULL;
  if(nlhs>1)                    { mexErrMsgIdAndTxt("MATLAB:multi_mex:maxlhs","Too many output arguments."); }
  if(!(IS_STRT(nrhs,prhs,N0)))  { MATLAB_ERROR("multi_power: The 1st-arg should be Struct."); }
  if(!(IS_STRT(nrhs,prhs,N0+1))){ MATLAB_ERROR("multi_power: The 2nd-arg should be Struct."); }
  // allocate by clone
  x=mxArray_to_array(prhs[N0]);
  y=mxArray_to_array(prhs[N0+1]);
  if(!(_M(x)==_N(x) && _L(y)==1))        { MATLAB_ERROR("multi_mpower: z=x^y: x should be square matrix or scalar."); }
  if(!(_M(y)==1 && _N(y)==1 && _L(y)==1)){ MATLAB_ERROR("multi_mpower: z=x^y: y should be scalar."); }
  // operation
  if(_M(x)==1 && _N(x)==1 && _L(x)==1){ // if x is scalar
         if(_T(x)=='r' && _T(y)=='r'){ z=multi_allocate('r',1,1,1); rpow_rr(_R(z)[0],_R(x)[0],_R(y)[0]); }
    else if(_T(x)=='r' && _T(y)=='c'){ z=multi_allocate('c',1,1,1); cpow_rc(_C(z)[0],_R(x)[0],_C(y)[0]); }
    else if(_T(x)=='c' && _T(y)=='r'){ z=multi_allocate('c',1,1,1); cpow_cr(_C(z)[0],_C(x)[0],_R(y)[0]); }
    else if(_T(x)=='c' && _T(y)=='c'){ z=multi_allocate('c',1,1,1); cpow_cc(_C(z)[0],_C(x)[0],_C(y)[0]); }
    else{ MATLAB_ERROR("multi_mpower: Unkown type"); }
  }else if(_M(y)==1 && _N(y)==1 && _L(y)==1){ // if x is square matrix
    if(!(_T(y)=='r' && ris_integer(_R(y)[0]))) { MATLAB_ERROR("multi_mpower: z=x^y: y should be integer."); }
         if(_T(x)=='r'){ z=multi_allocate('r',_M(x),_M(x),1); rmat_power(_M(x),_R(z),_LD1(z),_R(x),_LD1(x),(int)rget_d(_R(y)[0])); }
    else if(_T(x)=='c'){ z=multi_allocate('c',_M(x),_M(x),1); cmat_power(_M(x),_C(z),_LD1(z),_C(x),_LD1(x),(int)rget_d(_R(y)[0])); }
    else{ MATLAB_ERROR("multi_mpower: Unkown type"); }    
  }else{ MATLAB_ERROR("multi_mpower: z=x^y: x should be square matrix."); } 
  // done
  plhs[0]=array_to_mxArray(z);  
  x=array_free(x);
  y=array_free(y);
  z=array_free(z);
  return;
}

//EOF
