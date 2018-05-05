/**
 * @brief lambda=eig(A) or [V,D]=eig(A)
 */
void multi_eig(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  int debug=0;
  array *A=NULL,*rlambda=NULL,*clambda=NULL,*V=NULL,*D=NULL;  
  if(nlhs>=3)                 { mexErrMsgIdAndTxt("MATLAB:multi_mex:maxlhs","Too many output arguments."); }
  if(nrhs>N0+2)               { mexErrMsgIdAndTxt("MATLAB:multi_mex:maxrhs","Too many input arguments."); }
  if(!(IS_STRT(nrhs,prhs,N0))){ MATLAB_ERROR("multi_eig: The 1st-arg should be Struct."); }
  if(IS_NUMR(nrhs,prhs,N0+1)){ debug=GET_DOUBLE(prhs[N0+1])[0]; }
  // allocate by clone
  A=mxArray_to_array(prhs[N0]);
  // check size
  if(!(_M(A)==_N(A) && _L(A)==1)){ MATLAB_ERROR("multi_eig: eig(A): A should be square matrix."); }
  // output=1
  if(nlhs<=1){
    // allocate
    clambda=multi_allocate('c',_M(A),1,1);
    // operations
    if     (_T(A)=='r'){ reig_hqr(_M(A),_C(clambda),_R(A),_LD1(A),debug); }
    else if(_T(A)=='c'){ ceig_hqr(_M(A),_C(clambda),_C(A),_LD1(A),debug); }
    else{ MATLAB_ERROR("multi_eig: Unkown type"); }
    // done
    if(cvec_is_real(_M(A),_C(clambda))){
      rlambda=multi_allocate('r',_M(A),1,1);
      cvec_real(_M(A),_R(rlambda),_C(clambda));
      plhs[0]=array_to_mxArray(rlambda);
    }else{
      plhs[0]=array_to_mxArray(clambda);
    }
  }else if(nlhs==2){
    MATLAB_ERROR("multi_eig: Not Implemented Yet.");
  }else{ MATLAB_ERROR("multi_eig: Too many output arguments."); }
  // free
  A=array_free(A);
  rlambda=array_free(rlambda);
  clambda=array_free(clambda);
  V=array_free(V);
  D=array_free(D);
  return;
}

//EOF
