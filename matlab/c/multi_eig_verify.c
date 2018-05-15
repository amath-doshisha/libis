/**
 * @brief [V,D,EV,ED]=eig_verify(A,prec_verify)
 * @detail 
 * D=eig_verify(A,prec_verify)
 * [V,D]=eig_verify(A,prec_verify)
 * [V,D,EV,ED]=eig_verify(A,prec_verify)
 */
void multi_eig_verify(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  int m=0,prec_verify=128,prec=-1,kprec=-1,LDA,debug=0;
  array *A=NULL,*cA=NULL,*D=NULL,*rD=NULL,*V=NULL,*rV=NULL,*ED=NULL,*rED=NULL,*EV=NULL,*rEV=NULL;
  if(nlhs>=5)                 { mexErrMsgIdAndTxt("MATLAB:multi_mex:maxlhs","Too many output arguments."); }
  if(nrhs>N0+3)               { mexErrMsgIdAndTxt("MATLAB:multi_mex:maxrhs","Too many input arguments."); }
  if(!(IS_STRT(nrhs,prhs,N0))){ MATLAB_ERROR("multi_eig_verify: The 1st-arg should be Struct."); }
  if(IS_NUMR(nrhs,prhs,N0+1)){ prec_verify=GET_DOUBLE(prhs[N0+1])[0]; }
  if(IS_NUMR(nrhs,prhs,N0+2)){ debug=GET_DOUBLE(prhs[N0+2])[0]; }
  // allocate by clone
  A=mxArray_to_array(prhs[N0]);
  // check size
  if(!(_M(A)==_N(A) && _L(A)==1)){ MATLAB_ERROR("multi_eig: eig(A): A should be square matrix."); }
  m=_M(A);
  // allocate
  D=multi_allocate('c',m,1,1);
  ED=multi_allocate('c',m,1,1);
  V=multi_allocate('c',m,m,1);
  EV=multi_allocate('c',m,m,1);
  // compute
  if(_T(A)=='r'){    
    cA=multi_allocate('c',m,m,1);
    cmat_copy_rmat(m,m,_C(cA),_LD1(cA),_R(A),_LD1(A)); // cA=A
    chpeig_verify(m,_C(V),_LD1(V),_C(D),_C(EV),_LD1(EV),_C(ED),_C(cA),_LD1(cA),prec_verify,&prec,&kprec,debug);
  }else if(_T(A)=='c'){
    chpeig_verify(m,_C(V),_LD1(V),_C(D),_C(EV),_LD1(EV),_C(ED),_C(A),_LD1(A),prec_verify,&prec,&kprec,debug);
  }
  else{ MATLAB_ERROR("multi_eig_verify: Unkown type"); }
  // real?
  if(cvec_is_real(m,_C(D)) && cvec_is_real(m,_C(ED))){
    rD=multi_allocate('r',m,1,1);
    rED=multi_allocate('r',m,1,1);
    rvec_real_cvec(m,_R(rD),_C(D));
    rvec_real_cvec(m,_R(rED),_C(ED));
  }
  if(cmat_is_real(_M(V),_N(V),_C(V),_LD1(V)) && cmat_is_real(_M(EV),_N(EV),_C(EV),_LD1(EV))){
    rV=multi_allocate('r',m,m,1);
    rEV=multi_allocate('r',m,m,1);
    cmat_real(m,m,_R(rV),_LD1(rV),_C(V),_LD1(V));
    cmat_real(m,m,_R(rEV),_LD1(rEV),_C(EV),_LD1(EV));
  }
  // done
  if(nlhs==1){ // D=eig_verify(A,prec_verify)
    if(rD!=NULL){ plhs[0]=array_to_mxArray(rD); }
    else        { plhs[0]=array_to_mxArray(D); }
  }else if(nlhs==2){ // [V,D]=eig_verify(A,prec_verify)
    plhs[0]=array_to_mxArray(V);
    if(rD!=NULL){ plhs[1]=array_to_mxArray(rD); }
    else        { plhs[1]=array_to_mxArray(D); }
  }else if(nlhs==3){
    MATLAB_ERROR("multi_eig_verify: 3 output arguments is illeagal.");
  }else if(nlhs==4){ // [V,D,EV,ED]=eig_verify(A,prec_verify)
    if(rV!=NULL){ plhs[0]=array_to_mxArray(rV); }
    else        { plhs[0]=array_to_mxArray(V); }    
    if(rD!=NULL){ plhs[1]=array_to_mxArray(rD); }
    else        { plhs[1]=array_to_mxArray(D); }
    if(rEV!=NULL){ plhs[2]=array_to_mxArray(rEV); }
    else         { plhs[2]=array_to_mxArray(EV); }
    if(rED!=NULL){ plhs[3]=array_to_mxArray(rED); }
    else         { plhs[3]=array_to_mxArray(ED); }
  }else{ MATLAB_ERROR("multi_eig_verify: Too many output arguments."); }
  // free
  A=array_free(A);
  cA=array_free(cA);
  V=array_free(V);
  rV=array_free(rV);
  D=array_free(D);
  rD=array_free(rD);
  EV=array_free(EV);
  rEV=array_free(rEV);
  ED=array_free(ED);
  rED=array_free(rED);
  return;
}

//EOF
