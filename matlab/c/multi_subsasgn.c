/**
 * @breif x=subsasgn(x,s,y), i.e. x(s)=y,
 */
void multi_subsasgn(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  multi *x=NULL,*y=NULL,*z=NULL;
  subs_index_t *s=NULL;
  if(nlhs>1){ mexErrMsgIdAndTxt("MATLAB:multi_mex:maxlhs","Too many output arguments."); }
  if(!(IS_STRT(nrhs,prhs,N0)))  { MATLAB_ERROR("multi_subsasgn: The 1st-arg should be Struct."); }
  if(!(IS_STRT(nrhs,prhs,N0+1))){ MATLAB_ERROR("multi_subsasgn: The 2nd-arg should be Struct."); }
  if(!(IS_STRT(nrhs,prhs,N0+2))){ MATLAB_ERROR("multi_subsasgn: The 3rd-arg should be Struct."); }
  // load the 1st-arg
  x=multi_allocate_mxArray(prhs[N0]);
  // load the 2nd-arg
  s=subs_index_allocate(prhs[N0+1]);
  // load the 3rd-arg
  y=multi_allocate_mxArray(prhs[N0+2]);
  // asignment
  if(s->ndim==1 && s->dims[0]==0 && s->index[0]==NULL){ // A(:)
    if(_M(x)*_N(x)*_L(x)!=_M(y)*_N(y)*_L(y)){ MATLAB_ERROR("multi_subsasgn: The index is out of the dimensions of the matrix."); }
    // z=x
    z=multi_allocate(_T(x),_M(x),_N(x),_L(x));
         if(_T(z)=='r'){ rmat3_copy(_M(z),_N(z),_L(z),_R(z),_LD1(z),_LD2(z),_R(x),_LD1(x),_LD2(x)); }
    else if(_T(z)=='c'){ cmat3_copy(_M(z),_N(z),_L(z),_C(z),_LD1(z),_LD2(z),_C(x),_LD1(x),_LD2(x)); }
    else{ MATLAB_ERROR("multi_subsasgn: Unkown type"); }
    // z(:)=y
         if(_T(z)=='r'){ rvec_copy(_M(z)*_N(z)*_L(z),_R(z),_R(y)); }
    else if(_T(z)=='c'){ cvec_copy(_M(z)*_N(z)*_L(z),_C(z),_C(y)); }
    else{ MATLAB_ERROR("multi_subsasgn: Unkown type"); }
  }else if(s->ndim==1){ // A(*)
    if(s->dims[0]!=_M(y)*_N(y)*_L(y)){ MATLAB_ERROR("multi_subsasgn: The index is out of the dimensions of the matrix."); }
    // z=x
    z=multi_allocate(_T(x),_M(x),_N(x),_L(x));
         if(_T(z)=='r'){ rmat3_copy(_M(z),_N(z),_L(z),_R(z),_LD1(z),_LD2(z),_R(x),_LD1(x),_LD2(x)); }
    else if(_T(z)=='c'){ cmat3_copy(_M(z),_N(z),_L(z),_C(z),_LD1(z),_LD2(z),_C(x),_LD1(x),_LD2(x)); }
    else{ MATLAB_ERROR("multi_subsasgn: Unkown type"); }
    // z(s)=y
         if(_T(z)=='r'){ rvec_index_copy(s->dims[0],_R(z),_R(y),s->index[0]); }
    else if(_T(z)=='c'){ cvec_index_copy(s->dims[0],_C(z),_C(y),s->index[0]); }
    else{ MATLAB_ERROR("multi_subsasgn: Unkown type"); }
  }else{ // otherwise
    MATLAB_ERROR("multi_subsasgn: NOT SUPPORTED!");
  }  
  // done
  plhs[0]=mxCreateStructMulti(z);
  subs_index_free(s);
  x=multi_free(x);
  y=multi_free(y);
  z=multi_free(z);
  return;
}

//EOF
