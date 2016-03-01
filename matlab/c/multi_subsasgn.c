/**
 * @breif x=subsasgn(x,s,y), i.e. x(s)=y,
 */
void multi_subsasgn(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  int m,n,l;
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
    if(_M(x)*_N(x)*_L(x)!=_M(y)*_N(y)*_L(y)){ MATLAB_ERROR("multi_subsasgn: The dimensions of both sides are not equal."); }
    // allocate
    if(_T(x)=='c' || _T(y)=='c'){ z=multi_allocate('c',_M(x),_N(x),_L(x)); }
    else                        { z=multi_allocate('r',_M(x),_N(x),_L(x)); }
    // z=x
         if(_T(z)=='r' && _T(x)=='r'){ rmat3_copy      (_M(x),_N(x),_L(x),_R(z),_LD1(z),_LD2(z),_R(x),_LD1(x),_LD2(x)); }
    else if(_T(z)=='c' && _T(x)=='r'){ cmat3_copy_rmat3(_M(x),_N(x),_L(x),_C(z),_LD1(z),_LD2(z),_R(x),_LD1(x),_LD2(x)); }
    else if(_T(z)=='c' && _T(x)=='c'){ cmat3_copy      (_M(x),_N(x),_L(x),_C(z),_LD1(z),_LD2(z),_C(x),_LD1(x),_LD2(x)); }
    else{ MATLAB_ERROR("multi_subsasgn: Unkown type"); }
    // z(:)=y
         if(_T(z)=='r' && _T(y)=='r'){ rvec_copy     (_M(z)*_N(z)*_L(z),_R(z),_R(y)); }
    else if(_T(z)=='c' && _T(y)=='r'){ cvec_copy_rvec(_M(z)*_N(z)*_L(z),_C(z),_R(y)); }
    else if(_T(z)=='c' && _T(y)=='c'){ cvec_copy     (_M(z)*_N(z)*_L(z),_C(z),_C(y)); }
    else{ MATLAB_ERROR("multi_subsasgn: Unkown type"); }
  }else if(s->ndim==1){ // A(*)
    if(s->dims[0]==0 && s->index[0]==NULL){ s->dims[0]=_M(x)*_N(x)*_L(x); s->index[0]=ivec_allocate(s->dims[0]); ivec_grid(s->dims[0],s->index[0]); }
    if(s->dims[0]!=_M(y)*_N(y)*_L(y)){ MATLAB_ERROR("multi_subsasgn: The dimensions of both sides are not equal."); }
    if(ivec_min(s->dims[0],s->index[0])<0){ MATLAB_ERROR("multi_subsasgn: The index is out of the dimensions of the matrix."); }
    // check size
    m=ivec_max(s->dims[0],s->index[0])+1;
         if(_M(x)==0 && _N(x)==0 && _L(x)==0)           { n=m;     m=1;     l=1;     }
    else if(_M(x)==1 && _L(x)==1 && m>_M(x)*_N(x)*_L(x)){ n=m;     m=1;     l=1;     }
    else if(_M(x)==1 && _L(x)==1)                       { n=_N(x); m=_M(x); l=_L(x); }
    else if(_N(x)==1 && _L(x)==1 && m>_M(x)*_N(x)*_L(x)){ m=m;     n=1;     l=1;     }
    else if(_N(x)==1 && _L(x)==1)                       { m=_M(x); n=_N(x); l=_L(x); }
    else if(m>_M(x)*_N(x)*_L(x))                        { MATLAB_ERROR("multi_subsasgn: The index is out of the dimensions of the matrix."); }
    else                                                { m=_M(x); n=_N(x); l=_L(x); }
    // allocate
    if(_T(x)=='c' || _T(y)=='c'){ z=multi_allocate('c',m,n,l); }
    else                        { z=multi_allocate('r',m,n,l); }
    // z=0
         if(_T(z)=='r'){ rmat3_set_zeros(_M(z),_N(z),_L(z),_R(z),_LD1(z),_LD2(z)); }
    else if(_T(z)=='c'){ cmat3_set_zeros(_M(z),_N(z),_L(z),_C(z),_LD1(z),_LD2(z)); }
    else{ MATLAB_ERROR("multi_subsasgn: Unkown type"); }
    // z=x
         if(_T(z)=='r' && _T(x)=='r'){ rmat3_copy      (_M(x),_N(x),_L(x),_R(z),_LD1(z),_LD2(z),_R(x),_LD1(x),_LD2(x)); }
    else if(_T(z)=='c' && _T(x)=='r'){ cmat3_copy_rmat3(_M(x),_N(x),_L(x),_C(z),_LD1(z),_LD2(z),_R(x),_LD1(x),_LD2(x)); }
    else if(_T(z)=='c' && _T(x)=='c'){ cmat3_copy      (_M(x),_N(x),_L(x),_C(z),_LD1(z),_LD2(z),_C(x),_LD1(x),_LD2(x)); }
    else{ MATLAB_ERROR("multi_subsasgn: Unkown type"); }
    // z(s)=y
         if(_T(z)=='r' && _T(y)=='r'){ rvec_index_copy     (s->dims[0],_R(z),_R(y),s->index[0]); }
    else if(_T(z)=='c' && _T(y)=='r'){ cvec_index_copy_rvec(s->dims[0],_C(z),_R(y),s->index[0]); }
    else if(_T(z)=='c' && _T(y)=='c'){ cvec_index_copy     (s->dims[0],_C(z),_C(y),s->index[0]); }
    else{ MATLAB_ERROR("multi_subsasgn: Unkown type"); }
  }else if(s->ndim==2){ // A(*,*)
    if(s->dims[0]==0 && s->index[0]==NULL){ s->dims[0]=_M(x);       s->index[0]=ivec_allocate(s->dims[0]); ivec_grid(s->dims[0],s->index[0]); }
    if(s->dims[1]==0 && s->index[1]==NULL){ s->dims[1]=_N(x)*_L(x); s->index[1]=ivec_allocate(s->dims[1]); ivec_grid(s->dims[1],s->index[1]); }
    if(!(s->dims[0]==_M(y) && s->dims[1]==_N(y)*_L(y))){ MATLAB_ERROR("multi_subsasgn: The dimensions of both sides are not equal."); }    
    if(ivec_min(s->dims[0],s->index[0])<0){ MATLAB_ERROR("multi_subsasgn: The index is out of the dimensions of the matrix."); }
    if(ivec_min(s->dims[1],s->index[1])<0){ MATLAB_ERROR("multi_subsasgn: The index is out of the dimensions of the matrix."); }    
    // check size
    m=ivec_max(s->dims[0],s->index[0])+1;
    n=ivec_max(s->dims[1],s->index[1])+1;
         if(_M(x)==0 && _N(x)==0 && _L(x)==0){ m=m; n=n; l=1; }
    else if(_L(x)==1)                        { m=MAX2(_M(x),m); n=MAX2(_N(x),n); l=_L(x); }    
    else if(_L(x)>1 && n>_N(x)*_L(x))        { MATLAB_ERROR("multi_subsasgn: The index is out of the dimensions of the matrix."); }
    else                                     { m=MAX2(_M(x),m); n=_N(x); l=_L(x); }
    // allocate
    if(_T(x)=='c' || _T(y)=='c'){ z=multi_allocate('c',m,n,l); }
    else                        { z=multi_allocate('r',m,n,l); }
    // z=0
         if(_T(z)=='r'){ rmat3_set_zeros(_M(z),_N(z),_L(z),_R(z),_LD1(z),_LD2(z)); }
    else if(_T(z)=='c'){ cmat3_set_zeros(_M(z),_N(z),_L(z),_C(z),_LD1(z),_LD2(z)); }
    else{ MATLAB_ERROR("multi_subsasgn: Unkown type"); }
    // z=x
         if(_T(z)=='r' && _T(x)=='r'){ rmat3_copy      (_M(x),_N(x),_L(x),_R(z),_LD1(z),_LD2(z),_R(x),_LD1(x),_LD2(x)); }
    else if(_T(z)=='c' && _T(x)=='r'){ cmat3_copy_rmat3(_M(x),_N(x),_L(x),_C(z),_LD1(z),_LD2(z),_R(x),_LD1(x),_LD2(x)); }
    else if(_T(z)=='c' && _T(x)=='c'){ cmat3_copy      (_M(x),_N(x),_L(x),_C(z),_LD1(z),_LD2(z),_C(x),_LD1(x),_LD2(x)); }
    else{ MATLAB_ERROR("multi_subsasgn: Unkown type"); }
    // z(s)=y
         if(_T(z)=='r' && _T(y)=='r'){ rmat_index_copy     (s->dims[0],s->dims[1],_R(z),_LD1(z),_R(y),_LD1(y),s->index[0],s->index[1]); }
    else if(_T(z)=='c' && _T(y)=='r'){ cmat_index_copy_rmat(s->dims[0],s->dims[1],_C(z),_LD1(z),_R(y),_LD1(y),s->index[0],s->index[1]); }
    else if(_T(z)=='c' && _T(y)=='c'){ cmat_index_copy     (s->dims[0],s->dims[1],_C(z),_LD1(z),_C(y),_LD1(y),s->index[0],s->index[1]); }
    else{ MATLAB_ERROR("multi_subsasgn: Unkown type"); }
  }else if(s->ndim==3){ // A(*,*,*)
    if(s->dims[0]==0 && s->index[0]==NULL){ s->dims[0]=_M(x); s->index[0]=ivec_allocate(s->dims[0]); ivec_grid(s->dims[0],s->index[0]); }
    if(s->dims[1]==0 && s->index[1]==NULL){ s->dims[1]=_N(x); s->index[1]=ivec_allocate(s->dims[1]); ivec_grid(s->dims[1],s->index[1]); }
    if(s->dims[2]==0 && s->index[2]==NULL){ s->dims[2]=_L(x); s->index[2]=ivec_allocate(s->dims[2]); ivec_grid(s->dims[2],s->index[2]); }
    if(!(s->dims[0]==_M(y) && s->dims[1]==_N(y) && s->dims[2]==_L(y))){ MATLAB_ERROR("multi_subsasgn: The dimensions of both sides are not equal."); }
    if(ivec_min(s->dims[0],s->index[0])<0){ MATLAB_ERROR("multi_subsasgn: The index is out of the dimensions of the matrix."); }
    if(ivec_min(s->dims[1],s->index[1])<0){ MATLAB_ERROR("multi_subsasgn: The index is out of the dimensions of the matrix."); }    
    if(ivec_min(s->dims[2],s->index[2])<0){ MATLAB_ERROR("multi_subsasgn: The index is out of the dimensions of the matrix."); }
    // check size
    m=ivec_max(s->dims[0],s->index[0])+1;
    n=ivec_max(s->dims[1],s->index[1])+1;
    l=ivec_max(s->dims[2],s->index[2])+1;
    m=MAX2(_M(x),m);
    n=MAX2(_N(x),n);
    l=MAX2(_L(x),l);
    // allocate
    if(_T(x)=='c' || _T(y)=='c'){ z=multi_allocate('c',m,n,l); }
    else                        { z=multi_allocate('r',m,n,l); }
    // z=0
         if(_T(z)=='r'){ rmat3_set_zeros(_M(z),_N(z),_L(z),_R(z),_LD1(z),_LD2(z)); }
    else if(_T(z)=='c'){ cmat3_set_zeros(_M(z),_N(z),_L(z),_C(z),_LD1(z),_LD2(z)); }
    else{ MATLAB_ERROR("multi_subsasgn: Unkown type"); }
    // z=x
         if(_T(z)=='r' && _T(x)=='r'){ rmat3_copy      (_M(x),_N(x),_L(x),_R(z),_LD1(z),_LD2(z),_R(x),_LD1(x),_LD2(x)); }
    else if(_T(z)=='c' && _T(x)=='r'){ cmat3_copy_rmat3(_M(x),_N(x),_L(x),_C(z),_LD1(z),_LD2(z),_R(x),_LD1(x),_LD2(x)); }
    else if(_T(z)=='c' && _T(x)=='c'){ cmat3_copy      (_M(x),_N(x),_L(x),_C(z),_LD1(z),_LD2(z),_C(x),_LD1(x),_LD2(x)); }
    else{ MATLAB_ERROR("multi_subsasgn: Unkown type"); }
    // z(s)=y
         if(_T(z)=='r' && _T(y)=='r'){ rmat3_index_copy      (s->dims[0],s->dims[1],s->dims[2],_R(z),_LD1(z),_LD2(z),_R(y),_LD1(y),_LD2(y),s->index[0],s->index[1],s->index[2]); }
    else if(_T(z)=='c' && _T(y)=='r'){ cmat3_index_copy_rmat3(s->dims[0],s->dims[1],s->dims[2],_C(z),_LD1(z),_LD2(z),_R(y),_LD1(y),_LD2(y),s->index[0],s->index[1],s->index[2]); }
    else if(_T(z)=='c' && _T(y)=='c'){ cmat3_index_copy      (s->dims[0],s->dims[1],s->dims[2],_C(z),_LD1(z),_LD2(z),_C(y),_LD1(y),_LD2(y),s->index[0],s->index[1],s->index[2]); }
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
