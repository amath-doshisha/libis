/**
 * @breif y=x(s)
 */
void multi_subsref(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  multi *x=NULL,*y=NULL;
  subs_index_t *s=NULL;
  if(nlhs>1){ mexErrMsgIdAndTxt("MATLAB:multi_mex:maxlhs","Too many output arguments."); }
  if(!(IS_STRT(nrhs,prhs,N0))){ MATLAB_ERROR("multi_subsref: The 1st-arg should be Struct."); }
  if(!(IS_STRT(nrhs,prhs,N0+1))){ MATLAB_ERROR("multi_subsref: The 2nd-arg should be Struct."); }
  // load the 1st-arg
  x=multi_allocate_mxArray(prhs[N0]);
  // load the 2nd-arg
  s=subs_index_allocate(prhs[N0+1]);
  // copy by indecies
  if(s->ndim==1 && s->dims[0]==0 && s->index[0]==NULL){ // A(:)
    // convert to column vector
    y=multi_allocate(_T(x),_LD1(x)*_LD2(x)*_L(x),1,1);
         if(_T(x)=='r'){ rmat3_copy_rvec(_LD1(x),_LD2(x),_L(x),_R(y),_R(x),_LD1(x),_LD2(x)); }
    else if(_T(x)=='c'){ cmat3_copy_rvec(_LD1(x),_LD2(x),_L(x),_C(y),_C(x),_LD1(x),_LD2(x)); }
    else{ MATLAB_ERROR("multi_subsref: Unkown type"); }
  }else if(s->ndim==1){ // A(*)
    if(s->dims[0]==0 && s->index[0]==NULL){ s->dims[0]=_M(x)*_N(x)*_L(x); s->index[0]=ivec_allocate(s->dims[0]); ivec_grid(s->dims[0],s->index[0]); }
    if(ivec_min(s->dims[0],s->index[0])<0)                 { MATLAB_ERROR("multi_subsref: The index is out of the dimensions of the matrix."); }
    if(ivec_max(s->dims[0],s->index[0])>=_M(x)*_N(x)*_L(x)){ MATLAB_ERROR("multi_subsref: The index is out of the dimensions of the matrix."); }
    if(_N(x)==1 && _L(x)==1){ y=multi_allocate(_T(x),s->dims[0],1,1); }
    else                    { y=multi_allocate(_T(x),1,s->dims[0],1); }
         if(_T(x)=='r'){ rmat3_copy_rvec_index(_R(y),_R(x),_LD1(x),_LD2(x),s->dims[0],s->index[0]); }
    else if(_T(x)=='c'){ cmat3_copy_rvec_index(_C(y),_C(x),_LD1(x),_LD2(x),s->dims[0],s->index[0]); }
    else{ MATLAB_ERROR("multi_subsref: Unkown type"); }	       
  }else if(s->ndim==2){ // A(*,*)
    if(s->dims[0]==0 && s->index[0]==NULL){ s->dims[0]=_M(x);       s->index[0]=ivec_allocate(s->dims[0]); ivec_grid(s->dims[0],s->index[0]); }
    if(s->dims[1]==0 && s->index[1]==NULL){ s->dims[1]=_N(x)*_L(x); s->index[1]=ivec_allocate(s->dims[1]); ivec_grid(s->dims[1],s->index[1]); }
    if(ivec_min(s->dims[0],s->index[0])<0)           { MATLAB_ERROR("multi_subsref: The index is out of the dimensions of the matrix."); }
    if(ivec_min(s->dims[1],s->index[1])<0)           { MATLAB_ERROR("multi_subsref: The index is out of the dimensions of the matrix."); }
    if(ivec_max(s->dims[0],s->index[0])>=_M(x))      { MATLAB_ERROR("multi_subsref: The index is out of the dimensions of the matrix."); }
    if(ivec_max(s->dims[1],s->index[1])>=_N(x)*_L(x)){ MATLAB_ERROR("multi_subsref: The index is out of the dimensions of the matrix."); }
    y=multi_allocate(_T(x),s->dims[0],s->dims[1],1);
         if(_T(x)=='r'){ rmat3_copy_rmat_index(_R(y),_LD1(y),_R(x),_LD1(x),_LD2(x),s->dims[0],s->index[0],s->dims[1],s->index[1]); }
    else if(_T(x)=='c'){ cmat3_copy_cmat_index(_C(y),_LD1(y),_C(x),_LD1(x),_LD2(x),s->dims[0],s->index[0],s->dims[1],s->index[1]); }
    else{ MATLAB_ERROR("multi_subsref: Unkown type"); }	           
  }else if(s->ndim==3){ // A(*,*,*)
    if(s->dims[0]==0 && s->index[0]==NULL){ s->dims[0]=_M(x); s->index[0]=ivec_allocate(s->dims[0]); ivec_grid(s->dims[0],s->index[0]); }
    if(s->dims[1]==0 && s->index[1]==NULL){ s->dims[1]=_N(x); s->index[1]=ivec_allocate(s->dims[1]); ivec_grid(s->dims[1],s->index[1]); }
    if(s->dims[2]==0 && s->index[2]==NULL){ s->dims[2]=_L(x); s->index[2]=ivec_allocate(s->dims[2]); ivec_grid(s->dims[2],s->index[2]); }    
    if(ivec_min(s->dims[0],s->index[0])<0)     { MATLAB_ERROR("multi_subsref: The index is out of the dimensions of the matrix."); }
    if(ivec_min(s->dims[1],s->index[1])<0)     { MATLAB_ERROR("multi_subsref: The index is out of the dimensions of the matrix."); }
    if(ivec_min(s->dims[2],s->index[2])<0)     { MATLAB_ERROR("multi_subsref: The index is out of the dimensions of the matrix."); }
    if(ivec_max(s->dims[0],s->index[0])>=_M(x)){ MATLAB_ERROR("multi_subsref: The index is out of the dimensions of the matrix."); }
    if(ivec_max(s->dims[1],s->index[1])>=_N(x)){ MATLAB_ERROR("multi_subsref: The index is out of the dimensions of the matrix."); }
    if(ivec_max(s->dims[2],s->index[2])>=_L(x)){ MATLAB_ERROR("multi_subsref: The index is out of the dimensions of the matrix."); }
    y=multi_allocate(_T(x),s->dims[0],s->dims[1],s->dims[2]);
         if(_T(x)=='r'){ rmat3_copy_index(_M(y),_N(y),_L(y),_R(y),_LD1(y),_LD2(y),_R(x),_LD1(x),_LD2(x),s->index[0],s->index[1],s->index[2]); }
    else if(_T(x)=='c'){ cmat3_copy_index(_M(y),_N(y),_L(y),_C(y),_LD1(y),_LD2(y),_C(x),_LD1(x),_LD2(x),s->index[0],s->index[1],s->index[2]); }
    else{ MATLAB_ERROR("multi_subsref: NOT SUPPORTED!"); }
  }else{ // otherwise
    MATLAB_ERROR("multi_subsref: NOT SUPPORTED!");
  }
  // done
  plhs[0]=mxCreateStructMulti(y);  
  subs_index_free(s);
  x=multi_free(x);
  y=multi_free(y);
  return;
}

//EOF
