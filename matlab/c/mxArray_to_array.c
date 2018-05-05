/**
 * @brief MATLABの構造体から多倍長型のメモリの割当
 */

array *mxArray_to_array(const mxArray *x)
{
  char type='r',*buf=NULL;
  int ndim,*dim=NULL,i,m,n,l;
  array *A=NULL;  
  // get type
  if(mxIsStruct(x)){
         if(mxGetFieldNameByNumber(x,0)!=NULL && STR_EQ(mxGetFieldNameByNumber(x,0),r_field_names[0])){ type='r'; }
    else if(mxGetFieldNameByNumber(x,0)!=NULL && STR_EQ(mxGetFieldNameByNumber(x,0),c_field_names[0])){ type='c'; }
    else if(mxGetFieldNameByNumber(x,0)!=NULL && STR_EQ(mxGetFieldNameByNumber(x,0),R_field_names[0])){ type='R'; }
    else if(mxGetFieldNameByNumber(x,0)!=NULL && STR_EQ(mxGetFieldNameByNumber(x,0),C_field_names[0])){ type='C'; }
    else { MATLAB_ERROR("mxArray_to_array: The argument should be Struct of multi."); }
  }else if(mxIsDouble(x)){
    if(mxIsComplex(x)){ type='z'; }
    else              { type='d'; }
  }else if(mxIsCell(x)){
    type='s';
  }else{ MATLAB_ERROR("mxArray_to_array: Unsupported type.");}
  // get size
  ndim=mxGetNumberOfDimensions(x);
  dim=ivec_allocate(ndim);
  for(i=0; i<ndim; i++){ dim[i]=mxGetDimensions(x)[i]; }
  // allocate
  A=array_allocate(type,ndim,dim);
  // get data
  for(i=0; i<ARRAY_SIZE(A); i++){
    if(ARRAY_TYPE(A)=='r'){
      mxArray_to_rmulti(ARRAY_RVEC(A,i),x,i,r_field_names);
    }else if(ARRAY_TYPE(A)=='c'){
      mxArray_to_rmulti(C_R(ARRAY_CVEC(A,i)),x,i,cr_field_names);
      mxArray_to_rmulti(C_I(ARRAY_CVEC(A,i)),x,i,ci_field_names);
    }else if(ARRAY_TYPE(A)=='R'){
      mxArray_to_rmulti(ARRAY_RVEC0(A,i),x,i,R0_field_names);
      mxArray_to_rmulti(ARRAY_RVEC1(A,i),x,i,R1_field_names);
    }else if(ARRAY_TYPE(A)=='C'){
      mxArray_to_rmulti(C_R(ARRAY_CVEC0(A,i)),x,i,C0r_field_names);
      mxArray_to_rmulti(C_I(ARRAY_CVEC0(A,i)),x,i,C0i_field_names);
      mxArray_to_rmulti(C_R(ARRAY_CVEC1(A,i)),x,i,C1r_field_names);
      mxArray_to_rmulti(C_I(ARRAY_CVEC1(A,i)),x,i,C1i_field_names);
    }else if(ARRAY_TYPE(A)=='d'){
      ARRAY_DVEC(A,i)=mxGetPr(x)[i];
    }else if(_T(A)=='z'){
      Z_SET(ARRAY_ZVEC(A,i),mxGetPr(x)[i],mxGetPi(x)[i]);
    }else if(ARRAY_TYPE(A)=='s'){
      if(mxGetCell(x,i)!=NULL && mxIsChar(mxGetCell(x,i))){
	buf=mxArrayToString(mxGetCell(x,i));
	ARRAY_SVEC(A,i)=char_new(buf,NULL);
	mxFree(buf); buf=NULL;
      }else{ MATLAB_ERROR("mxArray_to_array: The entries of cell should be Char."); }
    }
    else{ MATLAB_ERROR("mxArray_to_array: Unsupported type.");}
  }
  dim=ivec_free(dim);
  return A;
}

//EOF
