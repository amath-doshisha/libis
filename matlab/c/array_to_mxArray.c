/**
 * @brief 多倍長型配列からmxArrayの構造体配列への変換
 */

mxArray *array_to_mxArray(array *A)
{
  mwSize ndim,*dim=NULL;
  mxArray *ret=NULL;
  int i;
  if(A==NULL || (A!=NULL && ARRAY_NDIM(A)<=0)){
    ndim=2;
    dim=(mwSize*)mxMalloc(sizeof(mwSize)*ndim);
    for(i=0; i<ndim; i++){ dim[i]=0; }
    ret=mxCreateStructArray(ndim,dim,4,r_field_names);
  }else{
    // get size
    ndim=ARRAY_NDIM(A);
    dim=(mwSize*)mxMalloc(sizeof(mwSize)*ndim);
    for(i=0; i<ndim; i++){ dim[i]=ARRAY_DIM(A,i); }
    // convert
         if(_T(A)=='r'){ ret=mxCreateStructArray(ndim,dim,4,r_field_names); }
    else if(_T(A)=='c'){ ret=mxCreateStructArray(ndim,dim,8,c_field_names); }
    else if(_T(A)=='R'){ ret=mxCreateStructArray(ndim,dim,8,R_field_names); }       
    else if(_T(A)=='C'){ ret=mxCreateStructArray(ndim,dim,16,C_field_names); }
    else if(_T(A)=='d'){ ret=mxCreateNumericArray(ndim,dim,mxDOUBLE_CLASS,mxREAL); }
    else if(_T(A)=='z'){ ret=mxCreateNumericArray(ndim,dim,mxDOUBLE_CLASS,mxCOMPLEX); }
    else if(_T(A)=='i'){ ret=mxCreateNumericArray(ndim,dim,mxINT32_CLASS,mxREAL); }
    else{	MATLAB_ERROR("array_to_mxArray: Unknown type."); }
    for(i=0; i<ARRAY_SIZE(A); i++){
      if(_T(A)=='r'){
	rmulti_to_mxArray(ret,i,ARRAY_RVEC(A,i),r_field_names);
      }else if(_T(A)=='c'){
	rmulti_to_mxArray(ret,i,C_R(ARRAY_CVEC(A,i)),cr_field_names);
	rmulti_to_mxArray(ret,i,C_I(ARRAY_CVEC(A,i)),ci_field_names);
      }else if(_T(A)=='R'){
	rmulti_to_mxArray(ret,i,ARRAY_RVEC0(A,i),R0_field_names);
	rmulti_to_mxArray(ret,i,ARRAY_RVEC1(A,i),R1_field_names);
      }else if(_T(A)=='C'){
	rmulti_to_mxArray(ret,i,C_R(ARRAY_CVEC0(A,i)),C0r_field_names);
	rmulti_to_mxArray(ret,i,C_I(ARRAY_CVEC0(A,i)),C0i_field_names);
	rmulti_to_mxArray(ret,i,C_R(ARRAY_CVEC1(A,i)),C1r_field_names);
	rmulti_to_mxArray(ret,i,C_I(ARRAY_CVEC1(A,i)),C1i_field_names);
      }else if(_T(A)=='d'){
	mxGetPr(ret)[i]=ARRAY_DVEC(A,i);
      }else if(_T(A)=='z'){
	mxGetPr(ret)[i]=Z_R(ARRAY_ZVEC(A,i));
	mxGetPi(ret)[i]=Z_I(ARRAY_ZVEC(A,i));
      }else if(_T(A)=='i'){
	((int*)mxGetData(ret))[i]=ARRAY_IVEC(A,i);
      }else{
	MATLAB_ERROR("array_to_mxArray: Unknown type.");
      }
    }
  }
  // done
  mxFree(dim);
  return ret;
}

//EOF
