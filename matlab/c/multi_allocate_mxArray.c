/**
 * @brief MATLABの構造体から多倍長型のメモリの割当
 */
multi *multi_allocate_mxArray(const mxArray *x)
{
  char type='r';
  int m=1,n=1,l=1,i,j,k,t;
  multi *A=NULL;
  // get type
  if(mxIsStruct(x)){
         if(mxGetFieldNameByNumber(x,0)!=NULL && STR_EQ(mxGetFieldNameByNumber(x,0),r_field_names[0])){ type='r'; }
    else if(mxGetFieldNameByNumber(x,0)!=NULL && STR_EQ(mxGetFieldNameByNumber(x,0),c_field_names[0])){ type='c'; }
    else if(mxGetFieldNameByNumber(x,0)!=NULL && STR_EQ(mxGetFieldNameByNumber(x,0),R_field_names[0])){ type='R'; }
    else if(mxGetFieldNameByNumber(x,0)!=NULL && STR_EQ(mxGetFieldNameByNumber(x,0),C_field_names[0])){ type='C'; }
    else { MATLAB_ERROR("multi_allocate_mxArray: The argument should be Struct of multi."); }
  }else if(mxIsDouble(x)){
    if(mxIsComplex(x)){ type='z'; }
    else              { type='d'; }
  }else if(mxIsCell(x)){
    type='s';
  }else{ MATLAB_ERROR("multi_allocate_mxArray: Unsupported type.");}
  // get size
  if(mxGetNumberOfDimensions(x)>=1){ m=mxGetDimensions(x)[0]; }
  if(mxGetNumberOfDimensions(x)>=2){ n=mxGetDimensions(x)[1]; }
  if(mxGetNumberOfDimensions(x)>=3){ l=mxGetDimensions(x)[2]; }
  if(mxGetNumberOfDimensions(x)>=4){ MATLAB_ERROR("multi_allocate_mxArray: Number of dimensions is too big."); }
  if(m==0 && n==0){ l=0; }
  // allocate
  A=multi_allocate(type,m,n,l);
  // get data
  for(k=0; k<_L(A); k++){
    for(j=0; j<_N(A); j++){
      for(i=0; i<_M(A); i++){
	t=k*_M(A)*_N(A)+j*_M(A)+i;
             if(_T(A)=='r'){ mxArrayToMulti(    MAT3(_R(A), i,j,k,_LD1(A),_LD2(A)), x,t,r_field_names); }
        else if(_T(A)=='c'){ mxArrayToMulti(C_R(MAT3(_C(A), i,j,k,_LD1(A),_LD2(A))),x,t,cr_field_names);
                             mxArrayToMulti(C_I(MAT3(_C(A), i,j,k,_LD1(A),_LD2(A))),x,t,ci_field_names); }
	else if(_T(A)=='R'){ mxArrayToMulti(    MAT3(_R0(A),i,j,k,_LD1(A),_LD2(A)), x,t,R0_field_names);
	                     mxArrayToMulti(    MAT3(_R1(A),i,j,k,_LD1(A),_LD2(A)), x,t,R1_field_names); }
	else if(_T(A)=='C'){ mxArrayToMulti(C_R(MAT3(_C0(A),i,j,k,_LD1(A),_LD2(A))),x,t,C0r_field_names);
	                     mxArrayToMulti(C_I(MAT3(_C0(A),i,j,k,_LD1(A),_LD2(A))),x,t,C0i_field_names);
                             mxArrayToMulti(C_R(MAT3(_C1(A),i,j,k,_LD1(A),_LD2(A))),x,t,C1r_field_names);
                             mxArrayToMulti(C_I(MAT3(_C1(A),i,j,k,_LD1(A),_LD2(A))),x,t,C1i_field_names); }
	else if(_T(A)=='d'){ MAT3(_D(A),i,j,k,_LD1(A),_LD2(A))=mxGetPr(x)[t]; }
	else if(_T(A)=='z'){ Z_SET(MAT3(_Z(A),i,j,k,_LD1(A),_LD2(A)),mxGetPr(x)[t],mxGetPi(x)[t]); }
	else if(_T(A)=='s'){
	  if(mxGetCell(x,t)!=NULL && mxIsChar(mxGetCell(x,t))){ MAT3(_S(A),i,j,k,_LD1(A),_LD2(A))=mxArrayToString(mxGetCell(x,t)); }
	  else{ MATLAB_ERROR("multi_allocate_mxArray: The entries of cell should be Char."); }
	}
	else{ MATLAB_ERROR("multi_allocate_mxArray: Unsupported type.");}
      }
    }
  }
  return A;
}

//EOF
