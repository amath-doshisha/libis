/**
 * @breif 多倍長型配列からmxArrayの構造体配列への変換
 */
mxArray *mxCreateStructMulti(multi *A)
{
  mwSize dims[3]={1,1,1};
  mxArray *ret=NULL;
  int i,j,k,t;
  dims[0]=_M(A);
  dims[1]=_N(A);
  dims[2]=_L(A);
       if(_T(A)=='r'){ ret=mxCreateStructArray(3,dims,4,r_field_names); }
  else if(_T(A)=='c'){ ret=mxCreateStructArray(3,dims,8,c_field_names); }
  else if(_T(A)=='R'){ ret=mxCreateStructArray(3,dims,8,R_field_names); }       
  else if(_T(A)=='C'){ ret=mxCreateStructArray(3,dims,16,C_field_names); }
  for(k=0; k<_L(A); k++){
    for(j=0; j<_N(A); j++){
      for(i=0; i<_M(A); i++){
	t=k*_M(A)*_N(A)+j*_M(A)+i;
	if(_T(A)=='r'){
	  mxCreateStructMultiValue(ret,t,MAT3(_R(A),i,j,k,_LD1(A),_LD2(A)),r_field_names);
	}else if(_T(A)=='c'){
	  mxCreateStructMultiValue(ret,t,C_R(MAT3(_C(A),i,j,k,_LD1(A),_LD2(A))),cr_field_names);
	  mxCreateStructMultiValue(ret,t,C_I(MAT3(_C(A),i,j,k,_LD1(A),_LD2(A))),ci_field_names);
	}else if(_T(A)=='R'){
	  mxCreateStructMultiValue(ret,t,MAT3(_R0(A),i,j,k,_LD1(A),_LD2(A)),R0_field_names);
	  mxCreateStructMultiValue(ret,t,MAT3(_R1(A),i,j,k,_LD1(A),_LD2(A)),R1_field_names);
	}else if(_T(A)=='C'){
	  mxCreateStructMultiValue(ret,t,C_R(MAT3(_C0(A),i,j,k,_LD1(A),_LD2(A))),C0r_field_names);
	  mxCreateStructMultiValue(ret,t,C_I(MAT3(_C0(A),i,j,k,_LD1(A),_LD2(A))),C0i_field_names);
	  mxCreateStructMultiValue(ret,t,C_R(MAT3(_C1(A),i,j,k,_LD1(A),_LD2(A))),C1r_field_names);
	  mxCreateStructMultiValue(ret,t,C_I(MAT3(_C1(A),i,j,k,_LD1(A),_LD2(A))),C1i_field_names);
	}else{
	  MATLAB_ERROR("mxCreateStructMulti: Unknown type.");
	}	
      }
    }
  }
  return ret;
}

//EOF
