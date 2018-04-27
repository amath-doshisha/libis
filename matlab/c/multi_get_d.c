/**
 * @brief double型へのキャスト y=double(x)
 */
void multi_get_d(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  mwSize ndim=3,dims[]={1,1,1};
  int m=1,n=1,l=1,i,j,k,t;
  multi *A=NULL;
  dcomplex z,z0,z1;
  if(nlhs>1){ mexErrMsgIdAndTxt("MATLAB:multi_mex:maxlhs","Too many output arguments."); }
  if(!(IS_STRT(nrhs,prhs,N0))){ MATLAB_ERROR("multi_get_d: The 1st-arg should be Struct."); }
  A=multi_allocate_mxArray(prhs[N0]);
  dims[0]=_M(A);  dims[1]=_N(A);  dims[2]=_L(A);  
  if(dims[2]==1){ ndim=2; }
       if(_T(A)=='r' || _T(A)=='R'){ plhs[0]=mxCreateNumericArray(ndim,dims,mxDOUBLE_CLASS,mxREAL); }
  else if(_T(A)=='c' || _T(A)=='C'){ plhs[0]=mxCreateNumericArray(ndim,dims,mxDOUBLE_CLASS,mxCOMPLEX); }
  else                             { MATLAB_ERROR("multi_get_d: Unknow type."); }
  for(k=0; k<_L(A); k++){
    for(j=0; j<_N(A); j++){
      for(i=0; i<_M(A); i++){
	t=i+j*_M(A)+k*_M(A)*_N(A);
	if(_T(A)=='r'){
	  mxGetPr(plhs[0])[t]=rget_d(MAT3(_R(A),i,j,k,_LD1(A),_LD2(A)));
	}else if(_T(A)=='R'){
	  mxGetPr(plhs[0])[t]=0.5*(rget_d(MAT3(_R0(A),i,j,k,_LD1(A),_LD2(A)))+rget_d(MAT3(_R1(A),i,j,k,_LD1(A),_LD2(A))));
	}else if(_T(A)=='c'){
	  z=cget_z(MAT3(_C(A),i,j,k,_LD1(A),_LD2(A)));
	  mxGetPr(plhs[0])[t]=Z_R(z);
	  mxGetPi(plhs[0])[t]=Z_I(z);
	}else if(_T(A)=='C'){
	  z0=cget_z(MAT3(_C0(A),i,j,k,_LD1(A),_LD2(A)));
	  z1=cget_z(MAT3(_C1(A),i,j,k,_LD1(A),_LD2(A)));
	  mxGetPr(plhs[0])[t]=0.5*(Z_R(z0)+Z_R(z1));
	  mxGetPi(plhs[0])[t]=0.5*(Z_I(z0)+Z_I(z1));
	}else{
	  MATLAB_ERROR("multi_get_d: Unknow type.");
	}
      }
    }
  } 
  A=multi_free(A);
}

//EOF
