/**
 * @breif double型へのキャスト y=double(x)
 */
void multi_get_d(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  mwSize ndim=3,dims[]={1,1,1};
  int m=1,n=1,l=1,i,j,k,t,s;
  multi *A=NULL;
  dcomplex z;
  if(nlhs>1){ mexErrMsgIdAndTxt("MATLAB:multi_mex:maxlhs","Too many output arguments."); }
  if(!(IS_STRT(nrhs,prhs,N0))){ MATLAB_ERROR("multi_get_d: The 1st-arg should be Struct."); }
  A=multi_allocate_mxArray(prhs[N0]);
  dims[0]=_M(A);  dims[1]=_N(A);  dims[2]=_L(A);  
  if(_T(A)=='R' || _T(A)=='C'){ dims[1]*=2; }
  if(dims[2]==1){ ndim=2; }
       if(_T(A)=='r' || _T(A)=='R'){ plhs[0]=mxCreateNumericArray(ndim,dims,mxDOUBLE_CLASS,mxREAL); }
  else if(_T(A)=='c' || _T(A)=='C'){ plhs[0]=mxCreateNumericArray(ndim,dims,mxDOUBLE_CLASS,mxCOMPLEX); }
  else                             { MATLAB_ERROR("multi_get_d: Unknow type."); }
  for(k=0; k<_L(A); k++){
    for(j=0; j<_N(A); j++){
      for(i=0; i<_M(A); i++){
	if(_T(A)=='r'||_T(A)=='c'){ t=i+j*_M(A)+k*_M(A)*_N(A); }
	else if(_T(A)=='R'||_T(A)=='C'){ t=i+2*j*_M(A)+k*_M(A)*_N(A)*2; s=i+(2*j+1)*_M(A)+k*_M(A)*_N(A)*2; }
	if(_T(A)=='r'){
	  mxGetPr(plhs[0])[t]=rget_d(MAT3(_R(A),i,j,k,_LD1(A),_LD2(A)));
	}else if(_T(A)=='c'){
	  z=cget_z(MAT3(_C(A),i,j,k,_LD1(A),_LD2(A)));
	  mxGetPr(plhs[0])[t]=Z_R(z); mxGetPi(plhs[0])[t]=Z_I(z);
	}else if(_T(A)=='R'){
	  mxGetPr(plhs[0])[t]=rget_d(MAT3(_R0(A),i,j,k,_LD1(A),_LD2(A)));
	  mxGetPr(plhs[0])[s]=rget_d(MAT3(_R1(A),i,j,k,_LD1(A),_LD2(A)));
	}else if(_T(A)=='C'){
	  z=cget_z(MAT3(_C0(A),i,j,k,_LD1(A),_LD2(A)));
	  mxGetPr(plhs[0])[t]=Z_R(z);
	  mxGetPi(plhs[0])[t]=Z_I(z);
	  z=cget_z(MAT3(_C1(A),i,j,k,_LD1(A),_LD2(A)));
	  mxGetPr(plhs[0])[s]=Z_R(z);
	  mxGetPi(plhs[0])[s]=Z_I(z);
	}else{
	  MATLAB_ERROR("multi_get_d: Unknow type.");
	}
      }
    }
  } 
  A=multi_free(A);
}

//EOF
