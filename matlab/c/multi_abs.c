
/**
 @brief rmulti型の3次元配列の絶対値 C=abs(A).
*/
int rmat3_abs(int m, int n, int l, rmulti **C, int LDC1, int LDC2, rmulti **A, int LDA1, int LDA2){
  int i,j,k,e=0;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	e+=rabs(MAT3(C,i,j,k,LDC1,LDC2),MAT3(A,i,j,k,LDA1,LDA2));
      }
    }
  }
  return e;
}

/////////////////////////////////////////////////////////////////////////////////////////////////

/**
 @brief cmulti型の3次元配列の絶対値 C=abs(A).
*/
int cmat3_abs(int m, int n, int l, rmulti **C, int LDC1, int LDC2, cmulti **A, int LDA1, int LDA2){
  int i,j,k,e=0;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	e+=cabsv(MAT3(C,i,j,k,LDC1,LDC2),MAT3(A,i,j,k,LDA1,LDA2));
      }
    }
  }
  return e;
}

/////////////////////////////////////////////////////////////////////////////////////////////////

/**
 * @breif z=abs(x,y)
 */
void multi_abs(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  multi *x=NULL,*y=NULL,*z=NULL,*u=NULL,*v=NULL;
  int info=-1;
  if(nlhs>1){ mexErrMsgIdAndTxt("MATLAB:multi_mex:maxlhs","Too many output arguments."); }
  if(!(IS_STRT(nrhs,prhs,N0))){ MATLAB_ERROR("multi_abs: The 1st-arg should be Struct."); }
  // allocate by clone
  x=multi_allocate_mxArray(prhs[N0]);
  // operation
       if(_T(x)=='r'){ z=multi_allocate('r',_M(x),_N(x),_L(x)); rmat3_abs(_M(z),_N(z),_L(z),_R(z),_LD1(z),_LD2(z),_R(x),_LD1(x),_LD2(x)); }
  else if(_T(x)=='c'){ z=multi_allocate('r',_M(x),_N(x),_L(x)); cmat3_abs(_M(z),_N(z),_L(z),_R(z),_LD1(z),_LD2(z),_C(x),_LD1(x),_LD2(x)); }
  else{ MATLAB_ERROR("multi_abs: Unkown type"); }
  plhs[0]=mxCreateStructMulti(z);
  x=multi_free(x);
  z=multi_free(z);
  return;
}

//EOF


