
/**
 @brief rmulti型の3次元配列の列方向の最大値 B=rmat3_maxm(A).
*/
int rmat3_maxm(int m, int n, int l, rmulti **B, int LDB1, int LDB2, rmulti **A, int LDA1, int LDA2)
{
  int i,j,k,e=0;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      e+=rcopy(MAT3(B,0,j,k,LDB1,LDB2),MAT3(A,0,j,k,LDA1,LDA2));
      for(i=1; i<m; i++){
	if(rcmp(MAT3(B,0,j,k,LDB1,LDB2),MAT3(A,i,j,k,LDA1,LDA2))<0){ e+=rcopy(MAT3(B,0,j,k,LDB1,LDB2),MAT3(A,i,j,k,LDA1,LDA2)); }
      }
    }
  }
  return e;
}

///////////////////////////////////////////////////////////////////////////////////////

/**
 @brief rmulti型の3次元配列の行方向の最大値 B=rmat3_maxn(A).
*/
int rmat3_maxn(int m, int n, int l, rmulti **B, int LDB1, int LDB2, rmulti **A, int LDA1, int LDA2)
{
  int j,k,e=0;
  for(k=0; k<l; k++){
    e+=rcopy(MAT3(B,0,0,k,LDB1,LDB2),MAT3(A,0,0,k,LDA1,LDA2));
    for(j=1; j<n; j++){
      if(rcmp(MAT3(B,0,0,k,LDB1,LDB2),MAT3(A,0,j,k,LDA1,LDA2))<0){ e+=rcopy(MAT3(B,0,0,k,LDB1,LDB2),MAT3(A,0,j,k,LDA1,LDA2)); }
    }
  }
  return e;
}

///////////////////////////////////////////////////////////////////////////////////////

/**
 @brief rmulti型の3次元配列のL方向の最大値 B=rmat3_maxl(A).
*/
int rmat3_maxl(int m, int n, int l, rmulti **B, int LDB1, int LDB2, rmulti **A, int LDA1, int LDA2)
{
  int k,e=0;
  e+=rcopy(MAT3(B,0,0,0,LDB1,LDB2),MAT3(A,0,0,0,LDA1,LDA2));
  for(k=0; k<l; k++){
    if(rcmp(MAT3(B,0,0,0,LDB1,LDB2),MAT3(A,0,0,k,LDA1,LDA2))<0){ e+=rcopy(MAT3(B,0,0,0,LDB1,LDB2),MAT3(A,0,0,k,LDA1,LDA2)); }
  }
  return e;
}

///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////

/**
 @brief cmulti型の3次元配列の列方向の最大値 B=maxm(A).
*/
int cmat3_maxm(int m, int n, int l, cmulti **B, int LDB1, int LDB2, cmulti **A, int LDA1, int LDA2)
{
  int i,j,k,e=0,buf,prec;
  prec=cmat3_get_prec_max(m,n,l,A,LDA1,LDA2);
  rmulti* angle1=rallocate_prec(prec);
  rmulti* angle2=rallocate_prec(prec);
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      e+=ccopy(MAT3(B,0,j,k,LDB1,LDB2),MAT3(A,0,j,k,LDA1,LDA2));
      for(i=1; i<m; i++){
	buf=cabs_cmp(MAT3(B,0,j,k,LDB1,LDB2),MAT3(A,i,j,k,LDA1,LDA2));
	if(buf<0){ e+=ccopy(MAT3(B,0,j,k,LDB1,LDB2),MAT3(A,i,j,k,LDA1,LDA2)); }
	else if(buf==0){
	  //位相角が大きい方を選ぶ
	  e+=ratan2(angle1,C_I(MAT3(B,0,j,k,LDB1,LDB2)),C_R(MAT3(B,0,j,k,LDB1,LDB2)));
	  e+=ratan2(angle2,C_I(MAT3(A,i,j,k,LDA1,LDA2)),C_R(MAT3(A,i,j,k,LDA1,LDA2)));
	  if(rcmp(angle1,angle2)<0){ e+=ccopy(MAT3(B,0,j,k,LDB1,LDB2),MAT3(A,i,j,k,LDA1,LDA2)); }
	}
      }
    }
  }
  rfree(angle1);
  rfree(angle2);
  return e;
}

///////////////////////////////////////////////////////////////////////////////////////

/**
 @brief cmulti型の3次元配列の行方向の最大値 B=maxn(A).
*/
int cmat3_maxn(int m, int n, int l, cmulti **B, int LDB1, int LDB2, cmulti **A, int LDA1, int LDA2)
{
  int j,k,e=0,buf,prec;
  prec=cmat3_get_prec_max(m,n,l,A,LDA1,LDA2);
  rmulti* angle1=rallocate_prec(prec);
  rmulti* angle2=rallocate_prec(prec);
  for(k=0; k<l; k++){
    e+=ccopy(MAT3(B,0,0,k,LDB1,LDB2),MAT3(A,0,0,k,LDA1,LDA2));
    for(j=1; j<n; j++){
      buf=cabs_cmp(MAT3(B,0,0,k,LDB1,LDB2),MAT3(A,0,j,k,LDA1,LDA2));
      if(buf<0){ e+=ccopy(MAT3(B,0,0,k,LDB1,LDB2),MAT3(A,0,j,k,LDA1,LDA2)); }
      else if(buf==0){
	//位相角が大きい方を選ぶ
	e+=ratan2(angle1,C_I(MAT3(B,0,0,k,LDB1,LDB2)),C_R(MAT3(B,0,0,k,LDB1,LDB2)));
	e+=ratan2(angle2,C_I(MAT3(A,0,j,k,LDA1,LDA2)),C_R(MAT3(A,0,j,k,LDA1,LDA2)));
	if(rcmp(angle1,angle2)<0){ e+=ccopy(MAT3(B,0,0,k,LDB1,LDB2),MAT3(A,0,j,k,LDA1,LDA2)); }
      }
    }
  }
  rfree(angle1);
  rfree(angle2);
  return e;
}

///////////////////////////////////////////////////////////////////////////////////////

/**
 @brief cmulti型の3次元配列のL方向の最大値 B=maxl(A).
*/
int cmat3_maxl(int m, int n, int l, cmulti **B, int LDB1, int LDB2, cmulti **A, int LDA1, int LDA2)
{
  int k,e=0,buf,prec;
  prec=cmat3_get_prec_max(m,n,l,A,LDA1,LDA2);
  rmulti* angle1=rallocate_prec(prec);
  rmulti* angle2=rallocate_prec(prec);
  e+=ccopy(MAT3(B,0,0,0,LDB1,LDB2),MAT3(A,0,0,0,LDA1,LDA2));
  for(k=0; k<l; k++){
    buf=cabs_cmp(MAT3(B,0,0,0,LDB1,LDB2),MAT3(A,0,0,k,LDA1,LDA2));
    if(buf<0){ e+=ccopy(MAT3(B,0,0,0,LDB1,LDB2),MAT3(A,0,0,k,LDA1,LDA2)); }
    else if(buf==0){
      //位相角が大きい方を選ぶ
      e+=ratan2(angle1,C_I(MAT3(B,0,0,0,LDB1,LDB2)),C_R(MAT3(B,0,0,0,LDB1,LDB2)));
      e+=ratan2(angle2,C_I(MAT3(A,0,0,k,LDA1,LDA2)),C_R(MAT3(A,0,0,k,LDA1,LDA2)));
      if(rcmp(angle1,angle2)<0){ e+=ccopy(MAT3(B,0,0,0,LDB1,LDB2),MAT3(A,0,0,k,LDA1,LDA2)); }
    }
  }
  rfree(angle1);
  rfree(angle2);
  return e;
}

///////////////////////////////////////////////////////////////////////////////////////

/**
 * @breif z=max(x)
 */
void multi_max(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  multi *x=NULL,*y=NULL,*z=NULL,*u=NULL,*v=NULL;
  int info=-1;
  if(nlhs>1){ mexErrMsgIdAndTxt("MATLAB:multi_mex:maxlhs","Too many output arguments."); }
  if(!(IS_STRT(nrhs,prhs,N0))){ MATLAB_ERROR("multi_max: The 1st-arg should be Struct."); }
  // allocate by clone
  x=multi_allocate_mxArray(prhs[N0]);
  // operation
  if(_M(x)!=1){
    //press m
         if(_T(x)=='r'){ z=multi_allocate('r',1,_N(x),_L(x)); rmat3_maxm(_M(x),_N(x),_L(x),_R(z),_LD1(z),_LD2(z),_R(x),_LD1(x),_LD2(x)); }
    else if(_T(x)=='c'){ z=multi_allocate('c',1,_N(x),_L(x)); cmat3_maxm(_M(x),_N(x),_L(x),_C(z),_LD1(z),_LD2(z),_C(x),_LD1(x),_LD2(x)); }
  }else if(_N(x)!=1){
    //press n
         if(_T(x)=='r'){ z=multi_allocate('r',1,1,_L(x)); rmat3_maxn(_M(x),_N(x),_L(x),_R(z),_LD1(z),_LD2(z),_R(x),_LD1(x),_LD2(x)); }
    else if(_T(x)=='c'){ z=multi_allocate('c',1,1,_L(x)); cmat3_maxn(_M(x),_N(x),_L(x),_C(z),_LD1(z),_LD2(z),_C(x),_LD1(x),_LD2(x)); }
  }else if(_L(x)!=1){
    //press l
         if(_T(x)=='r'){ z=multi_allocate('r',1,1,1); rmat3_maxl(_M(x),_N(x),_L(x),_R(z),_LD1(z),_LD2(z),_R(x),_LD1(x),_LD2(x)); }
    else if(_T(x)=='c'){ z=multi_allocate('c',1,1,1); cmat3_maxl(_M(x),_N(x),_L(x),_C(z),_LD1(z),_LD2(z),_C(x),_LD1(x),_LD2(x)); }
  }else if(_M(x)==1 && _N(x)==1 && _L(x)==1){
    //copy z=x
         if(_T(x)=='r'){ z=multi_allocate('r',_M(x),_N(x),_L(x)); rmat3_copy(_M(z),_N(z),_L(z),_R(z),_M(z),_N(z),_R(x),_M(x),_N(x)); }
    else if(_T(x)=='c'){ z=multi_allocate('c',_M(x),_N(x),_L(x)); cmat3_copy(_M(z),_N(z),_L(z),_C(z),_M(z),_N(z),_C(x),_M(x),_N(x)); }
  }else{
    MATLAB_ERROR("multi_max: z=max(x): This function is not implemented");
  }

  plhs[0]=mxCreateStructMulti(z);
  x=multi_free(x);
  z=multi_free(z);
  return;
}

//EOF


