//////////////////////////////////////////////////////////////

//[multi_allocate_clone.c]の追加

//////////////////////////////////////////////////////////////

//[rmat3.c/rmat3.h]に追加。

/**
 @brief rmulti型の行列2つを水平方向に結合する C=[A B].
*/

int rmat3_join_vert(int m, int n, int l, rmulti **C, int LDC1, int LDC2, rmulti **A, int LDA1, int LDA2, rmulti **B, int LDB1, int LDB2)
{
  int i,j,k,e=0;
  rmulti **Z=NULL;
  Z=rmat3_allocate(LDC1,LDC2,l);
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<LDA1; i++){
	e+=rcopy(MAT3(Z,i,j,k,m,n),MAT3(A,i,j,k,LDA1,LDA2));
      }
      for(i=0; i<LDB1; i++){
	e+=rcopy(MAT3(Z,LDA1+i,j,k,m,n),MAT3(B,i,j,k,LDB1,LDB2));
      }
    }
  }
  e+=rmat3_copy(m,n,l,C,LDC1,LDC2,Z,m,n); // C=Z
  Z=rmat3_free(m,n,l,Z);
  return e;
}

//////////////////////////////////////////////////////////////

//[cmat3.c/cmat3.h]に追加。

/**
 @brief 行列2つを水平方向に結合する C=[A B].
*/

int cmat3_join_vert_r1(int m, int n, int l, cmulti **C, int LDC1, int LDC2, rmulti **A, int LDA1, int LDA2, cmulti **B, int LDB1, int LDB2)
{
  int i,j,k,e=0;
  cmulti **Z=NULL;
  Z=cmat3_allocate(LDC1,LDC2,l);
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<LDA1; i++){
	e+=ccopy_r(MAT3(Z,i,j,k,m,n),MAT3(A,i,j,k,LDA1,LDA2));
      }
      for(i=0; i<LDB1; i++){
	e+=ccopy(MAT3(Z,LDA1+i,j,k,m,n),MAT3(B,i,j,k,LDB1,LDB2));
      }
    }
  }
  e+=cmat3_copy(m,n,l,C,LDC1,LDC2,Z,m,n); // C=Z
  Z=cmat3_free(m,n,l,Z);
  return e;
}

//////////////////////////////////////////////////////////////

//[cmat3.c/cmat3.h]に追加。

/**
 @brief 行列2つを水平方向に結合する C=[A B].
*/

int cmat3_join_vert_r2(int m, int n, int l, cmulti **C, int LDC1, int LDC2, cmulti **A, int LDA1, int LDA2, rmulti **B, int LDB1, int LDB2)
{
  int i,j,k,e=0;
  cmulti **Z=NULL;
  Z=cmat3_allocate(LDC1,LDC2,l);
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<LDA1; i++){
	e+=ccopy(MAT3(Z,i,j,k,m,n),MAT3(A,i,j,k,LDA1,LDA2));
      }
      for(i=0; i<LDB1; i++){
	e+=ccopy_r(MAT3(Z,LDA1+i,j,k,m,n),MAT3(B,i,j,k,LDB1,LDB2));
      }
    }
  }
  e+=cmat3_copy(m,n,l,C,LDC1,LDC2,Z,m,n); // C=Z
  Z=cmat3_free(m,n,l,Z);
  return e;
}

//////////////////////////////////////////////////////////////

//[cmat3.c/cmat3.h]に追加。

/**
 @brief cmulti型の行列2つを水平方向に結合する C=[A B].
*/

int cmat3_join_vert(int m, int n, int l, cmulti **C, int LDC1, int LDC2, cmulti **A, int LDA1, int LDA2, cmulti **B, int LDB1, int LDB2)
{
  int i,j,k,e=0;
  cmulti **Z=NULL;
  Z=cmat3_allocate(LDC1,LDC2,l);
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<LDA1; i++){
	e+=ccopy(MAT3(Z,i,j,k,m,n),MAT3(A,i,j,k,LDA1,LDA2));
      }
      for(i=0; i<LDB1; i++){
	e+=ccopy(MAT3(Z,LDA1+i,j,k,m,n),MAT3(B,i,j,k,LDB1,LDB2));
      }
    }
  }
  e+=cmat3_copy(m,n,l,C,LDC1,LDC2,Z,m,n); // C=Z
  Z=cmat3_free(m,n,l,Z);
  return e;
}

//////////////////////////////////////////////////////////////

/**
 * @breif y=[x1 x2 ...]
 */
void multi_vertcat(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  int i=0;
  multi **x=NULL,*y=NULL,*z=NULL;
  
  x=malloc((nrhs-N0)*sizeof(multi**));
  for(i=0; i<nrhs-N0; i++){
    x[i]=multi_allocate_mxArray(prhs[N0+i]);
  }
  
  z=multi_allocate_clone(x[0]);
  for(i=1; i<nrhs-N0; i++){

    //error check
    if(_N(z)!=_N(x[i])){MATLAB_ERROR("Dimension Y is not same.");}
    if(_L(z)!=_L(x[i])){MATLAB_ERROR("Dimension Z is not same.");}
    
         if(_T(z)=='r' && _T(x[i])=='r'){ y=multi_allocate('r',_M(x[i])+_M(z),_N(x[i]),_L(x[i])); rmat3_join_vert(_M(y),_N(y),_L(y),_R(y),_M(y),_N(y),_R(z),_M(z),_N(z),_R(x[i]),_M(x[i]),_N(x[i])); }
    else if(_T(z)=='r' && _T(x[i])=='c'){ y=multi_allocate('c',_M(x[i])+_M(z),_N(x[i]),_L(x[i])); cmat3_join_vert_r1(_M(y),_N(y),_L(y),_C(y),_M(y),_N(y),_R(z),_M(z),_N(z),_C(x[i]),_M(x[i]),_N(x[i])); }
    else if(_T(z)=='c' && _T(x[i])=='r'){ y=multi_allocate('c',_M(x[i])+_M(z),_N(x[i]),_L(x[i])); cmat3_join_vert_r2(_M(y),_N(y),_L(y),_C(y),_M(y),_N(y),_C(z),_M(z),_N(z),_R(x[i]),_M(x[i]),_N(x[i])); }
    else if(_T(z)=='c' && _T(x[i])=='c'){ y=multi_allocate('c',_M(x[i])+_M(z),_N(x[i]),_L(x[i])); cmat3_join_vert(_M(y),_N(y),_L(y),_C(y),_M(y),_N(y),_C(z),_M(z),_N(z),_C(x[i]),_M(x[i]),_N(x[i])); }
    else{ MATLAB_ERROR("multi_vertcat: Unkown type"); }

    multi_free(z);
    z=multi_allocate_clone(y);
    multi_free(y);
  }
  
  // done
  plhs[0]=mxCreateStructMulti(z);
  
  for(i=0; i<nrhs-N0; i++){
    x[i]=multi_free(x[i]);
  }
  free(x);
  multi_free(z);
  return;
}

//EOF
