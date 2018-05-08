#include<stdio.h>
#include<stdlib.h>
#include<isys.h>

/**
 @file  riep_dht.c
 @brief 多倍長精度実数型rmultiにおいて，dhToda方程式の逆固有値問題により指定した固有値をもつ行列の生成.
 @details rmulti型の関数に関する定義は@link rmulti.c@endlinkを参照のこと.
          rmulti型のベクトルの関数に関する定義は@link rvec.c@endlinkを参照のこと.
          rmulti型の行列の関数に関する定義は@link rmat.c@endlinkを参照のこと.
          rmulti型の固有値問題に関する定義は@link reig.c@endlinkを参照のこと.
 */

/**
 @brief QEより行列Aを生成
*/
void riep_dhToda_QE_to_A(int m, int M, rmulti **A, int LDA, rmulti **Q, int LDQ, rmulti **E, int debug)
{
  // init
  int k=0,n=0,l=0,LDR=0,prec=0;
  rmulti **R=NULL;

  // precision
  prec=rmat_get_prec_max(m,m,A,LDA);

  //allocate
  LDR=m; R=rmat_allocate_prec(LDR,m,prec);

  // set A=L
  for(n=0; n<m; n++){
    for(k=0; k<m; k++){
      if     (n==k)    { rset_one(MAT(A,n,k,LDA)); }      // A[n][k]=1
      else if((n-k)==1){ rcopy(MAT(A,n,k,LDA),E[k]); }    // A[n][k]=E[k][0]
      else             { rset_zero(MAT(A,n,k,LDA)); }     // A[n][k]=0
    }
  }
  // set A=A*R
  for(l=M-1;l>=0;l--){
    for(n=0;n<m;n++){
      for(k=0;k<m;k++){
	if     ((k-n)==1){ rset_one(MAT(R,n,k,LDR)); }             // R[n][k]=1
	else if(n==k)    { rcopy(MAT(R,n,k,LDR),MAT(Q,k,l,LDQ)); } // R[n][k]=Q[k][l]
	else             { rset_zero(MAT(R,n,k,LDR)); }            // R[n][k]=0
      }
    }
    rmat_prod(m,m,m,A,LDA,A,LDA,R,LDR); // A=A*R
  }
  // debug
  if(debug>0){

  }
  //done
  R=rmat_free(LDR,m,R);

}

/**
 @brief 指定固有値よりdhTodaによりTN行列を生成．漸化式は方程式．
 @param[in]      m        行列のサイズ，
 @param[in]      M        M>=1の整数．帯幅は min(M+2,m)．
 @param[in,out]  A        [in]サイズが(m,m)の行列．出力結果の格納用に初期化済みのA，またはA=NULLのときは出力はされない．[out]計算結果．
 @param[in]      LDA      Aの第1次元．
 @param[in,out]  Q        サイズが(m,M)の行列．出力結果の格納用に初期化済みのQ，またはQ=NULLのときは出力はされない．[out]計算結果．
 @param[in]      LDQ      Qの第1次元．
 @param[in,out]  E        サイズがmのベクトル．出力結果の格納用に初期化済みのE，またはE=NULLのときは出力はされない．[out]計算結果．
 @param[in]      lambda   サイズが(m,m)の行列．
 @param[in]      c        サイズが(m,m)の行列．
 @param[in]      debug    デバグレベル．
*/

void riep_dhToda_TN(int m, int M, rmulti **A, int LDA, rmulti **Q, int LDQ, rmulti **E, rmulti **lambda, rmulti **c, int debug)
{
  int prec=0,k=0,n=0,f_size=0,n_size=0,*Q0_size=NULL,*E0_size=NULL,LDQ1;
  rmulti *a=NULL,**f=NULL,***Q0=NULL,***E0=NULL,**sigma=NULL,**Q1=NULL,**E1=NULL;

  // init
  n_size=(M+1)*(m-1)+2*M;
  // precision
  prec=rmat_get_prec_max(m,m,A,LDA);
  
  // allocate
  if(Q==NULL){ LDQ1=m; Q1=rmat_allocate_prec(m,M,prec); }else{ LDQ1=LDQ; Q1=Q; }
  if(E==NULL){ E1=rvec_allocate_prec(m,prec); }else{ E1=E; }
  sigma=rvec_allocate_prec(m,prec);
  a=rallocate_prec(prec);
  Q0_size=ivec_allocate(m);
  Q0=malloc(m*sizeof(rmulti**));
  for(k=0; k<m; k++){
    Q0_size[k]=n_size-k*(M+1);
    Q0[k]=rvec_allocate_prec(Q0_size[k],prec);
  }
  E0_size=ivec_allocate(m);
  E0=malloc(m*sizeof(rmulti**));
  for(k=0; k<m; k++){
    E0_size[k]=n_size-(k+1)*(M+1)+1;
    E0[k]=rvec_allocate_prec(E0_size[k],prec);
  }
  f_size=Q0_size[0]+1; f=rvec_allocate_prec(f_size,prec);

  // generate sigma[n]
  rinv_d(a,M);
  rvec_pow_r(m,sigma,lambda,a);
  // generate f[n]
  for(n=0; n<f_size; n++){
    rset_d(f[n],0);
    for(k=0; k<m; k++){
      rpow_si(a,sigma[k],n); // a=(sigma[i])^n
      radd_mul(f[n],c[k],a); // f[i]=f[i]+c[i]*(sigma[i])^n
    }
  }
  // Q[n][0]=f[n+1]/f[n]
  for(n=0; n<Q0_size[0]; n++){
    if(n+1<f_size){ rdiv(Q0[0][n],f[n+1],f[n]); }
    else          { rset_nan(Q0[0][n]); }
  }
  // E[0][n]=Q[0][n+M]-Q[0][n];
  k=0;
  for(n=0; n<E0_size[k]; n++){
    if(n+M<Q0_size[k] && n<Q0_size[k]){ rsub(E0[k][n],Q0[k][n+M],Q0[k][n]); }
    else                              { rset_nan(E0[k][n]); } 
  }  
  // loop for QE-table
  for(k=1; k<m; k++){
    // Q[k][n]=(E[k-1][n+1]*Q[k-1][n+M])/E[k-1][n];
    for(n=0; n<Q0_size[k]; n++){ rdiv(a,E0[k-1][n+1],E0[k-1][n]); rmul(Q0[k][n],a,Q0[k-1][n+M]); }    
    // E[k][n]=Q[k][n+M]-Q[k][n]+E[k-1][n+1]
    for(n=0; n<E0_size[k]; n++){ rsub(a,Q0[k][n+M],Q0[k][n]); radd(E0[k][n],a,E0[k-1][n+1]); }
  }

  // debug
  if(debug>0){
    printf("Q=\n");
    for(n=0; n<n_size; n++){ for(k=0; k<m; k++){ if(n<Q0_size[k]){ mpfr_printf("%.3Re ",Q0[k][n]); } } printf("\n"); }
    printf("E=\n");
    for(n=0; n<n_size; n++){ for(k=0; k<m; k++){ if(n<E0_size[k]){ mpfr_printf("%.3Re ",E0[k][n]); } } printf("\n"); }
  }

  // generate vector E
  for(k=0; k<m; k++){ rcopy(E1[k],E0[k][0]); }

  // generate matrix Q
  for(n=0;n<M;n++){
    for(k=0;k<m;k++){
      rcopy(MAT(Q1,k,n,LDQ1),Q0[k][n]);
    }
  }


  // genrate matrix A
  if(A!=NULL){ riep_dhToda_QE_to_A(m,M,A,LDA,Q1,LDQ1,E1,debug); }

  // done
  if(Q==NULL){ Q1=rmat_free(LDQ1,M,Q1); }else{ Q1=NULL; }
  if(E==NULL){ E1=rvec_free(m,E1); }else{ E1=NULL; }
  a=rfree(a);
  f=rvec_free(f_size,f);
  sigma=rvec_free(m,sigma);
  for(k=0; k<m; k++){ Q0[k]=rvec_free(Q0_size[k],Q0[k]); } free(Q0); Q0=NULL;
  for(k=0; k<m; k++){ E0[k]=rvec_free(E0_size[k],E0[k]); } free(E0); E0=NULL;
  Q0_size=ivec_free(Q0_size);
  E0_size=ivec_free(E0_size);  
  return;
}


/**
 @brief QEより行列Aを生成
*/
void riep_EXTdhToda_QE_to_A(int m, int N, int M, rmulti **A, int LDA, rmulti **Q, int LDQ, rmulti **E, int LDE, int debug)
{
  // init
  int k=0,n=0,l=0,LDR=0,LDL=0,prec=0;
  rmulti **R=NULL,**L=NULL;
  // precision
  prec=rmat_get_prec_max(m,m,A,LDA);
  //allocate
  LDR=m; R=rmat_allocate_prec(LDR,m,prec);
  LDL=m; L=rmat_allocate_prec(LDL,m,prec);
  // A=I
  rmat_set_eye(m,m,A,LDA);
  // set A=A*L
  for(l=0; l<N; l++){
    for(n=0; n<m; n++){
      for(k=0; k<m; k++){
	if     (n==k)    { rset_one(MAT(L,n,k,LDL)); }              // L[n][k]=1
	else if((n-k)==1){ rcopy(MAT(L,n,k,LDA),MAT(E,k,l,LDE)); }  // L[n][k]=E[k][0]
	else             { rset_zero(MAT(L,n,k,LDL)); }             // L[n][k]=0
      }
    }
    rmat_prod(m,m,m,A,LDA,A,LDA,L,LDL); // A=A*L
  }   
  // set A=A*R
  for(l=M-1;l>=0;l--){
    for(n=0;n<m;n++){
      for(k=0;k<m;k++){
	if     ((k-n)==1){ rset_one(MAT(R,n,k,LDR)); }             // R[n][k]=1
	else if(n==k)    { rcopy(MAT(R,n,k,LDR),MAT(Q,k,l,LDQ)); } // R[n][k]=Q[k][l]
	else             { rset_zero(MAT(R,n,k,LDR)); }            // R[n][k]=0
      }
    }
    rmat_prod(m,m,m,A,LDA,A,LDA,R,LDR); // A=A*R
  }
  // debug
  if(debug>0){

  }
  //done
  R=rmat_free(LDR,m,R);
}

/**
 @brief 指定固有値よりdhTodaによりTN行列を生成．漸化式は方程式．
 @param[in]  m      行列のサイズ，
 @param[in]  M      M>=1の整数．帯幅は min(M+2,m)．
 @param[in]  A      サイズが(m,m)の行列．出力結果の格納用に初期化済みのA，またはA=NULLのときは出力はされない．
 @param[in]  LDA    Aの第1次元．
 @param[in]  Q      サイズが(m,M)の行列．出力結果の格納用に初期化済みのQ，またはQ=NULLのときは出力はされない．
 @param[in]  LDQ    Qの第1次元．
 @param[in]  E      サイズがmのベクトル．出力結果の格納用に初期化済みのE，またはE=NULLのときは出力はされない．
 @param[in]  lambda サイズが(m,m)の行列．
 @param[in]  c      サイズが(m,m)の行列．
 @param[in]  debug  デバグレベル．
 @param[out] A      計算結果．
 @param[out] Q      計算結果．
 @param[out] E      計算結果．
 @return            なし．
*/

void riep_EXTdhToda_TN(int m, int N, int M, rmulti **A, int LDA, rmulti **Q, int LDQ, rmulti **E, int LDE, rmulti **lambda, rmulti **c, int debug)
{
  int prec=0,k=0,n=0,f_size=0,n_size=0,*Q0_size=NULL,*E0_size=NULL,LDQ1,LDE1;
  rmulti *a=NULL,**f=NULL,***Q0=NULL,***E0=NULL,**sigma=NULL,**Q1=NULL,**E1=NULL;
  // init
  n_size=(M+N)*(m-1)+N*M+1;
  // precision
  prec=rmat_get_prec_max(m,m,A,LDA);
  // allocate
  if(Q==NULL){ LDQ1=m; Q1=rmat_allocate_prec(m,M,prec); }else{ LDQ1=LDQ; Q1=Q; }
  if(E==NULL){ LDE1=m; E1=rmat_allocate_prec(m,N,prec); }else{ LDE1=LDE; E1=E; }
  sigma=rvec_allocate_prec(m,prec);
  a=rallocate_prec(prec);
  Q0_size=ivec_allocate(m);
  Q0=malloc(m*sizeof(rmulti**));
  for(k=0; k<m; k++){
    Q0_size[k]=n_size-k*(M+N)-N;
    Q0[k]=rvec_allocate_prec(Q0_size[k],prec);
  }
  E0_size=ivec_allocate(m);
  E0=malloc(m*sizeof(rmulti**));
  for(k=0; k<m; k++){
    E0_size[k]=n_size-(k+1)*(M+N);
    E0[k]=rvec_allocate_prec(E0_size[k],prec);
  }
  f_size=n_size;
  f=rvec_allocate_prec(f_size,prec);
  // generate sigma[n]
  rinv_d(a,M*N);
  rvec_pow_r(m,sigma,lambda,a);
  // generate f[n]
  for(n=0; n<f_size; n++){
    rset_d(f[n],0);
    for(k=0; k<m; k++){
      rpow_si(a,sigma[k],n); // a=(sigma[i])^n
      radd_mul(f[n],c[k],a); // f[i]=f[i]+c[i]*(sigma[i])^n
    }
  }
  // Q[n][0]=f[n+1]/f[n]
  for(n=0; n<Q0_size[0]; n++){
    if(n+N<f_size){ rdiv(Q0[0][n],f[n+N],f[n]); }
    else          { rset_nan(Q0[0][n]); }
  }
  // E[0][n]=Q[0][n+M]-Q[0][n];
  k=0;
  for(n=0; n<E0_size[k]; n++){
    if(n+M<Q0_size[k] && n<Q0_size[k]){ rsub(E0[k][n],Q0[k][n+M],Q0[k][n]); }
    else                              { rset_nan(E0[k][n]); } 
  }
  // loop for QE-table
  for(k=1; k<m; k++){
    // Q[k][n]=(E[k-1][n+N]*Q[k-1][n+M])/E[k-1][n];
    for(n=0; n<Q0_size[k]; n++){ rdiv(a,E0[k-1][n+N],E0[k-1][n]); rmul(Q0[k][n],a,Q0[k-1][n+M]); }    
    // E[k][n]=Q[k][n+M]-Q[k][n]+E[k-1][n+N]
    for(n=0; n<E0_size[k]; n++){ rsub(a,Q0[k][n+M],Q0[k][n]); radd(E0[k][n],a,E0[k-1][n+N]); }
  }
  // debug
  if(debug>0){
    printf("Q=\n");
    for(n=0; n<n_size; n++){ for(k=0; k<m; k++){ if(n<Q0_size[k]){ mpfr_printf("%.3Re ",Q0[k][n]); } } printf("\n"); }
    printf("E=\n");
    for(n=0; n<n_size; n++){ for(k=0; k<m; k++){ if(n<E0_size[k]){ mpfr_printf("%.3Re ",E0[k][n]); } } printf("\n"); }
  }
  // generate vector E
  for(n=0; n<N; n++){
    for(k=0; k<m-1; k++){
      rcopy(MAT(E1,k,n,LDE1),E0[k][M*n]);
    }
  }
  // generate matrix Q
  for(n=0; n<M; n++){
    for(k=0; k<m; k++){
      rcopy(MAT(Q1,k,n,LDQ1),Q0[k][N*n]);
    }
  }
  // genrate matrix A
  if(A!=NULL){ riep_EXTdhToda_QE_to_A(m,N,M,A,LDA,Q1,LDQ1,E1,LDE1,debug); }
  // done
  if(Q==NULL){ Q1=rmat_free(LDQ1,M,Q1); }else{ Q1=NULL; }
  if(E==NULL){ E1=rmat_free(LDE1,N,E1); }else{ E1=NULL; }
  a=rfree(a);
  f=rvec_free(f_size,f);
  sigma=rvec_free(m,sigma);
  for(k=0; k<m; k++){ Q0[k]=rvec_free(Q0_size[k],Q0[k]); } free(Q0); Q0=NULL;
  for(k=0; k<m; k++){ E0[k]=rvec_free(E0_size[k],E0[k]); } free(E0); E0=NULL;
  Q0_size=ivec_free(Q0_size);
  E0_size=ivec_free(E0_size);  
  return;
}

