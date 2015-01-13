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

void riep_dhToda_matgen(int m, int M, rmulti **A, int LDA, rmulti ***Q, rmulti ***E, int debug)
{
  // init
  int k=0,n=0,l=0,LDR=0;

  //allocate
  LDR=m; R=rmat_allocate_prec(LDR,m,prec);

  // set A=L
  for(n=0; n<m; n++){
    for(k=0; k<m; k++){
      if     (n==k)    { rset_one(MAT(A,n,k,LDA)); }      // A[n][k]=1;
      else if((n-k)==1){ rcopy(MAT(A,n,k,LDA),V[k][0]); } // A[n][k]=E[k][0];
      else             { rset_zero(MAT(A,n,k,LDA)); }     // A[n][k]=0;
    }
  }
  // set A=A*R
  for(l=M-1;l>=0;l--){
    for(n=0;n<m;n++){
      for(k=0;k<m;k++){
	if     ((k-n)==1){ rset_one(MAT(R,n,k,LDR)); }      // R[n][k]=1;
	else if(n==k)    { rcopy(MAT(R,n,k,LDR),U[k][l]); } // R[n][k]=Q[k][l];
	else             { rset_zero(MAT(R,n,k,LDR)); }     // R[n][k]=0;
      }
    }
    rmat_prod(m,m,m,A,LDA,A,LDA,R,LDR); // A=A*R
  }
  //done
  R=rmat_free(LDR,m,R);

}

void rmat_gen_dhToda(int m, int M, rmulti **A, int LDA, rmulti **Q, int LDQ, rmulti **E, rmulti **lambda, rmulti **c, int debug)
{
  int prec=0,k=0,n=0,l=0,f_size=0,n_size=0,*U_size=NULL,*V_size=NULL,LDR=0;
  rmulti *a=NULL,**f=NULL,***U=NULL,***V=NULL,**sigma=NULL,**R=NULL;

  // init
  n_size=(M+1)*(m-1)+2*M;
  // precision
  prec=rmat_get_prec_max(m,m,A,LDA);
  
  // allocate
  sigma=rvec_allocate_prec(m,prec);
  a=rallocate_prec(prec);
  U_size=ivec_allocate(m);
  U=malloc(m*sizeof(rmulti**));
  for(k=0; k<m; k++){
    U_size[k]=n_size-k*(M+1);
    U[k]=rvec_allocate_prec(U_size[k],prec);
  }
  V_size=ivec_allocate(m);
  V=malloc(m*sizeof(rmulti**));
  for(k=0; k<m; k++){
    V_size[k]=n_size-(k+1)*(M+1)+1;
    V[k]=rvec_allocate_prec(V_size[k],prec);
  }
  f_size=U_size[0]+1; f=rvec_allocate_prec(f_size,prec);
  // set sigma
  rinv_d(a,M);
  rvec_pow_r(m,sigma,lambda,a);

  // set f
  for(n=0; n<f_size; n++){
    rset_d(f[n],0);
    for(k=0; k<m; k++){
      rpow_si(a,sigma[k],n); // a=(sigma[i])^n
      radd_mul(f[n],c[k],a); // f[i]=f[i]+c[i]*(sigma[i])^n
    }
  }
  // set Q[0]
  for(n=0; n<U_size[0]; n++){
    if(n+1<f_size){
      rdiv(U[0][n],f[n+1],f[n]); // Q[n][1]=f[n+1]/f[n]
    }else{ rset_nan(U[0][n]); }
  }
  // set E[0]
  k=0;
  for(n=0; n<V_size[k]; n++){
    if(n+M<U_size[k] && n<U_size[k]){
      rsub(V[k][n],U[k][n+M],U[k][n]);       //    E[k][n]=Q[k][n+M]-Q[k][n];
    }else{ rset_nan(E[k][n]); }
  }  
  // set QE-table
  for(k=1; k<m; k++){
    for(n=0; n<U_size[k]; n++){
      rdiv(a,V[k-1][n+1],V[k-1][n]);
      rmul(U[k][n],a,U[k-1][n+M]);     // Q[k][n]=(E[k-1][n+1]*Q[k-1][n+M])/E[k-1][n];
    }    
    for(n=0; n<V_size[k]; n++){
      //radd(a,Q[k][n+M],E[k-1][n+1]);
      //rsub(E[k][n],a,Q[k][n]);
      rsub(a,U[k][n+M],U[k][n]);
      radd(V[k][n],a,V[k-1][n+1]);     // E[k][n]=Q[k][n+M]-Q[k][n]+E[k-1][n+1];
    }
  }

  //genrate matrix
  riep_dhToda_matgen(m,M,A,LDA,V,U,debug);

  //set E_vec

  //set Q_mat
  
  // debug
  if(debug>0){
    printf("Q=\n");
    for(n=0; n<n_size; n++){ for(k=0; k<m; k++){ if(n<U_size[k]){ mpfr_printf("%.3Re ",U[k][n]); } } printf("\n"); }
    printf("E=\n");
    for(n=0; n<n_size; n++){ for(k=0; k<m; k++){ if(n<V_size[k]){ mpfr_printf("%.3Re ",V[k][n]); } } printf("\n"); }
  }
  // done
  a=rfree(a);
  f=rvec_free(f_size,f);
  sigma=rvec_free(m,sigma);
  for(k=0; k<m; k++){ U[k]=rvec_free(U_size[k],U[k]); } free(U); U=NULL;
  for(k=0; k<m; k++){ V[k]=rvec_free(V_size[k],V[k]); } free(V); V=NULL;
  U_size=ivec_free(U_size);
  V_size=ivec_free(V_size);
  
  return;
}
