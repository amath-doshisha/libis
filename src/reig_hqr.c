#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<isys.h>

/**
 @file  reig_hqr.c
 @brief 多倍長精度実数型rmultiによるHessenberg型行列に前処理をした後にQR法で固有値を計算.
 @details rmulti型の関数に関する定義は@link rmulti.c@endlinkを参照のこと.
          rmulti型のベクトルの関数に関する定義は@link rvec.c@endlinkを参照のこと.
          rmulti型の行列の関数に関する定義は@link rmat.c@endlinkを参照のこと.
          rmulti型の固有値問題に関する定義は@link reig.c@endlinkを参照のこと.
 */

#define RA(X,P) ((X)=rallocate_prec(P));
#define RF(X) ((X)=rmfree(X))
#define NAME_HQR "reig_hqr"

/** @name 固有値計算に関する関数 */
/** @{ */

/**
 @brief   Hessenberg型行列に前処理をした後にQR法で固有値を計算.
 @param[in]  A      固有値計算をする(n,n)型行列.
 @param[in]  n      行列サイズ.
 @param[in]  LDA    Aの第1次元.
 @param[in]  debug  デバッグ.
 @param[out] lambda 計算された固有値を返す.
 */
void reig_hqr(int n, cmulti **lambda, rmulti **A, int LDA, int debug)
{
  int prec;
  rmulti **B=NULL;
  // allocate
  prec=cvec_get_prec_max(n,lambda);
  B=rmat_allocate_prec(n,n,prec);
  if(debug>=3){ printf("[%s] input matrix:\n",NAME_HQR); rmat_print(n,n,A,LDA,"A=",'f',2); }
  // convert to Hessenberg-type matrix
  if(debug>=1){ printf("[%s] converting to Hessenberg-type:\n",NAME_HQR); }
  if(debug>=3){ rmat_print(n,n,B,n,"A=",'f',2); }
  rhsnbrg_simtr(n,B,n,A,LDA);
  // compute eigenvalues
  if(debug>=1){ printf("[%s] computing eigenvalues:\n",NAME_HQR); }
  reig_hqr_main(n,lambda,B,n,debug);
  // sort
  cvec_sort(n,lambda,NULL);
  cvec_reverse(n,lambda);
  // done
  if(debug>=3){ printf("[%s] done:\n",NAME_HQR); rmat_print(n,n,B,n,"A=",'f',2); }
  B=rmat_free(n,n,B);
}

/** @} */

/** @name 固有値計算に関する関数 */
/** @{ */

/**
 @brief   Hessenberg型行列に前処理をした後にQR法で固有値を計算.
 @param[in]  A      固有値計算をする(n,n)型行列.
 @param[in]  n      行列サイズ.
 @param[in]  LDA    Aの第1次元.
 @param[in]  debug  デバッグ.
 @param[out] B      QR法の相似変換後の行列.
 @param[out] lambda 計算された固有値を返す.
 */
void reig_hqr_mt(int m, int n, rmulti **B, int LDB, cmulti **lambda, rmulti **A, int LDA, int debug)
{
  int prec=53;
  // allocate
  prec=cvec_get_prec_max(n,lambda);
  if(debug>=3){ printf("[%s] input matrix:\n",NAME_HQR); rmat_print(m,n,A,LDA,"A=",'f',2); }
  // convert to Hessenberg-type matrix
  if(debug>=1){ printf("[%s] converting to Hessenberg-type:\n",NAME_HQR); }
  if(debug>=3){ rmat_print(m,n,B,LDB,"A=",'f',2); }
  rhsnbrg_simtr(n,B,LDB,A,LDA);
  // compute eigenvalues
  if(debug>=1){ printf("[%s] computing eigenvalues:\n",NAME_HQR); }
  reig_hqr_main(n,lambda,B,LDB,debug);
  // sort
  cvec_sort(n,lambda,NULL);
  cvec_reverse(n,lambda);
  // done
  if(debug>=3){ printf("[%s] done:\n",NAME_HQR); rmat_print(m,n,B,LDB,"A=",'f',2); }
}

/** @} */

/////////////////////////////////////////////////////////////////////////////////////////////////////

/** @name 固有値計算に関する関数（ユーザが直接呼び出してはいけない関数） */
/** @{ */

/**
 @brief    多倍長QR法
 @details A はhessenberg前提
 @param[in]      n   行列サイズ
 @param[in,out]  A   [in]行列．[out]QR法後の行列
 */
void reig_hqr_main(int n, cmulti **lambda, rmulti **A, int LDA, int debug)
{
  int nn; // 残りの大きさ
  int prec,i,j,k,l,m,its,mmin;
  rmulti *anorm=NULL,*p=NULL,*q=NULL,*r=NULL,*s=NULL,*t=NULL,*u=NULL,*v=NULL,*w=NULL,*x=NULL,*y=NULL,*z=NULL;  
  //allocate
  prec=cvec_get_prec_max(n,lambda);
  anorm=rallocate_prec(prec);
  RA(p,prec); RA(q,prec); RA(r,prec); RA(s,prec); RA(t,prec);
  RA(u,prec); RA(v,prec); RA(w,prec); RA(x,prec); RA(y,prec); RA(z,prec);
  // anorm=......
  reig_hqr_calculate_norm(anorm,n,A,LDA);  
  // t=0
  rset_zero(t);
  // compute
  nn=n;
  while(nn>=1){
    // 反復回数をリセット
    its=0;
    do{
      if(debug>=2){ printf("[%s] size=%d itr=%d\n",NAME_HQR,nn,its); }
      // 収束しているか
      for(l=nn; l>=2; l--){
	// s=abs(A(l-2,l-2))+abs(A(l-1,l-1))
	radd_abs_abs_rr(s,MAT(A,l-2,l-2,LDA),MAT(A,l-1,l-1,LDA));
	// s=anorm if s==0
	if(ris_zero(s)){ rset_r(s,anorm); }
	// break if s+abs(A(l-1,l-2))==s
	if(reig_hqr_test_UDE(s,MAT(A,l-1,l-2,LDA))){ break; }
      }
      // x=A(n-1,n-1)
      rset_r(x,MAT(A,nn-1,nn-1,LDA));
      // ブロック1
      if(l==nn){
	// 実固有値に収束している場合
	// lambda(n-1)=x+t
	radd_rr(C_R(lambda[nn-1]),x,t); rset_zero(C_I(lambda[nn-1]));
	if(debug>=1){ mpfr_printf("[%s] lambda[%04d]=%+.15Re %+.15Re\n",NAME_HQR,nn-1,C_R(lambda[nn-1]),C_I(lambda[nn-1])); }
	// サイズを1だけ減次
	nn=nn-1;
      }else{
	// ブロック2以上
	// y=A(n-2,n-2)
	rset_r(y,MAT(A,nn-2,nn-2,LDA));
	// w=A(n-1,n-2)*A(n-2,n-1)
	rmul_rr(w,MAT(A,nn-1,nn-2,LDA),MAT(A,nn-2,nn-1,LDA));
	//
	if(l==(nn-1)){
	  // p=(y-x)/2
	  rsub_rr(p,y,x); rmul_rd(p,p,0.5);
	  // q=(y-x)^2/4+w
	  rmul_rr(q,p,p); radd_rr(q,q,w);
	  // z=sqrt(abs((y-x)^2/4+w))
	  rabs_r(z,q); rsqrt_r(z,z);
	  // x=x+t
	  radd_rr(x,x,t);
	  // if q>=0 then
	  if(ris_positive(q) || ris_zero(q)){
	    // 2つの実固有値に収束している場合
	    // z=abs(z)*sgn(p)
	    rabs_r(z,z); rmul_rd(z,z,rget_sgn_not0(p));
	    // z+=p
	    radd_rr(z,z,p);
	    // if z==0, then
	    if(ris_zero(z)){
	      // lambda(n-1)=x+z
	      radd_rr(C_R(lambda[nn-1]),x,z); rset_zero(C_I(lambda[nn-1]));
	    }else{
	      // lambda(n-1)=x-w/z
	      rdiv_rr(C_R(lambda[nn-1]),w,z); rsub_rr(C_R(lambda[nn-1]),x,C_R(lambda[nn-1])); rset_zero(C_I(lambda[nn-1]));
	    }
	    // lambda(n-2)=x+z
	    radd_rr(C_R(lambda[nn-2]),x,z); rset_zero(C_I(lambda[nn-2]));
	  }else{
	    // 複素固有値に収束している場合
	    // lambda(n-1)=(x+p)+I*z
	    radd_rr(C_R(lambda[nn-1]),x,p); rset_r(C_I(lambda[nn-1]),z);
	    // lambda(n-2)=(x+p)-I*z
	    radd_rr(C_R(lambda[nn-2]),x,p); rneg_r(C_I(lambda[nn-2]),z);
	  }
	  if(debug>=1){
	    mpfr_printf("[%s] lambda[%04d]=%+.15Re %+.15Re\n",NAME_HQR,nn-1,C_R(lambda[nn-1]),C_I(lambda[nn-1]));
	    mpfr_printf("[%s] lambda[%04d]=%+.15Re %+.15Re\n",NAME_HQR,nn-2,C_R(lambda[nn-2]),C_I(lambda[nn-2]));
	  }
	  // サイズを2だけ減次
	  nn=nn-2;
	}else{
	  // 収束していない場合
	  if(its==100) { ERROR_AT; exit(0); }
	  if((its%10)==0){
	    // t+=x
	    radd_rr(t,t,x);
	    // 原点シフト
	    for(i=0; i<nn; i++){ rsub_rr(MAT(A,i,i,LDA),MAT(A,i,i,LDA),x); }
	    // s=abs(A(n-1,n-2))+abs(A(n-2,n-3))
	    radd_abs_abs_rr(s,MAT(A,nn-1,nn-2,LDA),MAT(A,nn-2,nn-3,LDA));
	    // y=s*0.75
	    rmul_rd(y,s,0.750);
	    // x=s*0.75
	    rmul_rd(x,s,0.750);
	    // w=s^2
	    rmul_rr(w,s,s);
	    // w=w*(-0.4375)
	    rmul_rd(w,w,-0.4375);
	  }
	  // 反復回数を増やす
	  ++its;
	  for(m=(nn-2); m>=l; m--){
	    rset_r(z,MAT(A,m-1,m-1,LDA));
	    rsub_rr(r,x,z);
	    rsub_rr(s,y,z);	    
	    reig_hqr_calculate_formula_11_6_23_p(p,r,s,w,MAT(A,m,m-1,LDA),MAT(A,m-1,m,LDA));
	    reig_hqr_calculate_formula_11_6_23_q(q,MAT(A,m,m,LDA),r,s,z);
	    rset_r(r,MAT(A,m+1,m,LDA));
	    radd_abs_abs_abs(s,p,q,r);
	    reig_hqr_scale(p,q,r,s);
	    if(m==l){ break; }
	    reig_hqr_calculate_formula_11_6_26_u(u,MAT(A,m-1,m-2,LDA),q,r);
	    reig_hqr_calculate_formula_11_6_26_v(v,MAT(A,m-2,m-2,LDA),MAT(A,m,m,LDA),p,z);
	    if(reig_hqr_test_formula_11_6_26(u,v)){ break; }
	  }
	  for(i=m+2; i<=nn; i++){
	    rset_d(MAT(A,i-1,i-3,LDA),0.0);
	    if(i!=(m+2)){ rset_d(MAT(A,i-1,i-4,LDA),0.0); }
	  }
	  for(k=m; k<=nn-1; k++){
	    if(k!=m){
	      rset_r(p,MAT(A,k-1,k-2,LDA));
	      rset_r(q,MAT(A,k,k-2,LDA));
	      rset_zero(r);
	      if(k!=(nn-1)){ rset_r(r,MAT(A,k+1,k-2,LDA)); }
	      radd_abs_abs_abs(x,p,q,r);
	      if(ris_zero(x)!=1){ reig_hqr_scale(p,q,r,x); }
	    }
	    reig_hqr_calculate_formula_11_6_25(s,p,q,r);
	    if(ris_zero(s)!=1){
	      if(k==m){
		if(l!=m) { rneg_r(MAT(A,k-1,k-2,LDA),MAT(A,k-1,k-2,LDA)); }
	      }else{ 
		rmul_rr(MAT(A,k-1,k-2,LDA),s,x); 
		rneg_r(MAT(A,k-1,k-2,LDA),MAT(A,k-1,k-2,LDA)); 
	      }
	      radd_rr(p,p,s);
	      rdiv_rr(x,p,s);
	      rdiv_rr(y,q,s);
	      rdiv_rr(z,r,s);
	      rdiv_rr(q,q,p);
	      rdiv_rr(r,r,p);
	      for(j=k;j<=nn;j++){
		rmul_rr(p,q,MAT(A,k,j-1,LDA));
		radd_rr(p,p,MAT(A,k-1,j-1,LDA));
		if(k!=(nn-1)){
		  radd_mul_rr(p,r,MAT(A,k+1,j-1,LDA));
		  rsub_mul_rr(MAT(A,k+1,j-1,LDA),p,z);
		}
		rsub_mul_rr(MAT(A,k,j-1,LDA),p,y);
		rsub_mul_rr(MAT(A,k-1,j-1,LDA),p,x);
	      }
	      if(nn<(k+3)) { mmin=nn;  }
	      else         { mmin=k+3; }
	      for(i=l;i<=mmin;i++){
		rmul_rr(p,x,MAT(A,i-1,k-1,LDA));
		radd_mul_rr(p,y,MAT(A,i-1,k,LDA));
		if(k!=(nn-1)){
		  radd_mul_rr(p,z,MAT(A,i-1,k+1,LDA));
		  rsub_mul_rr(MAT(A,i-1,k+1,LDA),p,r);
		}
		rsub_mul_rr(MAT(A,i-1,k,LDA),p,q);
		rsub_rr(MAT(A,i-1,k-1,LDA),MAT(A,i-1,k-1,LDA),p);
	      }
	    }
	  }
	}
      }
    }while(l<nn-1);
  }
  if(debug>=1){ printf("[%s] done.\n",NAME_HQR); }
  // free
  RF(anorm); RF(p); RF(q); RF(r); RF(s); RF(t); RF(u); RF(v); RF(w); RF(x); RF(y); RF(z);  
}

/////////////////////////////////////////////////////////////////////////////////////////////////

/**
 @brief 下副対角要素を探すためのノルム計算(ユーザーは直接使用してはならない)
*/
void reig_hqr_calculate_norm(rmulti *anorm, int n, rmulti **A, int LDA)
{
  int i,j;
  rabs_r(anorm,MAT(A,0,0,LDA));
  for(i=1; i<n; i++){
    for(j=(i-1); j<n; j++){
      radd_abs_r(anorm,MAT(A,i,j,LDA));
    }
  }
}

/**
 @brief p=(r*s-w)/A[m+1,m] + A[m,m+1] (ユーザーは直接使用してはならない)
*/
void reig_hqr_calculate_formula_11_6_23_p(rmulti *p, rmulti *r, rmulti *s, rmulti *w, rmulti *a1, rmulti *a2)
{
  rmul_rr(p,r,s);  // p=r*s
  rsub_rr(p,p,w);  // p=r*s-w
  rdiv_rr(p,p,a1); // p=(r*s-w)/a1
  radd_rr(p,p,a2); // p=(r*s-w)/a1+a2
}

/**
 @brief q=A[m+1,m+1]-z-r-s (ユーザーは直接使用してはならない)
*/
void reig_hqr_calculate_formula_11_6_23_q(rmulti *q, rmulti *a, rmulti *r, rmulti *s, rmulti *z)
{
  rsub_rr(q,a,z); // q=a-z
  rsub_rr(q,q,r); // q=a-(z+r)
  rsub_rr(q,q,s); // q=a-(z+r+s)
}

/**
 @brief s=SIGN(sqrt(p*p+q*q+r*r),p) (ユーザーは直接使用してはならない)
*/
void reig_hqr_calculate_formula_11_6_25(rmulti *s, rmulti *p, rmulti *q, rmulti *r)
{
  int prec;
  rmulti *X=NULL,*Y=NULL;
  // init
  prec=rget_prec(s);
  RA(X,prec); RA(Y,prec);
  // compute
  rmul_rr(s,p,p); // s=p^2
  rmul_rr(X,q,q); // X=q^2
  rmul_rr(Y,r,r); // Y=r^2
  radd_rr(s,s,X); // s=p^2+q^2
  radd_rr(s,s,Y); // s=p^2+q^2+r-2
  rsqrt_r(s,s);  // s=sqrt(p^2+q^2+r-2)
  rmul_rd(s,s,rget_sgn_not0(p)); // s=sqrt(p^2+q^2+r-2)*sgn(p)
  // done
  RF(X); RF(Y);
}

/**
 @brief u=abs(A(m,m-1))*(abs(q)+abs(r)) (ユーザーは直接使用してはならない)
*/
void reig_hqr_calculate_formula_11_6_26_u(rmulti *u, rmulti *a, rmulti *q, rmulti *r)
{
  int prec;
  rmulti *X=NULL;  
  // init
  prec=rget_prec(u);
  RA(X,prec);
  // compute
  radd_abs_abs_rr(X,q,r); // X=abs(q)+abs(r)
  rabs_r(u,a);          // u=abs(a)
  rmul_rr(u,u,X);        // u=abs(a)*(abs(q)+abs(r))
  // done
  RF(X);
}

/**
 @brief v=abs(p)*(abs(A(m-1,m-1))+abs(z)+abs(A(m+1,m+1))) (ユーザーは直接使用してはならない)
*/
void reig_hqr_calculate_formula_11_6_26_v(rmulti *v, rmulti *a1, rmulti *a2, rmulti *p, rmulti *z)
{
  int prec;
  rmulti *X=NULL;
  // init
  prec=rget_prec(v);
  RA(X,prec);
  // compute
  radd_abs_abs_abs(X,a1,z,a2); // X=abs(a1)+abs(z)+abs(a2)
  rabs_r(v,p);              // v=abs(p)
  rmul_rr(v,v,X);            // v=abs(p)*(abs(a1)+abs(z)+abs(a2))
  // done
  RF(X);
}

/**
 @brief u+v==vの判定(ユーザーは直接使用してはならない)
*/
int reig_hqr_test_formula_11_6_26(rmulti *u, rmulti *v)
{
  int k,prec;
  rmulti *X=NULL;
  // init
  prec=rget_prec(v);
  RA(X,prec);
  // compute
  radd_rr(X,u,v);  // X=u+v
  k=eq_rr(X,v);   // u+v==v
  // done
  RF(X);
  return k;
}

/**
 @brief (abs(a)+s==s)の判定(ユーザーは直接使用してはならない)
*/
int reig_hqr_test_UDE(rmulti *s, rmulti *a)
{
  int prec,k=0;
  rmulti *X=NULL;
  // init
  prec=rget_prec(s);
  RA(X,prec);
  // compute
  rabs_r(X,a);        // X=abs(a)
  radd_rr(X,X,s);      // X=abs(a)+s
  k=eq_rr(X,s);       // (abs(a)+s==s)
  // done
  RF(X);
  return k;
}

/**
 @brief オーバーフローやアンダーフローを避けるためのスケール動作(ユーザーは直接使用してはならない)
*/
void reig_hqr_scale(rmulti *p, rmulti *q, rmulti *r, rmulti *s)
{
  rdiv_rr(p,p,s); // p/=s
  rdiv_rr(q,q,s); // q/=s
  rdiv_rr(r,r,s); // r/=s
}

/** @} */

//EOF
