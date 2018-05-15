#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"is_macros.h"
#include"is_ivec.h"
#include"is_zvec.h"
#include"is_zmat.h"
#include"is_zhshldr.h"


// size(x)=n
// h[0]=0; ...; h[k-1]=0; h[k]=-s*xi; h[k+1]=x[k+1]; ...; h[n-1]=x[n-1];
void zhouseholder_vec(int n, int k, dcomplex *h, double *alpha, dcomplex *x)
{
  // function [h,alpha]=householder_vec(x,k)
  // y=H*x, H=I-2*h*h'/(h'*h), h=x-y
  // i<k   y(i)=x(i)   h(i)=0
  // i=k   y(k)=s*eta  h(k)=x(k)-y(k)=x(k)-s*eta=-s*xi
  // i>k   y(i)=0      h(i)=x(i)
  // eta=norm2(x(k:end))
  // s=sgn(x(k))
  // xi=-h(k)/s=-(x(k)-s*eta)/s=eta-x(k)/s=eta-|x(k)|
  //   =(eta^2-|x(k)|^2)/(eta+|x(k|)=(|x(k+1)|^2+..+|x(n)|^2)/(eta+|x(k)|)
  // h'*h=2*xi*eta
  // alpha=2/(h'*h)=1/(xi*eta)
  double eta,s,xi,axk;
  //----------- norm
  xi=dsum_pow2_abs_zvec(n-k-1,&x[k+1]); // xi=sum(abs(x((k+1):end)).^2);
  axk=Z_ABS2(x[k]);
  eta=sqrt(axk+xi);
  axk=sqrt(axk);
  if(eta==0) xi=eta-axk;        // eta=eta-|x(k)|
  else       xi=xi/(eta+axk);   // eta=xi/(eta+|x(k)|)
  //----------- h
  zvec_set_zeros(k,h);
  zvec_copy_zvec(n-k-1,&h[k+1],&x[k+1]);    // h((k+1):end)=x((k+1):end);
  if((Z_R(x[k])*Z_I(x[k]))==0.0){
    Z_SET(h[k],-xi,0);
  }else{
    s=-xi/axk;
    Z_R(h[k])=s*Z_R(x[k]);
    Z_I(h[k])=s*Z_I(x[k]);
  }
  //----------- alpha
  if(xi==0 || eta==0) (*alpha)=0;
  else                (*alpha)=1.0/(xi*eta);
}

void zhouseholder(int n, int k0, int nH, int k, dcomplex *h, double *alpha, dcomplex *H, int LDH, double *Alpha, dcomplex *x)
{
  int j,l;
  dcomplex value,*R=NULL; 
  R=zvec_allocate(n);
  zvec_copy_zvec(n,R,x);
  // R=(I-alpha*h*h')*R=R-alpha*h*(h'*R)    
  for(j=0; j<nH; j++){
    // R=R-alpha[j]*(H(:,j)'*R)*H(:,j)
    l=k0+j;
    value=zvec_dot(n-l,&MAT(H,l,j,LDH),&R[l]);
    Z_R(value)*=-Alpha[j];
    Z_I(value)*=-Alpha[j];      
    zvec_add_scaled(n-l,&R[l],value,&MAT(H,l,j,LDH));
  }
  zhouseholder_vec(n,k,h,alpha,R);
  R=zvec_free(R);
}

// B=A*H, H=I-alpha*h*h'
void zhouseholder_right(int m, int n, dcomplex *A, int LDA, int k, dcomplex *h, double alpha)
{
  dcomplex zalpha,*p=NULL;
  p=zvec_allocate(m);
  // p=A*h
  zvec_lintr(m,n-k,p,&COL(A,k,LDA),LDA,&h[k]);
  // B=A*H=A-alpha*(A*h)*h'=A-alpha*p*h'
  Z_SET(zalpha,-alpha,0);
  zmat_rank1op(m,n-k,&COL(A,k,LDA),LDA,zalpha,p,&h[k]);
  // done
  p=zvec_free(p);
}

// B=H*A, H=I-alpha*h*h'
void zhouseholder_left(int m, int n, dcomplex *A, int LDA, int k, dcomplex *h, double alpha)
{
  dcomplex zalpha,*p=NULL;
  p=zvec_allocate(m);
  // p=A'*h
  zvec_lintr_ct(n,m-k,p,&MAT(A,k,0,LDA),LDA,&h[k]);
  // B=H*A=A-alpha*h*(A'*h)'=A-alpha*h*p'
  Z_SET(zalpha,-alpha,0);
  zmat_rank1op(m-k,n,&MAT(A,k,0,LDA),LDA,zalpha,&h[k],p);
  // done
  p=zvec_free(p);
}

//EOF
