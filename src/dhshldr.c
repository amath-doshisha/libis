#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"is_macros.h"
#include"is_ivec.h"
#include"is_dvec.h"
#include"is_dmat.h"
#include"is_dhshldr.h"

// size(x)=n
// h[0]=0; ...; h[k-1]=0; h[k]=-s*xi; h[k+1]=x[k+1]; ...; h[n-1]=x[n-1];
void dhouseholder_vec(int n, int k, double *h, double *alpha, double *x)
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
  double xi,eta,s=1;
  //----------- norm
  xi=dvec_dot(n-k-1,&x[k+1],&x[k+1]); //  xi=sum(abs(x((k+1):end)).^2);
  eta=sqrt(x[k]*x[k]+xi);
  if(eta==0) xi=eta-fabs(x[k]);
  else       xi=xi/(fabs(x[k])+eta);
  if(x[k]<0) s=-1;
  //----------- h
  dvec_set_zeros(k,h);
  dvec_copy(n-k-1,&h[k+1],&x[k+1]); // h((k+1):end)=x((k+1):end);
  h[k]=-s*xi;
  //----------- alpha
  if(xi==0 || eta==0) (*alpha)=0;
  else                (*alpha)=1.0/(xi*eta);
}

void dhouseholder(int n, int k0, int nH, int k, double *h, double *alpha, double *H, int LDH, double *Alpha, double *x)
{
  int j,l;
  double value,*R=NULL;
  R=dvec_allocate(n);
  dvec_copy(n,R,x);
  // R=(I-alpha*h*h')*R=R-alpha*h*(h'*R)
  for(j=0; j<nH; j++){
    // R=R-alpha[j]*(H(:,j)'*R)*H(:,j)
    l=k0+j;
    value=dvec_dot(n-l,&MAT(H,l,j,LDH),&R[l]);
    dvec_add_scaled(n-l,&R[l],-Alpha[j]*value,&MAT(H,l,j,LDH));
  }
  dhouseholder_vec(n,k,h,alpha,R);
  dvec_free(R);
}

// B=A*H, H=I-alpha*h*h'
void dhouseholder_right(int m, int n, double *A, int LDA, int k, double *h, double alpha)
{
  double *p=NULL;
  p=dvec_allocate(m);
  // p=A*h
  dvec_lintr(m,n-k,p,&COL(A,k,LDA),LDA,&h[k]);
  // B=A*H=A-alpha*(A*h)*h'=A-alpha*p*h'
  dmat_rank1op(m,n-k,&COL(A,k,LDA),LDA,-alpha,p,&h[k]);
  // done
  dvec_free(p);
}

// B=H*A, H=I-alpha*h*h'
void dhouseholder_left(int m, int n, double *A, int LDA, int k, double *h, double alpha)
{
  double *p=NULL;
  p=dvec_allocate(m);
  // p=A'*h
  dvec_lintr_t(n,m-k,p,&MAT(A,k,0,LDA),LDA,&h[k]);
  // B=H*A=A-alpha*h*(A'*h)'=A-alpha*h*p'
  dmat_rank1op(m-k,n,&MAT(A,k,0,LDA),LDA,-alpha,&h[k],p);
  // done
  dvec_free(p);
}

//EOF
