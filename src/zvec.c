#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<stdarg.h>

#include"mt19937ar.h"
#include"is_macros.h"
#include"is_strings.h"
#include"is_dcomplex.h"
#include"is_rmulti.h"
#include"is_cmulti.h"
#include"is_irmulti.h"
#include"is_icmulti.h"
#include"is_svec.h"
#include"is_ivec.h"
#include"is_zvec.h"
#include"is_dvec.h"

/////////////////////////////////

#define FILE_NAME_LENGTH_MAX 100000

/////////////////////////////////


dcomplex* zvec_allocate(int n)
{
  return (dcomplex*)malloc(sizeof(dcomplex)*n);
}

dcomplex* zvec_free(dcomplex *x)
{
  if(x==NULL) return NULL;
  free(x);
  return x=NULL;
}


// x=nan(n,1)
void zvec_set_nan(int n, dcomplex *x)
{
  int i;
  for(i=0; i<n; i++){ Z_SET(x[i],NAN,NAN); }
}

// x=inf(n,1)
void zvec_set_inf(int n, dcomplex *x, int rsgn, int isgn)
{
  double zero=0;
  int i;
  for(i=0; i<n; i++){ Z_SET(x[i],rsgn/zero,isgn/zero); }
}


// x=zeros(n,1)
void zvec_set_zeros(int n, dcomplex *x)
{
  int i;
  for(i=0; i<n; i++){ Z_SET(x[i],0,0); }
}

// x=ones(n,1)
void zvec_set_ones(int n, dcomplex *x)
{
  int i;
  for(i=0; i<n; i++){ Z_SET(x[i],1,0); }
}

// x=ones(n,1)*a
void zvec_set_all(int n, dcomplex *x, dcomplex a)
{
  int i;
  for(i=0; i<n; i++){ x[i]=a; }
}

// x=ones(n,1)*a
void zvec_set_all_d(int n, dcomplex *x, double a)
{
  int i;
  for(i=0; i<n; i++){ Z_SET(x[i],a,0); }
}

// x=ones(n,1)*a
void zvec_set_all_dd(int n, dcomplex *x, double real, double imag)
{
  int i;
  for(i=0; i<n; i++){ Z_SET(x[i],real,imag); }
}

// x=zeros(n,1); x[k]=1
void zvec_set_unit(int n, dcomplex *x, int k)
{
  zvec_set_zeros(n,x);
  Z_SET(x[k],1,0);
}

// x[0]=0; x[1]=1; x[2]=2; ...; x[n-1]=n-1
void zvec_set_grid(int n, dcomplex *x)
{
  int i;
  for(i=0; i<n; i++){ Z_SET(x[i],i,0); }
}

// x=rand(n,1)*a+b
void zvec_set_rand(int n, dcomplex *x, double a, double b){
  double xr,xi;
  int i;
  for(i=0; i<n; i++){
    xr=genrand_int32();
    xi=genrand_int32();
    xr/=0xffffffff;
    xi/=0xffffffff;
    xr=xr*a+b;
    xi=xi*a+b;
    Z_SET(x[i],xr,xi);
  }
}

// y=x
void zvec_set_s(int n, dcomplex *y, char **x)
{
  strings *list=NULL;
  int i;
  for(i=0; i<n; i++){
    list=strings_split_number(x[i]);
    if(list==NULL){ Z_SET(y[i],NAN,NAN); }
    else{
      Z_SET(y[i],0,0);
      if(strings_size(list)>=1 && strings_at(list,0)!=NULL){ Z_R(y[i])=atof(strings_at(list,0)); }
      if(strings_size(list)>=2 && strings_at(list,1)!=NULL){ Z_I(y[i])=atof(strings_at(list,1)); }
    }    
    list=strings_del(list);
  }
}

// x[0]=x[n-1]; x[1]=x[n-2]; x[2]=x[n-3]; ...
void zvec_reverse(int n, dcomplex *x)
{
  int i;
  dcomplex a;
  for(i=0; i<n/2; i++){
    a=x[i];
    x[i]=x[n-i-1];
    x[n-i-1]=a;
  }
}

/////////////////////////////////////////////////////////////////////////////////////

// y=real(x)
void zvec_real(int n, double *y, dcomplex *x)
{
  int i;
  for(i=0; i<n; i++){ y[i]=Z_R(x[i]); }
}

/////////////////////////////////////////////////////////////////////////////////////

// x[i] <=> x[j]
void zvec_swap_at(dcomplex *x, int i, int j)
{
  dcomplex a;
  a=x[i];
  x[i]=x[j];
  x[j]=a;
}

// x[i] <=> x[I[i]]
void zvec_swap_index(int n, dcomplex *x, const int *I)
{
  dcomplex *y=NULL;
  y=zvec_allocate(n);
  zvec_copy_index(n,y,x,I);
  zvec_copy(n,x,y);
  y=zvec_free(y);
}

/////////////////////////////////////////////////////////////////////////////////////

// sort x as x[0] <= x[1] <= x[2] <= ... <= x[n-1]
// if I==NULL, then I is not ussed.
// if I!=NULL, then I is stored with sorted indexes
void zvec_sort(int n, dcomplex *x, int *I)
{
  if(I!=NULL){ ivec_set_grid(n,I); }
  zvec_quick_sort(n,x,I,0,n-1);
}

// Don't call this function directly!
void zvec_quick_sort(int n, dcomplex *x, int *I, int left, int right)
{
  int i,last;
  if(left>=right) return;
  zvec_swap_at(x,left,(left+right)/2);
  if(I!=NULL) ivec_swap_at(I,left,(left+right)/2);
  last=left;
  for(i=left+1; i<=right; i++){
    if(zlt(x[i],x[left])){
      ++last;
      zvec_swap_at(x,last,i);
      if(I!=NULL) ivec_swap_at(I,last,i);
    }
  }
  zvec_swap_at(x,left,last);
  if(I!=NULL) ivec_swap_at(I,left,last);
  zvec_quick_sort(n,x,I,left,last-1);
  zvec_quick_sort(n,x,I,last+1,right);
}

//////////////////////////////////////////////////////////////////////////////////////////

// store the list of indexes, x is not destroyed
void zvec_sort_index(int *I, int n, dcomplex *X)
{
  dcomplex *Y=NULL;
  Y=zvec_allocate(n);
  zvec_copy(n,Y,X);
  zvec_sort(n,Y,I);
  Y=zvec_free(Y);
}


/////////////////////////////////////////////////////////////////////////////////////

// y=x
void zvec_set_z(int n, dcomplex *y, const dcomplex *x)
{
  int i;
  for(i=0; i<n; i++){ y[i]=x[i]; }
}

// y=x
void zvec_set_d(int n, dcomplex *y, const double *x)
{
  int i;
  for(i=0; i<n; i++){ Z_SET(y[i],x[i],0); }
}

// y=x
void zvec_set_r(int n, dcomplex *y, rmulti **x)
{
  double a;
  int i;  
  for(i=0; i<n; i++){
    a=rget_d(x[i]);
    Z_SET(y[i],a,0);
  }
}

// y=x
void zvec_set_c(int n, dcomplex *y, cmulti **x)
{
  int i;
  for(i=0; i<n; i++){ y[i]=cget_z(x[i]); }
}


// y=x
void zvec_set_ir(int n, dcomplex *y, rmulti **x0, rmulti **x1)
{
  double a;
  int i;  
  for(i=0; i<n; i++){
    a=0.5*(rget_d(x0[i])+rget_d(x1[i]));
    Z_SET(y[i],a,0);
  }
}

// y=x
void zvec_set_ic(int n, dcomplex *y, cmulti **x0, cmulti **x1)
{
  dcomplex z0,z1;
  int i;
  for(i=0; i<n; i++){
    z0=cget_z(x0[i]);
    z1=cget_z(x1[i]);
    Z_ADD(z0,z1);
    Z_SCALE(z0,0.5);
    y[i]=z0;
  }
}

/////////////////////////////////////////////////////////////////////////////////////


// y=x
void zvec_copy(int n, dcomplex *y, const dcomplex *x)
{
  int i;
  for(i=0; i<n; i++) y[i]=x[i];
}

// y=x
void zvec_copy_d(int n, dcomplex *y, const double *x)
{
  int i;
  for(i=0; i<n; i++){
    Z_SET(y[i],x[i],0);
  }
}

// y=x
void zvec_copy_dd(int n, dcomplex *y, const double *x_r, const double *x_i)
{
  int i;
  for(i=0; i<n; i++){
    Z_SET(y[i],x_r[i],x_i[i]);
  }
}



// Y[i]=X[I[i]], 0<=i<n
void zvec_copy_index(int n, dcomplex *Y, const dcomplex *X, const int *I)
{
  int i;
  for(i=0; i<n; i++) Y[i]=X[I[i]];
}

/////////////////////////////////////////////////////////////////////////////////////

// y=conj(x)
void zvec_conj(int n, dcomplex *y, const dcomplex *x)
{
  int i;
  for(i=0; i<n; i++){ Z_SET_CONJ(y[i],x[i]); }
}

/////////////////////////////////////////////////////////////////////////////////////


// x=x+a
void zvec_add_scalar(int n, dcomplex *x, dcomplex a)
{
  int i;
  for(i=0; i<n; i++){
    Z_ADD(x[i],a);
  }
}

// x=x-a
void zvec_sub_scalar(int n, dcomplex *x, dcomplex a)
{
  int i;
  for(i=0; i<n; i++){
    Z_SUB(x[i],a);
  }
}

// z=x+y
void zvec_add(int n, dcomplex *z, dcomplex *x, dcomplex *y)
{
  int i;
  for(i=0; i<n; i++){
    Z_SET_PLUS(z[i],x[i],y[i]);
  }
}

// z=x+y
void zvec_add_d(int n, dcomplex *z, dcomplex *x, double *y)
{
  int i;
  for(i=0; i<n; i++){
    Z_R(z[i])=Z_R(x[i])+y[i];
    Z_I(z[i])=Z_I(x[i]);
  }
}

// z=x+y
void zvec_add_dz(int n, dcomplex *z, double *x, dcomplex *y)
{
  return zvec_add_d(n,z,y,x);
}

// z=x+y
void zvec_add_zd(int n, dcomplex *z, dcomplex *x, double *y)
{
  return zvec_add_d(n,z,x,y);
}

// z=x-y
void zvec_sub(int n, dcomplex *z, dcomplex *x, dcomplex *y)
{
  int i;
  for(i=0; i<n; i++){
    Z_SET_MINUS(z[i],x[i],y[i]);
  }
}

// y=y+a*x
void zvec_add_scaled(int n, dcomplex *y, dcomplex a, const dcomplex *x)
{
  int i;
  for(i=0; i<n; i++){
    Z_ADD_TIMES(y[i],a,x[i]);
  }
}

// y=y-a*x
void zvec_sub_scaled(int n, dcomplex *y, dcomplex a, const dcomplex *x)
{
  int i;
  for(i=0; i<n; i++){
    Z_SUB_TIMES(y[i],a,x[i]);
  }
}

// y=A*x
void zvec_lintr(int m, int n, dcomplex *y, const dcomplex *A, int LDA, const dcomplex *x)
{
  int i,j;
  for(i=0; i<m; i++){
    Z_SET(y[i],0,0);
    for(j=0; j<n; j++){
      Z_ADD_TIMES(y[i],MAT(A,i,j,LDA),x[j]);
    }
  }
}

// y=y+A*x
void zvec_add_lintr(int m, int n, dcomplex *y, const dcomplex *A, int LDA, const dcomplex *x)
{
  int i,j;
  for(i=0; i<m; i++){
    for(j=0; j<n; j++){
      Z_ADD_TIMES(y[i],MAT(A,i,j,LDA),x[j]);
    }
  }
}

// y=y-A*x
void zvec_sub_lintr(int m, int n, dcomplex *y, const dcomplex *A, int LDA, const dcomplex *x)
{
  int i,j;
  for(i=0; i<m; i++){
    for(j=0; j<n; j++){
      Z_SUB_TIMES(y[i],MAT(A,i,j,LDA),x[j]);
    }
  }
}

// y=A^T*x
void zvec_lintr_t(int m, int n, dcomplex *y, const dcomplex *A, int LDA, const dcomplex *x)
{
  int i,j;
  for(j=0; j<n; j++){
    Z_SET(y[j],0,0);
    for(i=0; i<m; i++){
      Z_ADD_TIMES(y[j],MAT(A,i,j,LDA),x[i]);
    }
  }
}

// y=y+A^T*x
void zvec_add_lintr_t(int m, int n, dcomplex *y, const dcomplex *A, int LDA, const dcomplex *x)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      Z_ADD_TIMES(y[j],MAT(A,i,j,LDA),x[i]);
    }
  }
}

// y=y-A^T*x
void zvec_sub_lintr_t(int m, int n, dcomplex *y, const dcomplex *A, int LDA, const dcomplex *x)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      Z_SUB_TIMES(y[j],MAT(A,i,j,LDA),x[i]);
    }
  }
}

// y=A'*x
void zvec_lintr_ct(int m, int n, dcomplex *y, const dcomplex *A, int LDA, const dcomplex *x)
{
  int i,j;
  for(j=0; j<n; j++){
    Z_SET(y[j],0,0);
    for(i=0; i<m; i++){
      Z_ADD_DOT(y[j],MAT(A,i,j,LDA),x[i]);
    }
  }
}


// y=y+A'*x
void zvec_add_lintr_ct(int m, int n, dcomplex *y, const dcomplex *A, int LDA, const dcomplex *x)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      Z_ADD_DOT(y[j],MAT(A,i,j,LDA),x[i]);
    }
  }
}

// y=y-A'*x
void zvec_sub_lintr_ct(int m, int n, dcomplex *y, const dcomplex *A, int LDA, const dcomplex *x)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      Z_SUB_DOT(y[j],MAT(A,i,j,LDA),x[i]);
    }
  }
}

// x=-x
void zvec_neg(int n, dcomplex *x)
{
  int i;
  for(i=0; i<n; i++){
    Z_NEG(x[i]);
  }
}


// x=a*x
void zvec_scale(int n, dcomplex *x, dcomplex a)
{
  int i;
  dcomplex value;
  for(i=0; i<n; i++){
    Z_SET_TIMES(value,a,x[i]);
    x[i]=value;
  }
}

// x=a*x
void zvec_scale_d(int n, dcomplex *x, double a)
{
  int i;
  for(i=0; i<n; i++){
    Z_SCALE(x[i],a);
  }
}

// x=(-a)*x
void zvec_scale_neg(int n, dcomplex *x, dcomplex a)
{
  dcomplex b;
  b=a;
  Z_NEG(b);
  zvec_scale(n,x,b);
}

// x=x/a
void zvec_scale_div(int n, dcomplex *x, dcomplex a)
{
  dcomplex b;
  Z_SET_INV(b,a);
  zvec_scale(n,x,b);
}

// x=x/sqrt(x'*x)
void zvec_normalize(int n, dcomplex *x)
{
  double norm=0;
  norm=zvec_norm2(n,x);
  zvec_scale_d(n,x,(1.0/norm));
}

// x=x/sqrt(x'*x) where x[i]/|x[i]|=1, abs(x[i]) is max
void zvec_normalize_sgn(int n, dcomplex *x)
{
  int k;
  double norm=0;
  dcomplex a;
  norm=zvec_norm2(n,x);
  zvec_max_abs_index(n,x,&k);
  a=znormalize(x[k]);
  Z_R(a)=Z_CONJ_R(a);
  Z_I(a)=Z_CONJ_I(a);
  Z_R(a)/=norm;
  Z_I(a)/=norm;
  zvec_scale(n,x,a);
}

// y=y-(x'*y)*x where x'*x=1
void zvec_orthogonalize(int n, dcomplex *y, const dcomplex *x)
{
  dcomplex dot;
  dot=zvec_dot(n,x,y);
  Z_R(dot)=-Z_R(dot);
  Z_I(dot)=-Z_I(dot);
  zvec_add_scaled(n,y,dot,x);
}

// x <=> y
void zvec_swap(int n, dcomplex *x, dcomplex *y)
{
  int i;
  dcomplex value;
  for(i=0; i<n; i++){
    value=x[i];
    x[i]=y[i];
    y[i]=value;
  }
}

// x'*y
dcomplex zvec_dot(int n, const dcomplex *x, const dcomplex *y)
{
  int i;
  dcomplex value;
  Z_SET(value,0,0);
  for(i=0; i<n; i++){
    Z_ADD_DOT(value,x[i],y[i]);
  }
  return value;
}


// abs(x'*y)/\sqrt(x'*x)/sqrt(y'*y)
double zvec_dcos_abs(int n, const dcomplex *x, const dcomplex *y)
{
  double norm_x,norm_y,dcos;
  dcomplex dot;
  norm_x=zvec_norm2(n,x);
  norm_y=zvec_norm2(n,y);
  dot=zvec_dot(n,x,y);
  dcos=Z_ABS(dot);
  dcos/=norm_x*norm_y;
  return dcos;
}

// acos(abs(x'*y)/\sqrt(x'*x)/sqrt(y'*y))
double zvec_angle(int n, const dcomplex *x, const dcomplex *y)
{
  double theta;
  theta=zvec_dcos_abs(n,x,y);
  if(theta>1) theta=1;
  theta=acos(theta);
  return theta;
}

// (180/PI)*acos(abs(x'*y)/\sqrt(x'*x)/sqrt(y'*y))
double zvec_angle_deg(int n, const dcomplex *x, const dcomplex *y)
{
  return zvec_angle(n,x,y)*(180.0/M_PI);
}

// value+=|x[i]-y[i]| (0<=i<n)
double zvec_dist_norm1(int n, const dcomplex *x, const dcomplex *y)
{
  int i;
  double value;
  value=0.0;
  for(i=0; i<n; i++) value+=Z_DIST(x[i],y[i]);
  return value;
}

// max(|x-y|)
double zvec_dist_norm_max(int n, const dcomplex *x, const dcomplex *y)
{
  int i;
  double value,a;
  value=Z_DIST2(x[0],y[0]);
  for(i=1; i<n; i++){
    a=Z_DIST2(x[i],y[i]);
    if(a>value) value=a;
  }
  value=sqrt(value);
  return value;
}

//////////////////////////////////////////////////////////////////////////

// z=abs(x-y)
void zvec_abs_sub(int n, double *z, const dcomplex *x, const dcomplex *y)
{
  int i;
  for(i=0; i<n; i++){ z[i]=Z_DIST(x[i],y[i]); }
}

// z=abs(x-y)
void zvec_abs_sub_scalar(int n, double *z, const dcomplex *x, dcomplex y)
{
  int i;
  for(i=0; i<n; i++){ z[i]=Z_DIST(x[i],y); }
}


//////////////////////////////////////////////////////////////////////

// sum(abs(x))
double zvec_norm1(int n, const dcomplex *x)
{
  int i;
  double value=0;
  for(i=0; i<n; i++){
    value+=Z_ABS(x[i]);
  }
  return value;
}

// sqrt(sum(abs(x).^2))
double zvec_norm2(int n, const dcomplex *x)
{
  int i;
  double value=0;
  for(i=0; i<n; i++){
    value+=Z_ABS2(x[i]);
  }
  value=sqrt(value);
  return value;
}

// max(abs(x))
double zvec_norm_max(int n, const dcomplex *x)
{
  int i;
  double value=0,foo;
  for(i=0; i<n; i++){
    foo=Z_ABS2(x[i]);
    if(foo>value) value=foo;
  }
  value=sqrt(value);
  return value;
}

//////////////////////////////////////////////////

// sum(abs(x).^2)
double zvec_sum_pow2_abs(int n, const dcomplex *x)
{
  double value=0;
  int i;
  for(i=0; i<n; i++){
    value+=Z_ABS2(x[i]);
  }
  return value;
}


//////////////////////////////////////////////////////////////

// max(abs(x))
double zvec_max_abs(int n, const dcomplex *x)
{
  int i;
  double value=0,foo;
  for(i=0; i<n; i++){
    foo=Z_ABS2(x[i]);
    if(foo>value){
      value=foo;
    }
  }
  value=sqrt(value);
  return value;
}

// min(abs(x))
double zvec_min_abs(int n, const dcomplex *x)
{
  int i;
  double value,foo;
  value=Z_ABS2(x[0]);
  for(i=1; i<n; i++){
    foo=Z_ABS2(x[i]);
    if(foo<value){
      value=foo;
    }
  }
  value=sqrt(value);
  return value;
}

//////////////////////////////////////////////////////////////

// max(abs(x))
// (*I)=k where k=max{ k | abs(x[k]) } unless I==NULL
double zvec_max_abs_index(int n, const dcomplex *x, int *I)
{
  int i;
  double value=0,foo;
  value=Z_ABS2(x[0]);
  if(I!=NULL) (*I)=0;
  for(i=1; i<n; i++){
    foo=Z_ABS2(x[i]);
    if(foo>value){
      value=foo;
      if(I!=NULL) (*I)=i;
    }
  }
  value=sqrt(value);
  return value;
}

// min(abs(x))
// (*I)=k where k=min{ k | abs(x[k]) } unless I==NULL
double zvec_min_abs_index(int n, const dcomplex *x, int *I)
{
  int i;
  double value,foo;
  value=Z_ABS2(x[0]);
  if(I!=NULL) (*I)=0;
  for(i=1; i<n; i++){
    foo=Z_ABS2(x[i]);
    if(foo<value){
      value=foo;
      if(I!=NULL) (*I)=i;
    }
  }
  value=sqrt(value);
  return value;
}

//////////////////////////////////////////////////////////////


// max(x)
dcomplex zvec_max(int n, const dcomplex *x)
{
  int i;
  dcomplex value;
  value=x[0];
  for(i=1; i<n; i++){
    if(zgt(x[i],value)) value=x[i];
  }
  return value;
}

// min(x)
dcomplex zvec_min(int n, const dcomplex *x)
{
  int i;
  dcomplex value=x[0];
  for(i=1; i<n; i++){
    if(zlt(x[i],value)) value=x[i];
  }
  return value;
}

///////////////////////////////////////////////////////////////////

// sum(x)
dcomplex zvec_sum(int n, const dcomplex *x)
{
  dcomplex value;
  int i;
  Z_SET(value,0,0);
  for(i=0; i<n; i++){
    Z_ADD(value,x[i]);
  }
  return value;
}

// sum(x)/n
dcomplex zvec_avarage(int n, const dcomplex *x)
{
  dcomplex value;
  value=zvec_sum(n,x);
  Z_R(value)/=n;
  Z_I(value)/=n;
  return value;
}

///////////////////////////////////////////////////////////////////

// y=int(x)
void zvec_get_si(int n, int *y, dcomplex *x)
{
  int i;
  for(i=0; i<n; i++){ y[i]=Z_R(x[i]); }
}

void zvec_get_s(int n, char **y, dcomplex *x, char format, int digits)
{
  char f[1024];
  int i;
  sprintf(f,"%%-.%d%c%%+.%d%ci",digits,format,digits,format);
  for(i=0; i<n; i++){ y[i]=char_renew_sprintf(y[i],NULL,f,Z_R(x[i]),Z_I(x[i])); }
}

///////////////////////////////////////////////////////////////////

void zvec_print(int n, dcomplex *x, char *name, char format, int digits)
{
  char **s=NULL;
  if(x==NULL){
    if(name!=NULL){ printf("%s=NULL\n",name); }
    else          { printf("NULL\n"); }
    return;
  }
  s=svec_allocate(n);
  zvec_get_s(n,s,x,format,digits);
  svec_print(n,s,name);
  s=svec_free(n,s);
  /*  
  int i;
  char format[128];
  if(STR_EQ(f,"e")){ sprintf(format,"(%%%d.%d%s, %%%d.%d%s)\n",digits+7,digits,f,digits+7,digits,f); }
  else             { sprintf(format,"%%-.%d%s%%+.%d%si\n",digits,f,digits,f); }
  if(x==NULL){
    if(name!=NULL){ printf("%s=NULL\n",name); }
    else          { printf("NULL\n"); }
    return;
  }
  if(name!=NULL){ printf("%s\n",name); }
  if(x==NULL) return;
  for(i=0; i<n; i++){
    printf(format,Z_R(x[i]),Z_I(x[i]));
  }
  */
}

///////////////////////////////////////////////////////////////////

void zvec_save(int n, dcomplex *x, int offset, char* fmt, ...)
{
  int i;
  char fname[FILE_NAME_LENGTH_MAX+1];
  va_list argp;
  FILE *fid;
  va_start(argp,fmt); vsprintf(fname,fmt,argp);   // file name
  if((fid=fopen(fname,"w"))==0){ ERROR_AT; printf("Cant't open file: %s\n",fname); exit(0); } // open
  for(i=0; i<n; i++){ fprintf(fid,"%d\t%+.20e\t%+.20e\n",i+offset,Z_R(x[i]),Z_I(x[i])); } // write
  fclose(fid); // close
}

///////////////////////////////////////////////////////////////////

void zvec_save_cplane(int n, dcomplex *x, char* fmt, ...)
{
  int i;
  char fname[FILE_NAME_LENGTH_MAX+1];
  va_list argp;
  FILE *fid;
  va_start(argp,fmt); vsprintf(fname,fmt,argp);   // file name
  if((fid=fopen(fname,"w"))==0){ ERROR_AT; printf("Cant't open file: %s\n",fname); exit(0); } // open
  for(i=0; i<n; i++){ fprintf(fid,"%+.20e\t%+.20e\n",Z_R(x[i]),Z_I(x[i])); } // write
  fclose(fid); // close
}

void zvec_save_cplane_label(int n, dcomplex *x, int offset, int *label, char* fmt, ...)
{
  int i;
  char fname[FILE_NAME_LENGTH_MAX+1];
  va_list argp;
  FILE *fid;
  va_start(argp,fmt); vsprintf(fname,fmt,argp); // file name
  if((fid=fopen(fname,"w"))==0){ ERROR_AT; printf("Cant't open file: %s\n",fname); exit(0); } // open
  for(i=0; i<n; i++){ fprintf(fid,"%+.20e\t%+.20e\t%d\n",Z_R(x[i]),Z_I(x[i]),label[i]+offset); } // write
  fclose(fid); // close
}

///////////////////////////////////////////////////////////////////

void zvec_bin_save(int n, dcomplex *x, char* fmt, ...)
{
  char fname[FILE_NAME_LENGTH_MAX+1];
  va_list argp;
  FILE *fid;
  va_start(argp,fmt); vsprintf(fname,fmt,argp); // file name
  if((fid=fopen(fname,"w"))==0){ ERROR_AT; printf("Cant't open file: %s\n",fname); exit(0); } // open
  fwrite("zvec",sizeof(char),strlen("zvec"),fid); // write header
  fwrite(&n,sizeof(int),1,fid); // write size
  fwrite(x,sizeof(dcomplex),n,fid); // write data
  fclose(fid); // close
}

dcomplex *zvec_bin_load(int *n, char* fmt, ...)
{
  int l;
  size_t k;
  dcomplex *zx=NULL;
  double *dx=NULL;
  char fname[FILE_NAME_LENGTH_MAX+1],*buf=NULL;
  va_list argp;
  FILE *fid;
   // file name
  va_start(argp,fmt); vsprintf(fname,fmt,argp);
  // open file
  if((fid=fopen(fname,"r"))==0){ fclose(fid); zx=NULL; (*n)=0; }
  else{
    // read header
    l=strlen("zvec");
    buf=malloc(sizeof(char)*(l+1));
    k=fread(buf,sizeof(char),l,fid);
    if(k==(size_t)l && strncmp(buf,"zvec",l)==0){ /* zvec */
      // read size
      k=fread(n,sizeof(int),1,fid);
      if(k!=1 || (*n)<=0){ ERROR_AT; printf("Failed to load the size from the file '%s'.\n",fname); exit(0); }
      // allocate
      zx=zvec_allocate((*n));
      // read data
      k=fread(zx,sizeof(dcomplex),(*n),fid);
      if(k!=(size_t)(*n)){ ERROR_AT; printf("Failed to load the data from the file '%s'.\n",fname); exit(0); }
      // close
      fclose(fid);
    }else if(k==(size_t)l && strncmp(buf,"dvec",l)==0){ /* dvec */
      fclose(fid);               // close
      dx=dvec_bin_load(n,fname); // load
      zx=zvec_allocate((*n));    // allocate
      zvec_copy_d((*n),zx,dx);   // copy
      dx=dvec_free(dx);          // free
    }else{ fclose(fid); zx=NULL; (*n)=0; }
  }
  // done
  free(buf);
  return zx;
}





/////////////////////////////////////////////////////////////////////////////

dcomplex zpow_si(dcomplex x, int n)
{
  int i;
  dcomplex y={0,0},a;
  if(n==0){
    Z_SET(y,1,0);
  }else if(n>0){
    Z_LET(y,x);
    for(i=1 ;i<n ;i++){
      Z_SET_TIMES(a,y,x); Z_LET(y,a);
    }
  }else if(n<0){
    Z_LET(y,x);
    for(i=1; i<(-n); i++){
      Z_SET_TIMES(a,y,x); Z_LET(y,a);
    }
    Z_SET_INV(a,y); Z_LET(y,a);
  }else{
    // error
  }
  return y;
}

//EOF
