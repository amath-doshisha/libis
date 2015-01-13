#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<stdarg.h>
#include"is_macros.h"
#include"is_ivec.h"
#include"is_dvec.h"
#include"is_strings.h"
#include"mt19937ar.h"

#define FILE_NAME_LENGTH_MAX 100000

/////////////////////////////////////////////////////////////

// x=[ n x 1 ]
double *dvec_allocate(int n)
{
  return (double*)malloc(sizeof(double)*n);
}

double *dvec_free(double *x)
{
  if(x==NULL) return NULL;
  free(x);
  return x=NULL;
}

/////////////////////////////////////////////////////////////

double *dvec_allocate_s(char *str, int *n)
{
  int i;
  strings *line=NULL;
  double *x=NULL;
  line=strings_split(str,",; \t\n",NULL,NULL,NULL);
  if((line->n)>=0){
    (*n)=line->n;
    x=dvec_allocate((*n));
    for(i=0; i<line->n; i++){ x[i]=atof(line->str[i]); }
  }else{ (*n)=0; x=NULL; }
  line=strings_del(line);
  return x;
}

/////////////////////////////////////////////////////////////

// x=zeros(n,1)
void dvec_zeros(int n, double *x)
{
  int i;
  for(i=0; i<n; i++) x[i]=0;
}

// x=ones(n,1)
void dvec_ones(int n, double *x)
{
  int i;
  for(i=0; i<n; i++) x[i]=1;
}

// x=ones(n,1)*a
void dvec_set(int n, double *x, double a)
{
  int i;
  for(i=0; i<n; i++) x[i]=a;
}

// x=zeros(n,1); x[k]=1
void dvec_unit(int n, double *x, int k)
{
  dvec_zeros(n,x);
  x[k]=1;
}

// x[0]=0; x[1]=1; x[2]=2; ...; x[n-1]=n-1
void dvec_set_grid(int n, double *x)
{
  int i;
  for(i=0; i<n; i++) x[i]=i;
}

// x=rand(n,1)*a+b
void dvec_rand(int n, double *x, double a, double b){
  int i;
  for(i=0; i<n; i++){
    x[i]=genrand_real3()*a+b;
  }
}

// x[0]=x[n-1]; x[1]=x[n-2]; x[2]=x[n-3]; ...
void dvec_reverse(int n, double *x)
{
  int i;
  double a;
  for(i=0; i<n/2; i++){
    a=x[i];
    x[i]=x[n-i-1];
    x[n-i-1]=a;
  }
}

///////////////////////////////////////////////////////////////////////////////

// x[i] <=> x[j]
void dvec_swap_at(double *x, int i, int j)
{
  double foo;
  foo=x[i];
  x[i]=x[j];
  x[j]=foo;
}

// x[i] <=> x[I[i]]
void dvec_swap_index(int n, double *x, const int *I)
{
  double *y=NULL;
  y=dvec_allocate(n);
  dvec_copy_index(n,y,x,I);
  dvec_copy(n,x,y);
  y=dvec_free(y);
}

///////////////////////////////////////////////////////////////////////////////

// Don't call this function directly!
void dvec_quick_sort(int n, double *x, int *I, int left, int right)
{
  int i,last;
  if(left>=right) return;
  dvec_swap_at(x,left,(left+right)/2);
  if(I!=NULL) ivec_swap_at(I,left,(left+right)/2);
  last=left;
  for(i=left+1; i<=right; i++){
    if(x[i]<x[left]){
      ++last;
      dvec_swap_at(x,last,i);
      if(I!=NULL) ivec_swap_at(I,last,i);      
    }
  }
  dvec_swap_at(x,left,last);
  if(I!=NULL) ivec_swap_at(I,left,last);
  dvec_quick_sort(n,x,I,left,last-1);
  dvec_quick_sort(n,x,I,last+1,right);
}

// sort x as x[0] <= x[1] <= x[2] <= ... <= x[n-1]
// if I==NULL, then I is not ussed.
// if I!=NULL, then I is stored with sorted indexes
void dvec_sort(int n, double *x, int *I)
{
  if(I!=NULL) ivec_grid(n,I);
  dvec_quick_sort(n,x,I,0,n-1);
}


//////////////////////////////////////////////////////////////////////////////////////////

// store the list of indexes, x is not destroyed
void dvec_sort_index(int *I, int n, double *X)
{
  double *Y=NULL;
  Y=dvec_allocate(n);
  dvec_copy(n,Y,X);
  dvec_sort(n,Y,I);
  Y=dvec_free(Y);
}

//////////////////////////////////////////////////////

// y=x
void dvec_copy(int n, double *y, const double *x)
{
  int i;
  for(i=0; i<n; i++){ y[i]=x[i]; }
}

// Y[i]=X[I[i]], 0<=i<n
void dvec_copy_index(int n, double *Y, const double *X, const int *I)
{
  int i;
  for(i=0; i<n; i++){ Y[i]=X[I[i]]; }
}

//////////////////////////////////////////////////////

// z=x+y
void dvec_add_d(int n, double *z, const double *x, double y)
{
  int i;
  for(i=0; i<n; i++){ z[i]=x[i]+y; }
}

// z=x/y
void dvec_div_d2(int n, double *z, const double *x, double y)
{
  int i;
  for(i=0; i<n; i++){ z[i]=x[i]/y; }
}

//////////////////////////////////////////////////////

// x=x+a
void dvec_add_scalar(int n, double *x, double a)
{
  int i;
  for(i=0; i<n; i++){
    x[i]+=a;
  }
}

// x=x-a
void dvec_sub_scalar(int n, double *x, double a)
{
  int i;
  for(i=0; i<n; i++){
    x[i]-=a;
  }
}

//////////////////////////////////////////////////////

// y=y+x
void dvec_add(int n, double *y, const double *x)
{
  int i;
  for(i=0; i<n; i++){
    y[i]+=x[i];
  }
}

// y=y-x
void dvec_sub(int n, double *y, const double *x)
{
  int i;
  for(i=0; i<n; i++){
    y[i]-=x[i];
  }
}

//////////////////////////////////////////////////////

// y=y+a*x
void dvec_add_scaled(int n, double *y, double a, const double *x)
{
  int i;
  for(i=0; i<n; i++){
    y[i]+=a*x[i];
  }
}

// y=y-a*x
void dvec_sub_scaled(int n, double *y, double a, const double *x)
{
  int i;
  for(i=0; i<n; i++){
    y[i]-=a*x[i];
  }
}

//////////////////////////////////////////////////////

// y=A*x
void dvec_lintr(int m, int n, double *y, const double *A, int LDA, const double *x){
  int i,j;
  for(i=0; i<m; i++){
    y[i]=0;
    for(j=0; j<n; j++){
      y[i]+=MAT(A,i,j,LDA)*x[j];
    }
  }
}

// y=y+A*x
void dvec_add_lintr(int m, int n, double *y, const double *A, int LDA, const double *x){
  int i,j;
  for(i=0; i<m; i++){
    for(j=0; j<n; j++){
      y[i]+=MAT(A,i,j,LDA)*x[j];
    }
  }
}

// y=y-A*x
void dvec_sub_lintr(int m, int n, double *y, const double *A, int LDA, const double *x){
  int i,j;
  for(i=0; i<m; i++){
    for(j=0; j<n; j++){
      y[i]-=MAT(A,i,j,LDA)*x[j];
    }
  }
}

// y=A'*x
void dvec_lintr_t(int m, int n, double *y, const double *A, int LDA, const double *x){
  int i,j;
  for(j=0; j<n; j++){
    y[j]=0;
    for(i=0; i<m; i++){
      y[j]+=MAT(A,i,j,LDA)*x[i];
    }
  }
}

// y=y+A'*x
void dvec_add_lintr_t(int m, int n, double *y, const double *A, int LDA, const double *x){
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      y[j]+=MAT(A,i,j,LDA)*x[i];
    }
  }
}

// y=y-A'*x
void dvec_sub_lintr_t(int m, int n, double *y, const double *A, int LDA, const double *x){
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      y[j]-=MAT(A,i,j,LDA)*x[i];
    }
  }
}

// x=-x
void dvec_neg(int n, double *x)
{
  int i;
  for(i=0; i<n; i++){
    x[i]=-x[i];
  }
}

// x=a*x
void dvec_scale(int n, double *x, double a)
{
  int i;
  for(i=0; i<n; i++){
    x[i]*=a;
  }
}

// x=x/sqrt(x'*x)
void dvec_normalize(int n, double *x)
{
  double norm=0;
  norm=dvec_norm2(n,x);
  dvec_scale(n,x,(1.0/norm));
}

// x=x/sqrt(x'*x) where x[i]>0, abs(x[i]) is max
void dvec_normalize_sgn(int n, double *x)
{
  int k;
  double a=0;
  a=dvec_norm2(n,x);
  a=1.0/a;
  dvec_max_abs_index(n,x,&k);
  if(x[k]<0) a=-a;
  dvec_scale(n,x,a);
}

// y=y-(x'*y)*x where x'*x=1
void dvec_orthogonalize(int n, double *y, const double *x)
{
  double dot;
  dot=dvec_dot(n,x,y);
  dvec_sub_scaled(n,y,dot,x);
}

/**
 @brief double型のべき乗 z=x.^y
*/
void dvec_pow(int n, double *z, double *x, double *y) 
{
  int i;
  for(i=0; i<n; i++){ z[i]=pow(x[i],y[i]); }
}

double pow_si(double x, int n)
{
  int i;
  double y=0;
  if(n==0){ y=1; }
  else if(n>0){
    y=x;
    for(i=1 ;i<n ;i++){ y=y*x; }
  }
  else if(n<0){
    y=x;
    for(i=1; i<(-n); i++){ y=y*x; }
    y=1/y;
  }
  return y;
}

/**
 @brief double型のべき乗 z=x^y
*/
void dvec_pow_si(int n, double *z, double *x, long y)
{
  int i;
  for(i=0; i<n; i++){ z[i]=pow_si(x[i],y); }
}

/**
 @brief double型のべき乗 z=x^y
*/
void dvec_pow_d(int n, double *z, double *x, double y)
{
  int i;
  for(i=0; i<n; i++){ z[i]=pow(x[i],y); }
}


// y=x.^2
void dvec_pow2(int n, double *y, const double *x){
  int i;
  for(i=0; i<n; i++) y[i]=x[i]*x[i];
}

// y=log2(abs(x))
void dvec_log2_abs(int n, double *y, const double *x)
{
  int i;
  for(i=0; i<n; i++) y[i]=log(fabs(x[i]))/log(2);
}

// x <=> y
void dvec_swap(int n, double *x, double *y)
{
  int i;
  double value;
  for(i=0; i<n; i++){
    value=x[i];
    x[i]=y[i];
    y[i]=value;
  }
}

// x'*y
double dvec_dot(int n, const double *x, const double *y)
{
  int i;
  double value=0;
  for(i=0; i<n; i++){
    value+=x[i]*y[i];
  }
  return value;
  /*
  int i;
  double a1=0,a2=0,*Z=NULL;
  Z=dvec_allocate(n);
  for(i=0; i<n; i++){
    Z[i]=x[i]*y[i];
  }
  for(i=0; i<n; i++){
    if(Z[i]<0) a2+=Z[i];
    else a1+=Z[i];
  }
  FREE(Z);
  return a1+a2;
  */
}

// abs(x'*y)/\sqrt(x'*x)/sqrt(y'*y)
double dvec_dcos_abs(int n, const double *x, const double *y){
  double norm_x,norm_y,dcos,dot;
  norm_x=dvec_norm2(n,x);
  norm_y=dvec_norm2(n,y);
  dot=dvec_dot(n,x,y);
  dcos=fabs(dot)/(norm_x*norm_y);
  return dcos;
}

// acos(abs(x'*y)/\sqrt(x'*x)/sqrt(y'*y))
double dvec_angle(int n, const double *x, const double *y)
{
  double theta;
  theta=dvec_dcos_abs(n,x,y);
  if(theta>1) theta=1;
  theta=acos(theta);
  return theta;
}

// (180/PI)*acos(abs(x'*y)/\sqrt(x'*x)/sqrt(y'*y))
double dvec_angle_deg(int n, const double *x, const double *y)
{
  return dvec_angle(n,x,y)*(180.0/M_PI);
}

// value=sum(abs(x-y))
double dvec_dist_norm1(int n, const double *x, const double *y)
{
  int i;
  double value;
  value=0.0;
  for(i=0; i<n; i++) value+=fabs(x[i]-y[i]);
  return value;
}

// value=max(abs(x-y))
double dvec_dist_norm_max(int n, const double *x, const double *y)
{
  int i;
  double value,a;
  value=fabs(x[0]-y[0]);
  for(i=1; i<n; i++){
    a=fabs(x[i]-y[i]);
    if(a>value) value=a;
  }
  return value;
}

//////////////////////////////////////////////////////////////////////////

// z=abs(x-y)
void dvec_abs_sub(int n, double *z, const double *x, const double *y)
{
  int i;
  for(i=0; i<n; i++){ z[i]=fabs(x[i]-y[i]); }
}

// z=abs(x-y)
void dvec_abs_sub_scalar(int n, double *z, const double *x, double y)
{
  int i;
  for(i=0; i<n; i++){ z[i]=fabs(x[i]-y); }
}

//////////////////////////////////////////////////////////////////////////

// sum(abs(x))
double dvec_norm1(int n, const double *x)
{
  int i;
  double value=0;
  for(i=0; i<n; i++){
    value+=fabs(x[i]);
  }
  return value;
}

// sqrt(sum(x.^2))
double dvec_norm2(int n, const double *x)
{
  int i;
  double value=0;
  for(i=0; i<n; i++){
    value+=x[i]*x[i];
  }
  value=sqrt(value);
  return value;
}

// max(abs(x))
double dvec_norm_max(int n, const double *x)
{
  int i;
  double value=0;
  for(i=0; i<n; i++){
    if(fabs(x[i])>value) value=fabs(x[i]);
  }
  return value;
}

// sum(x)
double dvec_sum(int n, const double *x)
{
  double value=0;
  int i;
  for(i=0; i<n; i++) value+=x[i];
  return value;
}

// sum(x)/n
double dvec_average(int n, const double *x)
{
  return dvec_sum(n,x)/n;
}

// sum(x.^2)
double dvec_sum_pow2(int n, const double *x)
{
  double value=0;
  int i;
  for(i=0; i<n; i++){
    value+=x[i]*x[i];
  }
  return value;
}

///////////////////////////////////////////////////////////////

// max(abs(x))
double dvec_max_abs(int n, const double *x)
{
  int i;
  double value=0;
  value=fabs(x[0]);
  for(i=1; i<n; i++){
    if(fabs(x[i])>value){
      value=fabs(x[i]);
    }
  }
  return value;
}

// min(abs(x))
double dvec_min_abs(int n, const double *x)
{
  int i;
  double value;
  value=fabs(x[0]);
  for(i=1; i<n; i++){
    if(fabs(x[i])<value){
      value=fabs(x[i]);
    }
  }
  return value;
}


///////////////////////////////////////////////////////////////

double dvec_max_abs_index(int n, const double *x, int *I)
{
  int i;
  double value=0;
  value=fabs(x[0]);
  if(I!=NULL) (*I)=0;
  for(i=1; i<n; i++){
    if(fabs(x[i])>value){
      value=fabs(x[i]);
      if(I!=NULL) (*I)=i;
    }
  }
  return value;
}

// min(abs(x))
// (*I)=k where k=min{ k | abs(x[k]) } unless I==NULL
double dvec_min_abs_index(int n, const double *x, int *I)
{
  int i;
  double value;
  value=fabs(x[0]);
  if(I!=NULL) (*I)=0;
  for(i=1; i<n; i++){
    if(fabs(x[i])<value){
      value=fabs(x[i]);
      if(I!=NULL) (*I)=i;
    }
  }
  return value;
}

///////////////////////////////////////////////////////////////

// max(x)
double dvec_max(int n, const double *x)
{
  int i;
  double value;
  value=x[0];
  for(i=1; i<n; i++){
    if(x[i]>value) value=x[i];
  }
  return value;
}

// min(x)
double dvec_min(int n, const double *x)
{
  int i;
  double value=x[0];
  for(i=1; i<n; i++){
    if(x[i]<value) value=x[i];
  }
  return value;
}

void dvec_print(int n, const double *x, const char *name, const char *f, int digits)
{
  int i;
  char format[128];
  if(STR_EQ(f,"f")){ sprintf(format,"%%%d.%d%s\n",digits+3,digits,f); }
  else             { sprintf(format,"%%%d.%d%s\n",digits+9,digits,f); }
  if(x==NULL){
    if(name!=NULL){ printf("%s=NULL\n",name); }
    else          { printf("NULL\n"); }
    return;
  }
  if(name!=NULL){ printf("%s\n",name); }
  if(x==NULL) return;
  for(i=0; i<n; i++){
    printf((const char*)format,x[i]);
  }
}

///////////////////////////////////////////////////////////////////////////////

void dvec_save(int n, double *x, int offset, char* fmt, ...)
{
  int i;
  char fname[FILE_NAME_LENGTH_MAX+1];
  va_list argp;
  FILE *fid;
  va_start(argp,fmt); vsprintf(fname,fmt,argp); // file name
  if((fid=fopen(fname,"w"))==0){ ERROR_AT; printf("Cant't open file: %s\n",fname); exit(0); } // open
  for(i=0; i<n; i++){ fprintf(fid,"%d\t%+.20e\n",i+offset,x[i]); } // write
  fclose(fid); // close
}

void dvec_save2(int n, double *x, double *y, char* fmt, ...)
{
  int i;
  char fname[FILE_NAME_LENGTH_MAX+1];
  va_list argp;
  FILE *fid;
  va_start(argp,fmt); vsprintf(fname,fmt,argp); // file name
  if((fid=fopen(fname,"w"))==0){ ERROR_AT; printf("Cant't open file: %s\n",fname); exit(0); } // open
  for(i=0; i<n; i++){ fprintf(fid,"%+.20e\t%+.20e\n",x[i],y[i]); } // write
  fclose(fid); // close
}

void dvec_save_log2_abs(int n, double *x, int offset, char* fmt, ...)
{
  int i;
  char fname[FILE_NAME_LENGTH_MAX+1];
  double *value=NULL;
  va_list argp;
  FILE *fid;
  value=dvec_allocate(n);
  // file name
  va_start(argp,fmt); 
  vsprintf(fname,fmt,argp); 
  // open
  if((fid=fopen(fname,"w"))==0){ ERROR_AT; printf("Cant't open file: %s\n",fname); exit(0); } // open
  // write
  for(i=0; i<n; i++){ 
    dvec_log2_abs(n,value,x);
    fprintf(fid,"%d\t%+.20e\n",i+offset,value[i]);
  } // write
  fclose(fid); // close
}
///////////////////////////////////////////////////////////////////////////////

void dvec_bin_save(int n, double *x, char* fmt, ...)
{
  char fname[FILE_NAME_LENGTH_MAX+1];
  va_list argp;
  FILE *fid;
  va_start(argp,fmt); vsprintf(fname,fmt,argp); // file name
  if((fid=fopen(fname,"w"))==0){ ERROR_AT; printf("Cant't open file: %s\n",fname); exit(0); } // open
  fwrite("dvec",sizeof(char),strlen("dvec"),fid); // write header
  fwrite(&n,sizeof(int),1,fid); // write size
  fwrite(x,sizeof(double),n,fid); // write data
  fclose(fid); // close
}

double *dvec_bin_load(int *n, char* fmt, ...)
{
  int l;
  size_t k;
  double *x=NULL;
  char fname[FILE_NAME_LENGTH_MAX+1],*buf=NULL;
  va_list argp;
  FILE *fid;
  // file name
  va_start(argp,fmt); vsprintf(fname,fmt,argp);
  // open
  if((fid=fopen(fname,"r"))==0){ x=NULL; (*n)=0; }
  else{
    // read header
    l=strlen("dvec");
    buf=malloc(sizeof(char)*(l+1));
    k=fread(buf,sizeof(char),l,fid);
    if(k!=(size_t)l || strncmp(buf,"dvec",l)!=0){ x=NULL; (*n)=0; }
    else{
      // read size
      k=fread(n,sizeof(int),1,fid);
      if(k!=1 || (*n)<=0){ ERROR_AT; printf("Failed to load the size from the file '%s'.\n",fname); exit(0); }
      // allocate
      x=dvec_allocate((*n));
      // read data
      k=fread(x,sizeof(double),(*n),fid);
      if(k!=(size_t)(*n)){ ERROR_AT; printf("Failed to load the data from the file '%s'.\n",fname); exit(0); }
    }
  }
  // close
  fclose(fid);
  // done
  free(buf);
  return x;
}

//////////////////////////////////////////////////////////////////////////////


//EOF
