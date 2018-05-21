#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<stdarg.h>
#include"is_macros.h"
#include"is_ivec.h"
#include"is_dvec.h"
#include"is_strings.h"
#include"is_svec.h"
#include"is_rmulti.h"
#include"mt19937ar.h"

/**
 @file  dvec.c
 @brief 倍精度型doubleのベクトルに関する関数の定義
 */


#define FILE_NAME_LENGTH_MAX 100000

/////////////////////////////////////////////////////////////

/** @name allocation */
/** @{ */

/**
 @brief double型ベクトルのメモリ割当
 */
double *dvec_allocate(int n)
{
  return (double*)malloc(sizeof(double)*n);
}

/**
 @brief double型ベクトルのメモリ割当とコピー
 */
double *dvec_allocate_clone(int n, double *x)
{
  double *y=NULL;
  y=dvec_allocate(n);
  dvec_set_dvec(n,y,x);
  return y;
}

/**
 @brief double型ベクトルのメモリの解放
 */
double *dvec_free(double *x)
{
  if(x==NULL) return NULL;
  free(x);
  return x=NULL;
}

/** @} */

/////////////////////////////////////////////////////////////

/** @name setting */
/** @{ */

/**
 @brief x=nan(n,1)
 */
void dvec_set_nan(int n, double *x)
{
  int i;
  for(i=0; i<n; i++){ x[i]=NAN; }
}

/**
 @brief x=inf(n,1)
 */
void dvec_set_inf(int n, double *x, int sgn)
{
  double zero=0;
  int i;
  for(i=0; i<n; i++){ x[i]=sgn/zero; }
}

/**
 @brief x=zeros(n,1)
 */
void dvec_set_zeros(int n, double *x)
{
  int i;
  for(i=0; i<n; i++){ x[i]=0; }
}

/**
 @brief x=ones(n,1)
 */
void dvec_set_ones(int n, double *x)
{
  int i;
  for(i=0; i<n; i++){ x[i]=1; }
}

/**
 @brief x=ones(n,1)*a
 */
void dvec_set_all_d(int n, double *x, double a)
{
  int i;
  for(i=0; i<n; i++){ x[i]=a; }
}

/**
 @brief x=zeros(n,1); x[k]=1
 */
void dvec_set_unit(int n, double *x, int k)
{
  dvec_set_zeros(n,x);
  x[k]=1;
}

/**
 @brief x[0]=0; x[1]=1; x[2]=2; ...; x[n-1]=n-1
 */
void dvec_set_grid(int n, double *x)
{
  int i;
  for(i=0; i<n; i++) x[i]=i;
}

/**
 @brief x=rand(n,1)*a+b
 */
void dvec_set_rand(int n, double *x, double a, double b){
  int i;
  for(i=0; i<n; i++){
    x[i]=genrand_int32();
    x[i]/=0xffffffff;
    x[i]=x[i]*a+b;
  }
}

/** @} */

/////////////////////////////////////////////////////////////

/** @name double型のベクトルの値の設定 */
/** @{ */

/**
 @brief y=x
 */
void dvec_set_dvec(int n, double *y, double *x)
{
  int i;
  for(i=0; i<n; i++){ y[i]=x[i]; }
}

/**
 @brief int ->  double
 */
void dvec_set_ivec(int n, double *y, int *x)
{
  int i;
  for(i=0; i<n; i++){ y[i]=x[i]; }
}

/**
 @brief double -> int
 */
void dvec_get_ivec(int n, int *y, double *x)
{
  int i;
  for(i=0; i<n; i++){ y[i]=x[i]; }
}

/**
 @brief double -> int
 */
void ivec_set_dvec(int n, int *y, double *x)
{
  int i;
  for(i=0; i<n; i++){ y[i]=x[i]; }
}

/**
 @brief char -> double
 */
void dvec_set_svec(int n, double *y, char **x)
{
  int i;
  rmulti *a=NULL;
  a=rallocate();
  for(i=0; i<n; i++){
    rset_s(a,x[i]);
    y[i]=rget_d(a);    
  }
  a=rfree(a);
}

/**
 @brief double -> char
 */
void dvec_get_svec(int n, char **y, double *x, char format, int digits)
{
  char f[1024];
  int i;
  sprintf(f,"%%-.%d%c",digits,format);
  for(i=0; i<n; i++){ y[i]=char_renew_sprintf(y[i],NULL,f,x[i]); }
}

/** @} */

/////////////////////////////////////////////////////////////

/** @name double型ベクトルの要素の並び替え */
/** @{ */

/**
 @brief x <=> y
 */
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

/**
 @brief x[0]=x[n-1]; x[1]=x[n-2]; x[2]=x[n-3]; ...
 */
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

/**
 @brief x[i] <=> x[j]
 */
void dvec_swap_at(double *x, int i, int j)
{
  double foo;
  foo=x[i];
  x[i]=x[j];
  x[j]=foo;
}

/**
 @brief x[i] <=> x[I[i]]
 */
void dvec_swap_index(int n, double *x, int *I)
{
  double *y=NULL;
  y=dvec_allocate(n);
  dvec_copy_dvec_index(n,y,x,I);
  dvec_copy_dvec(n,x,y);
  y=dvec_free(y);
}

/**
 @brief Don't call this function directly!
 */
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

/**
 @brief sort x as x[0] <= x[1] <= x[2] <= ... <= x[n-1]
 @detail if I==NULL, then I is not ussed. if I!=NULL, then I is stored with sorted indexes
 */
void dvec_sort(int n, double *x, int *I)
{
  if(I!=NULL){ ivec_set_grid(n,I); }
  dvec_quick_sort(n,x,I,0,n-1);
}

/**
 @brief store the list of indexes, x is not destroyed
 */
void dvec_sort_index(int *I, int n, double *X)
{
  double *Y=NULL;
  Y=dvec_allocate(n);
  dvec_copy_dvec(n,Y,X);
  dvec_sort(n,Y,I);
  Y=dvec_free(Y);
}

/** @} */

/////////////////////////////////////////////////////////////

/** @name double型ベクトルの関数 y=f(x) */
/** @{ */

/**
 @brief y=x
 */
void dvec_copy_dvec(int n, double *y, double *x)
{
  int i;
  for(i=0; i<n; i++){ y[i]=x[i]; }
}

// Y[i]=X[I[i]], 0<=i<n
void dvec_copy_dvec_index(int n, double *Y, double *X, int *I)
{
  int i;
  for(i=0; i<n; i++){ Y[i]=X[I[i]]; }
}

/**
 @brief y=-x
 */
void dvec_neg_dvec(int n, double *y, double *x)
{
  int i;
  for(i=0; i<n; i++){ y[i]=-x[i]; }
}

/**
 @brief y=x.^2
 */
void dvec_pow2_dvec(int n, double *y, double *x){
  int i;
  for(i=0; i<n; i++){ y[i]=x[i]*x[i]; }
}

/**
 @brief y=log2(abs(x))
 */
void dvec_log2_abs_dvec(int n, double *y, double *x)
{
  int i;
  for(i=0; i<n; i++) y[i]=log(fabs(x[i]))/log(2);
}

/**
 @brief y=x/sqrt(x'*x)
 */
void dvec_normalize_dvec(int n, double *y, double *x)
{
  double norm=0;
  norm=dnorm2_dvec(n,x);
  dvec_mul_dvec_dscalar(n,y,x,(1.0/norm));
}

/**
 @brief x=x/sqrt(x'*x) where x[i]>0, abs(x[i]) is max
 */
void dvec_normalize_sgn_dvec(int n, double *y, double *x)
{
  int k;
  double a=0;
  a=dnorm2_dvec(n,x);
  a=1.0/a;
  dmax_abs_dvec_index(n,x,&k);
  if(x[k]<0) a=-a;
  dvec_mul_dvec_dscalar(n,y,x,a);
}

/**
 @brief y=y-(x'*y)*x where x'*x=1
 */
void dvec_orthogonalize(int n, double *y, double *x)
{
  double dot;
  dot=dinnprod_dvec_dvec(n,x,y);
  dvec_sub_scaled(n,y,dot,x);
}

/**
 @brief y=sum(abs(x))
 */
double dnorm1_dvec(int n, double *x)
{
  int i;
  double value=0;
  for(i=0; i<n; i++){
    value+=fabs(x[i]);
  }
  return value;
}

/**
 @brief y=sqrt(sum(x.^2))
 */
double dnorm2_dvec(int n, double *x)
{
  int i;
  double value=0;
  for(i=0; i<n; i++){
    value+=x[i]*x[i];
  }
  value=sqrt(value);
  return value;
}

/**
 @brief y=max(abs(x))
 */
double dnorm_max_dvec(int n, double *x)
{
  int i;
  double value=0;
  for(i=0; i<n; i++){
    if(fabs(x[i])>value) value=fabs(x[i]);
  }
  return value;
}

/**
 @brief y=sum(x)
 */
double dsum_dvec(int n, double *x)
{
  double value=0;
  int i;
  for(i=0; i<n; i++) value+=x[i];
  return value;
}

/**
 @brief y=sum(x)/n
 */
double daverage_dvec(int n, double *x)
{
  return dsum_dvec(n,x)/n;
}

/**
 @brief y=sum(x.^2)
 */
double dsum_pow2_abs_dvec(int n, double *x)
{
  double value=0;
  int i;
  for(i=0; i<n; i++){
    value+=x[i]*x[i];
  }
  return value;
}

/**
 @brief y=max(abs(x))
 */
double dmax_abs_dvec(int n, double *x)
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

/**
 @brief y=min(abs(x))
 */
double dmin_abs_dvec(int n, double *x)
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

/**
 @brief y=max(abs(x[I]))
 */
double dmax_abs_dvec_index(int n, double *x, int *I)
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

/**
 @brief y=min(abs(x)), (*I)=k where k=min{ k | abs(x[k]) } unless I==NULL
 */
double dmin_abs_dvec_index(int n, double *x, int *I)
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

/**
 @brief y=max(x)
 */
double dmax_dvec(int n, double *x)
{
  int i;
  double value;
  value=x[0];
  for(i=1; i<n; i++){
    if(x[i]>value) value=x[i];
  }
  return value;
}

/**
 @brief y=min(x)
 */
double dmin_dvec(int n, double *x)
{
  int i;
  double value=x[0];
  for(i=1; i<n; i++){
    if(x[i]<value) value=x[i];
  }
  return value;
}

/** @} */

/////////////////////////////////////////////////////////////

/** @name double型ベクトルの関数 z=f(x,y) */
/** @{ */

/**
 @brief z=x+y
 */
void dvec_add_dvec_dvec(int n, double *z, double *x, double *y)
{
  int i;
  for(i=0; i<n; i++){ z[i]=x[i]+y[i]; }
}

/**
 @brief z=x+y
 */
void dvec_add_dvec_dscalar(int n, double *z, double *x, double y)
{
  int i;
  for(i=0; i<n; i++){ z[i]=x[i]+y; }
}

/**
 @brief z=x+y
 */
void dvec_add_dscalar_dvec(int n, double *z, double x, double *y)
{
  int i;
  for(i=0; i<n; i++){ z[i]=x+y[i]; }
}

////////

/**
 @brief z=x-y
 */
void dvec_sub_dvec_dvec(int n, double *z, double *x, double *y)
{
  int i;
  for(i=0; i<n; i++){ z[i]=x[i]-y[i]; }
}

/**
 @brief z=x-y
 */
void dvec_sub_dvec_dscalar(int n, double *z, double *x, double y)
{
  int i;
  for(i=0; i<n; i++){ z[i]=x[i]-y; }
}

/**
 @brief z=x-y
 */
void dvec_sub_dscalar_dvec(int n, double *z, double x, double *y)
{
  int i;
  for(i=0; i<n; i++){ z[i]=x-y[i]; }
}

////////

/**
 @brief z=x*y
 */
void dvec_mul_dvec_dvec(int n, double *z, double *x, double *y)
{
  int i;
  for(i=0; i<n; i++){ z[i]=x[i]*y[i]; }
}

/**
 @brief z=x*y
 */
void dvec_mul_dvec_dscalar(int n, double *z, double *x, double y)
{
  int i;
  for(i=0; i<n; i++){ z[i]=x[i]*y; }
}

/**
 @brief z=x*y
 */
void dvec_mul_dscalar_dvec(int n, double *z, double x, double *y)
{
  int i;
  for(i=0; i<n; i++){ z[i]=x*y[i]; }
}

////////

/**
 @brief z=x/y
 */
void dvec_div_dvec_dvec(int n, double *z, double *x, double *y)
{
  int i;
  for(i=0; i<n; i++){ z[i]=x[i]/y[i]; }
}

/**
 @brief z=x/y
 */
void dvec_div_dvec_dscalar(int n, double *z, double *x, double y)
{
  int i;
  for(i=0; i<n; i++){ z[i]=x[i]/y; }
}

/**
 @brief z=x/y
 */
void dvec_div_dscalar_dvec(int n, double *z, double x, double *y)
{
  int i;
  for(i=0; i<n; i++){ z[i]=x/y[i]; }
}

////////

/**
 @brief y=x^n
 */
double pow_i(double x, int n)
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
 @brief z=x^y
 */
void dvec_pow_dvec_iscalar(int n, double *z, double *x, int y)
{
  int i;
  for(i=0; i<n; i++){ z[i]=pow_i(x[i],(int)y); }
}

/**
 @brief double型のべき乗 z=x^y
 */
void dvec_pow_dvec_dscalar(int n, double *z, double *x, double y)
{
  int i;
  for(i=0; i<n; i++){ z[i]=pow(x[i],y); }
}

////////

/**
 @brief z=x'*y
 */
double dinnprod_dvec_dvec(int n, double *x, double *y)
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

/**
 @brief z=abs(x'*y)/\sqrt(x'*x)/sqrt(y'*y)
 */
double ddcos_abs_dvec_dvec(int n, double *x, double *y){
  double norm_x,norm_y,dcos,dot;
  norm_x=dnorm2_dvec(n,x);
  norm_y=dnorm2_dvec(n,y);
  dot=dinnprod_dvec_dvec(n,x,y);
  dcos=fabs(dot)/(norm_x*norm_y);
  return dcos;
}

/**
 @brief z=acos(abs(x'*y)/\sqrt(x'*x)/sqrt(y'*y))
 */
double dangle_dvec_dvec(int n, double *x, double *y)
{
  double theta;
  theta=ddcos_abs_dvec_dvec(n,x,y);
  if(theta>1) theta=1;
  theta=acos(theta);
  return theta;
}

/**
 @brief z=(180/PI)*acos(abs(x'*y)/\sqrt(x'*x)/sqrt(y'*y))
 */
double ddeg_angle_dvec_dvec(int n, double *x, double *y)
{
  return dangle_dvec_dvec(n,x,y)*(180.0/M_PI);
}

/**
 @brief z=sum(abs(x-y))
 */
double dsum_abs_sub_dvec_dvec(int n, double *x, double *y)
{
  int i;
  double value;
  value=0.0;
  for(i=0; i<n; i++) value+=fabs(x[i]-y[i]);
  return value;
}

/**
 @brief z=max(abs(x-y))
 */
double dmax_abs_sub_dvec_dvec(int n, double *x, double *y)
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

/**
 @brief z=abs(x-y)
 */
void dvec_abs_sub_dvec_dvec(int n, double *z, double *x, double *y)
{
  int i;
  for(i=0; i<n; i++){ z[i]=fabs(x[i]-y[i]); }
}

/**
 @brief z=abs(x-y)
 */
void dvec_abs_sub_dvec_dscalar(int n, double *z, double *x, double y)
{
  int i;
  for(i=0; i<n; i++){ z[i]=fabs(x[i]-y); }
}


/** @} */

/////////////////////////////////////////////////////////////

/** @name double型ベクトルの関数 y=y+f(x) */
/** @{ */

/**
 @brief y=y+a*x
 */
void dvec_add_scaled(int n, double *y, double a, double *x)
{
  int i;
  for(i=0; i<n; i++){ y[i]+=a*x[i]; }
}

/**
 @brief y=y-a*x
 */
void dvec_sub_scaled(int n, double *y, double a, double *x)
{
  int i;
  for(i=0; i<n; i++){ y[i]-=a*x[i]; }
}

/** @} */

/////////////////////////////////////////////////////////////

/** @name I/O */
/** @{ */

/**
 @brief load text formated data
 */
void dvec_load(int n, double *x, char* fmt, ...)
{
  int i;
  char fname[FILE_NAME_LENGTH_MAX];
  va_list argp;
  FILE *fid;
  // file name
  va_start(argp,fmt);
  vsprintf(fname,fmt,argp);
  // open file
  if((fid=fopen(fname,"r"))==0){ ERROR_AT; printf("Cant't open file: %s\n",fname); exit(0); }
  // real
  for(i=0;i<n;i++){ if(fscanf(fid,"%lf",&x[i])==EOF){ ERROR_EXIT("Error in dvec_load(), fscanf(), fname=%s\n",fname); } }
  // close
  fclose(fid);
}

/**
 @brief save data as text format
 */
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

/**
 @brief save 2D-data (x,y) as text format
 */
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

/**
 @brief save data log2(abs(x)) as text format
 */
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
    dvec_log2_abs_dvec(n,value,x);
    fprintf(fid,"%d\t%+.20e\n",i+offset,value[i]);
  } // write
  fclose(fid); // close
}

/**
 @brief save data as binary format
 */
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

/**
 @brief load binary formated data
 */
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

/** @} */

/////////////////////////////////////////////////////////////




void dvec_print(int n, double *x, char *name, char format, int digits)
{
  char **s=NULL;
  if(x==NULL){
    if(name!=NULL){ printf("%s=NULL\n",name); }
    else          { printf("NULL\n"); }
    return;
  }
  s=svec_allocate(n);
  dvec_get_svec(n,s,x,format,digits);
  svec_print(n,s,name);
  s=svec_free(n,s);
}

///////////////////////////////////////////////////////////////////////////////





//////////////////////////////////////////////////////////////////////////////


//EOF
