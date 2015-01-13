#include<stdio.h>
#include<stdlib.h>
#include<stdarg.h>
#include"is_macros.h"
#include"is_ivec.h"

void ivec_scale(int n, int *x, int a)
{
  int i;
  for(i=0; i<n; i++) x[i]*=a;
}

// X!=Y
int ivec_ne(int n, const int *X, const int *Y)
{
  int i;
  if(n<=0) return 0;
  for(i=0; i<n; i++){
    if(X[i]!=Y[i]) return 1;
  }
  return 0;
}

// X==Y
int ivec_eq(int n, const int *X, const int *Y)
{
  int i;
  if(n<=0) return 0;
  for(i=0; i<n; i++){
    if(X[i]!=Y[i]) return 0;
  }
  return 1;
}

// X<=>Y
int ivec_cmp(int n, const int *X, int m, const int *Y)
{
  int i;
  for(i=0; i<MIN2(n,m); i++){
    if     (X[i]<Y[i]){ return -1; } // X<Y
    else if(X[i]>Y[i]){ return +1; } // X>Y
  }
  if     (n<m){ return -1; } // X<Y
  else if(n>m){ return +1; } // X>Y
  else        { return  0; } // X=Y
}

int ivec_first_index_if(int n, const int *X, int value)
{
  int i;
  for(i=0; i<n; i++){
    if(X[i]==value) return i;
  }
  return -1;
}

int ivec_first_index_if_not(int n, const int *X, int value)
{
  int i;
  for(i=0; i<n; i++){
    if(X[i]!=value) return i;
  }
  return -1;
}

int* ivec_allocate_index_if(int n, const int *X, int value, int *pm)
{
  int i,j,m,*I=NULL;
  m=ivec_count_if(n,X,value); (*pm)=m;
  if(m<=0){ I=NULL; }
  else{
    I=ivec_allocate(m);
    for(i=0,j=0; i<n; i++){
      if(X[i]==value) I[j++]=i;
    }
  }
  return I;
}

int* ivec_allocate_index_if_not(int n, const int *X, int value, int *pm)
{
  int i,j,m,*I=NULL;
  m=ivec_count_if_not(n,X,value); (*pm)=m;
  if(m<=0){ I=NULL; }
  else{
    I=ivec_allocate(m);
    for(i=0,j=0; i<n; i++){
      if(X[i]!=value) I[j++]=i;
    }
  }
  return I;
}

int ivec_count_if(int n, const int *X, int value)
{
  int i,count=0;
  for(i=0; i<n; i++){
    if(X[i]==value) count++;
  }
  return count;
}

int ivec_count_if_not(int n, const int *X, int value)
{
  int i,count=0;
  for(i=0; i<n; i++){
    if(X[i]!=value) count++;
  }
  return count;
}

//////////////////////////////////////////////////////////////////////////////////////////

// max(X)
int ivec_max(int n, const int *X)
{
  int value,i;
  value=X[0];
  for(i=1; i<n; i++){
    if(X[i]>value) value=X[i];
  }
  return value;
}

// min(X)
int ivec_min(int n, const int *X)
{
  int value,i;
  value=X[0];
  for(i=1; i<n; i++){
    if(X[i]<value) value=X[i];
  }
  return value;
}

//////////////////////////////////////////////////////////////////////////////////////////

// sum(x)
int ivec_sum(int n, const int *x)
{
  int value=0,i;
  for(i=0; i<n; i++) value+=x[i];
  return value;
}

// sum(x)/n
double ivec_average(int n, const int *x)
{
  double value;
  value=ivec_sum(n,x);
  value/=n;
  return value;
}


//////////////////////////////////////////////////////////////////////////////////////////

// sort X as X[0] <= X[1] <= X[2] <= ... <= X[n-1]
// if I==NULL, then I is not ussed.
// if I!=NULL, then I is stored with sorted indexes
void ivec_sort(int n, int *X, int *I)
{
  if(I!=NULL) ivec_grid(n,I);
  ivec_quick_sort(n,X,I,0,n-1);
}

void ivec_quick_sort(int n, int *X, int *I, int left, int right)
{
  int i,last;
  if(left>=right) return;
  ivec_swap_at(X,left,(left+right)/2);
  if(I!=NULL) ivec_swap_at(I,left,(left+right)/2);
  last=left;
  for(i=left+1; i<=right; i++){
    if(X[i]<X[left]){
      ++last;
      ivec_swap_at(X,last,i);
      if(I!=NULL) ivec_swap_at(I,last,i);      
    }
  }
  ivec_swap_at(X,left,last);
  if(I!=NULL) ivec_swap_at(I,left,last);
  ivec_quick_sort(n,X,I,left,last-1);
  ivec_quick_sort(n,X,I,last+1,right);
}

//////////////////////////////////////////////////////////////////////////////////////////

// store the list of indexes, x is not destroyed
void ivec_sort_index(int *I, int n, int *X)
{
  int *Y=NULL;
  Y=ivec_allocate(n);
  ivec_copy(n,Y,X);
  ivec_sort(n,Y,I);
  Y=ivec_free(Y);
}

//////////////////////////////////////////////////////////////////////////////////////////

// X[i] <=> X[j]
void ivec_swap_at(int *X, int i, int j)
{
  int foo;
  foo=X[i];
  X[i]=X[j];
  X[j]=foo;
}

// x[i] <=> x[I[i]]
void ivec_swap_index(int n, int *x, const int *I)
{
  int *y=NULL;
  y=ivec_allocate(n);
  ivec_copy_index(n,y,x,I);
  ivec_copy(n,x,y);
  y=ivec_free(y);
}

///////////////////////////////////////////////////////////////////////

// X[0]=X[n-1]; X[1]=X[n-2]; X[2]=X[n-3]; ...
void ivec_reverse(int n, int *x)
{
  int i;
  int a;
  for(i=0; i<n/2; i++){
    a=x[i];
    x[i]=x[n-i-1];
    x[n-i-1]=a;
  }
}

void ivec_relocate(int n, int *X, int *I)
{
  int i,*Y=NULL;
  Y=ivec_allocate(n);
  ivec_copy(n,Y,X);
  for(i=0; i<n; i++) X[i]=Y[I[i]];
  Y=ivec_free(Y);
}

// Y=X
void ivec_copy(int n, int *Y, const int *X)
{
  int i;
  for(i=0; i<n; i++) Y[i]=X[i];
}

// Y[i]=X[I[i]], 0<=i<n
void ivec_copy_index(int n, int *Y, const int *X, const int *I)
{
  int i;
  for(i=0; i<n; i++) Y[i]=X[I[i]];
}

// x <=> y
void ivec_swap(int n, int *x, int *y)
{
  int i;
  int value;
  for(i=0; i<n; i++){
    value=x[i];
    x[i]=y[i];
    y[i]=value;
  }
}



// X[0]=0; X[1]=1; X[2]=2; ...; X[n-1]=n-1
void ivec_grid(int n, int *X)
{
  int i;
  for(i=0; i<n; i++) X[i]=i;
}

// X=ones(n,1)*a
void ivec_set(int n, int *X, int a)
{
  int i;
  for(i=0; i<n; i++) X[i]=a;
}

// X=ones(n,1)
void ivec_ones(int n, int *X)
{
  int i;
  for(i=0;i<n;i++) X[i]=1;
}

// X=zeros(n,1)
void ivec_zeros(int n, int *X)
{
  int i;
  for(i=0; i<n; i++) X[i]=0;
}

int* ivec_allocate(int n)
{
  return (int*)malloc(sizeof(int)*n);
}

int* ivec_free(int *x)
{
  if(x==NULL) return NULL;
  free(x);
  return x=NULL;
}

//////////////////////////////////////////////////////////////////

// y+=x
void ivec_add(int n, int *y, int *x)
{
  int i;
  for(i=0; i<n; i++) y[i]+=x[i];
}

// y+=x
void ivec_add_scalar(int n, int *y, int x)
{
  int i;
  for(i=0; i<n; i++) y[i]+=x;
}

//////////////////////////////////////////////////////////////////

// y-=x
void ivec_sub(int n, int *y, int *x)
{
  int i;
  for(i=0; i<n; i++) y[i]-=x[i];
}


//////////////////////////////////////////////////////////////////

void ivec_print(int n, const int *X, char *name)
{
  int i;
  if(X==NULL){
    if(name!=NULL){ printf("%s NULL\n",name); }
    else{ printf("NULL\n"); }
    return;
  }
  if(name!=NULL){ printf("%s",name); }
  if(X==NULL) return;
  for(i=0; i<n; i++){
    printf("%d ",X[i]);
  }
  printf("\n");
}

void ivec_save(int n, int *x, int offset, char *fmt, ...)
{
  int i;
  char fname[100000];
  va_list argp;
  FILE *fid;
   // file name
  va_start(argp,fmt);
  vsprintf(fname,fmt,argp);
  // open file
  if((fid=fopen(fname,"w"))==0){
    ERROR_EXIT("Cant't open file: %s\n",fname);
  }
  // loop
  for(i=0; i<n; i++){
    fprintf(fid,"%d\t%d\n",i+offset,x[i]);
  }
  // close
  fclose(fid);
  fprintf(stderr,"saved: %s\n",fname);
}



//EOF
