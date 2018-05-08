#include<stdio.h>
#include<stdlib.h>
#include<stdarg.h>
#include"is_macros.h"
#include"is_imat.h"
#include"is_ivec.h"
#include"is_strings.h"
#include"is_svec.h"
#include"mt19937ar.h"

// B+=A
void imat_add(int m, int n, int *B, int LDB, const int *A, int LDA)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      MAT(B,i,j,LDB)+=MAT(A,i,j,LDA);
    }
  }
}

// B+=a
void imat_add_scalar(int m, int n, int *B, int LDB, int a)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      MAT(B,i,j,LDB)+=a;
    }
  }
}

// B-=A
void imat_sub(int m, int n, int *B, int LDB, const int *A, int LDA)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      MAT(B,i,j,LDB)-=MAT(A,i,j,LDA);
    }
  }
}

// B*=a
void imat_scale(int m, int n, int *B, int LDB, int a)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      MAT(B,i,j,LDB)*=a;
    }
  }
}


///////////////////////////////////////

int* imat_allocate(int m, int n)
{
  return (int*)malloc(sizeof(int)*m*n);
}

int* imat_free(int *A)
{
  if(A==NULL) return NULL;
  free(A);
  return A=NULL;
}


// A=zeros(m,n)
void imat_zeros(int m, int n, int *A, int LDA)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      MAT(A,i,j,LDA)=0;
    }
  }
}

// A=ones(m,n)
void imat_ones(int m, int n, int *A, int LDA)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      MAT(A,i,j,LDA)=1;
    }
  }
}

// A=ones(m,n)*a
void imat_set(int m, int n, int *A, int LDA, int a)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      MAT(A,i,j,LDA)=a;
    }
  }
}

// A=eye(m,n)
void imat_eye(int m, int n, int *A, int LDA)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      if(i==j) MAT(A,i,j,LDA)=1;
      else  MAT(A,i,j,LDA)=0;
    }
  }
}

// A=rand(m,n)*a+b
void imat_rand(int m, int n, int *A, int LDA, int a, int b)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      MAT(A,i,j,LDA)=(int)genrand_int31()*a+b;
    }
  }
}

///////////////////////////////////////////////////////////////////////

// B=A
void imat_copy(int m, int n, int *B, int LDB, const int *A, int LDA)
{
  int j;
  for(j=0; j<n; j++){
    ivec_copy(m,&COL(B,j,LDB),&COL(A,j,LDA));
  }
}

///////////////////////////////////////////////////////////////////////

// B(:,j)=A(:,I(j)) for 0<=j<n
void imat_copy_col_index(int m, int n, int *B, int LDB, const int *A, int LDA, const int *I)
{
  int j;
  for(j=0; j<n; j++){
    ivec_copy(m,&COL(B,j,LDB),&COL(A,I[j],LDA));
  }
}

///////////////////////////////////////////////////////////////////////

// A(:,j) <-> A(:,I[j]) for 0<=j<n
void imat_swap_index(int m, int n, int *A, int LDA, const int *I)
{
  int LDB;
  int *B=NULL;
  LDB=m; B=imat_allocate(LDB,n);
  imat_copy_col_index(m,n,B,LDB,A,LDA,I);
  imat_copy(m,n,A,LDA,B,LDB);
  B=imat_free(B);
}


///////////////////////////////////////////////////////////////////////

// A(:,k) <-> A(:,l)
void imat_swap_columns(int m, int n, int *A, int LDA, int k, int l)
{
  if(k<n && l<n && k!=l){
    ivec_swap(m,&COL(A,k,LDA),&COL(A,l,LDA));
  }
}


// A(k,:) <-> A(l,:)
void imat_swap_rows(int m, int n, int *A, int LDA, int k, int l)
{
  int j,a;
  if(k<m && l<m){
    for(j=0; j<n; j++){
      a=MAT(A,k,j,LDA);
      MAT(A,k,j,LDA)=MAT(A,l,j,LDA);
      MAT(A,l,j,LDA)=a;
    }
  }
}

//////////////////////////////////////

int imat_row_count_if(int m, int n, int *A, int LDA, int i, int value)
{
  int j,count=0;
  if(i<m){
    for(j=0; j<n; j++){
      if(MAT(A,i,j,LDA)==value) count++;
    }
  }
  return count;
}

int imat_column_count_if(int m, int n, int *A, int LDA, int j, int value)
{
  int i,count=0;
  if(j<n){
    for(i=0; i<m; i++){
      if(MAT(A,i,j,LDA)==value) count++;
    }
  }
  return count;
}

//////////////////////////////////////

// y=char(x)
void imat_get_s(int m, int n, char **B, int LDB, int *A, int LDA)
{
  int i,j;
  if(A==NULL){ return; }
  for(i=0; i<m; i++){
    for(j=0; j<n; j++){
      MAT(B,i,j,LDB)=char_renew_sprintf(MAT(B,i,j,LDB),NULL,"%d",MAT(A,i,j,LDA));
    }
  }
}

//////////////////////////////////////

void imat_print(int m, int n, int *A, int LDA, char *name)
{
  char **s=NULL;
  if(A==NULL){
    if(name!=NULL){ printf("%s=NULL\n",name); }
    else          { printf("NULL\n"); }
    return;
  }  
  s=svec_allocate(m*n);
  imat_get_s(m,n,s,m,A,LDA);
  smat_print(m,n,s,m,name);
  s=svec_free(m*n,s);
}

void imat_save(int m, int n, int *A, int LDA, char* fmt, ...)
{
  int i,j;
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
  fprintf(fid,"%d %d\n",m,n);
  for(i=0; i<m; i++){
    for(j=0; j<n; j++){
      if(j!=0) fprintf(fid," ");
      fprintf(fid,"%d",MAT(A,i,j,LDA));
    }
    fprintf(fid,"\n");
  }
  // close
  fclose(fid);
  fprintf(stderr,"saved: %s\n",fname);
}

int* imat_load_allocate(int *pm, int *pn, char* fmt, ...)
{
  int i,j;
  char fname[100000];
  va_list argp;
  FILE *fid;
  int m,n,*A,LDA;
   // file name
  va_start(argp,fmt);
  vsprintf(fname,fmt,argp);
  // open file
  if((fid=fopen(fname,"r"))==0){
    ERROR_EXIT("Cant't open file: %s\n",fname);
  }
  // get size
  fscanf(fid,"%d",&m);
  fscanf(fid,"%d",&n);
  // allocate
  LDA=m; A=imat_allocate(LDA,n);
  // loop
  for(i=0; i<m; i++){
    for(j=0; j<n; j++){
      fscanf(fid,"%d",&MAT(A,i,j,LDA));
    }
  }
  // close
  fclose(fid);
  fprintf(stderr,"loaded: %s\n",fname);
  // done
  (*pm)=m;
  (*pn)=n;
  return A;
}

//EOF
