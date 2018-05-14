#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<stdarg.h>
#include"is_macros.h"
#include"is_dvec.h"
#include"is_dmat.h"
#include"is_strings.h"
#include"is_svec.h"
#include"mt19937ar.h"

/////////////////////////////////////////

#define FILE_NAME_LENGTH_MAX 100000

/////////////////////////////////////////

double* dmat_allocate(int m, int n)
{
  return (double*)malloc(sizeof(double)*m*n);
}

double* dmat_free(double *A)
{
  if(A==NULL) return NULL;
  free(A);
  return A=NULL;
}

// A=zeros(m,n)
void dmat_set_zeros(int m, int n, double *A, int LDA)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      MAT(A,i,j,LDA)=0;
    }
  }
}

// A=ones(m,n)
void dmat_set_ones(int m, int n, double *A, int LDA)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      MAT(A,i,j,LDA)=1;
    }
  }
}

// A=ones(m,n)*a
void dmat_set_all_d(int m, int n, double *A, int LDA, double a)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      MAT(A,i,j,LDA)=a;
    }
  }
}

// A=eye(m,n)
void dmat_set_eye(int m, int n, double *A, int LDA)
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
void dmat_set_rand(int m, int n, double *A, int LDA, double a, double b)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      MAT(A,i,j,LDA)=genrand_real3()*a+b;
    }
  }
}

//////////////////////////////////////////////////////////////////////////

// A(:,k) <-> A(:,l)
void dmat_swap_columns(int m, int n, double *A, int LDA, int k, int l)
{
  if(k<n && l<n && k!=l){
    dvec_swap(m,&COL(A,k,LDA),&COL(A,l,LDA));
  }
}

// A(k,:) <-> A(l,:)
void dmat_swap_rows(int m, int n, double *A, int LDA, int k, int l)
{
  int j;
  double a;
  if(k<m && l<m){
    for(j=0; j<n; j++){
      a=MAT(A,k,j,LDA);
      MAT(A,k,j,LDA)=MAT(A,l,j,LDA);
      MAT(A,l,j,LDA)=a;
    }
  }
}

//////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////

// B=A
void dmat_copy(int m, int n, double *B, int LDB, double *A, int LDA)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      MAT(B,i,j,LDB)=MAT(A,i,j,LDA);
    }
  }
}

// B=A'
void dmat_copy_t(int m, int n, double *B, int LDB, double *A, int LDA)
{
  int i,j;
  for(i=0; i<m; i++){
    for(j=0; j<n; j++){
      MAT(B,j,i,LDB)=MAT(A,i,j,LDA);
    }
  }
}

//////////////////////////////////////////////////////////////////////////

// B(:,j)=A(:,I(j)) for 0<=j<n
void dmat_copy_col_index(int m, int n, double *B, int LDB, double *A, int LDA, int *I)
{
  int j;
  for(j=0; j<n; j++){
    dvec_copy(m,&COL(B,j,LDB),&COL(A,I[j],LDA));
  }
}

//////////////////////////////////////////////////////////////////////////

// A(:,j) <-> A(:,I[j]) for 0<=j<n
void dmat_swap_index(int m, int n, double *A, int LDA, int *I)
{
  int LDB;
  double *B=NULL;
  LDB=m; B=dmat_allocate(LDB,n);
  dmat_copy_col_index(m,n,B,LDB,A,LDA,I);
  dmat_copy(m,n,A,LDA,B,LDB);
  B=dmat_free(B);
}

//////////////////////////////////////////////////////////////////////////


// A(k,k+offset)=a for k=0,1,..

void dmat_diag_set_scalar(int m, int n, double *A, int LDA, double a, int offset)
{
  int i,j;
  for(i=0; i<m; i++){
    j=i+offset;
    if(i>=0 && j>=0 && i<m && j<n){
      MAT(A,i,j,LDA)=a;
    }
  }
}

void dmat_diag_add_scalar(int n, double *A, int LDA, double a)
{
  int k;
  for(k=0; k<n; k++){
    MAT(A,k,k,LDA)+=a;
  }
}

// A(k,k)-=a for k=0,1,..,n-1
void dmat_diag_sub_scalar(int n, double *A, int LDA, double a)
{
  int k;
  for(k=0; k<n; k++){
    MAT(A,k,k,LDA)-=a;
  }
}

// B=B+A
void dmat_add(int m, int n, double *B, int LDB, double *A, int LDA)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      MAT(B,i,j,LDB)+=MAT(A,i,j,LDA);
    }
  }
}

// B=B-A
void dmat_sub(int m, int n, double *B, int LDB, double *A, int LDA)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      MAT(B,i,j,LDB)-=MAT(A,i,j,LDA);
    }
  }
}

// C=A*B
void dmat_prod(int l, int m, int n, double *C, int LDC, double *A, int LDA, double *B, int LDB)
{
  int i,j,k;
  for(i=0; i<l; i++){
    for(j=0; j<n; j++){
      MAT(C,i,j,LDC)=0;
      for(k=0; k<m; k++){
	MAT(C,i,j,LDC)+=MAT(A,i,k,LDA)*MAT(B,k,j,LDB);
      }
    }
  }
}

// C=C+A*B
void dmat_add_prod(int l, int m, int n, double *C, int LDC, double *A, int LDA, double *B, int LDB)
{
  int i,j,k;
  for(i=0; i<l; i++){
    for(j=0; j<n; j++){
      for(k=0; k<m; k++){
	MAT(C,i,j,LDC)+=MAT(A,i,k,LDA)*MAT(B,k,j,LDB);
      }
    }
  }
}

// C=C-A*B
void dmat_sub_prod(int l, int m, int n, double *C, int LDC, double *A, int LDA, double *B, int LDB)
{
  int i,j,k;
  for(i=0; i<l; i++){
    for(j=0; j<n; j++){
      for(k=0; k<m; k++){
	MAT(C,i,j,LDC)-=MAT(A,i,k,LDA)*MAT(B,k,j,LDB);
      }
    }
  }
}

// A=A+a*x*y'
void dmat_rank1op(int m, int n, double *A, int LDA, double a, double *x, double *y)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      MAT(A,i,j,LDA)+=a*x[i]*y[j];
    }
  }
}

///////////////////////////////////////////////////////////////////

void dmat_normalize(int m, int n, double *A, int LDA)
{
  int j;
  for(j=0; j<n; j++){
    dvec_normalize(m,&COL(A,j,LDA),&COL(A,j,LDA));
  }
}

void dmat_normalize_sgn(int m, int n, double *A, int LDA)
{
  int j;
  for(j=0; j<n; j++){
    dvec_normalize_sgn(m,&COL(A,j,LDA),&COL(A,j,LDA));
  }
}

//////////////////////////////////////

// x(j)=max(abs(A(:,j)-B(:,j))) for j=1,2,..,n
void dmat_sub_norm_max(double *x, int m, int n, double *A, int LDA, double *B, int LDB)
{
  int j;
  for(j=0; j<n; j++){
    x[j]=dvec_dist_norm_max(m,&COL(A,j,LDA),&COL(B,j,LDB));
  }
}

//////////////////////////////////////

// y(j)=max(abs(A(:,j)-x)) for j=1,2,..,n
void dmat_dist_norm_max(double *y, int m, int n, double *A, int LDA, double *x)
{
  int j;
  for(j=0; j<n; j++){
    y[j]=dvec_dist_norm_max(m,&COL(A,j,LDA),x);
  }
}

//////////////////////////////////////

// y=char(x)
void dmat_get_s(int m, int n, char **B, int LDB, double *A, int LDA, char format, int digits)
{
  char f[1024];
  int i,j;
  sprintf(f,"%%-.%d%c",digits,format);
  if(A==NULL){ return; }
  for(i=0; i<m; i++){
    for(j=0; j<n; j++){
      MAT(B,i,j,LDB)=char_renew_sprintf(MAT(B,i,j,LDB),NULL,f,MAT(A,i,j,LDA));
    }
  }
}

//////////////////////////////////////

void dmat_print(int m, int n, double *A, int LDA, char *name, char format, int digits)
{
  char **s=NULL;
  if(A==NULL){
    if(name!=NULL){ printf("%s=NULL\n",name); }
    else          { printf("NULL\n"); }
    return;
  }
  s=svec_allocate(m*n);
  dmat_get_s(m,n,s,m,A,LDA,format,digits);
  smat_print(m,n,s,m,name);
  s=svec_free(m*n,s);
}

/////////////////////////////////////////////////////////////////////////////////////

void dmat_save(int m, int n, double *A, int LDA, char* fmt, ...)
{
  int i,j;
  char fname[FILE_NAME_LENGTH_MAX];
  va_list argp;
  FILE *fid;
  va_start(argp,fmt); vsprintf(fname,fmt,argp); // file name
  if((fid=fopen(fname,"w"))==0){ ERROR_AT; printf("Cant't open file: %s\n",fname); exit(0); } // open
  // write
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      fprintf(fid,"%d\t%d\t%+.20e\n",i,j,MAT(A,i,j,LDA));
    }
  }
  // close
  fclose(fid);
}

//////////////////////////////////////////////////////////////////////////////////

void dmat_bin_save(int m, int n, double *A, int LDA, char* fmt, ...)
{
  int j;
  char fname[FILE_NAME_LENGTH_MAX];
  va_list argp;
  FILE *fid;
  va_start(argp,fmt); vsprintf(fname,fmt,argp); // file name
  if((fid=fopen(fname,"w"))==0){ ERROR_AT; printf("Cant't open file: %s\n",fname); exit(0); } // open
  fwrite("dmat",sizeof(char),strlen("dmat"),fid); // write header
  fwrite(&m,sizeof(int),1,fid); // write row-size
  fwrite(&n,sizeof(int),1,fid); // write col-size
  for(j=0; j<n; j++){ fwrite(&COL(A,j,LDA),sizeof(double),m,fid); } // write data
  fclose(fid); // close
}

double *dmat_bin_load(int *m, int *n, char* fmt, ...)
{
  double *A=NULL;
  int l,j;
  size_t k;
  char fname[FILE_NAME_LENGTH_MAX],*buf=NULL;
  va_list argp;
  FILE *fid;
  // file name
  va_start(argp,fmt); vsprintf(fname,fmt,argp);
  // open
  if((fid=fopen(fname,"r"))==0){ fclose(fid); A=NULL; (*m)=0; (*n)=0; }
  else{
    // read header
    l=strlen("dmat");
    buf=malloc(sizeof(char)*(l+1));
    k=fread(buf,sizeof(char),l,fid);
    if(k!=(size_t)l || strncmp(buf,"dmat",l)!=0){ A=NULL; (*m)=0; (*n)=0; }
    else{
      // read size
      k=fread(m,sizeof(int),1,fid); if(k!=1 || (*m)<=0){ ERROR_AT; printf("Failed to load the row-size from the file '%s'.\n",fname); exit(0); }
      k=fread(n,sizeof(int),1,fid); if(k!=1 || (*n)<=0){ ERROR_AT; printf("Failed to load the col-size from the file '%s'.\n",fname); exit(0); }
      // allocate
      A=dmat_allocate((*m),(*n));
      // read data
      for(j=0; j<(*n); j++){
	k=fread(&COL(A,j,(*m)),sizeof(double),(*m),fid);
	if(k!=(size_t)(*m)){ ERROR_AT; printf("Failed to load the data from the file '%s'.\n",fname); exit(0); }
      }
    }
    fclose(fid); 
  }
  // done
  free(buf);
  return A;
}

//EOF
