#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<stdarg.h>

#include"is_macros.h"
#include"is_strings.h"
#include"is_svec.h"
#include"is_ivec.h"
#include"is_zvec.h"
#include"is_zmat.h"
#include"is_dmat.h"
#include"mt19937ar.h"

//////////////////////////////////////////

#define FILE_NAME_LENGTH_MAX 100000


//////////////////////////////////////////

dcomplex* zmat_allocate(int m, int n)
{
  return (dcomplex*)malloc(sizeof(dcomplex)*m*n);
}

dcomplex* zmat_free(dcomplex *A)
{
  if(A==NULL) return NULL;
  free(A);
  return A=NULL;
}

// A=zeros(m,n)
void zmat_set_zeros(int m, int n, dcomplex *A, int LDA)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      Z_SET(MAT(A,i,j,LDA),0,0);
    }
  }
}

// A=ones(m,n)
void zmat_set_ones(int m, int n, dcomplex *A, int LDA)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      Z_SET(MAT(A,i,j,LDA),1,0);
    }
  }
}

// A=ones(m,n)*a
void zmat_set_all_z(int m, int n, dcomplex *A, int LDA, dcomplex a)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      MAT(A,i,j,LDA)=a;
    }
  }
}

// A=ones(m,n)*a
void zmat_set_all_d(int m, int n, dcomplex *A, int LDA, double a)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      Z_SET(MAT(A,i,j,LDA),a,0);
    }
  }
}

// A=eye(m,n); % identiy matrix
void zmat_set_eye(int m, int n, dcomplex *A, int LDA)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      if(i==j){
	Z_SET(MAT(A,i,j,LDA),1,0);
      }
      else{
	Z_SET(MAT(A,i,j,LDA),0,0);
      }
    }
  }
}

// A=rand(m,n)*a+b
void zmat_set_rand(int m, int n, dcomplex *A, int LDA, double a, double b)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      Z_SET(MAT(A,i,j,LDA),genrand_real3()*a+b,genrand_real3()*a+b);
    }
  }
}

////////////////////////////////

/**
 @brief 実部 B=real(A)
 */
void zmat_real(int m, int n, double *B, int LDB, dcomplex *A, int LDA)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      MAT(B,i,j,LDB)=Z_R(MAT(A,i,j,LDA));
    }
  }
}


////////////////////////////////////////////////////////////////

// A(:,k) <-> A(:,l)
void zmat_swap_columns(int m, int n, dcomplex *A, int LDA, int k, int l)
{
  if(k<n && l<n && k!=l){
    zvec_swap(m,&COL(A,k,LDA),&COL(A,l,LDA));
  }
}

// A(k,:) <-> A(l,:)
void zmat_swap_rows(int m, int n, dcomplex *A, int LDA, int k, int l)
{
  int j;
  dcomplex a;
  if(k<m && l<m){
    for(j=0; j<n; j++){
      a=MAT(A,k,j,LDA);
      MAT(A,k,j,LDA)=MAT(A,l,j,LDA);
      MAT(A,l,j,LDA)=a;
    }
  }
}

////////////////////////////////////////////////////////////////

// B=A
void zmat_copy(int m, int n, dcomplex *B, int LDB, dcomplex *A, int LDA)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      MAT(B,i,j,LDB)=MAT(A,i,j,LDA);
    }
  }
}

// B=A
void zmat_copy_d(int m, int n, dcomplex *B, int LDB, double *A, int LDA)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      Z_SET(MAT(B,i,j,LDB),MAT(A,i,j,LDA),0);
    }
  }
}

// B=A^T
void zmat_copy_t(int m, int n, dcomplex *B, int LDB, dcomplex *A, int LDA)
{
  int i,j;
  for(i=0; i<m; i++){
    for(j=0; j<n; j++){
      MAT(B,j,i,LDB)=MAT(A,i,j,LDA);
    }
  }
}

// B=A'
void zmat_copy_ct(int m, int n, dcomplex *B, int LDB, dcomplex *A, int LDA)
{
  int i,j;
  for(i=0; i<m; i++){
    for(j=0; j<n; j++){
      Z_SET_CONJ(MAT(B,j,i,LDB),MAT(A,i,j,LDA));
    }
  }
}

////////////////////////////////////////////////////////////////

// B(:,j)=A(:,I(j)) for 0<=j<n
void zmat_copy_col_index(int m, int n, dcomplex *B, int LDB, dcomplex *A, int LDA, int *I)
{
  int j;
  for(j=0; j<n; j++){
    zvec_copy(m,&COL(B,j,LDB),&COL(A,I[j],LDA));
  }
}

////////////////////////////////////////////////////////////////

// A(:,j) <-> A(:,I[j]) for 0<=j<n
void zmat_swap_index(int m, int n, dcomplex *A, int LDA, int *I)
{
  int LDB;
  dcomplex *B=NULL;
  LDB=m; B=zmat_allocate(LDB,n);
  zmat_copy_col_index(m,n,B,LDB,A,LDA,I);
  zmat_copy(m,n,A,LDA,B,LDB);
  B=zmat_free(B);
}


////////////////////////////////////////////////////////////////


// A(k,k+offset)=a for k=0,1,..
void zmat_diag_set_scalar(int m, int n, dcomplex *A, int LDA, dcomplex a, int offset)
{
  int i,j;
  for(i=0; i<m; i++){
    j=i+offset;
    if(i>=0 && j>=0 && i<m && j<n){
      MAT(A,i,j,LDA)=a;
    }
  }
}

void zmat_diag_set_scalar_d(int m, int n, dcomplex *A, int LDA, double a, int offset)
{
  int i,j;
  for(i=0; i<m; i++){
    j=i+offset;
    if(i>=0 && j>=0 && i<m && j<n){
      Z_SET(MAT(A,i,j,LDA),a,0);
    }
  }
}

// A(k,k)+=a for k=0,1,..,n-1
void zmat_diag_add_scalar(int n, dcomplex *A, int LDA, dcomplex a)
{
  int k;
  for(k=0; k<n; k++){
    Z_ADD(MAT(A,k,k,LDA),a);
  }
}

// A(k,k)-=a for k=0,1,..,n-1
void zmat_diag_sub_scalar(int n, dcomplex *A, int LDA, dcomplex a)
{
  int k;
  for(k=0; k<n; k++){
    Z_SUB(MAT(A,k,k,LDA),a);
  }
}

// B=B+A
void zmat_add(int m, int n, dcomplex *B, int LDB, dcomplex *A, int LDA)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      Z_ADD(MAT(B,i,j,LDB),MAT(A,i,j,LDA));
    }
  }
}

// B=B-A
void zmat_sub(int m, int n, dcomplex *B, int LDB, dcomplex *A, int LDA)
{
  int i,j;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      Z_SUB(MAT(B,i,j,LDB),MAT(A,i,j,LDA));
    }
  }
}

// C=A*B
void zmat_prod(int l, int m, int n, dcomplex *C, int LDC, dcomplex *A, int LDA, dcomplex *B, int LDB)
{
  int i,j,k;
  for(i=0; i<l; i++){
    for(j=0; j<n; j++){
      Z_SET(MAT(C,i,j,LDC),0,0);
      for(k=0; k<m; k++){
	Z_ADD_TIMES(MAT(C,i,j,LDC),MAT(A,i,k,LDA),MAT(B,k,j,LDB));
      }
    }
  }
}

// C=C+A*B
void zmat_add_prod(int l, int m, int n, dcomplex *C, int LDC, dcomplex *A, int LDA, dcomplex *B, int LDB)
{
  int i,j,k;
  for(i=0; i<l; i++){
    for(j=0; j<n; j++){
      for(k=0; k<m; k++){
	Z_ADD_TIMES(MAT(C,i,j,LDC),MAT(A,i,k,LDA),MAT(B,k,j,LDB));
      }
    }
  }
}

// C=C-A*B
void zmat_sub_prod(int l, int m, int n, dcomplex *C, int LDC, dcomplex *A, int LDA, dcomplex *B, int LDB)
{
  int i,j,k;
  for(i=0; i<l; i++){
    for(j=0; j<n; j++){
      for(k=0; k<m; k++){
	Z_SUB_TIMES(MAT(C,i,j,LDC),MAT(A,i,k,LDA),MAT(B,k,j,LDB));
      }
    }
  }
}

// A=A+a*x*y'
void zmat_rank1op(int m, int n, dcomplex *A, int LDA, dcomplex a, dcomplex *x, dcomplex *y)
{
  int i,j;
  dcomplex value;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      Z_SET_DOT(value,y[j],x[i]);
      Z_ADD_TIMES(MAT(A,i,j,LDA),a,value);
    }
  }
}

///////////////////////////////////////////////////////////////////

void zmat_normalize(int m, int n, dcomplex *A, int LDA)
{
  int j;
  for(j=0; j<n; j++){
    zvec_normalize(m,&COL(A,j,LDA));
  }
}

void zmat_normalize_sgn(int m, int n, dcomplex *A, int LDA)
{
  int j;
  for(j=0; j<n; j++){
    zvec_normalize_sgn(m,&COL(A,j,LDA));
  }
}


/////////////////////////////////////////

// x(j)=max(abs(A(:,j)-B(:,j))) for j=1,2,..,n
void zmat_sub_norm_max(double *x, int m, int n, dcomplex *A, int LDA, dcomplex *B, int LDB)
{
  int j;
  for(j=0; j<n; j++){
    x[j]=zvec_dist_norm_max(m,&COL(A,j,LDA),&COL(B,j,LDB));
  }
}

/////////////////////////////////////////


// y(j)=max(abs(A(:,j)-x)) for j=1,2,..,n
void zmat_dist_norm_max(double *y, int m, int n, dcomplex *A, int LDA, dcomplex *x)
{
  int j;
  for(j=0; j<n; j++){
    y[j]=zvec_dist_norm_max(m,&COL(A,j,LDA),x);
  }
}

/////////////////////////////////////////

// y=char(x)
void zmat_get_s(int m, int n, char **B, int LDB, dcomplex *A, int LDA, char format, int digits)
{
  char f[1024];
  int i,j;
  sprintf(f,"%%-.%d%c%%+.%d%ci",digits,format,digits,format);
  if(A==NULL){ return; }
  for(i=0; i<m; i++){
    for(j=0; j<n; j++){
      MAT(B,i,j,LDB)=char_renew_sprintf(MAT(B,i,j,LDB),NULL,f,Z_R(MAT(A,i,j,LDA)),Z_I(MAT(A,i,j,LDA)));
    }
  }
}

/////////////////////////////////////////

void zmat_print(int m, int n, dcomplex *A, int LDA, char *name, char format, int digits)
{
  char **s=NULL;
  if(A==NULL){
    if(name!=NULL){ printf("%s=NULL\n",name); }
    else          { printf("NULL\n"); }
    return;
  }
  s=svec_allocate(m*n);
  zmat_get_s(m,n,s,m,A,LDA,format,digits);
  smat_print(m,n,s,m,name);
  s=svec_free(m*n,s);
}

/////////////////////////////////////////////////////////////////////////

void zmat_save(int m, int n, dcomplex *A, int LDA, char* fmt, ...)
{
  int i,j;
  char fname[FILE_NAME_LENGTH_MAX];
  va_list argp;
  FILE *fid;
   // file name
  va_start(argp,fmt);
  vsprintf(fname,fmt,argp);
  // open file
  if((fid=fopen(fname,"w"))==0){
    fprintf(stderr,"Cant't open file: %s\n",fname);
    exit(0);
  }
  // loop
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      fprintf(fid,"%d\t%d\t%+.20e\t%+.20e\n",i,j,MAT(A,i,j,LDA).r,MAT(A,i,j,LDA).i);
    }
  }
  // close
  fclose(fid);
  fprintf(stderr,"saved: %s\n",fname);
}

/////////////////////////////////////////////////////////////////////////

void zmat_bin_save(int m, int n, dcomplex *A, int LDA, char* fmt, ...)
{
  int j;
  char fname[FILE_NAME_LENGTH_MAX];
  va_list argp;
  FILE *fid;
  va_start(argp,fmt); vsprintf(fname,fmt,argp); // file name
  if((fid=fopen(fname,"w"))==0){ ERROR_AT; printf("Can't open file: %s\n",fname); exit(0); } // open
  fwrite("zmat",sizeof(char),strlen("zmat"),fid); // write header
  fwrite(&m,sizeof(int),1,fid); // write row-size
  fwrite(&n,sizeof(int),1,fid); // write col-size
  for(j=0; j<n; j++){ fwrite(&COL(A,j,LDA),sizeof(dcomplex),m,fid); } // write data
  fclose(fid); // close
}

dcomplex *zmat_bin_load(int *m, int *n, char* fmt, ...)
{
  dcomplex *zA=NULL;
  double *dA=NULL;
  int l,j;
  size_t k;
  char fname[FILE_NAME_LENGTH_MAX],*buf=NULL;
  va_list argp;
  FILE *fid;
   // file name
  va_start(argp,fmt); vsprintf(fname,fmt,argp);
  // open file
  if((fid=fopen(fname,"r"))==0){ fclose(fid); zA=NULL; (*m)=0; (*n)=0; }
  else{
    // read header
    l=strlen("zmat");
    buf=malloc(sizeof(char)*(l+1));
    k=fread(buf,sizeof(char),l,fid);
    if(k==(size_t)l && strncmp(buf,"zmat",l)==0){ /* zmat */
      // read size
      k=fread(m,sizeof(int),1,fid); if(k!=1 || (*m)<=0){ ERROR_AT; printf("Failed to load the row-size from the file '%s'.\n",fname); exit(0); }
      k=fread(n,sizeof(int),1,fid); if(k!=1 || (*n)<=0){ ERROR_AT; printf("Failed to load the col-size from the file '%s'.\n",fname); exit(0); }
      // allocate
      zA=zmat_allocate((*m),(*n));
      // read data
      for(j=0; j<(*n); j++){    
	k=fread(&COL(zA,j,(*m)),sizeof(dcomplex),(*m),fid);
	if(k!=(size_t)(*m)){ ERROR_AT; printf("Failed to load the data from the file '%s'.\n",fname); exit(0); }
      }
      // close
      fclose(fid);
    }else if(k==(size_t)l && strncmp(buf,"dmat",l)==0){ /* dmat */
      fclose(fid);                            // close
      dA=dmat_bin_load(m,n,fname);            // load
      zA=zmat_allocate((*m),(*n));            // allocate
      zmat_copy_d((*m),(*n),zA,(*m),dA,(*m)); // copy
      dA=dmat_free(dA);                       // free
    }else{ fclose(fid); zA=NULL; (*m)=0; (*n)=0; }
  }
  // done
  free(buf);
  return zA;
}

//EOF
