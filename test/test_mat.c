#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#define MAT(A,I,J,LDA) ((A)[(I)+(J)*(LDA)])

/**
 @brief allocation
 */
double *dmat_allocate(int m, int n)
{
  return (double*)malloc(sizeof(double)*m*n);
}

/**
 @brief free
 */
double *dmat_free(double *A)
{
  if(A==NULL) return NULL;
  free(A);
  return A=NULL;
}

/**
 @brief output
 */
void dmat_print(int m, int n, double *A, int LDA)
{
  int i,j;
  if(A==NULL){ printf("NULL\n"); return; }
  for(i=0; i<m; i++){
    for(j=0; j<n; j++){ printf("%+.1e ",MAT(A,i,j,LDA)); }
    printf("\n");
  }
}

/**
 @brief A=ones(m,n)*value
 */
void dmat_set_all_d(int m, int n, double *A, int LDA, double value)
{
  int i,j;
  if(A==NULL){ return; }
  for(i=0; i<m; i++){
    for(j=0; j<n; j++){ MAT(A,i,j,LDA)=value; }
  }
}

/**
 @brief C=A+B
 */
void dmat_add_dmat_dmat(int m, int n, double *C, int LDC, double *A, int LDA, double *B, int LDB)
{
  int i,j;
  if(A==NULL || B==NULL){ return; }
  for(i=0; i<m; i++){
    for(j=0; j<n; j++){ MAT(C,i,j,LDC)=MAT(A,i,j,LDA)+MAT(B,i,j,LDB); }
  }
}

int main(void)
{
  int m,n,LDA,LDB;
  double *A=NULL,*B=NULL;

  m=10; n=10; 
  LDA=m; A=dmat_allocate(LDA,n);
  LDB=m; B=dmat_allocate(LDB,n);

  dmat_set_all_d(m,n,A,LDA,0);  
  printf("A=zeros(10,10)=\n");
  dmat_print(m,n,A,LDA);

  dmat_set_all_d(5,5,A,LDA,-1);
  printf("A(1:5,1:5)=-1*ones(5,5);\n");
  dmat_set_all_d(5,5,&MAT(A,5,5,LDA),LDA,-2);
  printf("A(6:10,6:10)=-2*ones(5,5)\n");
  printf("A=\n");
  dmat_print(m,n,A,LDA);

  dmat_set_all_d(m,n,B,LDB,1);    
  printf("B=ones(10,10)=\n");
  dmat_print(m,n,B,LDA);
  
  dmat_add_dmat_dmat(m,n,B,LDB,A,LDA,B,LDB);
  printf("B=B+A=\n");
  dmat_print(m,n,B,LDB);

  dmat_add_dmat_dmat(5,5,&MAT(B,0,0,LDB),LDB,&MAT(B,0,0,LDB),LDB,&MAT(B,5,5,LDB),LDB);
  printf("B(1:5,1:5)=B(1:5,1:5)+B(6:10,6:10);\n");
  printf("B=\n");
  dmat_print(m,n,B,LDB);

  A=dmat_free(A); 
  B=dmat_free(B); 
  return 0;
}
