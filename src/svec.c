#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<strings.h>
#include"is_macros.h"
#include"is_strings.h"
#include"is_svec.h"
#include"is_ivec.h"

char **svec_allocate(int n)
{
  int i;
  char **x=NULL;
  x=(char**)malloc(n*sizeof(char*));
  for(i=0; i<n; i++){ x[i]=NULL; }
  return x;
}

char **svec_free(int n, char **x)
{
  int i;
  if(x==NULL){ return NULL; }
  for(i=0; i<n; i++){ x[i]=char_del(x[i]); }
  free(x);
  x=NULL;
  return x;
}

// y=x
void svec_copy(int n, char **y, char **x)
{
  int i;
  NULL_EXC2(y,x);
  for(i=0; i<n; i++){ y[i]=char_new(x[i],NULL); }
}

////////////////////////////////////////////////////////////////////////////////

// y=x
void svec_set_si(int n, char **y, int *x)
{
  int i;
  NULL_EXC1(y);
  for(i=0; i<n; i++){ y[i]=char_renew_sprintf(y[i],NULL,"%d",x[i]); }
}

// y=ones(n,1)*x
void svec_set_all_d(int n, char **y, double x, char *format)
{
  int i;
  NULL_EXC1(y);
  for(i=0; i<n; i++){ y[i]=char_renew_sprintf(y[i],NULL,format,x); }
}

// y=ones(n,1)*x
void svec_set_all(int n, char **y, char *x)
{
  int i;
  NULL_EXC1(y);
  for(i=0; i<n; i++){ y[i]=char_renew(y[i],x,NULL); }
}

////////////////////////////////////////////////////////////////////////////////

// y=max(length(x))
int svec_max_length(int n, char **x)
{
  int value=0,a,i;
  for(i=0; i<n; i++){
    if(x[i]==NULL){ a=0; }
    else          { a=strlen(x[i]); }
    if(a>value){ value=a; }
  }
  return value;
}

// y=max(length(A(k,:))
int smat_max_length_row(int n, char **A, int LDA, int k)
{
  int value=0,a,j;
  for(j=0; j<n; j++){
    if(MAT(A,k,j,LDA)==NULL){ a=0; }
    else                    { a=strlen(MAT(A,k,j,LDA)); }
    if(a>value){ value=a; }
  }
  return value;
}

// y=max(length(A(:,k))
int smat_max_length_column(int m, char **A, int LDA, int k)
{
  int value=0,a,i;
  for(i=0; i<m; i++){
    if(MAT(A,i,k,LDA)==NULL){ a=0; }
    else                    { a=strlen(MAT(A,i,k,LDA)); }
    if(a>value){ value=a; }
  }
  return value;
}

////////////////////////////////////////////////////////////////////////////////

int svec_has(int n, char **x, char *match)
{
  int i,ret=1;
  if(x==NULL){ return 0; }
  for(i=0; ret==1 && i<n; i++){
    if(!str_has_any_char(x[i],match)){ ret=0; };
  }
  return ret;
}

////////////////////////////////////////////////////////////////////////////////

void svec_put(int n, char **x, char *sep)
{
  int i;
  for(i=0; i<n; i++){
    if(i>0){ if(sep==NULL){ printf(" "); }else{ printf("%s",sep); } }
    if(x[i]==NULL){ printf("NULL"); }
    else          { printf("%s",x[i]); }
  }
}

void svec_print(int n, char **x, const char *name)
{
  int i,l,k,a;
  if(x==NULL){
    if(name!=NULL){ printf("%s=NULL\n",name); }
    else          { printf("NULL\n"); }
    return;
  }  
  if(name!=NULL){ printf("%s=\n",name); }
  l=svec_max_length(n,x);
  for(i=0; i<n; i++){
    a=strlen(x[i]);
    for(k=0; k<l-a; k++){ printf(" "); }
    if(x[i]==NULL){ printf("NULL\n"); }
    else          { printf("%s\n",x[i]); }
  }
}

void smat_print(int m, int n, char **A, int LDA, const char *name)
{
  int i,j,a,k,*l=NULL;
  if(A==NULL){
    if(name!=NULL){ printf("%s=NULL\n",name); }
    else          { printf("NULL\n"); }
    return;
  }
  if(name!=NULL){ printf("%s=\n",name); }
  l=ivec_allocate(n);
  for(j=0; j<n; j++){ l[j]=smat_max_length_column(m,A,LDA,j); }
  for(i=0; i<m; i++){
    for(j=0; j<n; j++){
      if(j>0){ printf(" "); }
      a=strlen(MAT(A,i,j,LDA));
      for(k=0; k<l[j]-a; k++){ printf(" "); }
      if(MAT(A,i,j,LDA)==NULL){ printf("NULL"); }
      else                    { printf("%s",MAT(A,i,j,LDA)); }
    }
    printf("\n");
  }
  l=ivec_free(l);
}


//EOF
