#include<stdio.h>
#include<stdlib.h>
#include"is_macros.h"
#include"is_ivec.h"
#include"is_strings.h"
#include"is_print.h"
#include"is_func.h"
#include"is_rmulti.h"
#include"is_cmulti.h"

#define FR(F)      (func_retain(F))
#define FGET_SP(f) (func_get_split_part(f))
#define FGET_NB(f) (func_get_number_part(f))

/////////////////////////////////////////////////////////////

void func_a_rm_op(func_t *f, func_is_t *fis)
{
  int i;
  if(f==NULL) return;
  for(i=0; i<func_asize(f); i++){ if(fis(func_aget(f,i))){ func_adel(f,i); } }
  func_a_rm_null(f);
}

void func_a_del_op(func_t *f, func_is_t *fis)
{
  int i;
  if(f==NULL) return;
  for(i=0; i<func_asize(f); i++){ if(fis(func_aget(f,i))){ func_adel(f,i); } }
}

/////////////////////////////////////////////////////////////

void func_a_rm_not_op(func_t *f, func_is_t *fis)
{
  int i;
  if(f==NULL) return;
  for(i=0; i<func_asize(f); i++){ if(!fis(func_aget(f,i))){ func_adel(f,i); } }
  func_a_rm_null(f);
}

void func_a_del_not_op(func_t *f, func_is_t *fis)
{
  int i;
  if(f==NULL) return;
  for(i=0; i<func_asize(f); i++){ if(!fis(func_aget(f,i))){ func_adel(f,i); } }
}

void func_a_del_not_op_op(func_t *f, func_is_t *fis1, func_is_t *fis2)
{
  int i;
  if(f==NULL) return;
  for(i=0; i<func_asize(f); i++){
    if(!fis1(func_aget(f,i)) && !fis2(func_aget(f,i))){ func_adel(f,i); }
  }
}

/////////////////////////////////////////////////////////////

int func_a_has_op_pow1(func_t *f, func_is_t *fis)
{
  int value=0,i;
  for(i=0; i<func_asize(f); i++){
    if(fis(func_aget(f,i)) && func_has_power(func_aget(f,i)) && func_power(func_aget(f,i))==1){ value=1; }
  }
  return value;
}

/////////////////////////////////////////////////////////////

int func_a_has_op(func_t *f, func_is_t *fis)
{
  int i;
  for(i=0; i<func_asize(f); i++){
    if(fis(func_aget(f,i))){ return 1; }
  }
  return 0;
}

/////////////////////////////////////////////////////////////

void func_a_append(func_t *f, func_t *g)
{
  func_a_resize(f,func_asize(f)+1);  
  func_aset(f,func_asize(f)-1,g);
}

void func_a_rm_end(func_t *f)
{
  func_a_resize(f,func_asize(f)-1);
}

void func_a_insert(func_t *f, int n, func_t *g)
{
  int i,k;
  func_t **a;
  if(f==NULL){ FUNC_ERROR_ARG1("func_a_insert",f); }
  a=(func_t**)malloc(sizeof(func_t*)*(func_asize(f)+1));
  for(k=0,i=0; k<func_asize(f)+1; ){
    if(i==n){ a[k++]=g; a[k++]=FR(f->a[i++]); }
    else    { a[k++]=FR(f->a[i++]); }
  }
  func_a_replace(f,func_asize(f)+1,a);
}

void func_a_rm(func_t *f, int i)
{
  if(f==NULL || i<0 || i>=func_asize(f)){ FUNC_ERROR_ARG1("func_a_rm",f); }
  f->a[i]=func_del(f->a[i]);
  func_a_rm_null(f);
}

void func_a_rm_null(func_t *f)
{
  int i,k,n;
  func_t **a;
  if(f==NULL){ FUNC_ERROR_ARG1("func_a_rm_null",f); }
  n=func_args_count_non_null(f);
  if(n<=0){ n=0; a=NULL; }
  else{
    a=(func_t**)malloc(sizeof(func_t*)*n);
    for(k=0,i=0; i<func_asize(f); i++){
      if((f->a[i])!=NULL){ a[k++]=FR(f->a[i]); }
    }
  }
  func_a_replace(f,n,a);
}

/////////////////////////////////////////////////////////////

func_t *func_aget(func_t *f, int i)
{
  if(f==NULL || f->a==NULL){ FUNC_ERROR_ARG1("func_aget",f); }
  if(i<0 || i>=func_asize(f)){ printf("Error in func_aget(%s,%d)\n",func_op(f),i); exit(0); }
  return f->a[i];
}

void func_aset(func_t *f, int i, func_t *g)
{
  if(f==NULL){ FUNC_ERROR_ARG2("func_aset",f,g); }
  if(i<0 || i>=func_asize(f)){ printf("Error in func_aset(%s,%d,%s)\n",func_op(f),i,func_op(g)); exit(0); }
  f->a[i]=func_del(f->a[i]);
  f->a[i]=g;
}

void func_adel(func_t *f, int i)
{
  if(f==NULL){ FUNC_ERROR_ARG1("func_adel",f); }
  if(i<0 || i>=func_asize(f)){ printf("Error in func_adel(%s,%d)\n",func_op(f),i); exit(0); }
  f->a[i]=func_del(f->a[i]);
}



/////////////////////////////////////////////////////////////

void func_a_resize(func_t *f, int n)
{
  int i;
  func_t **arg=NULL;
  if(f==NULL || n<func_find_amin(f) || (func_find_amax(f)!=FUNC_NOLIMIT && n>func_find_amax(f))){ FUNC_ERROR_ARG1("func_a_resize",f); }
  arg=(func_t**)malloc(sizeof(func_t*)*n);
  for(i=0; i<n; i++){
    if(i<func_asize(f)){ arg[i]=FR(f->a[i]); }
    else               { arg[i]=NULL; }
  }
  func_a_replace(f,n,arg);
}

/////////////////////////////////////////////////////////////

void func_a_swap(func_t *f, func_t *g)
{
  int n;
  func_t **a;
  if(f==NULL || g==NULL || func_asize(g)<func_find_amin(f) || func_asize(f)<func_find_amin(g) || (func_find_amax(f)!=FUNC_NOLIMIT && func_asize(g)>func_find_amax(f)) || (func_find_amax(g)!=FUNC_NOLIMIT && func_asize(f)>func_find_amax(g))){ FUNC_ERROR_ARG2("func_a_swap",f,g); }
  n=f->n;  f->n=g->n;  g->n=n;
  a=f->a;  f->a=g->a;  g->a=a;
}

void func_a_replace(func_t *f, int n, func_t **arg)
{
  func_a_del(f);
  if(f->n!=0 || f->a!=NULL){ FUNC_ERROR_ARG1("func_a_replace",f); }
  f->n=n;
  f->a=arg;
}

/////////////////////////////////////////////////////////////

int func_asize(func_t *f)
{
  if(f==NULL){ return 0; }
  return f->n;
}

/////////////////////////////////////////////////////////////

void func_a_new(func_t *f, int n)
{
  int i;
  if(f==NULL) return;
  if(n<func_find_amin(f) || (func_find_amax(f)!=FUNC_NOLIMIT && n>func_find_amax(f))){ FUNC_ERROR_ARG1("func_a_new",f); }
  if(f->a!=NULL){ func_a_del(f); }
  if(n>0){
    f->n=n;
    f->a=(func_t**)malloc(sizeof(func_t*)*n);
    for(i=0; i<n; i++){ f->a[i]=NULL; }
  }else{
    f->n=0;
    f->a=NULL;
  }
}

void func_a_del(func_t *f)
{
  int i;
  if(f==NULL) return;
  if(f->a==NULL) return;
  for(i=0; i<f->n; i++){ f->a[i]=func_del(f->a[i]); }
  free(f->a);
  f->a=NULL;
  f->n=0;
}

void func_a_clone(func_t *f, func_t *g)
{
  int n,i;
  if(f==NULL || g==NULL) return;
  if(f->a!=NULL){ func_a_del(f); }
  n=func_asize(g);
  if(n>0){
    f->n=n;
    f->a=(func_t**)malloc(sizeof(func_t*)*n);
    for(i=0; i<n; i++){
      f->a[i]=FR(g->a[i]);
    }
  }
}

/////////////////////////////////////////////////////////////


///private//////////////////////////////////////////////////


// g in f ?
int func_args_have(func_t *f, func_t *g)
{
  int i;
  for(i=0; i<func_asize(f); i++){
    if(func_cmp(f->a[i],g)==0){ return 1; }
  }
  return 0;
}

/////private////////////////////////////////////////

int func_args_count_op(func_t *f, char *op)
{
  int k,i;
  if(func_find_amax(f)<0){
    for(k=0, i=0; i<func_asize(f); i++){
      if(func_has_power(f)){ if(func_is(f->a[i],op) && func_power(f->a[i])==1){ k+=func_asize(f->a[i]); } else{ k++; } }
      else                 { if(func_is(f->a[i],op))      { k+=func_asize(f->a[i]); } else{ k++; } }
    }
  }else{ FUNC_ERROR_ARG1("func_args_count_op",f); }
  return k;
}

/////private////////////////////////////////////////

void func_args_arrange(func_t *f, int *I)
{
  int i,n;
  func_t **arg=NULL;
  n=func_asize(f);
  arg=(func_t**)malloc(sizeof(func_t*)*n);
  for(i=0; i<n; i++){ arg[i]=FR(f->a[I[i]]); }
  func_a_replace(f,n,arg);
}

void func_args_rm_op(func_t *f, char *op)
{
  int i;
  if(f==NULL) return;
  for(i=0; i<func_asize(f); i++){
    if(func_is(f->a[i],op)){
      f->a[i]=func_del(f->a[i]);
    }
  }
  func_a_rm_null(f);
}

void func_args_rm_add_pow1(func_t *f)
{
  int i;
  if(f==NULL) return;
  for(i=0; i<func_asize(f); i++){
    if(func_is_add(f->a[i]) && func_power(f->a[i])==1){
      f->a[i]=func_del(f->a[i]);
    }
  }
}

int func_args_count_non_null(func_t *f)
{
  int i,count=0;
  if(f==NULL) return 0;
  for(i=0; i<func_asize(f); i++){
    if((f->a[i])!=NULL) count++;
  }
  return count;
}

///private//////////////////////////////////////////////////

void func_args_swap(func_t *f, int i, int j)
{
  func_t *foo;
  foo=f->a[i];
  f->a[i]=f->a[j];
  f->a[j]=foo;
}

///private//////////////////////////////////////////////////


void func_args_quick_sort(func_t *f, int *I, int left, int right)
{
  int i,last;
  if(f==NULL) return;
  if(left>=right) return;
  func_args_swap(f,left,(left+right)/2);
  if(I!=NULL) ivec_swap_at(I,left,(left+right)/2);
  last=left;
  for(i=left+1; i<=right; i++){
    if(func_lt(f->a[i],f->a[left])){
      ++last;
      func_args_swap(f,last,i);
      if(I!=NULL) ivec_swap_at(I,last,i);
    }
  }
  func_args_swap(f,left,last);
  if(I!=NULL) ivec_swap_at(I,left,last);
  func_args_quick_sort(f,I,left,last-1);
  func_args_quick_sort(f,I,last+1,right);
}

void func_args_sort(func_t *f)
{
  if(f==NULL) return;
  func_args_quick_sort(f,NULL,0,(func_asize(f))-1);
}

void func_args_sort_index(func_t *f, int *I)
{
  if(f==NULL) return;
  ivec_set_grid(func_asize(f),I);
  func_args_quick_sort(f,I,0,(func_asize(f))-1);
}

/////////////////////////////////////////////////////////////

void func_args_reverse(func_t *f)
{
  int i,n;
  func_t *g;
  if(f==NULL) return;
  n=func_asize(f);
  for(i=0; i<n/2; i++){
    g=f->a[i];
    f->a[i]=f->a[n-i-1];
    f->a[n-i-1]=g;
  }
}


//EOF
