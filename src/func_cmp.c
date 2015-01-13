#include<stdio.h>
#include<stdlib.h>
#include"is_macros.h"
#include"is_func.h"

#define FR(F) (func_retain(F))

/////////////////////////////////////////////////////////////////////////////////////

int __func_mono_order=FUNC_MONO_GRLEX;

void func_set_mono_order(int order)
{
  __func_mono_order=order;
}

int func_get_mono_order(void)
{
  return __func_mono_order;
}

/////////////////////////////////////////////////////////////////////////////////////

// f <=> g
int func_cmp(func_t *f, func_t *g)
{
  func_cmp_t *cmp1=NULL,*cmp2=NULL;
  func_t *a=NULL,*b=NULL;
  int i,j,value;
  // definite
  if     (f==NULL && g==NULL){ return  0; } // f=g
  else if(f==NULL)           { return -1; } // f<g
  else if(g==NULL)           { return +1; } // f>g
  else if(func_in_bigint(f)  && func_in_bigint(g)) { return func_bigint_cmp(f,g); }
  else if(func_in_real(f)    && func_in_real(g))   { return func_real_cmp(f,g); }
  else if(func_in_complex(f) && func_in_complex(g)){ return func_complex_cmp(f,g); }
  else if(func_is_number(f)  && func_is_number(g)) { return func_number_cmp(f,g); }
  else if(func_is_var(f)     && func_is_var(g))    { return func_var_cmp(f,g); }
  // split mul
  a=func_get_split_part(f);
  b=func_get_split_part(g);
  if(a!=NULL && b!=NULL && (a!=f || b!=g)){
    value=func_cmp(a,b);
    if(value){ return value; }
    a=func_get_number_part(f); if(a==NULL){ a=func_one(); }else{ FR(a); }
    b=func_get_number_part(g); if(b==NULL){ b=func_one(); }else{ FR(b); }
    value=func_cmp(a,b);
    a=func_del(a);
    b=func_del(b);
    return value;
  }
  // op
  value=char_cmp(func_op(f),func_op(g));
  if(value){
    i=func_find_order(f);
    j=func_find_order(g);
    if(i<j){ return -1; } // f<g
    if(i>j){ return +1; } // f>g
    return value;         // f<g, f>g
  }
  // p
  cmp1=func_find_p_cmp(f);
  cmp2=func_find_p_cmp(g);
  if(cmp1!=cmp2){ FUNC_ERROR_ARG2("func_cmp",f,g); }
  if(cmp1!=NULL){
    value=cmp1(f,g);
    if(value){ return value; } // f<g, f>g
  }
  // args
  if(func_find_amax(f) && func_find_amax(g)){
    for(i=0; i<MIN2(func_asize(f),func_asize(g)); i++){
      value=func_cmp(func_aget(f,i),func_aget(g,i));
      if(value){ return value; } // f<g or f>g
    }
    if     (func_asize(f)<func_asize(g)) { return -1; } // f<g
    else if(func_asize(f)>func_asize(g)) { return +1; } // f>g
  }
  // f=g
  return 0;
}

///////////////////////////////////////////////////////////////////////////////////////

// f<g ?
int func_lt(func_t *f, func_t *g)
{
  return func_cmp(f,g)<0;
}

// f<=g ?
int func_le(func_t *f, func_t *g)
{
  return func_cmp(f,g)<=0;
}

// f<g ?
int func_gt(func_t *f, func_t *g)
{
  return func_cmp(f,g)>0;
}

// f<=g ?
int func_ge(func_t *f, func_t *g)
{
  return func_cmp(f,g)>=0;
}

//EOF
