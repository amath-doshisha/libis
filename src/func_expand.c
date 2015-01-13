#include<stdio.h>
#include<stdlib.h>
#include"is_macros.h"
#include"is_ivec.h"
#include"is_strings.h"
#include"is_print.h"
#include"is_func.h"
#include"is_rmulti.h"
#include"is_cmulti.h"

#define FR(F) (func_retain(F))

///private////////////////////////////////////////////////////

func_t *func_expand_mul_pow_n(func_t *f)
{
  int i;
  func_t *g=NULL;
  if(func_is_mul(f) && func_power(f)!=1){
    g=func_one();
    for(i=0; i<func_asize(f); i++){
      g=func_mul(g,func_expand(func_pow_n(FR(f->a[i]),func_power(f))));
    }
  }else{ FUNC_ERROR_ARG1("func_expand_mul_pow_n",f); }
  f=func_del(f);
  return g;
}

func_t *func_expand_add_pow_n(func_t *f)
{
  int i,j,k;
  func_t *g=NULL,*h=NULL;
  if(func_is_add(f) && func_power(f)>=2){
    g=func_clone(FR(f)); func_set_power(g,1);
    for(k=1; k<func_power(f); k++){
      h=func_zero();
      for(i=0; i<func_asize(g); i++){	
	for(j=0; j<func_asize(f); j++){
	  h=func_add(h,func_expand(func_mul(FR(g->a[i]),FR(f->a[j]))));
	}
      }
      g=func_del(g); g=h; h=NULL;
    }
  }else{ FUNC_ERROR_ARG1("func_expand_add_pow_n",f); }
  f=func_del(f);
  return g;
}

func_t *func_expand_mul1_add1(func_t *f)
{
  int i,j,k;
  func_t *h=NULL,*g=NULL;
  if(func_is_mul(f) && func_power(f)==1 && func_a_has_op_pow1(f,func_is_add)){
    // get terms without (add)^1
    g=func_clone(FR(f));
    func_args_rm_add_pow1(g);
    func_a_rm_null(g);
    g=func_mul_eval(g);
    // loop
    for(i=0; i<func_asize(f); i++){
      if(func_is_add(f->a[i]) && func_power(f->a[i])==1){
	if(func_is_add(g) && func_power(g)==1){
	  h=func_zero();
	  for(k=0; k<func_asize(g); k++){
	    for(j=0; j<func_asize(f->a[i]); j++){
	      h=func_add(h,func_expand(func_mul(FR(g->a[k]),FR(f->a[i]->a[j]))));
	    }
	  }
	}else{	  
	  h=func_zero();
	  for(j=0; j<func_asize(f->a[i]); j++){
	    h=func_add(h,func_expand(func_mul(FR(g),FR(f->a[i]->a[j]))));
	  }
	}
	g=func_del(g); g=h; h=NULL;
      }
    }
  }else{ FUNC_ERROR_ARG1("func_expand_mul1_add1",f); }
  f=func_del(f);
  return g;
}

func_t *func_expand_add_pow1(func_t *f)
{
  int i;
  func_t *g=NULL;
  if(func_is_add(f) && func_power(f)==1){
    g=func_zero();
    for(i=0; i<func_asize(f); i++){
      g=func_add(g,func_expand(FR(f->a[i])));
    }
  }else{ FUNC_ERROR_ARG1("func_expand_add_pow1",f); }
  f=func_del(f);
  return g;
}

func_t *func_expand_mul_pow1(func_t *f)
{
  int i;
  func_t *g=NULL,*h=NULL;
  if(func_is_mul(f) && func_power(f)==1){
    g=func_one();
    for(i=0; i<func_asize(f); i++){
      h=func_expand(FR(f->a[i]));
      if((func_is_add(g) && func_power(g)==1) || (func_is_add(h) && func_power(h)==1)){
	g=func_expand(func_mul(g,h));
      }else{
	g=func_mul(g,h);
      }
    }
  }else{ FUNC_ERROR_ARG1("func_expand_mul_pow1",f); }
  f=func_del(f);
  return g;
}

func_t *func_expand_list(func_t *f)
{
  int i;
  func_t *g=NULL;
  if(func_is_list(f)){
    g=func_list(func_asize(f));
    for(i=0; i<func_asize(f); i++){
      g->a[i]=func_expand(FR(f->a[i]));
    }
  }else{ FUNC_ERROR_ARG1("func_expand_list",f); }
  f=func_del(f);
  return g;
}

///////////////////////////////////////////////////////

func_t *func_expand(func_t *f)
{
  return func_expand_eval(func_arg1_new("expand",f));
}

func_t *func_expand_eval(func_t *g)
{
  func_t *f=NULL,*h=NULL;
  if(g==NULL || !func_is(g,"expand") || func_asize(g)!=1 || func_aget(g,0)==NULL){ FUNC_ERROR_ARG1("func_expand_eval",g); }
  f=func_aget(g,0);
  h=func_clone(FR(f));
  if(func_is_mul(h)){ h=func_flatten(h,"mul"); }   // (mul)^1
  if(func_is_add(h)){ h=func_flatten(h,"add"); }   // (add)^1
  if     (func_is_mul(h) && func_power(h)!=1)                              { h=func_expand(func_expand_mul_pow_n(h)); } // (mul)^n with n!=1
  else if(func_is_add(h) && func_power(h)>=2)                              { h=func_expand_add_pow_n(h); }              // (add)^n with n>=2
  else if(func_is_mul(h) && func_power(h)==1 && func_a_has_op_pow1(h,func_is_add)){ h=func_expand_mul1_add1(h); }              // (mul ... add^1 ...)^1
  else if(func_is_add(h) && func_power(h)==1)                              { h=func_expand_add_pow1(h); }               // (add)^1
  else if(func_is_mul(h) && func_power(h)==1)                              { h=func_expand_mul_pow1(h); }               // (mul)^1
  else if(func_is_list(h))                                                 { h=func_expand_list(h); }                   // {...}
  else                                                                     { }                                          // otherwise
  g=func_del(g);
  return h;
}

//EOF
