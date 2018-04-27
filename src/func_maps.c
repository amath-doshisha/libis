#include<stdio.h>
#include<stdlib.h>
#include"is_macros.h"
#include"is_ivec.h"
#include"is_strings.h"
#include"is_print.h"
#include"is_func.h"
#include"is_rmulti.h"
#include"is_cmulti.h"
#include"is_rvec.h"
#include"is_cvec.h"

#define FR(F) func_retain(F)

///////////////////////////////////////////////////////////////

func_t *func_maps_var(func_t *f, int n0, func_t *x)
{
  int i,k,p;
  func_t *g=NULL,*h=NULL;
  if(!func_is_var(f)){ FUNC_ERROR_ARG1("func_maps_var",f); }
  g=func_clone(FR(f));
  if(func_is_list(x)){
    h=func_one();
    for(i=0; i<func_asize(x); i++){
      k=func_var_get_index(g,n0+i);
      if(k>=0){
	p=g->p.var->pow[k];
	g->p.var->pow[k]=0;
	h=func_mul(h,func_pow_n(FR(x->a[i]),p));
      }
    }
    g=func_var_check_mul(g);
    g=func_mul(g,FR(h));
    h=func_del(h);
  }else if(func_is(x,"rvec")){
    h=func_one();
    for(i=0; i<func_rvec_size(x); i++){
      k=func_var_get_index(g,n0+i);
      if(k>=0){
	p=g->p.var->pow[k];
	g->p.var->pow[k]=0;
	h=func_mul(h,func_pow_n(func_rvec_get(x,i),p));
      }
    }
    g=func_var_check_mul(g);
    g=func_mul(g,FR(h));
    h=func_del(h);
  }else if(func_is(x,"cvec")){
    h=func_one();
    for(i=0; i<func_cvec_size(x); i++){
      k=func_var_get_index(g,n0+i);
      if(k>=0){
	p=g->p.var->pow[k];
	g->p.var->pow[k]=0;
	h=func_mul(h,func_pow_n(func_cvec_get(x,i),p));
      }
    }
    g=func_var_check_mul(g);
    g=func_mul(g,FR(h));
    h=func_del(h);
  }else{
    k=func_var_get_index(g,n0);
    if(k>=0){
      p=g->p.var->pow[k];
      g->p.var->pow[k]=0;
      g=func_var_check_mul(g);
      g=func_mul(g,func_pow_n(FR(x),p));
    }
  }
  f=func_del(f);
  x=func_del(x);
  return g;
}

// y=f(x)
func_t *func_maps(func_t *f, int n0, func_t *x)
{
  int i;
  func_t *g=NULL;
  func_arg1_t *eval=NULL;
  if(f==NULL)            { FUNC_ERROR_ARG1("func_maps",f); }
  else if(func_is_var(f)){ g=func_maps_var(FR(f),n0,FR(x)); }
  else{
    g=func_clone(FR(f));
    for(i=0; i<func_asize(f); i++){
      func_aset(g,i,func_maps(FR(func_aget(f,i)),n0,FR(x)));
    }
    eval=func_find_eval(g);
    if(eval!=NULL){ g=eval(g); }
  }
  f=func_del(f);
  x=func_del(x);
  return g;
}

//EOF
