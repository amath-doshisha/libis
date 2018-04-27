#include<stdio.h>
#include<stdlib.h>
#include"is_macros.h"
#include"is_func.h"
#include"is_rmulti.h"
#include"is_rvec.h"
#include"is_cvec.h"

#define FR(F) func_retain(F)

///////////////////////////////////////////////////

void func_print_rvec(func_t *f)
{
  int i;
  //printf("real[");
  printf("[");
  for(i=0; i<func_rvec_size(f); i++){
    mpfr_printf("%.5Rg",func_rvec_at(f,i));
    if(i<(func_rvec_size(f))-1) printf(" ");
  }
  printf("]");
}

///////////////////////////////////////////////////

void func_rvec_p_new(func_t *f)
{
  if(f->p.mem!=NULL){ func_rvec_p_del(f); }
  f->ptype=FUNC_P_RVEC;
  f->p.mem=malloc(sizeof(func_rvec_struct));
  f->p.rvec->n=0;
  f->p.rvec->x=NULL;
}

void func_rvec_p_del(func_t *f)
{
  if(func_ptype(f)!=FUNC_P_RVEC){ FUNC_ERROR_ARG1("func_rvec_p_del",f); }
  if(f->p.rvec!=NULL){
    f->p.rvec->x=rvec_free(func_rvec_size(f),f->p.rvec->x);
    free(f->p.rvec);
    f->p.rvec=NULL;
    f->ptype=FUNC_P_NULL;
  }
}

void func_rvec_p_clone(func_t *f, func_t *g)
{
  if(func_ptype(g)!=FUNC_P_RVEC){ FUNC_ERROR_ARG2("func_rvec_p_clone",f,g); }
  if(f->p.mem!=NULL){ func_rvec_p_del(f); }
  f->ptype=FUNC_P_RVEC;
  f->p.mem=malloc(sizeof(func_rvec_struct));
  f->p.rvec->n=g->p.rvec->n;
  f->p.rvec->x=rvec_allocate_clone(func_rvec_size(g),g->p.rvec->x);
}

///////////////////////////////////////////////////

void func_rvec_p_resize(func_t *f, int n)
{
  if(func_ptype(f)!=FUNC_P_RVEC){ FUNC_ERROR_ARG1("func_rvec_p_resize",f); }
  if(f->p.rvec==NULL){ func_rvec_p_new(f); }
  if(n>0){
    f->p.rvec->x=rvec_free(func_rvec_size(f),f->p.rvec->x);
    f->p.rvec->n=n;
    f->p.rvec->x=rvec_allocate(n);
  }
}

///////////////////////////////////////////////////

int func_rvec_size(func_t *f)
{
  if(f==NULL || func_ptype(f)!=FUNC_P_RVEC || f->p.rvec==NULL){ FUNC_ERROR_ARG1("func_rvec_size",f); }
  return f->p.rvec->n;
}

rmulti *func_rvec_at(func_t *f, int i)
{
  if(f==NULL || func_ptype(f)!=FUNC_P_RVEC || f->p.rvec==NULL || f->p.rvec->x==NULL || i<0 || i>=f->p.rvec->n){ FUNC_ERROR_ARG1("func_rvec_at",f); }  
  return f->p.rvec->x[i];
}

rmulti **func_rvec_p(func_t *f)
{
  if(f==NULL || func_ptype(f)!=FUNC_P_RVEC || f->p.rvec==NULL){ FUNC_ERROR_ARG1("func_rvec_p",f); }
  return f->p.rvec->x;
}

///////////////////////////////////////////////

func_t *func_rvec_get(func_t *f, int i)
{
  return func_real_rmulti(func_rvec_at(f,i));
}

void func_rvec_set(func_t *f, int i, func_t *g)
{
  func_t *a=NULL;
  if(f==NULL || func_ptype(f)!=FUNC_P_RVEC || f->p.rvec==NULL || i<0 || i>=f->p.rvec->n){ FUNC_ERROR_ARG2("func_rvec_set",f,g); }
  a=func_evalf(FR(g));
  if     (func_is_real(a))                          { rcopy(func_rvec_at(f,i),func_real_p(a)); }
  else if(func_is(a,"rvec") && func_rvec_size(a)==1){ rcopy(func_rvec_at(f,i),func_rvec_at(a,0)); }
  else{ FUNC_ERROR_ARG2("func_rvec_set",f,g); } 
  a=func_del(a);
}

///////////////////////////////////////////////

func_t *func_rvec(int n)
{
  func_t *f=NULL;
  f=func_new("rvec");
  func_rvec_p_resize(f,n);
  return f;
}

int func_rvec_cmp(func_t *f, func_t *g)
{
  if(!func_is(f,"rvec") || !func_is(g,"rvec")){ FUNC_ERROR_ARG2("func_rvec_cmp",f,g); }
  return rvec_cmp(func_rvec_size(f),f->p.rvec->x,func_rvec_size(g),g->p.rvec->x);
}

///////////////////////////////////////////////

func_t *func_rvec_get_cvec(func_t *g)
{
  func_t *f=NULL;
  if(g==NULL || !func_is(g,"rvec")){ FUNC_ERROR_ARG1("func_rvec_get_cvec",g); }
  f=func_cvec(func_rvec_size(g));
  cvec_copy_rvec(func_cvec_size(f),func_cvec_p(f),func_rvec_p(g));
  return f;
}

//EOF
