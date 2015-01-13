#include<stdio.h>
#include<stdlib.h>
#include"is_macros.h"
#include"is_func.h"
#include"is_cmulti.h"
#include"is_cvec.h"

#define FR(F) func_retain(F)

////////////////////////////////////////////////

void func_print_cvec(func_t *f)
{
  int i;
  //printf("complex[");
  printf("[");
  for(i=0; i<func_cvec_size(f); i++){
    if(cis_real(func_cvec_at(f,i))){
      mpfr_printf("%.5Rg",C_R(func_cvec_at(f,i)));
    }else if(cis_pure_imaginary(func_cvec_at(f,i))){
      mpfr_printf("%.5Rg*I",C_I(func_cvec_at(f,i)));
    }else{
      mpfr_printf("%.5Rg%+.5Rg*I",C_R(func_cvec_at(f,i)),C_I(func_cvec_at(f,i)));
    }
    if(i<(func_cvec_size(f))-1) printf(" ");
  }
  printf("]");
}

////////////////////////////////////////////////

void func_cvec_p_new(func_t *f)
{
  if(f->p.mem!=NULL){ func_cvec_p_del(f); }
  f->ptype=FUNC_P_CVEC;
  f->p.mem=malloc(sizeof(func_cvec_struct));
  f->p.cvec->n=0;
  f->p.cvec->x=NULL;
}

void func_cvec_p_del(func_t *f)
{
  if(func_ptype(f)!=FUNC_P_CVEC){ FUNC_ERROR_ARG1("func_cvec_p_del",f); }
  if(f->p.cvec!=NULL){
    f->p.cvec->x=cvec_free(f->p.cvec->n,f->p.cvec->x);
    free(f->p.cvec);
    f->p.cvec=NULL;
    f->ptype=FUNC_P_NULL;
  }
}

void func_cvec_p_clone(func_t *f, func_t *g)
{
  if(func_ptype(g)!=FUNC_P_CVEC){ FUNC_ERROR_ARG2("func_cvec_p_clone",f,g); }
  if(f->p.mem!=NULL){ func_cvec_p_del(f); }
  f->ptype=FUNC_P_CVEC;
  f->p.mem=malloc(sizeof(func_cvec_struct));
  f->p.cvec->n=g->p.cvec->n;
  f->p.cvec->x=cvec_allocate_clone(g->p.cvec->n,g->p.cvec->x);
}

////////////////////////////////////////////////

void func_cvec_p_resize(func_t *f, int n)
{
  if(func_ptype(f)!=FUNC_P_CVEC){ FUNC_ERROR_ARG1("func_cvec_p_resize",f); }
  if(f->p.cvec==NULL){ func_cvec_p_new(f); }
  if(n>0){
    f->p.cvec->x=cvec_free(f->p.cvec->n,f->p.cvec->x);
    f->p.cvec->n=n;
    f->p.cvec->x=(cmulti**)cvec_allocate(n);
  }
}

////////////////////////////////////////////////

func_t *func_cvec(int n)
{
  func_t *f=NULL;
  f=func_new("cvec");
  func_cvec_p_resize(f,n);
  return f;
}

int func_cvec_cmp(func_t *f, func_t *g)
{
  if(!func_is(f,"cvec") || !func_is(g,"cvec")){ FUNC_ERROR_ARG2("func_cvec_cmp",f,g); }
  return cvec_cmp(f->p.cvec->n,f->p.cvec->x,g->p.cvec->n,g->p.cvec->x);
}

////////////////////////////////////////////////

int func_cvec_size(func_t *f)
{
  if(f==NULL || func_ptype(f)!=FUNC_P_CVEC || f->p.cvec==NULL){ FUNC_ERROR_ARG1("func_cvec_size",f); }
  return f->p.cvec->n;
}

cmulti *func_cvec_at(func_t *f, int i)
{
  if(f==NULL || func_ptype(f)!=FUNC_P_CVEC || f->p.cvec==NULL || f->p.cvec->x==NULL || i<0 || i>=f->p.cvec->n){ FUNC_ERROR_ARG1("func_cvec_at",f); }  
  return f->p.cvec->x[i];
}

cmulti **func_cvec_p(func_t *f)
{
  if(f==NULL || func_ptype(f)!=FUNC_P_CVEC || f->p.cvec==NULL){ FUNC_ERROR_ARG1("func_cvec_p",f); }
  return f->p.cvec->x;
}

////////////////////////////////////////////////////////////

func_t *func_cvec_get(func_t *f, int i)
{
  return func_complex_cmulti(func_cvec_at(f,i));
}

void func_cvec_set(func_t *f, int i, func_t *g)
{
  func_t *a=NULL;
  if(f==NULL || func_ptype(f)!=FUNC_P_CVEC || f->p.rvec==NULL || i<0 || i>=f->p.rvec->n){ FUNC_ERROR_ARG2("func_cvec_set",f,g); }
  a=func_evalf(FR(g));
  if(func_is_real(a))        { ccopy_r(func_cvec_at(f,i),func_real_p(a)); }
  else if(func_is_complex(a)){ ccopy  (func_cvec_at(f,i),func_complex_p(a)); }
  else                       { FUNC_ERROR_ARG2("func_cvec_set",f,g); }
  a=func_del(a);
}

//EOF
