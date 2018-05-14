#include<stdio.h>
#include<stdlib.h>
#include"is_macros.h"
#include"is_func.h"
#include"is_rmulti.h"
#include"is_rmat.h"
#include"is_cmat.h"

#define FR(F) func_retain(F)

////////////////////////////////////////////////

void func_print_rmat(func_t *f)
{
  int i,j;
  //printf("real[");
  printf("[");
  for(i=0; i<func_rmat_rows(f); i++){
    for(j=0; j<func_rmat_cols(f); j++){
      mpfr_printf("%.5Rg",func_rmat_at(f,i,j));
      if(j<(func_rmat_cols(f))-1) printf(" ");
    }
    if(i<(func_rmat_rows(f))-1) printf("; ");
  }
  printf("]");
}

////////////////////////////////////////////////

void func_rmat_p_new(func_t *f)
{
  if(f->p.mem!=NULL){ func_rmat_p_del(f); }
  f->ptype=FUNC_P_RMAT;
  f->p.mem=malloc(sizeof(func_rmat_struct));
  f->p.rmat->LD=0;
  f->p.rmat->m=0;
  f->p.rmat->n=0;
  f->p.rmat->A=NULL;
}

void func_rmat_p_del(func_t *f)
{
  if(func_ptype(f)!=FUNC_P_RMAT){ FUNC_ERROR_ARG1("func_rmat_p_del",f); }
  if(f->p.rmat!=NULL){
    f->p.rmat->A=rmat_free(f->p.rmat->LD,f->p.rmat->n,f->p.rmat->A);
    free(f->p.rmat);
    f->p.rmat=NULL;
    f->ptype=FUNC_P_NULL;
  }
}

void func_rmat_p_clone(func_t *f, func_t *g)
{
  if(func_ptype(g)!=FUNC_P_RMAT){ FUNC_ERROR_ARG2("func_rmat_p_clone",f,g); }
  if(f->p.mem!=NULL){ func_rmat_p_del(f); }
  f->ptype=FUNC_P_RMAT;
  f->p.mem=malloc(sizeof(func_rmat_struct));
  f->p.rmat->LD=g->p.rmat->LD;
  f->p.rmat->m=g->p.rmat->m;
  f->p.rmat->n=g->p.rmat->n;
  f->p.rmat->A=rmat_allocate_clone(g->p.rmat->LD,g->p.rmat->n,g->p.rmat->A);
}

////////////////////////////////////////////////

void func_rmat_p_resize(func_t *f, int m, int n)
{
  if(func_ptype(f)!=FUNC_P_RMAT){ FUNC_ERROR_ARG1("func_rmat_p_resize",f); }
  if(f->p.rmat==NULL){ func_rmat_p_new(f); }
  if(m>0 && n>0){
    f->p.rmat->A=rmat_free(f->p.rmat->LD,f->p.rmat->n,f->p.rmat->A);
    f->p.rmat->LD=m;
    f->p.rmat->m=m;
    f->p.rmat->n=n;
    f->p.rmat->A=rmat_allocate(m,n);
  }
}

////////////////////////////////////////////////

func_t *func_rmat(int m, int n)
{
  func_t *f=NULL;
  f=func_new("rmat");
  func_rmat_p_resize(f,m,n);
  return f;
}

int func_rmat_cmp(func_t *f, func_t *g)
{
  if(!func_is(f,"rmat") || !func_is(g,"rmat")){ FUNC_ERROR_ARG2("func_rmat_cmp",f,g); }
  return rmat_cmp(f->p.rmat->m,f->p.rmat->n,f->p.rmat->A,f->p.rmat->LD,
		  g->p.rmat->m,g->p.rmat->n,g->p.rmat->A,g->p.rmat->LD);
}

////////////////////////////////////////////////

int func_rmat_rows(func_t *f)
{
  if(f==NULL || func_ptype(f)!=FUNC_P_RMAT || f->p.rmat==NULL){ FUNC_ERROR_ARG1("func_rmat_rows",f); }
  return f->p.rmat->m;
}

int func_rmat_cols(func_t *f)
{
  if(f==NULL || func_ptype(f)!=FUNC_P_RMAT || f->p.rmat==NULL){ FUNC_ERROR_ARG1("func_rmat_cols",f); }
  return f->p.rmat->n;
}

int func_rmat_ld(func_t *f)
{
  if(f==NULL || func_ptype(f)!=FUNC_P_RMAT || f->p.rmat==NULL){ FUNC_ERROR_ARG1("func_rmat_ld",f); }
  return f->p.rmat->LD;
}

rmulti *func_rmat_at(func_t *f, int i, int j)
{
  if(f==NULL || func_ptype(f)!=FUNC_P_RMAT || f->p.rmat==NULL || f->p.rmat->A==NULL || i<0 || j<0 || i>=f->p.rmat->m || j>=f->p.rmat->n){ FUNC_ERROR_ARG1("func_rmat_at",f); }
  return MAT(f->p.rmat->A,i,j,f->p.rmat->LD);
}

rmulti **func_rmat_p(func_t *f)
{
  if(f==NULL || func_ptype(f)!=FUNC_P_RMAT || f->p.rmat==NULL){ FUNC_ERROR_ARG1("func_rmat_p",f); }
  return f->p.rmat->A;
}

///////////////////////////////////////////////

void func_rmat_set_row(func_t *f, int i, func_t *g)
{
  int j;
  func_t *a=NULL;
  if(f==NULL || func_ptype(f)!=FUNC_P_RMAT || f->p.rmat==NULL || i<0 || i>=func_rmat_rows(f)){ FUNC_ERROR_ARG2("func_rmat_set_row",f,g); }
  a=func_evalf(FR(g));
  if(func_is_real(a)){
    if(func_rmat_cols(f)==1){ rset_r(func_rmat_at(f,i,0),func_real_p(a)); }
    else{ FUNC_ERROR_ARG2("func_rmat_set_row",f,g); }
  }else if(func_is(a,"rvec")){
    if(func_rmat_cols(f)==func_rvec_size(a)){
      for(j=0; j<func_rmat_cols(f); j++){ rset_r(MAT(func_rmat_p(f),i,j,func_rmat_ld(f)),func_rvec_at(a,j)); }
    }else{ FUNC_ERROR_ARG2("func_rmat_set_row",f,g); }
  }else{ FUNC_ERROR_ARG2("func_rmat_set_row",f,g); }
  a=func_del(a);
}

///////////////////////////////////////////////

func_t *func_rmat_get_cmat(func_t *g)
{
  func_t *f=NULL;
  if(g==NULL || !func_is(g,"rmat")){ FUNC_ERROR_ARG1("func_rmat_get_cmat",g); }
  f=func_cmat(func_rmat_rows(g),func_rmat_cols(g));
  cmat_copy_rmat(func_cmat_rows(f),func_cmat_cols(f),func_cmat_p(f),func_cmat_ld(f),func_rmat_p(g),func_rmat_ld(g));
  return f;
}

//EOF
