#include<stdio.h>
#include<stdlib.h>
#include"is_macros.h"
#include"is_func.h"
#include"is_cmulti.h"
#include"is_cmat.h"

#define FR(F) func_retain(F)

////////////////////////////////////////////////

void func_print_cmat(func_t *f)
{
  int i,j;
  //printf("complex[");
  printf("[");
  for(i=0; i<func_cmat_rows(f); i++){
    for(j=0; j<func_cmat_cols(f); j++){
      if(cis_real(func_cmat_at(f,i,j))){
	mpfr_printf("%.5Rg",C_R(func_cmat_at(f,i,j)));
      }else if(cis_pure_imaginary(func_cmat_at(f,i,j))){
	mpfr_printf("%.5Rg*I",C_I(func_cmat_at(f,i,j)));
      }else{
	mpfr_printf("%.5Rg%+.5Rg*I",C_R(func_cmat_at(f,i,j)),C_I(func_cmat_at(f,i,j)));
      }
      if(j<(func_cmat_cols(f))-1) printf(" ");
    }
    if(i<(func_cmat_rows(f))-1) printf("; ");
  }
  printf("]");
}

////////////////////////////////////////////////

void func_cmat_p_new(func_t *f)
{
  if(f->p.mem!=NULL){ func_cmat_p_del(f); }
  f->ptype=FUNC_P_CMAT;
  f->p.mem=malloc(sizeof(func_cmat_struct));
  f->p.cmat->LD=0;
  f->p.cmat->m=0;
  f->p.cmat->n=0;
  f->p.cmat->A=NULL;
}

void func_cmat_p_del(func_t *f)
{
  if(func_ptype(f)!=FUNC_P_CMAT){ FUNC_ERROR_ARG1("func_cmat_p_del",f); }
  if(f->p.cmat!=NULL){
    f->p.cmat->A=cmat_free(func_cmat_ld(f),func_cmat_cols(f),func_cmat_p(f));
    free(f->p.cmat);
    f->p.cmat=NULL;
    f->ptype=FUNC_P_NULL;
  }
}

void func_cmat_p_clone(func_t *f, func_t *g)
{
  if(func_ptype(g)!=FUNC_P_CMAT){ FUNC_ERROR_ARG2("func_cmat_p_clone",f,g); }
  if(f->p.mem!=NULL){ func_cmat_p_del(f); }
  f->ptype=FUNC_P_CMAT;
  f->p.mem=malloc(sizeof(func_cmat_struct));
  f->p.cmat->LD=func_cmat_ld(g);
  f->p.cmat->m=func_cmat_rows(g);
  f->p.cmat->n=func_cmat_cols(g);
  f->p.cmat->A=cmat_allocate_clone(func_cmat_ld(g),func_cmat_cols(g),func_cmat_p(g));
}

////////////////////////////////////////////////

void func_cmat_p_resize(func_t *f, int m, int n)
{
  if(func_ptype(f)!=FUNC_P_CMAT){ FUNC_ERROR_ARG1("func_cmat_p_resize",f); }
  if(f->p.cmat==NULL){ func_cmat_p_new(f); }
  if(m>0 && n>0){
    f->p.cmat->A=cmat_free(func_cmat_ld(f),func_cmat_cols(f),func_cmat_p(f));
    f->p.cmat->LD=m;
    f->p.cmat->m=m;
    f->p.cmat->n=n;
    f->p.cmat->A=cmat_allocate(m,n);
  }
}

////////////////////////////////////////////////


func_t *func_cmat(int m, int n)
{
  func_t *f=NULL;
  f=func_new("cmat");
  func_cmat_p_resize(f,m,n);
  return f;
}

int func_cmat_cmp(func_t *f, func_t *g)
{
  if(!func_is(f,"cmat") || !func_is(g,"cmat")){ FUNC_ERROR_ARG2("func_cmat_cmp",f,g); }
  return cmat_cmp(func_cmat_rows(f),func_cmat_cols(f),func_cmat_p(f),func_cmat_ld(f),
		  func_cmat_rows(g),func_cmat_cols(g),func_cmat_p(g),func_cmat_ld(g));
}

///////////////////////////////////////////////////////////////

int func_cmat_rows(func_t *f)
{
  if(f==NULL || func_ptype(f)!=FUNC_P_CMAT || f->p.cmat==NULL){ FUNC_ERROR_ARG1("func_cmat_rows",f); }
  return f->p.cmat->m;
}

int func_cmat_cols(func_t *f)
{
  if(f==NULL || func_ptype(f)!=FUNC_P_CMAT || f->p.cmat==NULL){ FUNC_ERROR_ARG1("func_cmat_cols",f); }
  return f->p.cmat->n;
}

int func_cmat_ld(func_t *f)
{
  if(f==NULL || func_ptype(f)!=FUNC_P_CMAT || f->p.cmat==NULL){ FUNC_ERROR_ARG1("func_cmat_ld",f); }
  return f->p.cmat->LD;
}

cmulti *func_cmat_at(func_t *f, int i, int j)
{
  if(f==NULL || func_ptype(f)!=FUNC_P_CMAT || f->p.cmat==NULL || f->p.cmat->A==NULL || i<0 || j<0 || i>=f->p.cmat->m || j>=f->p.cmat->n){ FUNC_ERROR_ARG1("func_cmat_at",f); }  
  return MAT(f->p.cmat->A,i,j,f->p.cmat->LD);
}

cmulti **func_cmat_p(func_t *f)
{
  if(f==NULL || func_ptype(f)!=FUNC_P_CMAT || f->p.cmat==NULL){ FUNC_ERROR_ARG1("func_cmat_p",f); }
  return f->p.cmat->A;
}

///////////////////////////////////////////////

void func_cmat_set_row(func_t *f, int i, func_t *g)
{
  int j;
  func_t *a=NULL;
  if(f==NULL || func_ptype(f)!=FUNC_P_CMAT || f->p.cmat==NULL || i<0 || i>=func_cmat_rows(f)){ FUNC_ERROR_ARG2("func_cmat_set_row",f,g); }
  a=func_evalf(FR(g));
  if     (func_is_real(a)    && func_cmat_cols(f)==1){ cset_r(func_cmat_at(f,i,0),func_real_p(a)); }
  else if(func_is_complex(a) && func_cmat_cols(f)==1){ cset_c  (func_cmat_at(f,i,0),func_complex_p(a)); }
  else if(func_is_rvec(a)){
    if(func_cmat_cols(f)==func_rvec_size(a)){
      for(j=0; j<func_cmat_cols(f); j++){ cset_r(MAT(func_cmat_p(f),i,j,func_cmat_ld(f)),func_rvec_at(a,j)); }
    }else{ FUNC_ERROR_ARG2("func_cmat_set_row",f,g); }
  }else if(func_is_cvec(a)){
    if(func_cmat_cols(f)==func_cvec_size(a)){
      for(j=0; j<func_cmat_cols(f); j++){ cset_c(MAT(func_cmat_p(f),i,j,func_cmat_ld(f)),func_cvec_at(a,j)); }
    }else{ FUNC_ERROR_ARG2("func_cmat_set_row",f,g); }
  }else{ FUNC_ERROR_ARG2("func_cmat_set_row",f,g); }
  a=func_del(a);
}

//EOF
