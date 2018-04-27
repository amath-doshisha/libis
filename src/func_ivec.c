#include<stdio.h>
#include<stdlib.h>
#include"is_macros.h"
#include"is_ivec.h"
#include"is_func.h"

#define FR(F) func_retain(F)

///////////////////////////////////////////////////

void func_print_ivec(func_t *f)
{
  int i;
  printf("int[");
  for(i=0; i<func_ivec_size(f); i++){
    printf("%d",func_ivec_at(f,i));
    if(i<(func_ivec_size(f))-1) printf(" ");
  }
  printf("]");
}

///////////////////////////////////////////////////

void func_ivec_p_new(func_t *f)
{
  if(f->p.mem!=NULL){ func_ivec_p_del(f); }
  f->ptype=FUNC_P_IVEC;
  f->p.mem=malloc(sizeof(func_ivec_struct));
  f->p.ivec->n=0;
  f->p.ivec->x=NULL;
}

void func_ivec_p_del(func_t *f)
{
  if(func_ptype(f)!=FUNC_P_IVEC){ FUNC_ERROR_ARG1("func_ivec_p_del",f); }
  if(f->p.ivec!=NULL){
    f->p.ivec->x=ivec_free(f->p.ivec->x);
    free(f->p.ivec);
    f->p.ivec=NULL;
    f->ptype=FUNC_P_NULL;
  }
}

void func_ivec_p_clone(func_t *f, func_t *g)
{
  int n;
  if(func_ptype(g)!=FUNC_P_IVEC){ FUNC_ERROR_ARG2("func_ivec_p_clone",f,g); }
  if(f->p.mem!=NULL){ func_ivec_p_del(f); }
  f->ptype=FUNC_P_IVEC;
  f->p.mem=malloc(sizeof(func_ivec_struct));
  n=g->p.ivec->n;
  f->p.ivec->n=n;
  f->p.ivec->x=ivec_allocate(n);
  ivec_copy(n,f->p.ivec->x,g->p.ivec->x);
}

/////////////////////////////////////////

void func_ivec_p_resize(func_t *f, int n)
{
  if(func_ptype(f)!=FUNC_P_IVEC){ FUNC_ERROR_ARG1("func_ivec_p_resize",f); }
  if(f->p.ivec==NULL){ func_ivec_p_new(f); }
  if(n>0){
    f->p.ivec->x=ivec_free(f->p.ivec->x);
    f->p.ivec->n=n;
    f->p.ivec->x=ivec_allocate(n);
  }
}

/////////////////////////////////////////

int func_ivec_size(func_t *f)
{
  if(f==NULL || func_ptype(f)!=FUNC_P_IVEC || f->p.ivec==NULL){ FUNC_ERROR_ARG1("func_ivec_size",f); }
  return f->p.ivec->n;
}

int func_ivec_at(func_t *f, int i)
{
  if(f==NULL || func_ptype(f)!=FUNC_P_IVEC || f->p.ivec==NULL || f->p.ivec->x==NULL || i<0 || i>=f->p.ivec->n){ FUNC_ERROR_ARG1("func_ivec_at",f); }
  return f->p.ivec->x[i];
}

int *func_ivec_p(func_t *f)
{
  if(f==NULL || func_ptype(f)!=FUNC_P_IVEC || f->p.ivec==NULL){ FUNC_ERROR_ARG1("func_ivec_p",f); }
  return f->p.ivec->x;
}

/////////////////////////////////////////

func_t *func_ivec(int n)
{
  func_t *f=NULL;
  f=func_new("ivec");
  func_ivec_p_resize(f,n);
  return f;
}

int func_ivec_cmp(func_t *f, func_t *g)
{
  if(!func_is(f,"ivec") || !func_is(g,"ivec")){ FUNC_ERROR_ARG2("func_ivec_cmp",f,g); }
  return ivec_cmp(func_ivec_size(f),f->p.ivec->x,func_ivec_size(g),g->p.ivec->x);
}

//EOF
