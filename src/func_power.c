#include<stdio.h>
#include<stdlib.h>
#include"is_macros.h"
#include"is_func.h"

#define FR(F) (func_retain(F))

//////////////////////////////////////////////////////////////////

void func_pow_p_new(func_t *f)
{
  if(f->p.mem!=NULL){ func_pow_p_del(f); }
  f->ptype=FUNC_P_POWER;
  f->p.mem=malloc(sizeof(int)*1);
  f->p.i[0]=1;
}

void func_pow_p_del(func_t *f)
{
  if(func_ptype(f)!=FUNC_P_POWER){ FUNC_ERROR_ARG1("func_pow_p_del",f); }
  if(f->p.mem!=NULL){
    free(f->p.mem);
    f->p.mem=NULL;
    f->ptype=FUNC_P_NULL;
  }
}

void func_pow_p_clone(func_t *f, func_t *g)
{
  if(func_ptype(g)!=FUNC_P_POWER){ FUNC_ERROR_ARG2("func_pow_p_clone",f,g); }
  if(f->p.mem!=NULL){ func_pow_p_del(f); }
  f->ptype=FUNC_P_POWER;
  f->p.mem=malloc(sizeof(int)*1);
  f->p.i[0]=g->p.i[0];
}

//////////////////////////////////////////////////////////////////

int func_has_power(func_t *f)
{
  return (func_ptype(f)==FUNC_P_POWER);
}

int func_power(func_t *f)
{
  if(f!=NULL && f->p.i!=NULL && func_has_power(f)){ return f->p.i[0]; }
  else{ FUNC_ERROR_ARG1("func_power",f); }
}

void func_set_power(func_t *f, int pow)
{
  if(f!=NULL && f->p.i!=NULL && func_has_power(f)){ f->p.i[0]=pow; }
  else{ FUNC_ERROR_ARG1("func_set_power",f); }
}

//////////////////////////////////////////////////////////////////////////

int func_power_cmp(func_t *f, func_t *g)
{
  if(f==NULL || g==NULL || func_ptype(f)!=FUNC_P_POWER || func_ptype(g)!=FUNC_P_POWER || f->p.i==NULL || g->p.i==NULL){ FUNC_ERROR_ARG2("func_power_cmp",f,g); }
  if((f->p.i[0]) < (g->p.i[0])){ return -1; } // f<g
  if((f->p.i[0]) > (g->p.i[0])){ return +1; } // f>g
  return 0;                                   // f==g
}

//EOF
