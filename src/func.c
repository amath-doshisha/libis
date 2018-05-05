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

///////////////////////////////////////////////

static long long int __func_struct_new=0;
static long long int __func_struct_del=0;

func_t *func_struct_new(void)
{
  func_t *f=NULL;
  f=(func_t*)malloc(sizeof(func_t)*1);
  if(f==NULL){ ERROR_AT; printf("Cannot allocate struct func_t.\n"); exit(0); }
  f->count=1;
  f->op=NULL;
  f->n=0;
  f->a=NULL;
  f->ptype=FUNC_P_NULL;
  f->p.mem=NULL;
  __func_struct_new++;
  return f;
}

func_t *func_struct_del(func_t *f)
{
  __func_struct_del++;
  free(f);
  return NULL;
}

int func_new_del_check_sum()
{
  return (int)(__func_struct_new-__func_struct_del);
}

////////////////////////////////////////////////

func_t *func_new(char *op)
{
  func_t *f=NULL;
  if(op==NULL || char_eq(op,"NULL") || char_eq(op,"null")){ return f; }
  f=func_struct_new();
  f->op=char_new(op," \t\n");
  func_p_new(f);
  func_a_new(f,func_find_amin(f));
  return f;
}

func_t *func_del(func_t *f)
{
  if(f==NULL) return NULL;
  f->count--;
  if(f->count<=0){
    func_a_del(f);         // 1st (The order is important!)
    func_p_del(f);         // 2nd
    f->op=char_del(f->op); // 3rd
    f=func_struct_del(f);  // finished
  }
  return NULL;
}

func_t *func_clone(func_t *g)
{
  func_t *f=NULL;
  if(g==NULL){ return f; }
  f=func_struct_new();
  f->op=char_new(g->op," \t\n");
  func_p_clone(f,g);
  func_a_clone(f,g);
  g=func_del(g);
  return f;
}

func_t *func_retain(func_t *f)
{
  if(f!=NULL){ f->count++; }
  return f;
}

func_t *func_replace(func_t *f, func_t *g)
{
  f=func_del(f);
  return g;
}

//////////////////////////////////////////////////////////////////////

void func_p_new(func_t *f)
{
  func_set_t *p_new=NULL;
  if(f==NULL) return;
  if(f->p.mem!=NULL){ func_p_del(f); }
  p_new=func_find_p_new(f);
  if(p_new!=NULL){
    p_new(f);
  }
}

void func_p_del(func_t *f)
{
  func_set_t *p_del=NULL;
  if(f==NULL) return;
  p_del=func_find_p_del(f);
  if(p_del!=NULL && f->p.mem!=NULL){ p_del(f); }
  if(f->p.mem!=NULL){ FUNC_ERROR_ARG1("func_p_del",f); }
}

void func_p_clone(func_t *f, func_t *g)
{
  func_set2_t *p_clone=NULL;
  if(f==NULL || g==NULL) return;
  if(f->p.mem!=NULL){ func_p_del(f); }
  p_clone=func_find_p_clone(f);
  if(p_clone!=NULL){ p_clone(f,g); }
}

//////////////////////////////////////////////////////////////////

int func_ptype(func_t *f)
{
  if(f==NULL){ return FUNC_P_NULL; }
  return f->ptype;
}

char *func_op(func_t *f)
{
  if(f==NULL){ return "NULL"; }
  if(f->op==NULL){ return "null"; }
  return f->op;
}

int func_is(func_t *f, char *op)
{
  if(op==NULL)   { return 0; }  
  if(f==NULL)    { if(char_eq(op,"NULL")) return 1; else return 0; }
  if(f->op==NULL){ if(char_eq(op,"null")) return 1; else return 0; }
  return char_eq(f->op,op);
}

///////////////////////////////////////////////////////////////////////

int func_is_coeff(func_t *f)
{
  if     (f==NULL)             { return 0; }
  else if(func_in_complex(f))  { return 1; }
  else                         { return 0; }
}

////////////////////////////////////////////////////////////

func_t *func_flatten(func_t *f, char *op)
{
  int i,j,k=0,redo;
  func_t *g=NULL;
  if(!func_is(f,op) || func_find_amax(f)!=FUNC_NOLIMIT){ FUNC_ERROR_ARG1("func_flatten",f); }
  k=func_args_count_op(f,op);
  if(k>func_asize(f)){    
    do{
      g=func_new(op);
      func_a_resize(g,k);
      if(func_has_power(g)){ func_set_power(g,func_power(f)); }
      for(k=0, i=0; i<func_asize(f); i++){
	if((func_has_power(f) && func_is(f->a[i],op) && func_power(f->a[i])==1) || (!func_has_power(f) && func_is(f->a[i],op))){
	  for(j=0; j<func_asize(f->a[i]); j++){
	    g->a[k++]=FR(f->a[i]->a[j]);
	  }
	}else{ g->a[k++]=FR(f->a[i]); }
      }
      k=func_args_count_op(g,op);
      if(k>func_asize(g)){ f=func_del(f); f=g; redo=1; }else{ redo=0; }
    }while(redo);
  }else{ g=FR(f); }
  f=func_del(f);
  return g;
}

///private//////////////////////////////////////////////////

func_t *func_arg1_new(char *op, func_t *g)
{
  func_t *f=NULL;
  f=func_new(op);
  func_a_resize(f,1);
  f->a[0]=g;
  return f;
}

func_t *func_arg2_new(char *op, func_t *g0, func_t *g1)
{
  func_t *f=NULL;
  f=func_new(op);
  func_a_resize(f,2);
  f->a[0]=g0;
  f->a[1]=g1;
  return f;
}

/////////////////////////////////////////////////////////////////

int func_size(func_t *f)
{
  if     (func_is_list(f))  { return func_asize(f); }
  else if(func_is(f,"rvec")){ return func_rvec_size(f); }
  else if(func_is(f,"cvec")){ return func_cvec_size(f); }
  else                      { return 1; }
}

//////////////////////////////////////////////

func_t *func_eval_eval(func_t *g)
{
  func_t *f=NULL;
  if(!func_is(g,"eval") || func_asize(g)!=1){ FUNC_ERROR_ARG1("func_eval_eval",g); }
  f=func_eval(FR(func_aget(g,0)));
  g=func_del(g);
  return f;
}

func_t *func_eval(func_t *f)
{
  int i;
  func_t *g=NULL;
  func_arg1_t *eval=NULL;
  if(f!=NULL){
    g=func_clone(FR(f));
    for(i=0; i<func_asize(f); i++){
      func_aset(g,i,func_eval(FR(f->a[i])));
    }
    eval=func_find_eval(g);
    if(eval!=NULL){ g=eval(g); }
  }
  f=func_del(f);
  return g;
}

func_t *func_op_eval_new(void)
{
  func_t *f=func_builtin_new("eval");
  func_builtin_p(f)->amin=1;
  func_builtin_p(f)->amax=1;
  func_builtin_p(f)->eval=func_eval_eval;
  return f;
}

///////////////////////////////////////////////////////////////////////////////


//EOF
