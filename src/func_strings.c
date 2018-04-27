#include<stdio.h>
#include<stdlib.h>
#include"is_macros.h"
#include"is_strings.h"
#include"is_func.h"

#define SKIP_SPC " \t\n"

////////////////////////////////

static const char *__func_strings="@S";

////////////////////////////////

func_t *func_op_strings_new(void)
{
  func_t *f=func_builtin_new(__func_strings);
  func_builtin_p(f)->order=FUNC_ORDER_STRINGS;
  func_builtin_p(f)->ptype=FUNC_P_STRINGS;
  func_builtin_p(f)->p_new=func_strings_p_new;
  func_builtin_p(f)->p_del=func_strings_p_del;
  func_builtin_p(f)->p_clone=func_strings_p_clone;
  func_builtin_p(f)->p_cmp=func_strings_p_cmp;
  func_builtin_p(f)->print=func_strings_print;
  return f;
}

/////////////////////////////////////

void func_strings_p_new(func_t *f)
{
  if(f->p.s!=NULL){ func_strings_p_del(f); }
  f->ptype=FUNC_P_STRINGS;
  f->p.s=strings_new(0);
}

void func_strings_p_del(func_t *f)
{
  if(func_ptype(f)!=FUNC_P_STRINGS){ FUNC_ERROR_ARG1("func_strings_p_del",f); }
  if(f->p.s!=NULL){
    f->ptype=FUNC_P_NULL;
    f->p.s=strings_del(f->p.s);
  }
}

void func_strings_p_clone(func_t *f, func_t *g)
{
  if(func_ptype(g)!=FUNC_P_STRINGS){ FUNC_ERROR_ARG2("func_strings_p_clone",f,g); }
  if(f->p.s!=NULL){ func_strings_p_del(f); }
  f->ptype=FUNC_P_STRINGS;
  f->p.s=strings_new_clone(g->p.s);
}

int func_strings_p_cmp(func_t *f, func_t *g)
{
  return strings_cmp(f->p.s,g->p.s);
}

/////////////////////////////////////

void func_strings_p_replace(func_t *f, strings *s)
{
  if(func_ptype(f)!=FUNC_P_STRINGS){ FUNC_ERROR_ARG1("func_strings_p_replace",f); }
  func_strings_p_del(f);
  f->ptype=FUNC_P_STRINGS;
  f->p.s=s;
}

/////////////////////////////////////

strings *func_strings_p(func_t *f)
{
  if(f==NULL || func_ptype(f)!=FUNC_P_STRINGS){ FUNC_ERROR_ARG1("func_strings_p",f); }
  return f->p.s;
}

int func_strings_size(func_t *f)
{
  if(f==NULL || func_ptype(f)!=FUNC_P_STRINGS || f->p.s==NULL){ FUNC_ERROR_ARG1("func_strings_size",f); }
  return f->p.s->n;
}

char *func_strings_at(func_t *f, int i)
{
  if(f==NULL || func_ptype(f)!=FUNC_P_STRINGS || f->p.s==NULL || f->p.s->str==NULL){ FUNC_ERROR_ARG1("func_strings_at",f); }
  return f->p.s->str[i];
}


void func_strings_set(func_t *f, int i, const char *str)
{
  if(f==NULL || func_ptype(f)!=FUNC_P_STRINGS || f->p.s==NULL){ FUNC_ERROR_ARG1("func_strings_set",f); }
  strings_item_set(f->p.s,i,str,SKIP_SPC);
  return;
}

/////////////////////////////////////

func_t *func_strings_split(const char *str, const char *sep, const char *mask_begin, const char *mask_end, const char *skip)
{
  func_t *f=NULL;
  strings *a=NULL;
  a=strings_split(str,sep,mask_begin,mask_end,skip);
  f=func_strings(0);
  func_strings_p_replace(f,a);
  return f;
}

/////////////////////////////////////

void func_strings_print(func_t *f)
{
  int i;
  printf("[");
  for(i=0; i<func_strings_size(f); i++){
    printf("'%s'",func_strings_at(f,i));
    if(i<func_strings_size(f)-1) printf(" ");
  }
  printf("]");
}

/////////////////////////////////////

int func_is_strings(func_t *f)
{
  return func_is(f,__func_strings);
}

func_t *func_strings(int n)
{
  func_t *f=NULL;
  f=func_new(__func_strings);
  func_strings_p_replace(f,strings_new(n));
  return f;
}

func_t *func_strings_str(const char *str[])
{
  func_t *f=NULL;
  f=func_new(__func_strings);
  func_strings_p_replace(f,strings_new_str(str,NULL));
  return f;
}

func_t *func_strings_strings(strings *str)
{
  func_t *f=NULL;
  f=func_new(__func_strings);
  func_strings_p_replace(f,strings_new_clone(str));
  return f;
}

func_t *func_strings_char(const char *str)
{
  func_t *f=NULL;
  f=func_new(__func_strings);
  func_strings_set(f,0,str);
  return f;
}

////////////////////////////////

//EOF
