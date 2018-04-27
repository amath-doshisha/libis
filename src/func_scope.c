#include<stdio.h>
#include<stdlib.h>
#include"is_macros.h"
#include"is_ivec.h"
#include"is_strings.h"
#include"is_print.h"
#include"is_func.h"
#include"is_strings.h"
#include"GeneralHashFunctions.h"

#define FR(f) func_retain(f)

///////////////////////////////////////////////

static func_t *__func_scope_current=NULL;

///////////////////////////////////////////////

func_t *func_scope_find(int n, const char *name)
{
  func_t *f=NULL,*scope=NULL,*table=NULL;
  func_init();
  if(!func_is_scope(__func_scope_current)){ ERROR_EXIT("Error in func_scope_find.\nThe global variable '__func_scope_current' is destroyed.\n%s",""); }
  scope=func_scope(n);
  do{
    if(!func_is_scope(scope)){ ERROR_EXIT("Error in func_scope_find.\nThe global variable '__func_scope_current' is destroyed.\n%s",""); }
    table=func_aget(scope,2);
    f=func_table_find(table,name);
    scope=func_aget(scope,0);
  }while(f==NULL && scope!=NULL);
  return f;
}

///////////////////////////////////////////////

void func_scope_set(int n, func_t *f)
{
  int i;
  func_t *scope=NULL,*vlist=NULL,*table=NULL;
  func_init();
  if(!func_is_scope(__func_scope_current)){ ERROR_EXIT("Error in func_scope_set.\nThe global variable '__func_scope_current' is destroyed.\n%s",""); }
  scope=func_scope(n);
  if(!func_is_scope(scope)){ ERROR_EXIT("Error in func_scope_set.\nThe global variable '__func_scope_current' is destroyed.\n%s",""); }
  vlist=func_aget(scope,1);
  table=func_aget(scope,2);
  if(func_key_name(f)==NULL){ FUNC_ERROR_ARG1("func_scope_set",f); }
  if((i=strings_index(func_strings_p(vlist),func_key_name(f)))>=0){
    printf("Error in func_scope_set(func_t *f) with f="); func_print(f); printf(".\n");
    printf("The f is already set as the independent variable '#%d'.\n",i);
    exit(0);
  }
  func_table_set(table,f);
}

///////////////////////////////////////////////

void func_init(void)
{
  if(__func_scope_current==NULL){
    __func_scope_current=func_scope_new(NULL,NULL);
    func_import_basic(0);
  }
}

void func_clear(void)
{
  func_t *previous=NULL;
  while(__func_scope_current!=NULL){
    if(!func_is_scope(__func_scope_current)){ ERROR_EXIT("Error in func_clear.\nThe global variable '__func_scope_end' is destroyed.\n%s",""); }
    // delete table except for bulit-in op
    func_a_del_not_op_op(__func_scope_current->a[2],func_is_builtin,func_is_def);
    func_a_del_not_op(__func_scope_current->a[2],func_is_builtin);
    // delete scope
    previous=func_aget(__func_scope_current,0);   
    __func_scope_current->a[0]=NULL;
    __func_scope_current=func_del(__func_scope_current);
    // go to previous
    __func_scope_current=previous;

  }
}

///////////////////////////////////////////////

void func_scope_init_vlist_to_table(int n)
{
  int i;
  func_t *scope=NULL,*table=NULL,*vlist=NULL;
  func_init();
  scope=func_scope(n);
  if(!func_is_scope(scope)){ ERROR_EXIT("Error in func_scope_init_vlist_to_table(int n).\nThe global variable '__func_scope_current' is destroyed.\n%s",""); }
  vlist=func_aget(scope,1);
  table=func_aget(scope,2);
  for(i=0; i<func_strings_size(vlist); i++){
    func_table_set(table,func_def(func_strings_at(vlist,i),func_var1(i,1),0,0));
  }
}

///////////////////////////////////////////////

func_t *func_scope_table(int n)
{
  func_t *scope=NULL;
  func_init();
  if(!func_is_scope(__func_scope_current)){ ERROR_EXIT("Error in func_scope_table(int n).\nThe global variable '__func_scope_current' is destroyed.\n%s",""); }
  scope=func_scope(n);
  if(!func_is_scope(scope)){ ERROR_EXIT("Error in func_scope_table(int n).\nThe global variable '__func_scope_current' is destroyed.\n%s",""); }
  return func_aget(scope,2);
}

func_t *func_scope_vlist(int n)
{
  func_t *scope=NULL;
  func_init();
  if(!func_is_scope(__func_scope_current)){ ERROR_EXIT("Error in func_scope_vlist(int n).\nThe global variable '__func_scope_current' is destroyed.\n%s",""); }
  scope=func_scope(n);
  if(!func_is_scope(scope)){ ERROR_EXIT("Error in func_scope_vlist(int n).\nThe global variable '__func_scope_current' is destroyed.\n%s",""); }
  return func_aget(scope,1);
}

func_t *func_scope_previous(int n)
{
  func_t *scope=NULL;
  func_init();
  if(!func_is_scope(__func_scope_current)){ ERROR_EXIT("Error in func_scope_previous(int n).\nThe global variable '__func_scope_current' is destroyed.\n%s",""); }
  scope=func_scope(n);
  if(!func_is_scope(scope)){ ERROR_EXIT("Error in func_scope_previous(int n).\nThe global variable '__func_scope_current' is destroyed.\n%s",""); }
  return func_aget(scope,0);
}

func_t *func_scope(int n)
{
  int i;
  func_t *f=NULL;
  func_init();
  if(!func_is_scope(__func_scope_current)){ ERROR_EXIT("Error in func_scope(int n).\nThe global variable '__func_scope_current' is destroyed.\n%s",""); }
  n=abs(n);
  f=__func_scope_current;
  for(i=0; i<n; i++){
    f=func_aget(f,0);
    if(f==NULL || !func_is_scope(f)){ ERROR_EXIT("Error in func_scope(int n).\nNo scope with n=%d.\n",n); }
  }
  return f;
}

///////////////////////////////////////////////

void func_scope_end(void)
{
  func_t *previous=NULL;
  func_init();
  if(!func_is_scope(__func_scope_current)){ ERROR_EXIT("Error in func_scope_end().\nThe global variable '__func_scope_end' is destroyed.\n%s",""); }
  previous=__func_scope_current->a[0];
  if(previous==NULL){ ERROR_EXIT("Error in func_scope_end().\nThe default scope cannot be cleared.\n%s",""); }
  // delete table except for bulit-in op
  func_a_del_not_op_op(__func_scope_current->a[2],func_is_builtin,func_is_def);
  func_a_del_not_op(__func_scope_current->a[2],func_is_builtin);
  // delete scope
  __func_scope_current->a[0]=NULL;
  __func_scope_current=func_del(__func_scope_current);
  // change current scope
  __func_scope_current=previous;
}

void func_scope_begin(func_t *vlist)
{
  func_init();
  if(!func_is_scope(__func_scope_current)){ ERROR_EXIT("Error in func_scope_begin(strings).\nThe global variable '__func_scope_current' is destroyed.\n%s",""); }
  if(!(vlist==NULL || func_is_strings(vlist))){ FUNC_ERROR_ARG1("func_scope_begin",vlist); }
  __func_scope_current=func_scope_new(__func_scope_current,vlist);
  func_scope_init_vlist_to_table(0);
}


///////////////////////////////////////////////

void func_print_scope(func_t *f)
{
  printf("scope_vlist="); func_print(func_aget(f,1)); printf("\n"); 
  printf("scope_table="); func_print(func_aget(f,2)); printf("\n");
  if(func_aget(f,0)==NULL){ printf("scope_previous=NULL\n"); }
  else{ printf("scope_previous=0x%lx\n",(unsigned long int)func_aget(f,0)); func_print(func_aget(f,0)); }
}

///////////////////////////////////////////////

static const char *__func_scope="@s";

///////////////////////////////////////////////

func_t *func_op_scope_new(void)
{
  func_t *f=func_builtin_new(__func_scope);
  func_builtin_p(f)->order=FUNC_ORDER_SCOPE;
  func_builtin_p(f)->amin=3;
  func_builtin_p(f)->amax=3;
  func_builtin_p(f)->print=func_print_scope;
  return f;
}

func_t *func_scope_new(func_t *previous, func_t *vlist)
{
  func_t *f=NULL;
  f=func_new(__func_scope);
  if(!(previous==NULL || func_is_scope(previous)) || !(vlist==NULL || func_is_strings(vlist))){ FUNC_ERROR_ARG2("func_scope_new",previous,vlist); }
  func_aset(f,0,previous);
  if(vlist==NULL){ func_aset(f,1,func_strings(0)); }else{ func_aset(f,1,vlist); }
  func_aset(f,2,func_table());
  return f;
}

int func_is_scope(func_t *f)
{
  return func_is(f,__func_scope);
}


//EOF
