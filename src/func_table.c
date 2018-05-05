#include<stdio.h>
#include<stdlib.h>
#include"is_macros.h"
#include"is_func.h"
#include"GeneralHashFunctions.h"

#define FR(f) func_retain(f)

#define TABLE_SIZE_LIST_NUM 12
int __func_table_size[]={ 97,241,397,499,743,997,1499,1999,3989,4999,7499,9973 };

/////////////////////////////////////////////////////////////////////////

static char *__func_table="@T";

int func_is_table(func_t *f)
{
  return func_is(f,__func_table);
}

/////////////////////////////////////////////////////////////////////////

func_t *func_table_find(func_t *f, char *name)
{
  int i;
  if(f==NULL || !func_is_table(f)){ FUNC_ERROR_ARG1("func_table_find",f); }
  i=func_table_key_index(f,name);
  if(i>=0 && i<func_asize(f)){ return func_aget(f,i); }
  return NULL;
}

void func_table_set(func_t *f, func_t *g)
{
  int i;
  if(f==NULL || !func_is_table(f) || func_key_name(g)==NULL){ FUNC_ERROR_ARG2("func_table_set",f,g); }
  i=func_table_key_index(f,func_key_name(g));
  if(i<0 || i>=func_asize(f)){
    func_table_grow(f);
    i=func_table_key_index(f,func_key_name(g));
  }
  func_aset(f,i,g);
}

/////////////////////////////////////////////////////////////////////////

void func_table_grow(func_t *f)
{
  int i,k;
  func_t *g=NULL;
  if(f==NULL || !func_is_table(f)){ FUNC_ERROR_ARG1("func_table_grow",f); }
  g=func_table_bigger(f);
  for(i=0; i<func_asize(f); i++){
    if(func_aget(f,i)!=NULL && func_key_name(func_aget(f,i))!=NULL){
      k=func_table_key_index(g,func_key_name(func_aget(f,i)));
      if(k<0 || k>=func_asize(g)){ FUNC_ERROR_ARG1("func_table_grow",f); }
      func_aset(g,k,FR(func_aget(f,i)));
    }
  }
  func_a_swap(f,g);
  g=func_del(g);
}

int func_table_bigger_size(func_t *f)
{
  int i;
  if(!func_is_table(f)){ FUNC_ERROR_ARG1("func_table_bigger_size",f); }
  for(i=0; i<TABLE_SIZE_LIST_NUM; i++){
    if(func_asize(f)<__func_table_size[i]){ return __func_table_size[i]; }
  }
  return 2*(func_asize(f))-1;
}

/////////////////////////////////////////////////////////////////////////

int func_table_key_index(func_t *f, char *name)
{
  int i,k,n;
  if(f==NULL || !func_is_table(f) || func_asize(f)<=0 || name==NULL){ FUNC_ERROR_ARG1("func_table_key_index",f); }
  n=PJWHash(name);
  i=(n % func_asize(f));
  if(f->a[i]==NULL || char_eq(func_key_name(f->a[i]),name)){ return i; }
  k=i;
  do{
    if(++i>=func_asize(f)){ i=0; }
    if(f->a[i]==NULL || char_eq(func_key_name(f->a[i]),name)){ return i; }
  }while(k!=i);
  return -1;
}

/////////////////////////////////////////////////////////////////////////

char *func_key_name(func_t *f)
{
  if     (f==NULL)           { return NULL; }
  else if(func_is_builtin(f)){ return func_builtin_name(f); }
  else if(func_is_def(f))    { return func_def_name(f); }
  else                       { return NULL; }
}

/////////////////////////////////////////////////////////////////////////

func_t *func_table(void)
{
  func_t *f=NULL;
  f=func_new(__func_table);
  func_a_resize(f,__func_table_size[0]);
  return f;
}

func_t *func_table_bigger(func_t *g)
{
  func_t *f=NULL;
  f=func_new(__func_table);
  func_a_resize(f,func_table_bigger_size(g));
  return f;
}


////////////////////////////////////////////////

void func_print_table(func_t *f)
{
  int i;
  printf("{\n");
  for(i=0; i<func_asize(f); i++){
    if(func_aget(f,i)!=NULL){      
      printf("[%d] ",i);
      //      if(func_key_name(func_aget(f,i))==NULL){ printf("NULL"); }
      //      else{ printf("%s",func_key_name(func_aget(f,i))); }
      //      printf(" => ");
      func_print(func_aget(f,i));
      printf("\n");
    }
  }
  printf("}");
}

////////////////////////////////////////////////

func_t *func_op_table_new(void)
{
  func_t *f=func_builtin_new(__func_table);
  func_builtin_p(f)->order=FUNC_ORDER_TABLE;
  func_builtin_p(f)->amax=FUNC_NOLIMIT;
  func_builtin_p(f)->print=func_print_table;
  return f;
}

//EOF
