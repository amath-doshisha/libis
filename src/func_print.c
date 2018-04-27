#include<stdio.h>
#include<stdlib.h>
#include"is_macros.h"
#include"is_strings.h"
#include"is_print.h"
#include"is_func.h"
#include"is_rmulti.h"
#include"is_cmulti.h"
#include"is_bigint.h"
#include"is_ivec.h"
#include"is_rvec.h"
#include"is_cvec.h"

#define FR(F) func_retain(F)

////////////////////////////////////////////////////////////////

void func_put(func_t *f)
{  
  int i;
  if(f==NULL){ printf("[NULL]"); return; }
  printf("[%s",func_op(f));
  if(func_has_power(f)){
    printf("^%d",func_power(f));
  }
  if((func_asize(f))>0){
    for(i=0; i<func_asize(f); i++){
      printf(" ");
      func_put(func_aget(f,i));
    }
  }
  printf("]");
}

void func_print_hat(int pow)
{
  if(pow<0){
    printf("^(%d)",pow);
  }else if(pow!=1){
    printf("^%d",pow);
  }
}


///////////////////////////////////////////////////////////


void func_print_args(func_t *f)
{
  int i;
  printf("%s",func_op(f));
  if(func_has_power(f)){ func_print_hat(func_power(f)); }
  if(func_asize(f)>0){
    printf("(");
    for(i=0; i<func_asize(f); i++){
      func_print(func_aget(f,i));
      if(i<func_asize(f)-1){ printf(","); }
    }
    printf(")");
  }
}

///////////////////////////////////////////////////////////

void func_print(func_t *f)
{
  func_set_t *print=NULL;
  if     (f==NULL)          { printf("NULL"); }
  else if(func_is(f,"null")){ printf("null"); }
  else{
    print=func_find_print(f);
    if(print!=NULL){ print(f); }
    else{
      func_print_args(f);
      if(func_ptype(f)!=FUNC_P_NULL){ printf("_%d",func_ptype(f)); }
    }
  }
}

//////////////////////////////////

func_t *func_print_eval(func_t *f)
{
  if(f==NULL || !func_is(f,"print") || func_asize(f)!=1){ FUNC_ERROR_ARG1("func_print_eval",f); }
  func_print(func_aget(f,0));
  f=func_del(f);
  return NULL;
}

////////////////////////////////////////////////

//EOF
