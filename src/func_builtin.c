#include<stdio.h>
#include<stdlib.h>
#include"is_macros.h"
#include"is_func.h"


////////////////////////////////////////

#define SPC " \t\n"
#define BLK0 "({["
#define BLK1 ")}]"

func_t *func_builtin_script(func_t *f)
{
  int i;
  strings *a=NULL,*b=NULL;
  func_t *g=NULL;
  if(func_is_strings(f) && func_strings_size(f)==1){
    if(str_match(func_strings_at(f,0),"^%A%I*$",BLK0,BLK1)){
      // f
      g=func_scope_find(0,func_strings_at(f,0));
      if(func_is_builtin(g) && func_builtin_name(g)!=NULL){
	g=func_new(func_builtin_name(g));
	if(func_find_amin(g)!=0){ g=func_del(g); }
      }else{ g=NULL; }
    }else if(str_match(func_strings_at(f,0),"^%A%I*(%m+)$",BLK0,BLK1)){
      // f(.....) 
      a=strings_split_mask(func_strings_at(f,0),BLK0,BLK1,SPC);
      if(strings_size(a)==2){
	g=func_scope_find(0,strings_at(a,0));
	b=strings_split(strings_at(a,1),",",BLK0,BLK1,SPC);
	if(func_is_builtin(g) && func_builtin_name(g)!=NULL && char_eq(func_builtin_name(g),"begin")){
	  g=func_begin(func_strings_strings(b));
	}else if(func_is_builtin(g) && func_builtin_name(g)!=NULL){
	  g=func_new(func_builtin_name(g));
	  if(func_find_amin(g)<=strings_size(b) && (strings_size(b)<=func_find_amax(g) || func_find_amax(g)!=FUNC_NOLIMIT)){
	    for(i=0; i<strings_size(b); i++){
	      if(func_asize(g)<=i+1){ func_a_resize(g,i+1); }
	      func_aset(g,i,func_script(strings_at(b,i)));
	    }
	  }else{ g=func_del(g); }
	}else{ g=NULL; }
      }
    }
  }
  if(g==NULL){ g=f; }else{ f=func_del(f); }
  a=strings_del(a);
  b=strings_del(b);
  return g;
}


/////////////////////////////////////////////////////////////

static const char *__func_builtin="@B";

/////////////////////////////////////////////////////////////

func_t *func_op_builtin_new(void)
{
  func_t *f=func_builtin_new(__func_builtin);
  func_builtin_p(f)->order=FUNC_ORDER_BUILTIN;
  func_builtin_p(f)->ptype=FUNC_P_BUILTIN;
  func_builtin_p(f)->p_new=func_builtin_p_new;
  func_builtin_p(f)->p_del=func_builtin_p_del;
  func_builtin_p(f)->p_clone=func_builtin_p_clone;
  func_builtin_p(f)->p_cmp=func_builtin_p_cmp;
  func_builtin_p(f)->print=func_builtin_print;
  return f;
}

/////////////////////////////////////////////////////////////

int func_is_builtin(func_t *f)
{
  return func_is(f,__func_builtin);
}

func_t *func_builtin_new(const char *name)
{
  func_t *f=NULL;
  f=func_new(__func_builtin);
  f->p.builtin->name=char_new(name," \t\n");
  return f;
}

void func_builtin_p_new(func_t *f)
{
  if(f->p.builtin!=NULL){ func_builtin_p_del(f); }
  f->ptype=FUNC_P_BUILTIN;
  f->p.builtin=func_builtin_struct_new();
}

void func_builtin_p_del(func_t *f)
{
  if(func_ptype(f)!=FUNC_P_BUILTIN){ FUNC_ERROR_ARG1("func_builtin_p_del",f); }
  if(f->p.builtin!=NULL){
    f->p.builtin=func_builtin_struct_del(f->p.builtin);
    f->ptype=FUNC_P_NULL;
  }
}

void func_builtin_p_clone(func_t *f, func_t *g)
{
  if(func_ptype(g)!=FUNC_P_BUILTIN){ FUNC_ERROR_ARG2("func_builtin_p_clone",f,g); }
  if(f->p.builtin!=NULL){ func_builtin_p_del(f); }
  f->ptype=FUNC_P_BUILTIN;
  f->p.builtin=func_builtin_struct_new();
  (*(f->p.builtin))=(*(g->p.builtin));
  f->p.builtin->name=char_new(g->p.builtin->name," \t\n");
}

//////////////////////////////////////////////////////////////////////////////////////////////

int func_builtin_p_cmp(func_t *f, func_t *g)
{
  if(f==NULL || g==NULL || func_ptype(f)!=FUNC_P_BUILTIN || func_ptype(g)!=FUNC_P_BUILTIN || f->p.builtin==NULL || g->p.builtin==NULL){ FUNC_ERROR_ARG2("func_builtin_cmp",f,g); }
  return char_cmp(f->p.builtin->name,g->p.builtin->name);
}
////////////////////////////////////////////////////////////////////////////////////////////

func_builtin_t *func_builtin_p(func_t *f)
{
  if(f==NULL || func_ptype(f)!=FUNC_P_BUILTIN){ FUNC_ERROR_ARG1("func_builtin_p",f); }
  return f->p.builtin;
}

char *func_builtin_name(func_t *f)
{
  if(f==NULL || func_ptype(f)!=FUNC_P_BUILTIN || f->p.builtin==NULL){ FUNC_ERROR_ARG1("func_builtin_name",f); }
  return f->p.builtin->name;
}

int func_builtin_order(func_t *f)
{
  if(f==NULL || func_ptype(f)!=FUNC_P_BUILTIN || f->p.builtin==NULL){ FUNC_ERROR_ARG1("func_builtin_order",f); }
  return f->p.builtin->order;
}

int func_builtin_amin(func_t *f)
{
  if(f==NULL || func_ptype(f)!=FUNC_P_BUILTIN || f->p.builtin==NULL){ FUNC_ERROR_ARG1("func_builtin_amin",f); }
  return f->p.builtin->amin;
}

int func_builtin_amax(func_t *f)
{
  if(f==NULL || func_ptype(f)!=FUNC_P_BUILTIN || f->p.builtin==NULL){ FUNC_ERROR_ARG1("func_builtin_amax",f); }
  return f->p.builtin->amax;
}

int func_builtin_type(func_t *f)
{
  if(f==NULL || func_ptype(f)!=FUNC_P_BUILTIN || f->p.builtin==NULL){ FUNC_ERROR_ARG1("func_builtin_type",f); }
  return f->p.builtin->type;
}

int func_builtin_ptype(func_t *f)
{
  if(f==NULL || func_ptype(f)!=FUNC_P_BUILTIN || f->p.builtin==NULL){ FUNC_ERROR_ARG1("func_builtin_p_p",f); }
  return f->p.builtin->ptype;
}

func_set_t *func_builtin_p_new_call(func_t *f)
{
  if(f==NULL || func_ptype(f)!=FUNC_P_BUILTIN || f->p.builtin==NULL){ FUNC_ERROR_ARG1("func_builtin_p_new_call",f); }
  return f->p.builtin->p_new;
}

func_set_t *func_builtin_p_del_call(func_t *f)
{
  if(f==NULL || func_ptype(f)!=FUNC_P_BUILTIN || f->p.builtin==NULL){ FUNC_ERROR_ARG1("func_builtin_p_del_call",f); }
  return f->p.builtin->p_del;
}

func_set2_t *func_builtin_p_clone_call(func_t *f)
{
  if(f==NULL || func_ptype(f)!=FUNC_P_BUILTIN || f->p.builtin==NULL){ FUNC_ERROR_ARG1("func_builtin_p_clone_call",f); }
  return f->p.builtin->p_clone;
}

func_arg1_t *func_builtin_eval_call(func_t *f)
{
  if(f==NULL || func_ptype(f)!=FUNC_P_BUILTIN || f->p.builtin==NULL){ FUNC_ERROR_ARG1("func_builtin_eval_call",f); }
  return f->p.builtin->eval;
}

func_set_t *func_builtin_print_call(func_t *f)
{
  if(f==NULL || func_ptype(f)!=FUNC_P_BUILTIN || f->p.builtin==NULL){ FUNC_ERROR_ARG1("func_builtin_print_call",f); }
  return f->p.builtin->print;
}

func_cmp_t *func_builtin_cmp_call(func_t *f)
{
  if(f==NULL || func_ptype(f)!=FUNC_P_BUILTIN || f->p.builtin==NULL){ FUNC_ERROR_ARG1("func_builtin_cmp_call",f); }
  return f->p.builtin->p_cmp;
}

/////////////////////////////////////////////////////////////

func_builtin_t *func_builtin_struct_new(void)
{
  func_builtin_t *builtin=NULL;
  builtin=(func_builtin_t*)malloc(sizeof(func_builtin_t)*1);
  builtin->name=NULL;
  builtin->order=FUNC_ORDER_DEFAULT;
  builtin->amin=0;
  builtin->amax=0;
  builtin->type=FUNC_SCALAR_NO;
  builtin->ptype=FUNC_P_NULL;
  builtin->p_new=NULL;
  builtin->p_del=NULL;
  builtin->p_clone=NULL;
  builtin->p_cmp=NULL;
  builtin->eval=NULL;
  builtin->print=NULL;
  return builtin;
}

func_builtin_t *func_builtin_struct_del(func_builtin_t *builtin)
{
  if(builtin!=NULL){
    builtin->name=char_del(builtin->name);
    free(builtin);
  }
  return NULL;
}


/////////////////////////////////////////////////////////////

void func_builtin_print(func_t *f)
{
  int i;
  void *p=NULL;
  char *ptype[]= { "FUNC_P_NULL", "FUNC_P_BUILTIN", "FUNC_P_DEF", "FUNC_P_POWER",
		   "FUNC_P_BIGINT", "FUNC_P_REAL", "FUNC_P_COMPLEX", "FUNC_P_STRINGS",
		   "FUNC_P_VAR", "FUNC_P_IVEC", "FUNC_P_RVEC", "FUNC_P_CVEC", "FUNC_P_RMAT", "FUNC_P_CMAT" };
  
  printf("%s",func_builtin_name(f));
  if(!(func_builtin_amin(f)==0 && func_builtin_amax(f)==0)){
    printf("(");
    for(i=0; i<MAX2(func_builtin_amin(f),func_builtin_amax(f)); i++){
      if(i!=0){ if(i==func_builtin_amin(f)){ printf(";"); }else{ printf(","); } }
      printf("#%d",i);
    }
    if(func_builtin_amax(f)<0){
      if(func_builtin_amin(f)>0){ printf(","); }
      printf("...");
    }
    printf(")");
  }
  printf("=built-in operator,");
  if(func_builtin_order(f)<0){ printf(" order=any"); }
  else                       { printf(" order=%d",func_builtin_order(f)); }  
  if     (func_builtin_type(f)==FUNC_SCALAR_NO)    { printf(" type=scalar_no"); }
  else if(func_builtin_type(f)==FUNC_SCALAR)       { printf(" type=scalar"); }
  else if(func_builtin_type(f)==FUNC_SCALAR_DEPEND){ printf(" type=scalar_depend"); }
  else if(func_builtin_type(f)==FUNC_COMMAND)      { printf(" type=command"); }
  else                                             { printf(" type=%d",func_builtin_type(f)); }
  if(func_builtin_ptype(f)<=FUNC_P_NULL || func_builtin_ptype(f)>=FUNC_P_SIZE){ printf(" p=no"); }
  else{ printf(" p=%s",ptype[func_builtin_ptype(f)]); }
  if(func_builtin_ptype(f)!=FUNC_P_NULL){
    p=(void*)func_builtin_p_new_call(f);   printf(" p_new=0x%lx",  (unsigned long int)p);
    p=(void*)func_builtin_p_del_call(f);   printf(" p_del=0x%lx",  (unsigned long int)p);
    p=(void*)func_builtin_p_clone_call(f); printf(" p_clone=0x%lx",(unsigned long int)p);
  }
  p=(void*)func_builtin_eval_call(f);    printf(" eval=0x%lx", (unsigned long int)p);
  p=(void*)func_builtin_print_call(f);   printf(" print=0x%lx",(unsigned long int)p);
}


//EOF
