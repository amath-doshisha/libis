#include<stdio.h>
#include<stdlib.h>
#include"is_macros.h"
#include"is_strings.h"
#include"is_func.h"

#define FR(f) func_retain(f)

///////////////////////////////////////////////

static const char *__func_def="@D";


////////////////////////////////////////

func_t *func_def_any_eval(func_t *f)
{
  int i;
  func_t *g=NULL,*x=NULL;
  g=func_scope_find(0,func_op(f));
  if(func_is_def(g) && func_aget(g,0)!=NULL){
    g=FR(func_aget(g,0));
    if(func_asize(f)>0){
      x=func_list(func_asize(f));
      for(i=0; i<func_asize(f); i++){ func_aset(x,i,FR(func_aget(f,i))); }
      g=func_maps(g,0,FR(x));
    }
    f=func_del(f);
  }else{ g=f; }
  x=func_del(x);
  return g;
}

////////////////////////////////////////

#define SPC " \t\n"
#define BLK0 "({["
#define BLK1 ")}]"

func_t *func_def_script(func_t *f)
{
  int i;
  strings *a=NULL,*b=NULL;
  func_t *g=NULL;
  if(func_is_strings(f) && func_strings_size(f)==1){
    if(str_match(func_strings_at(f,0),"^%A%I*$",BLK0,BLK1)){
      // f
      g=func_scope_find(0,func_strings_at(f,0));
      if(func_is_def(g) && func_def_name(g)!=NULL){
	g=func_new(func_def_name(g));
	if(func_find_amin(g)>0){ g=func_del(g); }
      }else{ g=NULL; }
    }else if(str_match(func_strings_at(f,0),"^%A%I*(%m+)$",BLK0,BLK1)){
      // f(.....) 
      a=strings_split_mask(func_strings_at(f,0),BLK0,BLK1,SPC);
      if(strings_size(a)==2){
	g=func_scope_find(0,strings_at(a,0));
	b=strings_split(strings_at(a,1),",",BLK0,BLK1,SPC);
	if(func_is_def(g) && func_def_name(g)!=NULL){
	  g=func_new(func_def_name(g));
	  if(func_find_amin(g)<=strings_size(b) && (strings_size(b)<=func_find_amax(g) || func_find_amax(g)==FUNC_NOLIMIT)){
	    func_a_resize(g,strings_size(b));
	    for(i=0; i<strings_size(b); i++){
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

////////////////////////////////////////

void func_def_print(func_t *f)
{
  func_t *vars=NULL;
  int i;
  vars=func_scope_vlist(0);
  printf("%s",func_def_name(f));
  if(!(func_def_amin(f)==0 && func_def_amax(f)==0)){
    printf("(");
    for(i=0; i<MAX2(func_def_amin(f),func_def_amax(f)); i++){
      if(i!=0){ if(i==func_def_amin(f)){ printf(";"); }else{ printf(","); } }
      if(vars!=NULL && i>=0 && i<(func_strings_size(vars))){
	printf("%s",func_strings_at(vars,i));
      }else{
	printf("#%d",i);
      }
    }
    if(func_def_amax(f)<0){
      if(func_def_amin(f)>0){ printf(","); }
      printf("...");
    }
    printf(")");
  }
  printf("=");
  func_print(func_aget(f,0));
}

////////////////////////////////////////

func_t *func_op_def_new(void)
{
  func_t *f=func_builtin_new(__func_def);
  func_builtin_p(f)->order=FUNC_ORDER_DEF;
  func_builtin_p(f)->amin=1;
  func_builtin_p(f)->amax=1;
  func_builtin_p(f)->ptype=FUNC_P_DEF;
  func_builtin_p(f)->p_new=func_def_p_new;
  func_builtin_p(f)->p_del=func_def_p_del;
  func_builtin_p(f)->p_clone=func_def_p_clone;
  func_builtin_p(f)->p_cmp=func_def_p_cmp;
  func_builtin_p(f)->print=func_def_print;
  return f;
}

////////////////////////////////////////

// def
func_t *func_def(const char *name, func_t *g, int amin, int amax)
{
  func_t *f=NULL;
  f=func_new(__func_def);
  f->p.def->name=char_new(name," \t\n");
  f->p.def->amin=amin;
  f->p.def->amax=amax;
  func_aset(f,0,g);
  return f;
}

//////////////////////

int func_is_def(func_t *f)
{
  return func_is(f,__func_def);
}


////////////////////////////////////////

void func_def_p_new(func_t *f)
{
  if(f->p.def!=NULL){ func_def_p_del(f); }
  f->ptype=FUNC_P_DEF;
  f->p.def=func_def_struct_new();
}

void func_def_p_del(func_t *f)
{
  if(func_ptype(f)!=FUNC_P_DEF){ FUNC_ERROR_ARG1("func_def_p_del",f); }
  if(f->p.def!=NULL){
    f->p.def=func_def_struct_del(f->p.def);
    f->ptype=FUNC_P_NULL;
  }
}

void func_def_p_clone(func_t *f, func_t *g)
{
  if(func_ptype(g)!=FUNC_P_DEF){ FUNC_ERROR_ARG2("func_def_p_clone",f,g); }
  if(f->p.def!=NULL){ func_def_p_del(f); }
  f->ptype=FUNC_P_DEF;
  f->p.def=func_def_struct_new();
  (*(f->p.def))=(*(g->p.def));
  f->p.def->name=char_new(g->p.def->name," \t\n");
}

int func_def_p_cmp(func_t *f, func_t *g)
{
  if(f==NULL || g==NULL || func_ptype(f)!=FUNC_P_DEF || func_ptype(g)!=FUNC_P_DEF || f->p.def==NULL || g->p.def==NULL){ FUNC_ERROR_ARG2("func_def_cmp",f,g); }
  return char_cmp(f->p.def->name,g->p.def->name);
}

///////////////////////////////////////////////////

func_def_t *func_def_p(func_t *f)
{
  if(f==NULL || func_ptype(f)!=FUNC_P_DEF){ FUNC_ERROR_ARG1("func_def_p",f); }
  return f->p.def;
}

char *func_def_name(func_t *f)
{
  if(f==NULL || func_ptype(f)!=FUNC_P_DEF || f->p.def==NULL){ FUNC_ERROR_ARG1("func_def_name",f); }
  return f->p.def->name;
}

int func_def_amin(func_t *f)
{
  if(f==NULL || func_ptype(f)!=FUNC_P_DEF || f->p.def==NULL){ FUNC_ERROR_ARG1("func_def_amin",f); }
  return f->p.def->amin;
}

int func_def_amax(func_t *f)
{
  if(f==NULL || func_ptype(f)!=FUNC_P_DEF || f->p.def==NULL){ FUNC_ERROR_ARG1("func_def_amax",f); }
  return f->p.def->amax;
}

///////////////////////////////////////////////////

func_def_t *func_def_struct_new(void)
{
  func_def_t *def=NULL;
  def=(func_def_t*)malloc(sizeof(func_def_t)*1);
  def->name=NULL;
  def->amin=0;
  def->amax=FUNC_NOLIMIT;
  return def;
}

func_def_t *func_def_struct_del(func_def_t *def)
{
  if(def!=NULL){
    def->name=char_del(def->name);
    free(def);
  }
  return NULL;
}

//EOF
