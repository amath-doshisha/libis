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
#define FGET_SP(f) (func_get_split_part(f))
#define FGET_NB(f) (func_get_number_part(f))

//////////////////////////////////////////////////////////

#define SPC " \t\n"
#define BLK0 "({["
#define BLK1 ")}]"

func_t *func_mul_script(func_t *f)
{
  int i;
  strings *a=NULL;
  func_t *g=NULL;
  if(func_is_strings(f) && func_strings_size(f)==1){
    if(str_match(func_strings_at(f,0),"%*",BLK0,BLK1)){
      a=strings_split(func_strings_at(f,0),"*",BLK0,BLK1,SPC);
      if(strings_size(a)>0){
	g=func_mul_new(strings_size(a));
	for(i=0; i<func_asize(g); i++){
	  func_aset(g,i,func_script(strings_at(a,i)));
	}
	g=func_mul_eval(g);
      }
    }
  }
  if(g==NULL){ g=f; }else{ f=func_del(f); }
  a=strings_del(a);
  return g;
}



/////////////////////////////////////

static char *__func_mul="mul";

///////////////////////////////////

func_t *func_op_mul_new(void)
{
  func_t *f=func_builtin_new(__func_mul);
  func_builtin_p(f)->order=FUNC_ORDER_MUL;
  func_builtin_p(f)->amax=FUNC_NOLIMIT;
  func_builtin_p(f)->type=FUNC_SCALAR_DEPEND;
  func_builtin_p(f)->ptype=FUNC_P_POWER;
  func_builtin_p(f)->p_new=func_pow_p_new;
  func_builtin_p(f)->p_del=func_pow_p_del;
  func_builtin_p(f)->p_clone=func_pow_p_clone;
  func_builtin_p(f)->p_cmp=func_power_cmp;
  func_builtin_p(f)->eval=func_mul_eval;
  func_builtin_p(f)->print=func_mul_print;
  return f;
}

///////////////////////////////////

void func_mul_print(func_t *f)
{
  int i;
  if(func_power(f)!=1) printf("(");
  for(i=0; i<func_asize(f); i++){
    if(func_is_one(func_aget(f,i))){} // do nothing
    else if(i==0 && func_is_bigint(func_aget(f,i)) && bigint_is_neg_one(func_bigint_p(func_aget(f,i)))){ printf("-"); }
    else{
      if(func_is_add(func_aget(f,i)) && func_power(func_aget(f,i))==1) printf("(");
      if(func_is_mul(func_aget(f,i))) printf("(");
      func_print(func_aget(f,i));
      if(func_is_add(func_aget(f,i)) && func_power(func_aget(f,i))==1) printf(")");
      if(func_is_mul(func_aget(f,i))) printf(")");
      if(i<(func_asize(f))-1) printf("*");
    }
  }
  if(func_power(f)!=1) printf(")");
  func_print_hat(func_power(f));
}

///////////////////////////////////


int func_is_mul(func_t *f)
{
  return func_is(f,__func_mul);
}

////////////////////////////////////////////////////

func_t *func_mul_split_pow_list(func_t *f)
{
  func_t *g=NULL;
  if(f==NULL){ FUNC_ERROR_ARG1("func_mul_split_pow_list",f); }
  g=func_list(2);
  if     (func_is(f,"pow")) { func_aset(g,0,FR(func_aget(f,0))); func_aset(g,1,FR(func_aget(f,1))); }
  else if(func_has_power(f)){ func_aset(g,0,func_clone(FR(f))); func_set_power(func_aget(g,0),1); func_aset(g,1,func_bigint_int(func_power(f),1)); }
  else                      { func_aset(g,0,FR(f)); func_aset(g,1,func_one()); }
  f=func_del(f);
  return g;
}

/////private////////////////////////////////////////

void func_mul_args_pow(func_t *f)
{
  int i,j;
  func_t *g=NULL;
  if(!func_is_mul(f)){ FUNC_ERROR_ARG1("func_mul_args_pow",f); }  
  // convert args to split list
  for(i=0; i<func_asize(f); i++){
    g=func_mul_split_pow_list(FR(func_aget(f,i)));
    func_aset(f,i,FR(g));
    g=func_del(g);
  }
  // add a arg to same arg
  for(i=0; i<func_asize(f); i++){
    for(j=i+1; func_aget(f,i)!=NULL && j<func_asize(f); j++){
      if(func_aget(f,j)!=NULL && func_cmp(func_aget(func_aget(f,i),0),func_aget(func_aget(f,j),0))==0){
	g=func_list(2);
	func_aset(g,0,FR(func_aget(func_aget(f,i),0)));
	func_aset(g,1,func_add(FR(func_aget(func_aget(f,i),1)),FR(func_aget(func_aget(f,j),1))));
	func_aset(f,i,FR(g));
	func_adel(f,j);
	g=func_del(g);
      }
    }
  }
  // convert list to pow
  for(i=0; i<func_asize(f); i++){
    if(func_aget(f,i)!=NULL){
      g=func_pow(FR(func_aget(func_aget(f,i),0)),FR(func_aget(func_aget(f,i),1)));
      func_aset(f,i,FR(g));
      g=func_del(g);
    }
  }
  func_a_rm_null(f);
}


////////////////////////////////////////////////////

func_t *func_mul_split_list(func_t *f, func_is_t *fis)
{
  int i;
  func_t *g=NULL;
  if(f==NULL){ FUNC_ERROR_ARG1("func_list_split_coeff",f); }
  g=func_list(2);
  if(func_is_mul(f) && func_power(f)==1){
    f=func_flatten(f,__func_mul);
    func_aset(g,0,func_clone(FR(f)));
    func_aset(g,1,func_clone(FR(f)));
    for(i=0; i<func_asize(f); i++){
      if(fis(func_aget(f,i))){ func_adel(func_aget(g,1),i); }
      else                   { func_adel(func_aget(g,0),i); }
    }
    func_aset(g,0,func_mul_eval_size(FR(func_aget(g,0))));
    func_aset(g,1,func_mul_eval_size(FR(func_aget(g,1))));
  }else if(fis(f)){
    func_aset(g,0,FR(f));
    func_aset(g,1,func_one());
  }else{
    func_aset(g,0,func_one());
    func_aset(g,1,FR(f));
  }
  f=func_del(f);
  return g;
}


///////////////////////////////////////////////////////////////////

void func_mul_args(func_t *f, func_is_t *fis, func_is_t *fin)
{
  int i,j;
  if(!func_is_mul(f)){ FUNC_ERROR_ARG1("func_mul_args",f); }
  for(i=0; i<func_asize(f); i++){
    for(j=0; fis(func_aget(f,i)) && j<func_asize(f); j++){
      if(i!=j && fin(func_aget(f,j))){
        func_aset(f,i,func_mul(FR(func_aget(f,i)),FR(func_aget(f,j))));
	func_adel(f,j);
      }
    }
  }
  func_a_rm_null(f);
}

///private/////////////////////////////////////////////////////////

func_t *func_mul_eval_size(func_t *f)
{
  func_t *g;
  if(!func_is_mul(f)){ FUNC_ERROR_ARG1("func_mul_eval_size",f); }
  func_a_rm_null(f);
  func_args_sort(f);
  if     (func_asize(f)<=0){ f=func_del(f); return func_one(); }
  else if(func_asize(f)==1){ g=func_pow(FR(func_aget(f,0)),func_bigint_int(func_power(f),1)); f=func_del(f); return g; }
  else                     { return f; }
}

////////////////////////////////////////////////////////////

func_t *func_mul_eval(func_t *f)
{
  if(!func_is_mul(f)) return f;
  f=func_flatten(f,__func_mul);
  func_mul_args(f,func_is_complex,func_is_number);
  func_mul_args(f,func_is_real,func_is_number);
  func_mul_args(f,func_in_bigint,func_in_bigint);
  func_mul_args(f,func_is_var,func_is_var);
  func_mul_args_pow(f);
  func_a_rm_op(f,func_is_one);
  f=func_mul_split(f);
  return func_mul_eval_size(f);
}

////////////////////////////////////////////////////////////

int func_is_split_mul(func_t *f)
{
  if(func_is_mul(f) && func_power(f)==1 &&
     func_asize(f)==2 &&
     func_is_number(func_aget(f,0)) &&
     func_is_split(func_aget(f,1))){ return 1; }
  else { return 0; }
}

int func_is_split(func_t *f)
{
  int i;
  if     (f==NULL)          { return 0; }
  else if(func_is_number(f)){ return 0; }
  else if(!(func_is_mul(f) && func_power(f)==1)){ return 1; }
  // mul^1
  for(i=0; i<func_asize(f); i++){
    if(func_is_number(func_aget(f,i))){ return 0; }
  }
  return 1;
}

func_t *func_get_number_part(func_t *f)
{
  if     (f==NULL)             { return NULL; }
  else if(func_is_split_mul(f)){ return func_aget(f,0); }
  else if(func_is_number(f))   { return f; }
  else if(func_is_split(f))    { return NULL; }
  else                         { return NULL; }
}

func_t *func_get_split_part(func_t *f)
{
  if     (f==NULL)             { return NULL; }
  else if(func_is_split_mul(f)){ return func_aget(f,1); }
  else if(func_is_number(f))   { return NULL; }
  else if(func_is_split(f))    { return f; }
  else                         { return NULL; }
}

////////////////////////////////////////////////////////////

func_t *func_mul_split(func_t *f)
{
  func_t *g=NULL,*h=NULL;
  if(!(func_is_mul(f) && func_power(f)==1)){ return f; }
  if(func_is_number(f))                    { return f; }
  if(func_is_split(f))                     { return f; }
  if(func_is_split_mul(f))                 { return f; }
  h=func_mul_split_list(FR(f),func_is_number);
  g=func_mul(FR(func_aget(h,0)),FR(func_aget(h,1)));
  f=func_del(f);
  h=func_del(h);
  return g;
}

////////////////////////////////////////////////////////////

func_t *func_mul(func_t *f1, func_t *f2)
{
  if     (f1==NULL || f2==NULL)                      { FUNC_ERROR_ARG2("func_mul",f1,f2); }
  else if(func_is_zero(f1))                          { f1=func_del(f1); f2=func_del(f2); return func_zero(); }
  else if(func_is_zero(f2))                          { f1=func_del(f1); f2=func_del(f2); return func_zero(); }
  else if(func_is_one(f1))                           { f1=func_del(f1); return f2; }
  else if(func_is_one(f2))                           { f2=func_del(f2); return f1; }
  else if(func_is_var(f1)     && func_is_var(f2))    { return func_var_mul(f1,f2); }
  else if(func_in_bigint(f1)  && func_in_bigint(f2)) { return func_bigint_mul(f1,f2); }
  else if(func_in_real(f1)    && func_in_real(f2))   { return func_real_mul(f1,f2); }
  else if(func_in_complex(f1) && func_in_complex(f2)){ return func_complex_mul(f1,f2); }
  else if(func_is_real(f1)    && func_is_number(f2)) { return func_mul(f1,func_evalf(f2)); }
  else if(func_is_complex(f1) && func_is_number(f2)) { return func_mul(f1,func_evalf(f2)); }
  else if(func_is_number(f1)  && func_is_real(f2))   { return func_mul(func_evalf(f1),f2); }
  else if(func_is_number(f1)  && func_is_complex(f2)){ return func_mul(func_evalf(f1),f2); }
  else if(func_is_number(f1)  && func_is_split(f2))  { return func_arg2_new(__func_mul,f1,f2); }
  else if(func_is_number(f2)  && func_is_split(f1))  { return func_arg2_new(__func_mul,f2,f1); }
  else                                               { return func_mul_eval(func_arg2_new(__func_mul,f1,f2)); }
}

///////////////////////////////////////////////////////

func_t *func_mul_pow_n(func_t *f, int pow)
{
  func_t *g=NULL,*c=NULL,*x=NULL;
  if(f!=NULL && pow==0){ g=func_one(); }
  else if(func_is_split_mul(f)){
    c=func_get_number_part(f);
    if(c==NULL){ c=func_one(); }else{ c=FR(c); }
    c=func_pow_n(c,pow);
    x=func_pow_n(FR(func_get_split_part(f)),pow);
    g=func_mul(FR(c),FR(x));
  }else if(func_is_mul(f)){
    g=func_clone(FR(f));
    func_set_power(g,func_power(g)*pow);
  }else{ FUNC_ERROR_ARG1("func_mul_pow_n",f); }
  c=func_del(c);
  x=func_del(x);
  f=func_del(f);
  return g;
}

/////////////////////////////////////////////////////

func_t *func_mul_new(int n)
{
  func_t *f=NULL;
  f=func_new(__func_mul);
  func_a_resize(f,n);
  return f;
}


//EOF
