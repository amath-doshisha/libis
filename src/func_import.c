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

#define SPC " \t\n"
#define BLK0 "({["
#define BLK1 ")}]"

func_t *func_set_script(func_t *f)
{
  strings *a=NULL,*b=NULL,*vlist=NULL;
  func_t *g=NULL;
  if(func_is_strings(f) && func_strings_size(f)==1){
    if(str_match(func_strings_at(f,0),"=",BLK0,BLK1)){
      a=strings_split(func_strings_at(f,0),"=",BLK0,BLK1,SPC);
      if(strings_size(a)==2 && str_match(strings_at(a,0),"^%A%I*$",BLK0,BLK1)){
	g=func_set(func_def(strings_at(a,0),func_script(strings_at(a,1)),0,FUNC_NOLIMIT));
      }else if(strings_size(a)==2 && str_match(strings_at(a,0),"^%A%I*(%m*)$",BLK0,BLK1)){
	b=strings_split_mask(strings_at(a,0),BLK0,BLK1,SPC);
	if(strings_size(b)==1){
	  g=func_set(func_def(strings_at(b,0),func_script(strings_at(a,1)),0,0));
	}else if(strings_size(b)>1){
	  vlist=strings_split(strings_at(b,1),",",BLK0,BLK1,SPC);
	  func_scope_begin(func_strings_strings(vlist));
	  g=func_set(func_def(strings_at(b,0),func_eval(func_script(strings_at(a,1))),strings_size(vlist),strings_size(vlist)));
	  func_scope_end();
	}
      }
    }
  }
  if(g==NULL){ g=f; }else{ f=func_del(f); }
  a=strings_del(a);
  b=strings_del(b);
  vlist=strings_del(vlist);
  return g;
}

func_t *func_op_set_new(void)
{
  func_t *f=func_builtin_new("set");
  func_builtin_p(f)->amin=1;
  func_builtin_p(f)->amax=1;
  func_builtin_p(f)->type=FUNC_COMMAND;
  func_builtin_p(f)->eval=func_set_eval;
  return f;
}

func_t *func_set(func_t *f)
{
  func_init();
  return func_arg1_new("set",f);
}

func_t *func_set_eval(func_t *f)
{
  func_init();
  if(f==NULL || !func_is(f,"set") || func_asize(f)!=1){ FUNC_ERROR_ARG1("func_set_eval",f); }
  func_scope_set(0,FR(func_aget(f,0)));
  f=func_del(f);
  return f;
}

////////////////////////////////////////////

func_t *func_op_end_new(void)
{
  func_t *f=func_builtin_new("end");
  func_builtin_p(f)->type=FUNC_COMMAND;
  func_builtin_p(f)->eval=func_end_eval;
  return f;
}

func_t *func_end(func_t *f)
{
  func_init();
  return func_arg1_new("end",f);
}

func_t *func_end_eval(func_t *f)
{
  func_init();
  if(f==NULL || !func_is(f,"end") || func_asize(f)!=0){ FUNC_ERROR_ARG1("func_end_eval",f); }
  func_scope_end();
  f=func_del(f);
  return NULL;
}

///////////////////////////////////////

func_t *func_op_begin_new(void)
{
  func_t *f=func_builtin_new("begin");
  func_builtin_p(f)->amin=0;
  func_builtin_p(f)->amax=1;
  func_builtin_p(f)->type=FUNC_COMMAND;
  func_builtin_p(f)->eval=func_begin_eval;
  return f;
}

func_t *func_begin(func_t *f)
{
  func_init();
  return func_arg1_new("begin",f);
}

func_t *func_begin_eval(func_t *f)
{
  func_init();
  if      (f==NULL || !func_is(f,"begin"))                    { FUNC_ERROR_ARG1("func_begin_eval",f); }
  else if(func_asize(f)==0)                                   { func_scope_begin(NULL);        f=func_del(f); }
  else if(func_asize(f)==1 && func_is_strings(func_aget(f,0))){ func_scope_begin(FR(f->a[0])); f=func_del(f); }
  return f;
}

//////////////////////////////////////////////////////

void func_import_basic(int n)
{
  func_init();
  // default op for debug only
  func_scope_set(n,func_op_builtin_new());
  func_scope_set(n,func_op_def_new());
  func_scope_set(n,func_op_list_new());
  func_scope_set(n,func_op_table_new());
  func_scope_set(n,func_op_strings_new());
  func_scope_set(n,func_op_scope_new());
  // basic
  func_scope_set(n,func_op_nan_new());
  func_scope_set(n,func_op_inf_new());
  func_scope_set(n,func_op_zero_new());
  func_scope_set(n,func_op_one_new());
  func_scope_set(n,func_op_bigint_new());
  func_scope_set(n,func_op_real_new());
  func_scope_set(n,func_op_complex_new());
  func_scope_set(n,func_op_var_new());
  func_scope_set(n,func_op_add_new());
  func_scope_set(n,func_op_mul_new());
  func_scope_set(n,func_op_sqrt_new());
  func_scope_set(n,func_op_exp_new());
  func_scope_set(n,func_op_log_new());
  func_scope_set(n,func_op_pow_new());
  func_scope_set(n,func_op_sin_new());
  func_scope_set(n,func_op_cos_new());
  func_scope_set(n,func_op_tan_new());
  func_scope_set(n,func_op_asin_new());
  func_scope_set(n,func_op_acos_new());
  func_scope_set(n,func_op_atan_new());
  func_scope_set(n,func_op_sinh_new());
  func_scope_set(n,func_op_cosh_new());
  func_scope_set(n,func_op_tanh_new());
  func_scope_set(n,func_op_asinh_new());
  func_scope_set(n,func_op_acosh_new());
  func_scope_set(n,func_op_atanh_new());
  func_scope_set(n,func_op_ivec_new());
  func_scope_set(n,func_op_rvec_new());
  func_scope_set(n,func_op_cvec_new());
  func_scope_set(n,func_op_rmat_new());
  func_scope_set(n,func_op_cmat_new());
  func_scope_set(n,func_op_begin_new());
  func_scope_set(n,func_op_end_new());
  func_scope_set(n,func_op_set_new());
  func_scope_set(n,func_op_print_new());
  func_scope_set(n,func_op_eval_new());
  func_scope_set(n,func_op_evalf_new());
  func_scope_set(n,func_op_expand_new());
  func_scope_set(n,func_op_diff_new());
  func_scope_set(n,func_op_grad_new());
  func_scope_set(n,func_op_gbasis_new());
  func_scope_set(n,func_def("I", func_complex_i(),0,0));
  func_scope_set(n,func_def("i", func_sqrt(func_bigint_int(-1,1)),0,0));
  func_scope_set(n,func_def("PI",func_mul(func_atan(func_one()),func_bigint_int(4,1)),0,0));
  func_scope_set(n,func_def("Pi",func_mul(func_atan(func_one()),func_bigint_int(4,1)),0,0));
  func_scope_set(n,func_def("pi",func_mul(func_atan(func_one()),func_bigint_int(4,1)),0,0));
}


//EOF
