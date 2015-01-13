#include<stdio.h>
#include<stdlib.h>
#include"is_macros.h"
#include"is_func.h"

///////////////////////////////////////////

int func_find_ptype(func_t *f) {
  func_t *g=NULL;
  if     (func_is_builtin(f)){ return FUNC_P_BUILTIN; }
  else if(func_is_list(f))   { return FUNC_P_NULL; }
  else if(func_is_table(f))  { return FUNC_P_NULL; }
  else if(func_is_strings(f)){ return FUNC_P_STRINGS; }
  else if(func_is_scope(f))  { return FUNC_P_NULL; }
  else if(func_is_def(f))    { return FUNC_P_DEF; }
  g=func_scope_find(0,func_op(f));
  if     (func_is_builtin(g)){ return func_builtin_ptype(g); }
  else if(func_is_def(g))    { return FUNC_P_NULL; }
  else                       { FUNC_ERROR_ARG1("func_find_ptype",f); }

}

func_set_t *func_find_p_new(func_t *f){
  func_t *g=NULL;
  if     (func_is_builtin(f)){ return func_builtin_p_new; }
  else if(func_is_list(f))   { return NULL; }
  else if(func_is_table(f))  { return NULL; }
  else if(func_is_strings(f)){ return func_strings_p_new; }
  else if(func_is_scope(f))  { return NULL; }
  else if(func_is_def(f))    { return func_def_p_new; }
  g=func_scope_find(0,func_op(f));
  if     (g==NULL)           { return NULL; }
  else if(func_is_builtin(g)){ return func_builtin_p_new_call(g); }
  else if(func_is_def(g))    { return NULL; }
  else                       { FUNC_ERROR_ARG1("func_find_p_new",f); }
}

func_set_t *func_find_p_del(func_t *f){
  func_t *g=NULL;
  if     (func_is_builtin(f)){ return func_builtin_p_del; }
  else if(func_is_list(f))   { return NULL; }
  else if(func_is_table(f))  { return NULL; }
  else if(func_is_strings(f)){ return func_strings_p_del; }
  else if(func_is_scope(f))  { return NULL; }
  else if(func_is_def(f))    { return func_def_p_del; }
  g=func_scope_find(0,func_op(f));
  if     (g==NULL)           { return NULL; }
  else if(func_is_builtin(g)){ return func_builtin_p_del_call(g); }
  else if(func_is_def(g))    { return NULL; }
  else                       { FUNC_ERROR_ARG2("func_find_p_del",f,g); }
}

func_set2_t *func_find_p_clone(func_t *f)
{
  func_t *g=NULL;
  if     (func_is_builtin(f)){ return func_builtin_p_clone; }
  else if(func_is_list(f))   { return NULL; }
  else if(func_is_table(f))  { return NULL; }
  else if(func_is_strings(f)){ return func_strings_p_clone; }
  else if(func_is_scope(f))  { return NULL; }
  else if(func_is_def(f))    { return func_def_p_clone; }
  g=func_scope_find(0,func_op(f));
  if     (g==NULL)           { return NULL; }
  else if(func_is_builtin(g)){ return func_builtin_p_clone_call(g); }
  else if(func_is_def(g))    { return NULL; }
  else                       { FUNC_ERROR_ARG1("func_find_p_clone",f); }
}

func_cmp_t *func_find_p_cmp(func_t *f)
{
  func_t *g=NULL;
  if     (func_is_builtin(f)){ return func_builtin_p_cmp; }
  else if(func_is_list(f))   { return NULL; }
  else if(func_is_table(f))  { return NULL; }
  else if(func_is_strings(f)){ return func_strings_p_cmp; }
  else if(func_is_scope(f))  { return NULL; }
  else if(func_is_def(f))    { return func_def_p_cmp; }
  g=func_scope_find(0,func_op(f));
  if     (func_is_builtin(g)){ return func_builtin_cmp_call(g); }
  else if(func_is_def(g))    { return NULL; }
  else                       { FUNC_ERROR_ARG1("func_find_p_cmp",f); }
}

func_arg1_t *func_find_eval(func_t *f)
{
  func_t *g=NULL;
  if     (func_is_builtin(f)){ return NULL; }
  else if(func_is_list(f))   { return func_list_eval; }
  else if(func_is_table(f))  { return NULL; }
  else if(func_is_strings(f)){ return NULL; }
  else if(func_is_scope(f))  { return NULL; }
  else if(func_is_def(f))    { return NULL; }
  g=func_scope_find(0,func_op(f));
  if     (g==NULL)           { return NULL; }
  else if(func_is_builtin(g)){ return func_builtin_eval_call(g); }
  else if(func_is_def(g))    { return func_def_any_eval; }
  else                       { FUNC_ERROR_ARG1("func_find_eval",f); }
}

func_set_t *func_find_print(func_t *f)
{
  func_t *g=NULL;
  if     (func_is_builtin(f)){ return func_builtin_print; }
  else if(func_is_list(f))   { return func_list_print; }
  else if(func_is_table(f))  { return func_print_table; }
  else if(func_is_strings(f)){ return func_strings_print; }
  else if(func_is_scope(f))  { return func_print_scope; }
  else if(func_is_def(f))    { return func_def_print; }
  g=func_scope_find(0,func_op(f));
  if     (g==NULL)           { return NULL; }
  else if(func_is_builtin(g)){ return func_builtin_print_call(g); }
  else if(func_is_def(g))    { return NULL; }
  else                       { FUNC_ERROR_ARG1("func_find_print",f); }
}

int func_find_amin(func_t *f)
{
  func_t *g=NULL;
  if     (func_is_builtin(f)){ return 0; }
  else if(func_is_list(f))   { return 0; }
  else if(func_is_table(f))  { return 0; }
  else if(func_is_strings(f)){ return 0; }
  else if(func_is_scope(f))  { return 3; }
  else if(func_is_def(f))    { return 1; }
  g=func_scope_find(0,func_op(f));
  if     (func_is_builtin(g)){ return func_builtin_amin(g); }
  else if(func_is_def(g))    { return func_def_amin(g); }
  else                       { FUNC_ERROR_ARG1("func_find_amin",f); }
}

int func_find_amax(func_t *f)
{
  func_t *g=NULL;
  if     (func_is_builtin(f)){ return 0; }
  else if(func_is_list(f))   { return FUNC_NOLIMIT; }
  else if(func_is_table(f))  { return FUNC_NOLIMIT; }
  else if(func_is_strings(f)){ return 0; }
  else if(func_is_scope(f))  { return 3; }
  else if(func_is_def(f))    { return 1; }
  g=func_scope_find(0,func_op(f));
  if     (func_is_builtin(g)){ return func_builtin_amax(g); }
  else if(func_is_def(g))    { return func_def_amax(g); }
  else                       { FUNC_ERROR_ARG1("func_find_amax",f); }
}

int func_find_order(func_t *f)
{
  func_t *g=NULL;
  if     (func_is_builtin(f)){ return FUNC_ORDER_BUILTIN; }
  else if(func_is_list(f))   { return FUNC_ORDER_LIST; }
  else if(func_is_table(f))  { return FUNC_ORDER_TABLE; }
  else if(func_is_strings(f)){ return FUNC_ORDER_STRINGS; }
  else if(func_is_scope(f))  { return FUNC_ORDER_SCOPE; }
  else if(func_is_def(f))    { return FUNC_ORDER_DEF; }
  g=func_scope_find(0,func_op(f));
  if     (func_is_builtin(g)){ return func_builtin_order(g); }
  else if(func_is_def(g))    { return FUNC_ORDER_DEFAULT; }
  else                       { FUNC_ERROR_ARG1("func_find_order",f); }
}

int func_find_type(func_t *f)
{
  func_t *g=NULL;
  if     (func_is_builtin(f)){ return FUNC_SCALAR_NO; }
  else if(func_is_list(f))   { return FUNC_SCALAR_NO; }
  else if(func_is_table(f))  { return FUNC_SCALAR_NO; }
  else if(func_is_strings(f)){ return FUNC_SCALAR_NO; }
  else if(func_is_scope(f))  { return FUNC_SCALAR_NO; }
  else if(func_is_def(f))    { return FUNC_SCALAR_NO; }
  g=func_scope_find(0,func_op(f));
  if     (func_is_builtin(g)){ return func_builtin_type(g); }
  else if(func_is_def(g))    { return FUNC_SCALAR_NO; }
  else                       { FUNC_ERROR_ARG1("func_find_type",f); }
}

//EOF
