#include<stdio.h>
#include<stdlib.h>
#include"is_macros.h"
#include"is_ivec.h"
#include"is_strings.h"
#include"is_print.h"
#include"is_func.h"
#include"is_rmulti.h"
#include"is_cmulti.h"

#define FR(F) func_retain(F)

int func_is_mono_not_num(func_t *f)
{
  if(f==NULL) return 0;
  if(func_is_var(f)) return 1;
  if(func_is_split_mul(f) && func_is_var(func_get_split_part(f))){ return 1; }
  return 0;
}

int func_is_mono(func_t *f)
{
  if(f==NULL) return 0;
  if(func_is_number(f)) return 1;
  return func_is_mono_not_num(f);
}


////////////////////////////////////////////////

int func_mono_n_var(func_t *f)
{
  if(!func_is_mono(f)) return -1;
  if(func_is_number(f)) return 0;
  return func_var_n_var(func_mono_get_var(f));
}

//////////////////////////////////////////////////////////////////////////////////

int func_mono_exist_varn(func_t *f, int n)
{
  if(!func_is_mono_not_num(f)) return 0;
  return func_var_exist_varn(func_mono_get_var(f),n);
}

//////////////////////////////////////////////////////////////////////////////////


int func_mono_var1n(func_t *f)
{
  if(!func_is_mono_not_num(f)) return -1;
  return func_var_var1n(func_mono_get_var(f));
}

///////////////////////////////////////////////////////////////////

int func_mono_can_div(func_t *f, func_t *g)
{
  return func_var_can_div(func_mono_get_var(f),func_mono_get_var(g));
}

int func_mono_degree(func_t *f)
{
  if(func_is_mono_not_num(f)){ return func_var_degree(func_mono_get_var(f)); }
  return 0; 
}

////////private///////////////////////////

func_t *func_mono_get_var(func_t *f)
{
  if(!func_is_mono_not_num(f)){ return NULL; }
  else if(func_is_var(f))     { return f; }
  else                        { return func_get_split_part(f); }
}

func_t *func_mono_get_coeff(func_t *f)
{
  if(!func_is_mono_not_num(f)){ return NULL; }
  else if(func_is_var(f))     { return NULL; }
  else                        { return func_get_number_part(f); }
}

//EOF
