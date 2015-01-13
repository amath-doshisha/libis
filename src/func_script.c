#include<stdio.h>
#include<stdlib.h>
#include<stdarg.h>
#include"is_macros.h"
#include"is_ivec.h"
#include"is_strings.h"
#include"is_print.h"
#include"is_func.h"
#include"is_rmulti.h"
#include"is_cmulti.h"

#define FR(F) func_retain(F)


/////////////////////////////////////////////

func_t *func_scriptf(char* fmt, ...)
{
  func_t *f=NULL;
  char *str=NULL;
  va_list argp;
  str=(char*)malloc(strlen(fmt)*2+1024);
  va_start(argp,fmt);
  vsprintf(str,fmt,argp);  
  f=func_script(str);
  free(str);
  return f;
}

func_t *func_scriptf_echo(char* fmt, ...)
{
  func_t *f=NULL;
  char *str=NULL;
  va_list argp;
  str=(char*)malloc(strlen(fmt)*2+1024);
  va_start(argp,fmt);
  vsprintf(str,fmt,argp);
  printf("%s",str);
  f=func_script(str);
  free(str);
  return f;
}

//////////////////////////////////////////////////////////

func_t *func_script(const char *str)
{
  func_t *f=NULL;
  func_init();
  f=func_strings_char(str);
  f=func_null_script(f);    // NULL
  f=func_nan_script(f);     // nan
  f=func_inf_script(f);     // inf
  f=func_zero_script(f);    // 0
  f=func_one_script(f);     // 1
  f=func_bigint_script(f);  // Z, Q
  f=func_real_script(f);    // R
  f=func_var_script(f);     // var: #0, #1, #2, ...
  f=func_set_script(f);     // =, def
  f=func_add_script(f);     // +
  f=func_sub_script(f);     // -
  f=func_mul_script(f);     // *
  f=func_div_script(f);     // /
  f=func_pow_script(f);     // ^
  f=func_bracket_script(f); // (......)
  f=func_list_script(f);    // list: {......}, [......]
  f=func_builtin_script(f); // built-in: f, f(...)
  f=func_def_script(f);     // def: f, f(...)
  if(func_is_strings(f)){
    printf("Error in func_script(const char *str)\n");
    printf("str='%s'\n",str);
    printf("f="); func_print(f); printf("\n");
    exit(0);
  }
  return f;
}

//////////////////////////////////////////////////////////

#define SPC " \t\n"
#define BLK0 "({["
#define BLK1 ")}]"

func_t *func_bracket_script(func_t *f)
{
  strings *a=NULL;
  func_t *g=NULL;
  if(func_is_strings(f) && func_strings_size(f)==1){
    if(str_match(func_strings_at(f,0),"^(%m+)$",BLK0,BLK1)){
      a=strings_split_mask(func_strings_at(f,0),BLK0,BLK1,SPC);
      if(strings_size(a)==1){	
	g=func_script(strings_at(a,0));
      }
    }
  }
  if(g==NULL){ g=f; }else{ f=func_del(f); }
  a=strings_del(a);
  return g;
}

func_t *func_null_script(func_t *f)
{
  func_t *g=NULL;
  if(func_is_strings(f) && func_strings_size(f)==1 && str_match(func_strings_at(f,0),"^NULL$",BLK0,BLK1)){
    g=NULL;
    f=func_del(f);
  }else{ g=f; }
  return g;
}




//EOF
