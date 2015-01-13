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

//////////////////////////////////////////////////////////

#define SPC " \t\n"

func_t *func_var_script(func_t *f)
{
  strings *a=NULL;
  func_t *g=NULL;
  if(func_is_strings(f) && func_strings_size(f)==1){
    if(str_match(func_strings_at(f,0),"^#%d+$",NULL,NULL)){ // variable
      a=strings_split(func_strings_at(f,0),"#",NULL,NULL,SPC);
      if(strings_size(a)==1){
	g=func_var1(atoi(strings_at(a,0)),1);
      }
    }else if(str_match(func_strings_at(f,0),"^#%d+%^%d+$",NULL,NULL)){ // variable with power
      a=strings_split(func_strings_at(f,0),"#^",NULL,NULL,SPC);
      if(strings_size(a)==2){
	g=func_var1(atoi(strings_at(a,0)),atoi(strings_at(a,1)));
      }
    }else if(str_match(func_strings_at(f,0),"^#%d+%^(%p?%d+)$",NULL,NULL)){ // variable with power
      a=strings_split(func_strings_at(f,0),"#^()",NULL,NULL,SPC);
      if(strings_size(a)==2){
	g=func_var1(atoi(strings_at(a,0)),atoi(strings_at(a,1)));
      }
    }
  }
  if(g==NULL){ g=f; }else{ f=func_del(f); }
  a=strings_del(a);
  return g;
}

/////////////////////////////////////////////////////

void func_var_print(func_t *f)
{
  func_t *vars=NULL;
  int i;
  vars=func_scope_vlist(0);
  for(i=0; i<func_var_size(f); i++){
    if(vars!=NULL && func_var_num(f,i)>=0 && func_var_num(f,i)<(func_strings_size(vars))){
      printf("%s",func_strings_at(vars,func_var_num(f,i)));
    }else{
      printf("#%d",func_var_num(f,i));
    }
    func_print_hat(func_var_pow(f,i));
    if(i<func_var_size(f)-1) printf("*");
  }
}

/////////////////////////////////////////////////////

void func_var_p_new(func_t *f)
{
  if(f->p.mem!=NULL){ func_var_p_del(f); }
  f->ptype=FUNC_P_VAR;
  f->p.mem=malloc(sizeof(func_var_struct));
  f->p.var->n=0;
  f->p.var->num=NULL;
  f->p.var->pow=NULL;
}

void func_var_p_del(func_t *f)
{
  if(func_ptype(f)!=FUNC_P_VAR){ FUNC_ERROR_ARG1("func_var_p_del",f); }
  if(f->p.var!=NULL){
    f->p.var->num=ivec_free(f->p.var->num);
    f->p.var->pow=ivec_free(f->p.var->pow);
    free(f->p.var);
    f->p.var=NULL;
    f->ptype=FUNC_P_NULL;
  }
}

void func_var_p_clone(func_t *f, func_t *g)
{
  int n;
  if(func_ptype(g)!=FUNC_P_VAR){ FUNC_ERROR_ARG2("func_var_p_clone",f,g); }
  if(f->p.mem!=NULL){ func_var_p_del(f); }
  f->ptype=FUNC_P_VAR;
  f->p.mem=malloc(sizeof(func_var_struct));
  n=g->p.var->n;
  f->p.var->n=n;
  f->p.var->num=ivec_allocate(n);
  f->p.var->pow=ivec_allocate(n);
  ivec_copy(n,f->p.var->num,g->p.var->num);
  ivec_copy(n,f->p.var->pow,g->p.var->pow);
}

/////////////////////////////////////////////////////////////////////////////

void func_var_p_resize(func_t *f, int n)
{
  if(func_ptype(f)!=FUNC_P_VAR){ FUNC_ERROR_ARG1("func_var_p_resize",f); }
  if(f->p.var==NULL){ func_var_p_new(f); }
  if(n>0){
    f->p.var->num=ivec_free(f->p.var->num);
    f->p.var->pow=ivec_free(f->p.var->pow);
    f->p.var->n=n;
    f->p.var->num=ivec_allocate(n);
    f->p.var->pow=ivec_allocate(n);
  }
}

/////////////////////////////////////////////////////

int func_var_size(func_t *f)
{
  if(f==NULL || func_ptype(f)!=FUNC_P_VAR){ FUNC_ERROR_ARG1("func_var_size",f); }
  return f->p.var->n;
}

int func_var_num(func_t *f, int i)
{
  if(f==NULL || func_ptype(f)!=FUNC_P_VAR || f->p.var==NULL || f->p.var->num==NULL || i<0 || i>=func_var_size(f)){ FUNC_ERROR_ARG1("func_var_num",f); }
  return f->p.var->num[i];
}

int func_var_pow(func_t *f, int i)
{
  if(f==NULL || func_ptype(f)!=FUNC_P_VAR || f->p.var==NULL || f->p.var->pow==NULL || i<0 || i>=func_var_size(f)){ FUNC_ERROR_ARG1("func_var_pow",f); }
  return f->p.var->pow[i];
}

/////////////////////////////////////////////////////

static const char *__func_var="@#";

func_t *func_var(int n)
{
  func_t *f=NULL;
  f=func_new(__func_var);
  func_var_p_resize(f,n);
  return f;
}

func_t *func_var1(int var, int pow)
{
  func_t *f=NULL;
  if(pow==0) return func_one();
  if(var>=0){
    f=func_var(1);
    f->p.var->num[0]=var;
    f->p.var->pow[0]=pow;
  }
  return f;
}

/////////////////////////////////////////////////////

int func_is_var(func_t *f)
{
  return func_is(f,__func_var);
}

int func_is_1var(func_t *f)
{
  return (f!=NULL && func_is_var(f) && func_var_size(f)==1 && func_var_pow(f,0)==1);
}

/////////////////////////////////////////////////////

func_t *func_op_var_new(void)
{
  func_t *f=func_builtin_new(__func_var);
  func_builtin_p(f)->order=FUNC_ORDER_VAR;
  func_builtin_p(f)->ptype=FUNC_P_VAR;
  func_builtin_p(f)->p_new=func_var_p_new;
  func_builtin_p(f)->p_del=func_var_p_del;
  func_builtin_p(f)->p_clone=func_var_p_clone;
  func_builtin_p(f)->p_cmp=func_var_cmp;
  func_builtin_p(f)->print=func_var_print;
  return f;
}


//private///////////////////////////////////////////////////

void func_var_sort(func_t *f)
{
  int *I=NULL;
  if(func_is_var(f)){
    I=ivec_allocate(func_var_size(f));
    ivec_sort(func_var_size(f),f->p.var->num,I);
    ivec_relocate(func_var_size(f),f->p.var->pow,I);
    I=ivec_free(I);
  }else{ FUNC_ERROR_ARG1("func_var_sort",f); }
}


void func_var_mv_pow_mul(func_t *f)
{
  int i,j,flag;
  if(func_is_var(f)){
    for(i=0; i<(func_var_size(f))-1; i++){
      flag=1;
      if(f->p.var->pow[i]==0) flag=0;
      for(j=i+1; flag && j<func_var_size(f); j++){
	if(f->p.var->num[i]==f->p.var->num[j]){
	  f->p.var->pow[j]+=f->p.var->pow[i];
	  f->p.var->pow[i]=0;
	  flag=0;
	}
      }
    }
  }else{ FUNC_ERROR_ARG1("func_var_mv_pow_mul",f); }
}

void func_var_mv_pow_lcm(func_t *f)
{
  int i,j,flag;
  if(func_is_var(f)){
    for(i=0; i<(func_var_size(f))-1; i++){
      flag=1;
      for(j=i+1; flag && j<func_var_size(f); j++){
	if(f->p.var->num[i]==f->p.var->num[j]){
	  if(f->p.var->pow[i] >= f->p.var->pow[j]){
	    f->p.var->pow[j]=f->p.var->pow[i];
	  }
	  f->p.var->pow[i]=0;
	  flag=0;
	}
      }
    }
  }else{ FUNC_ERROR_ARG1("func_var_mv_pow_lcm",f); }
}


void func_var_replace(func_t *f, int n, int *num, int *pow)
{
  f->p.var->num=ivec_free(f->p.var->num);
  f->p.var->pow=ivec_free(f->p.var->pow);
  if(n<=0){
    f->p.var->n=0;
    f->p.var->num=NULL;
    f->p.var->pow=NULL;
  }else{
    f->p.var->n=n;
    f->p.var->num=num;
    f->p.var->pow=pow;
  }
}

void func_var_rm_pow0(func_t *f)
{
  int n,*I=NULL,*pow=NULL,*num=NULL;
  if(func_is_var(f)){
    I=ivec_allocate_index_if_not(func_var_size(f),f->p.var->pow,0,&n);
    if(n>=func_var_size(f)){
      // do nothing
    }else if(n<=0){
      func_var_replace(f,0,NULL,NULL);
    }else{
      num=ivec_allocate(n);
      pow=ivec_allocate(n);
      ivec_copy_index(n,num,f->p.var->num,I);
      ivec_copy_index(n,pow,f->p.var->pow,I);
      func_var_replace(f,n,num,pow);
    }
    I=ivec_free(I);
  }else{ FUNC_ERROR_ARG1("func_var_rm_pow0",f); }
}

func_t *func_var_check_mul(func_t *f)
{
  if(func_is_var(f)){
    func_var_sort(f);
    func_var_mv_pow_mul(f);
    func_var_rm_pow0(f);
    if(func_var_size(f)<=0){ f=func_del(f); f=func_one(); }
  }else{ FUNC_ERROR_ARG1("func_var_check_mul",f); }
  return f;
}

func_t *func_var_check_lcm(func_t *f)
{
  if(func_is_var(f)){
    func_var_sort(f);
    func_var_mv_pow_lcm(f);
    func_var_rm_pow0(f);
    if(func_var_size(f)<=0){ f=func_del(f); f=func_one(); }
  }else{ FUNC_ERROR_ARG1("func_var_check_lcm",f); }
  return f;
}

//////// private ////////////////

func_t *func_var_mul(func_t *f1, func_t *f2)
{
  int n,m;
  func_t *f=NULL;
  if(func_is_var(f1) && func_is_var(f2)){
    n=func_var_size(f1);
    m=func_var_size(f2);
    f=func_var(n+m);
    ivec_copy(n,f->p.var->num,f1->p.var->num);
    ivec_copy(n,f->p.var->pow,f1->p.var->pow);
    ivec_copy(m,f->p.var->num+n,f2->p.var->num);
    ivec_copy(m,f->p.var->pow+n,f2->p.var->pow);
    f=func_var_check_mul(f);
  }else{ FUNC_ERROR_ARG2("func_var_mul",f1,f2); }
  f1=func_del(f1);
  f2=func_del(f2);
  return f;
}

func_t *func_var_lcm(func_t *f1, func_t *f2)
{
  int n,m;
  func_t *f=NULL;
  if(func_is_var(f1) && func_is_var(f2)){
    n=func_var_size(f1);
    m=func_var_size(f2);
    f=func_var(n+m);
    ivec_copy(n,f->p.var->num,f1->p.var->num);
    ivec_copy(n,f->p.var->pow,f1->p.var->pow);
    ivec_copy(m,f->p.var->num+n,f2->p.var->num);
    ivec_copy(m,f->p.var->pow+n,f2->p.var->pow);
    f=func_var_check_lcm(f);
  }else{ FUNC_ERROR_ARG2("func_var_lcm",f1,f2); }
  f1=func_del(f1);
  f2=func_del(f2);
  return f;
}

// f=f1/f2
func_t *func_var_div(func_t *f1, func_t *f2)
{
  int n,m;
  func_t *f=NULL;
  if(func_is_var(f1) && func_is_var(f2)){
    n=func_var_size(f1);
    m=func_var_size(f2);
    f=func_var(n+m);
    ivec_copy(n,f->p.var->num,f1->p.var->num);
    ivec_copy(n,f->p.var->pow,f1->p.var->pow);
    ivec_copy(m,f->p.var->num+n,f2->p.var->num);
    ivec_copy(m,f->p.var->pow+n,f2->p.var->pow);
    ivec_scale(m,f->p.var->pow+n,-1);
    f=func_var_check_mul(f);
  }else{ FUNC_ERROR_ARG2("func_var_div",f1,f2); }
  f1=func_del(f1);
  f2=func_del(f2);
  return f;
}


//////// private /////////////////////////////////////

func_t *func_var_pow_n(func_t *f, int p)
{
  int i;
  func_t *g=NULL;
  if(func_is_var(f)){
    if(p==0){ g=func_one(); }
    else if(p==1){ g=FR(f); }
    else{
      g=func_clone(FR(f));
      for(i=0; i<func_var_size(g); i++){ g->p.var->pow[i]*=p; }
    }
  }else{ FUNC_ERROR_ARG1("func_var_pow_n",f); }
  f=func_del(f);
  return g;
}

//private//////////////////////////////////////////////////


int func_var_degree(func_t *f)                              //次数の出力
{
  if(!func_is_var(f)) { FUNC_ERROR_ARG1("func_var_degree",f); }
  return ivec_sum(func_var_size(f),f->p.var->pow);
}


int func_var_cmp_grlex(func_t *f, func_t *g)
{
  int n,m;
  if(!func_is_var(f) || !func_is_var(g)) { FUNC_ERROR_ARG2("func_var_cmp_grlex",f,g); }
  n=func_var_degree(f);
  m=func_var_degree(g);
  if(n<m)  return -1;
  if(n>m)  return +1;
  return func_var_cmp_lex(f,g);
}

int func_var_cmp_lex(func_t *f, func_t *g)
{
  int i,j,a,b;
  i=0; j=0;
  if(!func_is_var(f) || !func_is_var(g)) { FUNC_ERROR_ARG2("func_var_cmp_lex",f,g); }
  while(i<func_var_size(f) || j<func_var_size(g)){
    if(i<func_var_size(f) && j<func_var_size(g) && f->p.var->num[i] == g->p.var->num[j]){
      a=f->p.var->pow[i++]; b=g->p.var->pow[j++];
    }else if(i<func_var_size(f) && j<func_var_size(g) && f->p.var->num[i] < g->p.var->num[j]){
      a=f->p.var->pow[i++]; b=0;
    }else if(i<func_var_size(f) && j<func_var_size(g) && f->p.var->num[i] > g->p.var->num[j]){
      a=0; b=g->p.var->pow[j++];
    }else if(i>=func_var_size(f) && j<func_var_size(g)){
      a=0; b=g->p.var->pow[j++];
    }else if(j>=func_var_size(g) && i<func_var_size(f)){
      a=f->p.var->pow[i++]; b=0;
    }else{
      printf("Error!\n");
      a=0; b=0;
    }
    if(a<b) return -1;
    if(a>b) return +1;
  }
  return 0;
}

int func_var_cmp_grevlex(func_t *f, func_t *g)
{
  int i,j,n,m,a,b;
  if(!func_is_var(f) || !func_is_var(g)) { FUNC_ERROR_ARG2("func_var_cmp_grevlex",f,g); }
  n=func_var_degree(f);
  m=func_var_degree(g);
  if(n<m)  return -1;
  if(n>m)  return +1;
  i=func_var_size(f)-1; j=func_var_size(g)-1;
  while(i>=0 || j>=0){
    if(i>=0 && j>=0 && f->p.var->num[i] == g->p.var->num[j]){
      a=f->p.var->pow[i--]; b=g->p.var->pow[j--];
    }else if(i>=0 && j>=0 && f->p.var->num[i] < g->p.var->num[j]){
      a=0; b=g->p.var->pow[j];
    }else if(i>=0 && j>=0 && f->p.var->num[i] > g->p.var->num[j]){
      a=f->p.var->pow[i]; b=0;
    }else if(i<0 && j>0){
      a=0; b=g->p.var->pow[j--];
    }else if(j<0 && i>0){
      a=f->p.var->pow[i--]; b=0;
    }else{
      printf("Error!\n");
      a=0; b=0;
    }
    if(a<b) return +1;
    if(a>b) return -1;
  }
  return 0;
}

int func_var_cmp(func_t *f, func_t *g)
{
  if(!func_is_var(f) || !func_is_var(g)) { FUNC_ERROR_ARG2("func_var_cmp",f,g); }
  if(func_get_mono_order()==FUNC_MONO_LEX)          { return func_var_cmp_lex(f,g);  }
  else if(func_get_mono_order()==FUNC_MONO_GRLEX)   { return func_var_cmp_grlex(f,g);  }
  else if(func_get_mono_order()==FUNC_MONO_GREVLEX) { return func_var_cmp_grevlex(f,g);  }
  else                                              { return func_var_cmp_grlex(f,g);  }
}

//////////////////////////////////////////////////////////////////////////////////

// 変数の個数を数える
int func_var_n_var(func_t *f)
{
  int i,count=1;
  if(func_is_number(f)) return 0;
  if(!func_is_var(f)) return -1;
  else{
    for(i=1; i<func_var_size(f); i++){
      if(f->p.var->num[0]!=f->p.var->num[i]) count=count+1;
    } 
  }
  return count;
}

///////////////////////////////////////////

//第n関数が存在するか
int func_var_exist_varn(func_t *f, int n)
{
  int i;
  if(func_is_number(f)) return 0;
  if(!func_is_var(f)) return 0;
  for(i=0; i<func_var_size(f); i++){
    if(f->p.var->num[i]==n){ return 1; }
  }
  return 0;
}

int func_var_get_index(func_t *f, int n)
{
  int i;
  if(!func_is_var(f)){ FUNC_ERROR_ARG1("func_var_get_index",f); }
  for(i=0; i<func_var_size(f); i++){
    if(f->p.var->num[i]==n){ return i; }
  }
  return -1;
}

///////////////////////////////////////////


//１変数多項式の判定 文字番号出力
int func_var1n(func_t *f)
{
  if(func_is_number(f))return -1;
  if(func_is_var(f))  return func_var_var1n(f);
  if(func_is_mono_not_num(f)) return func_mono_var1n(f);
  if(func_is_poly(f))         return func_poly_var1n(f);
  return -1;
}

int func_var_var1n(func_t *f)
{
  int i;
  if(!func_is_var(f)) return -1;
  for(i=1; i<func_var_size(f); i++){
    if(f->p.var->num[0]!=f->p.var->num[i]) return -1;
  }
  return f->p.var->num[0];
}

//varが割れるかどうか
int func_var_can_div(func_t *f, func_t *g)
{
  int i,j,a,b;
  if(!func_is_var(f) || !func_is_var(g)) { return 0; }
  i=0; j=0;
  while(i<func_var_size(f) || j<func_var_size(g)){
    if(i<func_var_size(f) && j<func_var_size(g) && f->p.var->num[i]==g->p.var->num[j]){
      a=f->p.var->pow[i++]; b=g->p.var->pow[j++];
      if(a<b)  return 0;                                             //割り算不可
    }else if(i<func_var_size(f) && j<func_var_size(g) && f->p.var->num[i] < g->p.var->num[j]){
      a=f->p.var->pow[i++]; b=0;
      if(a<b)  return 0;
    }else if(i<func_var_size(f) && j<func_var_size(g) && f->p.var->num[i] > g->p.var->num[j]){
      a=0; b=g->p.var->pow[j++];
      if(a<b)  return 0;
    }else if(i>=func_var_size(f) && j<func_var_size(g)){
      a=0; b=g->p.var->pow[j++];
      if(a<b)  return 0;
    }else if(j>=func_var_size(g) && i<func_var_size(f)){
      a=f->p.var->pow[i++]; b=0;
      if(a<b)  return 0;
    }else{
      printf("Error!\n");
      a=0; b=0;
      return 0;
    }
  }
  return 1;                                                          //割り算可能
}

//EOF
