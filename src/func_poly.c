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

//////////////////////////////////////////////////////////////////

int func_degree(func_t *f)
{
  if(f==NULL)                     { return 0; }
  else if(func_is_var(f))         { return func_var_degree(f); }
  else if(func_is_mono_not_num(f)){ return func_mono_degree(f); }
  else if(func_is_poly(f))        { return func_poly_degree_max(f); }
  else if(func_is_poly_list(f))   { return func_list_degree_max(f); }
  else                            { return 0; }
}

//////////////////////////////////////////////////////////////////

int func_is_poly(func_t *f)                                     //多項式であるか
{
  int i;
  if(f==NULL) return 0;
  if(func_is_mono_not_num(f)) return 1;
  if(func_is_number(f)) return 1;
  if(func_is_add(f) && func_power(f)==1){
    for(i=0; i<func_asize(f); i++){
      if(!func_is_mono_not_num(f->a[i]) && !func_is_number(f->a[i])) return 0;
    }
    return 1;
  }
  return 0;
}

int func_is_poly_list(func_t *f)                                 //listかつ多項式であるか
{
  int i;
  if(f==NULL) return 0;
  if(!func_is_list(f)) return 0;
  for(i=0;i<func_asize(f);i++){
    if(!func_is_poly(f->a[i])) return 0;
  }
  return 1;
}

/////////////////////////////////////////////////////////////////////

func_t *func_poly_get_lt(func_t *f)
{
  if     (!func_is_poly(f))                                   { return NULL; }
  else if(func_is_add(f) && func_power(f)==1 && func_asize(f)>0){ return f->a[0]; }
  return f;
}

func_t *func_poly_get_lc(func_t *f)
{
  return func_mono_get_coeff(func_poly_get_lt(f));
}

func_t *func_poly_get_lm(func_t *f)
{
  return func_mono_get_var(func_poly_get_lt(f));
}

/////////////////////////////////////////////////////////////////////

//n次の項をget x^n
func_t *func_poly_get_mono_ntarm(func_t *f,int n)
{ 
  int i,j;
  if(func_is_mono(f)){
    j=func_mono_degree(f);
    if(n==j){ return f; }else{ return NULL; }
  }
  else if(func_is_poly(f)){
    for(i=0; i<func_asize(f); i++){
      j=func_mono_degree(f->a[i]);
      if(n==j){ return f->a[i]; }
    }
  }
  return NULL;
}

func_t *func_poly_clone_coeff_ntarm(func_t *f,int n)
{
  func_t *a;
  if(!func_is_poly(f)){ return NULL; }
  a=func_poly_get_mono_ntarm(f,n);
  if(a==NULL){ return func_zero(); }
  if(func_is_number(a)){ return func_clone(FR(a)); }
  a=func_mono_get_coeff(a);
  if(a==NULL){ return func_one(); }
  return func_clone(FR(a));
}

/////////////////////////////////////////////////////////////////////

//変数の個数を数える
int func_poly_n_var(func_t *f)
{
  int i,j,k,a,b,count=0;
  func_t *g=NULL,*h=NULL;
  if(func_is_mono(f)) return func_mono_n_var(f);
  if(!func_is_poly(f)) return -1;
  for(i=0;i<func_asize(f);i++){
    a=func_mono_n_var(f->a[i]);
    count=count+a;
  }
  for(i=0;i<func_asize(f)-1;i++){
    g=func_mono_get_var(f->a[i]);
    h=func_mono_get_var(f->a[i+1]);
    if(g!=NULL && h!=NULL){
      for(j=0;j<func_var_size(g); j++){
	a=g->p.var->num[j];
	for(k=0;k<func_var_size(h); k++){
	  b=h->p.var->num[k];
	  if(a==b) count=count-1;
	}
      }
    }
    g=NULL;
    h=NULL;
  }
  g=func_del(g);
  h=func_del(h);    
  return count;
}

//第n関数が存在するか
int func_poly_exist_varn(func_t *f, int n)
{
  int i;
  if(func_is_mono(f)) return func_mono_exist_varn(f,n);
  if(!func_is_poly(f)) return 0;
  for(i=0; i<func_asize(f); i++){
    if(func_mono_exist_varn(f->a[i],n)){ return 1; }
  }
  return 0;
}

//////////////////////////////////////////////////////////////////////////

//１変数多項式の判定 文字番号出力
int func_poly_var1n(func_t *f)
{
  int i,a,b;
  if     (func_is_number(f))             { return -1; }
  else if(func_is_mono_not_num(f))       { return func_mono_var1n(f); }
  else if(!func_is_poly(f))              { return -1; }
  else if(!(func_is_add(f) && func_power(f)==1)){ return -1; }
  a=func_mono_var1n(f->a[0]);
  for(i=1; i<func_asize(f); i++){
    b=func_mono_var1n(f->a[i]);
    if(!func_is_number(f->a[i]) && a!=b) return -1;
  }
  return a;
}

//////////////////////////////////////////////////////////////////////////

func_t *func_poly_lm_lcm(func_t *f, func_t *g)
{
  func_t *lcm=NULL,*f_lm=NULL,*g_lm=NULL;
  if(func_is_poly(f) && func_is_poly(g)){
    f_lm=func_poly_get_lm(f);
    g_lm=func_poly_get_lm(g);
    if     (f_lm!=NULL && g_lm!=NULL) { lcm=func_var_lcm(FR(f_lm),FR(g_lm)); }
    else if(f_lm==NULL && g_lm!=NULL) { lcm=FR(g_lm); }
    else if(f_lm!=NULL && g_lm==NULL) { lcm=FR(f_lm); }
    else                              { lcm=func_one(); }
  }
  f_lm=NULL;
  g_lm=NULL;
  f=func_del(f);
  g=func_del(g);
  return lcm;
}

////////////////////////////////////////////////////////////////////

int func_poly_can_div(func_t *f, func_t *g)
{
  return func_mono_can_div(func_poly_get_lt(f),func_poly_get_lt(g));
}

////////////////////////////////////////////////////////////////////

int func_poly_degree_max(func_t *f)
{
  int a,b,i;
  if(func_is_number(f)) return 0;
  if(func_is_mono_not_num(f)) return func_mono_degree(f);
  if(func_is_poly(f)){
    a=func_mono_degree(f->a[0]);
    for(i=1; i<func_asize(f); i++){
      b=func_mono_degree(f->a[i]);
      if(a<b) a=b;
    }
    return a;
  }
  return 0;
}

////////////////////////////////////////////////////////

//h={q,r}
func_t *func_poly_div_r_and_q(func_t *f, func_t *g)
{
  func_t *h=NULL;
  h=func_list(2);
  h->a[0]=func_poly_div_q(FR(f),FR(g));
  h->a[1]=func_poly_div_r(FR(f),FR(g));
  f=func_del(f);
  g=func_del(g);
  return h;
}

// q=f/g, where f=g*q+r
func_t *func_poly_div_q(func_t *f, func_t *g)
{
  func_t *r=NULL,*neg=NULL,*q=NULL,*a;
  if(func_is_poly(f) && func_is_poly(g)){
    neg=func_bigint_int(-1,1);
    r=FR(f);                                    // r=f
    q=func_zero();
    while(func_is_poly(r) && func_poly_can_div(r,g)) {
      a=func_div(FR(func_poly_get_lt(r)),
		 FR(func_poly_get_lt(g))); // a=LT(r)/LT(g)
      q=func_add(q,FR(a));                      // q+=LT(r)/LT(g)
      a=func_mul(FR(neg),a);                    // a=(-1)*LT(r)/LT(g)
      a=func_expand(func_mul(a,FR(g)));         // a=((-1)*LT(r)/LT(g))*g
      r=func_add(r,a);                          // r-=(LT(r)/LT(g))*g
    }
  }
  f=func_del(f);
  g=func_del(g);
  neg=func_del(neg);
  r=func_del(r);
  return q;
}

// r=mod(f,g), where f=g*q+r
func_t *func_poly_div_r(func_t *f, func_t *g)
{
  func_t *r=NULL,*neg=NULL,*a;
  if(func_is_poly(f) && func_is_poly(g)){
    neg=func_bigint_int(-1,1);
    r=FR(f);                                    // r=f
    while(func_is_poly(r) && func_poly_can_div(r,g)) {
      a=func_div(FR(func_poly_get_lt(r)),
		 FR(func_poly_get_lt(g))); // a=LT(r)/LT(g)
      a=func_mul(FR(neg),a);                    // a=(-1)*LT(r)/LT(g)
      a=func_expand(func_mul(a,FR(g)));         // a=((-1)*LT(r)/LT(g))*g
      r=func_add(r,a);                                   // r-=(LT(r)/LT(g))*g
    }
  }
  f=func_del(f);
  g=func_del(g);
  neg=func_del(neg);
  return r;
}

//////////////////////////////////////////////////////////////

func_t *func_poly_rm_lt(func_t *f)
{
  func_t *g=NULL;
  if(func_is_poly(f)){
    if(func_is_add(f) && func_power(f)==1 && func_asize(f)>=2){
      g=func_clone(FR(f));
      g->a[0]=func_del(g->a[0]);
      g=func_add_eval(g);
    }else{
      g=func_zero();
    }
  }
  f=func_del(f);
  return g;
}

////////////////////////////////////////////////////////////////

func_t *func_poly_list_div_r(func_t *f, func_t *g)
{
  int i,flag;
  func_t *p=NULL,*r=NULL,*neg=NULL,*a=NULL;
  if(!(func_is_poly(f) && func_is_poly_list(g))){ FUNC_ERROR_ARG2("func_poly_list_div_r",f,g); }
  neg=func_bigint_int(-1,1);
  p=func_clone(FR(f));
  r=func_zero();
  while(!func_is_zero(p)){
    i=0;
    flag=0;
    while(i<func_asize(g) && flag==0){
      if(func_poly_can_div(p,g->a[i])){
	a=func_div(FR(func_poly_get_lt(p)),FR(func_poly_get_lt(g->a[i])));  // a:=LT(p)/LT(g[i])
	a=func_mul(FR(neg),a);                                                    // a:=-LT(p)/LT(g[i])
	a=func_expand(func_mul(a,func_poly_rm_lt(FR(g->a[i]))));            // a:=expand((-LT(p)/LT(g[i]))*(g[i]-LT(g[i])))
	p=func_poly_rm_lt(p);                                                     // p:=p-LT(p)
	p=func_add(p,FR(a));                                                      // p:=(p-LT(p))-(LT(p)/LT(g[i]))*(g[i]-LT(g[i])))=p-(LT(p)/LT(g[i]))*g[i]
	a=func_del(a);
	flag=1;
      }else{
	i=i+1;
      }
    }
    if(flag==0){
      r=func_add(r,FR(func_poly_get_lt(p)));
      p=func_poly_rm_lt(p);
    }
  }
  p=func_del(p);
  f=func_del(f);
  g=func_del(g);
  neg=func_del(neg);
  a=func_del(a);
  return r;
}

/////////////////////////////////////////////////////////////

func_t *func_poly_monic(func_t *f)
{
  func_t *a=NULL;
  if(!(func_is_poly(f))){ FUNC_ERROR_ARG1("func_poly_monic",f); }
  a=func_poly_get_lc(f);
  if(a==NULL || func_is_one(a)){ return f; }  
  return func_expand(func_mul(func_inv(FR(a)),f));
}

func_t *func_poly_list_monic(func_t *f)
{
  int i;
  func_t *g=NULL;
  if(!(func_is_poly_list(f))){ FUNC_ERROR_ARG1("func_poly_list_monic",f); }
  g=func_clone(FR(f));
  for(i=0; i<func_asize(g); i++){
    func_aset(g,i,func_poly_monic(FR(g->a[i])));
  }
  f=func_del(f);
  return g;
}



//EOF
