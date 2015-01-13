#include<stdio.h>
#include<stdlib.h>
#include"is_macros.h"
#include"is_ivec.h"
#include"is_strings.h"
#include"is_print.h"
#include"is_func.h"
#include"is_rmulti.h"
#include"is_cmulti.h"

#define FR(F)      func_retain(F)
#define FGET_SP(f) (func_get_split_part(f))
#define FGET_NB(f) (func_evalf_part(f))


//////////////////////////////////////////////////////////

#define SPC " \t\n"
#define BLK0 "({["
#define BLK1 ")}]"
#define N0 "["
#define N1 "]"

func_t *func_list_script(func_t *f)
{
  int i,j;
  strings *a=NULL,*b=NULL;
  func_t *g=NULL;
  if(func_is_strings(f) && func_strings_size(f)==1){
    if(str_match(func_strings_at(f,0),"^{%m+}$",BLK0,BLK1)){
      a=strings_split_mask(func_strings_at(f,0),BLK0,BLK1,SPC);
      if(strings_size(a)==1){
	b=strings_split(strings_at(a,0),",",BLK0,BLK1,SPC);
	g=func_list(strings_size(b));
	for(i=0; i<func_asize(g); i++){
	  func_aset(g,i,func_script(strings_at(b,i)));
	}
      }
    }else if(str_match(func_strings_at(f,0),"^[%m+]$",BLK0,BLK1)){
      a=strings_split_mask(func_strings_at(f,0),N0,N1,SPC);
      if(strings_size(a)==1){
	strings_item_replace(a,0,strings_split(strings_at(a,0),";",N0,N1,SPC));
	if(strings_size(a)==1){
	  b=strings_split(strings_at(a,0),", ",NULL,NULL,SPC);
	  g=func_list(strings_size(b));
	  for(i=0; i<func_asize(g); i++){
	    func_aset(g,i,func_script(strings_at(b,i)));
	  }
	}else if(strings_size(a)>1){
	  g=func_list(strings_size(a));
	  for(i=0; i<strings_size(a); i++){
	    b=strings_split(strings_at(a,i),", ",NULL,NULL,SPC);
	    func_aset(g,i,func_list(strings_size(b)));	    
	    for(j=0; j<strings_size(b); j++){
	      func_aset(func_aget(g,i),j,func_script(strings_at(b,j)));
	    }
	  }
	}
      }
    }
  }
  if(g==NULL){ g=f; }else{ f=func_del(f); }
  a=strings_del(a);
  b=strings_del(b);
  return g;
}


////////////////////////////////////////////////////

void func_list_print(func_t *f)
{
  int i;
  printf("{");
  for(i=0; i<func_asize(f); i++){
    func_print(func_aget(f,i));
    if(i<(func_asize(f))-1) printf(", ");
  }
  printf("}");
}

////////////////////////////////////////////////////

func_t *func_var1_list(int n)
{
  int i=0;
  func_t *f=NULL;
  f=func_list(n);
  for(i=0; i<n; i++){
    func_aset(f,i,func_var1(i,1));
  }
  return f;
}

////////////////////////////////////////////////////

int func_list_real_array_size(func_t *f)
{
  int i,c;
  func_t *g=NULL;
  if(f==NULL || !func_is_list(f)){ FUNC_ERROR_ARG1("func_list_real_array_size",f); }
  if(func_asize(f)<=0){ return 0; }
  g=func_aget(f,0);
  if     (func_in_real(g)){ c=1; }
  else if(func_is_rvec(g)){ c=func_rvec_size(g); }
  else if(func_is_list(g) && func_list_real_array_size(g)==1){ c=func_asize(g); }
  else                    { c=0; }
  for(i=1; i<func_asize(f); i++){
    g=func_aget(f,i);
    if     (func_in_real(g) && c==1){}
    else if(func_is_rvec(g) && func_rvec_size(g)==c){}
    else if(func_is_list(g) && func_list_real_array_size(g)==1 && func_asize(g)==c){}
    else{ return 0; }
  }
  return c;
}

int func_list_complex_array_size(func_t *f)
{
  int i,c;
  func_t *g=NULL;
  if(f==NULL || !func_is_list(f)){ FUNC_ERROR_ARG1("func_list_complex_array_size",f); }
  if(func_asize(f)<=0){ return 0; }
  g=func_aget(f,0);
  if     (func_in_complex(g)){ c=1; }
  else if(func_is_rvec(g))   { c=func_rvec_size(g); }
  else if(func_is_cvec(g))   { c=func_cvec_size(g); }
  else if(func_is_list(g) && func_list_complex_array_size(g)==1){ c=func_asize(g); }
  else                       { c=0; }
  for(i=1; i<func_asize(f); i++){
    g=func_aget(f,i);
    if     (func_in_complex(g) && c==1){}
    else if(func_is_rvec(g) && func_rvec_size(g)==c){}
    else if(func_is_cvec(g) && func_cvec_size(g)==c){}
    else if(func_is_list(g) && func_list_complex_array_size(g)==1 && func_asize(g)==c){}
    else{ return 0; }
  }
  return c;
}

////////////////////////////////////////////////////


func_t *func_list_eval(func_t *g)
{
  int i,k;
  func_t *f=NULL;
  if(func_a_has_op(g,func_is_real) || func_a_has_op(g,func_is_rvec)){
    k=func_list_real_array_size(g);
    if(k==1){
      f=func_rvec(func_asize(g));
      for(i=0; i<func_asize(g); i++){ func_rvec_set(f,i,func_aget(g,i)); }
    }else if(k>1){
      f=func_rmat(func_asize(g),k);
      for(i=0; i<func_asize(g); i++){ func_rmat_set_row(f,i,func_aget(g,i)); }
    }
  }
  if(f==NULL && (func_a_has_op(g,func_is_complex) || func_a_has_op(g,func_is_cvec))){
    k=func_list_complex_array_size(g);
    if(k==1){
      f=func_cvec(func_asize(g));
      for(i=0; i<func_asize(g); i++){ func_cvec_set(f,i,func_aget(g,i)); }
    }else if(k>1){
      f=func_cmat(func_asize(g),k);
      for(i=0; i<func_asize(g); i++){ func_cmat_set_row(f,i,func_aget(g,i)); }
    }
  }
  if(f==NULL){ f=FR(g); }
  g=func_del(g);
  return f;
}

////////////////////////////////////////////////////

func_t *func_list_zeros(int n)
{
  int i;
  func_t *f=NULL;
  f=func_list(n);
  for(i=0; i<n; i++){ func_aset(f,i,func_zero()); }
  return f;
}

func_t *func_list_zeros2(int m, int n)
{
  int i;
  func_t *f=NULL;
  f=func_list(m);
  for(i=0; i<m; i++){ func_aset(f,i,func_list_zeros(n)); }
  return f;
}

///////////////////////////////////////////////

static const char *__func_list="@L";

///////////////////////////////////////////////

int func_is_list(func_t *f)
{
  return func_is(f,__func_list);
}

func_t *func_list(int n)
{
  func_t *f=NULL;
  f=func_new(__func_list);
  func_a_resize(f,n);
  return f;
}

func_t *func_op_list_new(void)
{
  func_t *f=func_builtin_new(__func_list);
  func_builtin_p(f)->order=FUNC_ORDER_LIST;
  func_builtin_p(f)->amax=FUNC_NOLIMIT;
  func_builtin_p(f)->eval=func_list_eval;
  func_builtin_p(f)->print=func_list_print;
  return f;
}

////////////////////////////////////////////////////////////////////////

//１変数多項式を見つける 式出力
func_t *func_list_find_var1n(func_t *f)
{
  int i;
  func_t *a=NULL;
  if(func_is_poly_list(f)){
    for(i=0;i<func_asize(f);i++){
      if(func_poly_var1n(f->a[i]) != -1){
	a=FR(f->a[i]);
      }
    }
  }
  f=func_del(f);
  return a;
}

int func_list_degree_max(func_t *f)
{
  int a,b,i;
  if(func_is_poly_list(f)){
    a=func_poly_degree_max(f->a[0]);
    for(i=1; i<func_asize(f); i++){
      b=func_poly_degree_max(f->a[i]);
      if(a<b) a=b;
    }
    return a;
  }
  return 0;
}

////////////////////////////////////////////////////////////////////////

func_t *func_list_concat(func_t *f, func_t *g)
{
  int i,k;
  func_t *a=NULL;
  k=0;
  if(func_is_list(f)){ k+=func_asize(f); }else{ k+=1; }
  if(func_is_list(g)){ k+=func_asize(g); }else{ k+=1; }
  a=func_list(k);
  k=0;
  if(func_is_list(f)){
    for(i=0; i<func_asize(f); i++){ func_aset(a,k++,FR(f->a[i])); }
  }else{
    func_aset(a,k++,FR(f));
  }
  if(func_is_list(g)){
    for(i=0; i<func_asize(g); i++){ func_aset(a,k++,FR(g->a[i])); }
  }else{
    func_aset(a,k++,FR(g));
  }
  f=func_del(f);
  g=func_del(g);
  return a;
}

////////////////////////////////////////////////////////////////////////
 
func_t *func_list_sol_convert_tree(func_t *f)
{
  int i,j,k;
  func_t *a=NULL,*b=NULL,*c=NULL,*g=NULL,*h=NULL;
  if(f!=NULL){
    a=func_list(0);
    for(i=0; i<func_asize(f); i++){
      g=FR(f->a[i]->a[2]);
      if(g!=NULL){
	for(j=0; j<func_asize(g); j++){
	  h=FR(g->a[j]->a[2]);
	  if(h!=NULL){
	    b=func_list_sol_convert_tree(FR(g));
	    for(k=0; k<func_asize(b); k++){
	      c=func_list(2);
	      c->a[0]=func_list_concat(FR(b->a[k]->a[0]),FR(f->a[i]->a[0]) );
	      c->a[1]=func_list_concat(FR(b->a[k]->a[1]),FR(f->a[i]->a[1]) );
	      func_a_append(a,FR(c));
	      c=func_del(c);
	    }
	    b=func_del(b);
	    h=func_del(h);
	    break;
	  }else{
	    b=func_list(2);
	    b->a[0]=func_list_concat(FR(g->a[j]->a[0]),FR(f->a[i]->a[0]));
	    b->a[1]=func_list_concat(FR(g->a[j]->a[1]),FR(f->a[i]->a[1]));
	    func_a_append(a,FR(b));
	    b=func_del(b); 
	  }
	  h=func_del(h);
	}
      }
      g=func_del(g);
    }
  }
  f=func_del(f);
  b=func_del(b);
  c=func_del(c);
  g=func_del(g);
  h=func_del(h);
  return a;
}

//EOF
