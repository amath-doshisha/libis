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

//////////////////////////////////////////////////////////////
//グレブナ基底
//////////////////////////////////////////////////////////////

func_t *func_gbasis_eval(func_t *g)
{
  func_t *f=NULL;
  if(g==NULL || !func_is(g,"gbasis") || func_asize(g)!=1){ FUNC_ERROR_ARG1("func_gbasis_eval",g); }
  f=func_poly_list_reduced_groebner(FR(func_aget(g,0)),0);
  g=func_del(g);
  return f;
}

/////////////////////////////////////////////////////////////


func_t *func_poly_s(func_t *f, func_t *g)
{
  func_t *lcm=NULL,*a=NULL,*b=NULL,*neg=NULL,*s=NULL;
  if(func_is_poly(f) && func_is_poly(g)){
    neg=func_bigint_int(-1,1);
    lcm=func_poly_lm_lcm(FR(f),FR(g));
    a=func_div(FR(lcm),FR(func_poly_get_lt(f)));
    b=func_div(FR(lcm),FR(func_poly_get_lt(g)));
    f=func_poly_rm_lt(f);
    g=func_poly_rm_lt(g);
    a=func_expand(func_mul(a,FR(f)));
    b=func_mul(FR(neg),b);
    b=func_expand(func_mul(b,FR(g)));
    s=func_expand(func_add(FR(a),FR(b)));
  }
  f=func_del(f);
  g=func_del(g);
  a=func_del(a);
  b=func_del(b);
  neg=func_del(neg);
  lcm=func_del(lcm);
  return s;
}


func_t *func_poly_list_buchberger(func_t *f, int debug)
{
  int i,j,done;
  func_t *s=NULL,*g=NULL,*h=NULL;
  if(!func_is_poly_list(f)){ FUNC_ERROR_ARG1("func_poly_list_berger",f); }
  // G:=F
  g=func_clone(FR(f));
  if(debug==1){ printf("func_poly_list_buchberger: %d",func_asize(g)); fflush(stdout); }
  // repeat
  done=0;
  while(!done){
    // sort G in ascending order
    func_args_sort(g);
    // make G monic
    g=func_poly_list_monic(g);
    if(debug>=2){ printf("sorted="); func_print(g); printf("\n"); }
    // H:=G
    h=func_clone(FR(g));
    if(debug>=2){ printf("<size=%d>\n",func_asize(h)); fflush(stdout); }
    // for each pair {p,q}, p!=q in H
    for(i=0; !done && i<func_asize(h)-1; i++){
      if(debug>=2){ printf("%04d: ",i); fflush(stdout); }
      for(j=i+1; !done && j<func_asize(h); j++){
	// S:=S(p,q)
	s=func_poly_s(FR(h->a[i]),FR(h->a[j]));
	// S:=remainder S(p,q) of H
	s=func_poly_list_div_r(s,FR(h));
	s=func_poly_monic(s);
	if(func_is_zero(s)){
	  if(debug>=2){ printf("0"); fflush(stdout); }
	}else if(func_args_have(g,s)){
	  if(debug>=2){ printf("-"); fflush(stdout); }	
	}else if(func_is_number(s)){
	  if(debug>=2){ printf("n"); fflush(stdout); }	  
	}else{
	  if(debug>=2){ printf("+"); fflush(stdout); }
	  func_a_append(g,FR(s));
	}
	s=func_del(s);
      }
      if(debug>=2){ printf("\n"); fflush(stdout); }
    }
    if(debug==1){ printf("+%d",func_asize(g)-func_asize(h)); fflush(stdout); }
    if(func_asize(h)==func_asize(g)){ done=1; }
    h=func_del(h);
  }
  if(debug==1){ printf("=%d\n",func_asize(g)); fflush(stdout); }
  f=func_del(f);
  s=func_del(s);
  return g;
}

func_t *func_poly_list_buchberger_simple(func_t *f, int debug)
{
  int i,j,k;
  func_t *s=NULL,*g=NULL,*B=NULL;
  if(!func_is_poly_list(f)){ FUNC_ERROR_ARG1("func_poly_list_berger",f); }
  // G:=F, sort G in ascending order, and make G monic
  g=func_clone(FR(f));
  func_args_sort(g);
  g=func_poly_list_monic(g);
  // create index-list B
  k=func_asize(g);
  B=func_list(k*(k-1)/2);
  for(k=0, i=0; i<func_asize(g)-1; i++){
    for(j=i+1; j<func_asize(g); j++){
      func_aset(B,k,func_ivec(2));
      B->a[k]->p.ivec->x[0]=i;
      B->a[k]->p.ivec->x[1]=j;
      k++;
    }
  }
  if(debug){ printf("size=%d ",func_asize(g)); fflush(stdout); }
  // repeat them until B is empty
  while(func_asize(B)>0){
    // select (i,j) in B
    i=B->a[0]->p.ivec->x[0];
    j=B->a[0]->p.ivec->x[1];
    if(debug){  fflush(stdout); }
    // S:=S(p,q)
    s=func_poly_s(FR(g->a[i]),FR(g->a[j]));
    // S:=remainder S(p,q) of G
    s=func_poly_list_div_r(s,FR(g));
    s=func_poly_monic(s);
    if(func_is_zero(s)){
      if(debug){ printf("(%d,%d)",i,j); fflush(stdout); }
    }else if(func_args_have(g,s)){
      if(debug){ printf("[%d,%d]",i,j); fflush(stdout); }	
    }else if(func_is_number(s)){
      if(debug){ printf("[[%d,%d]]",i,j); fflush(stdout); }	  
    }else{
      if(debug){ printf("<%d,%d>",i,j); fflush(stdout); }
      func_a_append(g,FR(s));	
      for(k=0; k<func_asize(g)-1; k++){
	func_a_append(B,func_ivec(2));
	B->a[func_asize(B)-1]->p.ivec->x[0]=k;
	B->a[func_asize(B)-1]->p.ivec->x[1]=func_asize(g)-1;
      }
    if(debug){ printf("\nsize=%d ",func_asize(g)); fflush(stdout); }
    }
    s=func_del(s);
    // next
    func_a_rm(B,0);
  }
  if(debug){ printf("\n"); fflush(stdout); }
  g=func_poly_list_monic(g);
  func_args_sort(g);
  f=func_del(f);
  s=func_del(s);
  B=func_del(B);
  return g;
}

func_t *func_poly_list_buchberger_selected(func_t *f, int debug)
{
  int i,j,k,*I=NULL,i0,j0,flag;
  func_t *s=NULL,*g=NULL,*B=NULL,*Bdone=NULL,*lcm=NULL,*a=NULL,*b=NULL;
  if(!func_is_poly_list(f)){ FUNC_ERROR_ARG1("func_poly_list_berger",f); }
  // G:=F, sort G in ascending order, and make G monic
  g=func_clone(FR(f));
  func_args_sort(g);
  g=func_poly_list_monic(g);
  // create the list of index
  k=func_asize(g);
  B=func_list(k*(k-1)/2);
  for(k=0, i=0; i<func_asize(g)-1; i++){
    for(j=i+1; j<func_asize(g); j++){
      func_aset(B,k,func_ivec(2));
      B->a[k]->p.ivec->x[0]=i;
      B->a[k]->p.ivec->x[1]=j;
      k++;
    }
  }
  Bdone=func_list(0);
  // create the list of LCMs
  lcm=func_list(func_asize(B));
  for(k=0; k<func_asize(lcm); k++){
    i=B->a[k]->p.ivec->x[0];
    j=B->a[k]->p.ivec->x[1];
    func_aset(lcm,k,func_var_lcm(FR(func_poly_get_lm(g->a[i])),FR(func_poly_get_lm(g->a[j]))));
  }
  I=ivec_allocate(func_asize(B));
  func_args_sort_index(lcm,I);
  func_args_arrange(B,I);
  // repeat them until B is empty
  if(debug){ printf("func_poly_list_buchberger_selected: "); fflush(stdout); }
  while(func_asize(B)>0){
    flag=0;
    // select (i,j) in B
    i=B->a[0]->p.ivec->x[0];
    j=B->a[0]->p.ivec->x[1];
    // a=LCM(g[i],g[j])/(LM(g[i]*LM(g[j]))
    a=func_div(FR(lcm->a[0]),func_mul(FR(g->a[i]),FR(g->a[j])));
    if(func_is_one(a)){
      flag=1;
      if(debug){ printf("-"); fflush(stdout); }
    }
    a=func_del(a);
    // criterion(f[i],f[j],B)
    if(!flag){
      for(k=0; !flag && k<func_asize(g); k++){
	if(k!=i && k!=j){
	  a=func_ivec(2);
	  b=func_ivec(2);
	  if(i<k){ a->p.ivec->x[0]=i; a->p.ivec->x[1]=k; }else{ a->p.ivec->x[0]=k; a->p.ivec->x[1]=i; }
	  if(i<k){ b->p.ivec->x[0]=j; b->p.ivec->x[1]=k; }else{ b->p.ivec->x[0]=k; b->p.ivec->x[1]=j; }
	  if(func_args_have(Bdone,a) && func_args_have(Bdone,b)){
	    a=func_del(a);
	    a=func_div(FR(lcm->a[0]),FR(func_poly_get_lm(g->a[k])));
	    if(func_is_one(a)){
	      flag=1;
	      if(debug){ printf("."); }
	    }
	  }
	  a=func_del(a);
	  b=func_del(b);
	}
      }
    }
    if(!flag){
      // S:=S(p,q)
      s=func_poly_s(FR(g->a[i]),FR(g->a[j]));
      // S:=remainder S(p,q) of G
      s=func_poly_list_div_r(s,FR(g));
      s=func_poly_monic(s);
      if(func_is_zero(s)){
	if(debug){ printf("0"); fflush(stdout); }
      }else if(func_args_have(g,s)){
	if(debug){ printf("="); fflush(stdout); }
      }else if(func_is_number(s)){
	if(debug){ printf("[[%d,%d]]",i,j); fflush(stdout); }	  
      }else{
	if(debug){ printf("*"); fflush(stdout); }
	func_a_append(g,FR(s));
	for(k=0; k<func_asize(g)-1; k++){
	  func_a_append(B,func_ivec(2));
	  i0=k;
	  j0=func_asize(g)-1;
	  B->a[func_asize(B)-1]->p.ivec->x[0]=i0;
	  B->a[func_asize(B)-1]->p.ivec->x[1]=j0;
	  func_a_append(lcm,func_var_lcm(FR(func_poly_get_lm(g->a[i0])),FR(func_poly_get_lm(g->a[j0]))));
	}
      }
      s=func_del(s);
    }
    // next
    func_a_append(Bdone,FR(B->a[0]));
    func_a_rm(B,0);
    func_a_rm(lcm,0);
    I=ivec_free(I);
    I=ivec_allocate(func_asize(B));
    func_args_sort_index(lcm,I);
    func_args_arrange(B,I);
  }
  if(debug){ printf("\n"); fflush(stdout); }
  g=func_poly_list_monic(g);
  func_args_sort(g);
  f=func_del(f);
  s=func_del(s);
  B=func_del(B);
  Bdone=func_del(Bdone);
  lcm=func_del(lcm);
  I=ivec_free(I);
  return g;
}

func_t *func_poly_list_minimal_groebner(func_t *f, int debug)
{
  int i,j;
  func_t *g=NULL,*a=NULL;
  g=func_poly_list_buchberger_selected(FR(f),debug);
  if(g!=NULL){
    for(i=0; i<func_asize(g); i++){
      for(j=0; g->a[i]!=NULL && j<func_asize(g); j++){
	if(i!=j && g->a[j]!=NULL && func_poly_can_div(g->a[i],g->a[j])){
	  g->a[i]=func_del(g->a[i]);
	}
      }
    }
    func_a_rm_null(g);
  }
  func_args_sort(g);
  f=func_del(f);
  a=func_del(a);
  return g;
}

func_t *func_poly_list_reduced(func_t *f)
{
  int i;
  func_t *g=NULL,*h=NULL;
  g=func_clone(FR(f));
  func_args_sort(g);
  for(i=0; i<func_asize(g); i++){
    h=func_clone(FR(g));
    func_a_rm(h,i);
    func_aset(g,i,func_poly_list_div_r(FR(g->a[i]),FR(h)));
    h=func_del(h);
    if(func_is_zero(g->a[i])){
      func_a_rm(g,i);
      func_a_rm_null(g);
      i--;
    }
  }
  func_args_sort(g);
  f=func_del(f);
  h=func_del(h);
  return g;
}

func_t *func_poly_list_reduced_groebner(func_t *f, int debug)
{
  f=func_poly_list_minimal_groebner(f,debug);
  f=func_poly_list_reduced(f);
  return f;
}

//EOF
