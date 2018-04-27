#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"is_macros.h"
#include"is_rmulti.h"
#include"is_bigint.h"
#include"is_strings.h"

#define BIGINT_EXIT_ROUNDING(S) {if(e){printf("ROUDING in %s\n",(S));exit(-1);};}

////////////////////////////////////////

// x=(bigint)value
int rset_bi(rmulti *x, bigint *value)
{
  return rdiv(x,BIGINT_NUM(value),BIGINT_DEN(value));
}

////////////////////////////////////////

bigint *bigint_allocate(void)
{
  bigint *x=NULL;
  x=(bigint*)malloc(sizeof(bigint));
  BIGINT_NUM(x)=rallocate_prec(BIGINT_DEFAULT_PREC);
  BIGINT_DEN(x)=rallocate_prec(BIGINT_DEFAULT_PREC);
  return x;
}

bigint *bigint_allocate_int(int num, int den)
{
  bigint *x=NULL;
  x=bigint_allocate();
  bigint_set_int(x,num,den);
  return x;
}

bigint *bigint_allocate_str(const char *num, const char *den)
{
  bigint *x=NULL;
  x=bigint_allocate();
  bigint_set_str(x,num,den);
  return x;
}

bigint *bigint_allocate_script(const char *str)
{
  bigint *x=NULL;
  x=bigint_allocate();
  bigint_set_script(x,str);
  return x;
}

// y=x
bigint *bigint_allocate_clone(bigint *x)
{
  bigint *y=NULL;
  y=(bigint*)malloc(sizeof(bigint));
  BIGINT_NUM(y)=rallocate_clone(BIGINT_NUM(x));
  BIGINT_DEN(y)=rallocate_clone(BIGINT_DEN(x));
  return y;
}

void bigint_check(bigint *x)
{
  int e=0,s=1;
  rmulti *gcd=NULL;
  // special cases
  if     (ris_nan(BIGINT_NUM(x))  || ris_nan(BIGINT_DEN(x)))  { bigint_set_nan(x);  return; }
  else if(ris_zero(BIGINT_NUM(x)) && ris_zero(BIGINT_DEN(x))) { bigint_set_nan(x);  return; }
  else if(ris_zero(BIGINT_DEN(x)))                            { bigint_set_inf(x);  return; }
  else if(ris_zero(BIGINT_NUM(x)))                            { bigint_set_zero(x); return; }
  else if(req(BIGINT_NUM(x),BIGINT_DEN(x)))                   { bigint_set_one(x);  return; }
  // check sgn
  if(rget_sgn(BIGINT_NUM(x))<0)                               { s*=-1; e+=rneg(BIGINT_NUM(x),BIGINT_NUM(x)); }
  if(rget_sgn(BIGINT_DEN(x))<0)                               { s*=-1; e+=rneg(BIGINT_DEN(x),BIGINT_DEN(x)); }  
  // check common factors
  gcd=rallocate_prec(MAX2(rget_prec(BIGINT_NUM(x)),rget_prec(BIGINT_DEN(x))));
  e+=bigint_gcd(gcd,BIGINT_NUM(x),BIGINT_DEN(x));
  if(!ris_zero(gcd) && !req_si(gcd,1)){
    e+=rdiv(BIGINT_NUM(x),BIGINT_NUM(x),gcd);
    e+=rdiv(BIGINT_DEN(x),BIGINT_DEN(x),gcd);
  }
  // set sgn
  if(s<0){ e+=rneg(BIGINT_NUM(x),BIGINT_NUM(x)); }
  // done
  BIGINT_EXIT_ROUNDING("bigint_check()");
  gcd=rfree(gcd);
}

int bigint_gcd(rmulti *gcd, rmulti *m0, rmulti *n0)
{
  int e=0,mode;
  rmulti *q=NULL,*r=NULL,*m=NULL,*n=NULL;
  if(ris_zero(n0))    { return e+=rset_si(gcd,0); }
  if(req_si(n0,1))    { return e+=rset_si(gcd,1); }
  if     (rlt(m0,n0)) { return bigint_gcd(gcd,n0,m0); }
  else if(req(m0,n0)) { e+=rcopy(gcd,n0); }
  else{
    mode=get_auto_prec_mode();
    set_auto_prec_enabled();
    r=rallocate_prec(rget_prec(m0));
    q=rallocate_prec(rget_prec(m0));
    m=rallocate_prec(rget_prec(m0));
    n=rallocate_prec(rget_prec(n0));
    e+=rcopy(m,m0); e+=rcopy(n,n0);
    do{
      rdiv(q,m,n); rfloor(q,q);          // q:=floor(m/n)
      e+=rcopy(r,m); e+=rsub_mul(r,q,n); // r:=m-q*n=m%n
      e+=rcopy(m,n); // m:=n
      e+=rcopy(n,r); // n:=r
    }while(!ris_zero(n));
    e+=rcopy(gcd,m); // gcd:=m
    set_auto_prec_mode(mode);
  }
  q=rfree(q);
  r=rfree(r);
  m=rfree(m);
  n=rfree(n);
  return e;
}

bigint *bigint_free(bigint *x)
{
  if(x==NULL) return NULL; 
  BIGINT_NUM(x)=rfree(BIGINT_NUM(x));
  BIGINT_DEN(x)=rfree(BIGINT_DEN(x));
  free(x);
  x=NULL;
  return x;
}

////////////////////////////////////////////

void bigint_set_nan(bigint *x)
{
  int e=0;
  rset_nan(BIGINT_NUM(x));
  e+=rset_si(BIGINT_DEN(x),1);
  BIGINT_EXIT_ROUNDING("bigint_set_nan()");
}

void bigint_set_inf(bigint *x)
{
  int e=0;
  rset_inf(BIGINT_NUM(x),rget_sgn(BIGINT_NUM(x)));
  e+=rset_si(BIGINT_DEN(x),1);
  BIGINT_EXIT_ROUNDING("bigint_set_inf()");
}

void bigint_set_zero(bigint *x)
{
  int e=0;
  e+=rset_si(BIGINT_NUM(x),0);
  e+=rset_si(BIGINT_DEN(x),1);
  BIGINT_EXIT_ROUNDING("bigint_set_zero()");
}

void bigint_set_one(bigint *x)
{
  int e=0;
  e+=rset_si(BIGINT_NUM(x),1);
  e+=rset_si(BIGINT_DEN(x),1);
  BIGINT_EXIT_ROUNDING("bigint_set_one()");
}

void bigint_set_int(bigint *x, int num, int den)
{
  int e=0;
  e+=rset_si(BIGINT_NUM(x),num);
  e+=rset_si(BIGINT_DEN(x),den);
  BIGINT_EXIT_ROUNDING("bigint_set_int()");
  bigint_check(x);
}

void bigint_set_str(bigint *x, const char *num, const char *den)
{
  int e=0,p_num,p_den;
  p_num=(((int)(((int)strlen(num))/(log10(2.0))))/BIGINT_DEFAULT_PREC+1)*BIGINT_DEFAULT_PREC;
  p_den=(((int)(((int)strlen(den))/(log10(2.0))))/BIGINT_DEFAULT_PREC+1)*BIGINT_DEFAULT_PREC;
  rround(BIGINT_NUM(x),p_num);
  rround(BIGINT_DEN(x),p_den);
  rset_s(BIGINT_NUM(x),num);
  rset_s(BIGINT_DEN(x),den);
  BIGINT_EXIT_ROUNDING("bigint_set_str()");
  bigint_check(x);
}

void bigint_set_script(bigint *x, const char *str)
{
  strings *list=NULL;
  list=strings_split(str,"/",NULL,NULL," \t\n");
  if(list!=NULL && list->n==1){
    bigint_set_str(x,list->str[0],"1");
  }else if(list!=NULL && list->n==2){
    bigint_set_str(x,list->str[0],list->str[1]);
  }else{
    printf("Error in void bigint_set_script(bigint *x, const char *str)\n");
    printf("str=%s",str);
    exit(0);
  }
  list=strings_del(list);
}

///////////////////////////////////////////////////////

// y=x
void bigint_copy(bigint *y, bigint *x)
{
  int e=0,mode;
  mode=get_auto_prec_mode();
  set_auto_prec_enabled();
  e+=rcopy(BIGINT_NUM(y),BIGINT_NUM(x));
  e+=rcopy(BIGINT_DEN(y),BIGINT_DEN(x));
  BIGINT_EXIT_ROUNDING("bigint_copy()");
  set_auto_prec_mode(mode);
}

// x<->y
void bigint_swap(bigint *x, bigint *y)
{
  rswap(BIGINT_NUM(x),BIGINT_NUM(y));
  rswap(BIGINT_DEN(x),BIGINT_DEN(y));
}

///////////////////////////////////////

// z=-x
void bigint_neg(bigint *z, bigint *x)
{
  int e=0,mode;
  mode=get_auto_prec_mode();
  set_auto_prec_enabled();
  e+=rneg(BIGINT_NUM(z),BIGINT_DEN(x));
  e+=rcopy(BIGINT_DEN(z),BIGINT_NUM(x));
  BIGINT_EXIT_ROUNDING("bigint_neg()");
  set_auto_prec_mode(mode);
  bigint_check(z);
}


// z=x*y
void bigint_mul(bigint *z, bigint *x, bigint *y)
{
  int e=0,mode;
  bigint *a=NULL;
  a=bigint_allocate();
  mode=get_auto_prec_mode();
  set_auto_prec_enabled();
  e+=rmul(BIGINT_NUM(a),BIGINT_NUM(x),BIGINT_NUM(y));
  e+=rmul(BIGINT_DEN(a),BIGINT_DEN(x),BIGINT_DEN(y));
  bigint_swap(z,a);
  BIGINT_EXIT_ROUNDING("bigint_mul()");
  set_auto_prec_mode(mode);
  a=bigint_free(a);
  bigint_check(z);
}

// z=x/y
void bigint_div(bigint *z, bigint *x, bigint *y)
{
  int e=0,mode;
  bigint *a=NULL;
  a=bigint_allocate();
  mode=get_auto_prec_mode();
  set_auto_prec_enabled();
  e+=rmul(BIGINT_NUM(a),BIGINT_NUM(x),BIGINT_DEN(y));
  e+=rmul(BIGINT_DEN(a),BIGINT_DEN(x),BIGINT_NUM(y));
  bigint_swap(z,a);
  BIGINT_EXIT_ROUNDING("bigint_div()");
  set_auto_prec_mode(mode);
  a=bigint_free(a);
  bigint_check(z);
}

// z=1/x
void bigint_inv(bigint *z, bigint *x)
{
  int e=0,mode;
  bigint *a=NULL;
  a=bigint_allocate();
  mode=get_auto_prec_mode();
  set_auto_prec_enabled();
  e+=rcopy(BIGINT_NUM(a),BIGINT_DEN(x));
  e+=rcopy(BIGINT_DEN(a),BIGINT_NUM(x));
  bigint_swap(z,a);
  BIGINT_EXIT_ROUNDING("bigint_inv()");
  set_auto_prec_mode(mode);
  a=bigint_free(a);
  bigint_check(z);
}

////////////////////////////////

// z=x+y
void bigint_add(bigint *z, bigint *x, bigint *y)
{
  int e=0,mode;
  bigint *a=NULL;
  a=bigint_allocate();
  mode=get_auto_prec_mode();
  set_auto_prec_enabled();
  e+=rmul    (BIGINT_NUM(a),BIGINT_NUM(x),BIGINT_DEN(y));
  e+=radd_mul(BIGINT_NUM(a),BIGINT_DEN(x),BIGINT_NUM(y));
  e+=rmul    (BIGINT_DEN(a),BIGINT_DEN(x),BIGINT_DEN(y));
  bigint_swap(z,a);
  BIGINT_EXIT_ROUNDING("bigint_add()");
  set_auto_prec_mode(mode);
  a=bigint_free(a);
  bigint_check(z);
}

// z=x-y
void bigint_sub(bigint *z, bigint *x, bigint *y)
{
  int e=0,mode;
  bigint *a=NULL;
  a=bigint_allocate();
  mode=get_auto_prec_mode();
  set_auto_prec_enabled();
  e+=rmul    (BIGINT_NUM(a),BIGINT_NUM(x),BIGINT_DEN(y));
  e+=rsub_mul(BIGINT_NUM(a),BIGINT_DEN(x),BIGINT_NUM(y));
  e+=rmul    (BIGINT_DEN(a),BIGINT_DEN(x),BIGINT_DEN(y));
  bigint_swap(z,a);
  BIGINT_EXIT_ROUNDING("bigint_sub()");
  set_auto_prec_mode(mode);
  a=bigint_free(a);
  bigint_check(z);
}

// z=x^p
void bigint_pow_n(bigint *z, bigint *x, int p)
{
  int e=0,i,mode;
  bigint *a=NULL;
  if     (p==0){ bigint_set_one(z); return; }
  else if(p==1){ bigint_copy(z,x); return; }
  a=bigint_allocate_int(1,1);
  mode=get_auto_prec_mode();
  set_auto_prec_enabled();
  if(p<0){
    for(i=0; i<(-p); i++){ e+=rmul(BIGINT_NUM(a),BIGINT_NUM(a),BIGINT_DEN(x)); e+=rmul(BIGINT_DEN(a),BIGINT_DEN(a),BIGINT_NUM(x)); }
  }else{
    for(i=0; i<p; i++)   { e+=rmul(BIGINT_NUM(a),BIGINT_NUM(a),BIGINT_NUM(x)); e+=rmul(BIGINT_DEN(a),BIGINT_DEN(a),BIGINT_DEN(x)); }
  }
  bigint_swap(z,a);
  BIGINT_EXIT_ROUNDING("bigint_pow_n()");
  set_auto_prec_mode(mode);
  a=bigint_free(a);
  bigint_check(z);
}

int bigint_cmp(bigint *x, bigint *y)
{
  int e=0,value,mode;
  rmulti *a=NULL,*b=NULL;
  mode=get_auto_prec_mode();
  set_auto_prec_enabled();
  a=rallocate();
  b=rallocate();
  e+=rmul(a,BIGINT_NUM(x),BIGINT_DEN(y));
  e+=rmul(b,BIGINT_DEN(x),BIGINT_NUM(y));
  value=rcmp(a,b);
  a=rfree(a);
  b=rfree(b);
  BIGINT_EXIT_ROUNDING("bigint_cmp()");
  set_auto_prec_mode(mode);
  return value;
}

/////////////////////////////////////////////

int bigint_is_integer(bigint *x)
{
  return req_si(BIGINT_DEN(x),1);
}

int bigint_is_nan(bigint *x)
{
  return (ris_nan(BIGINT_NUM(x)) ||
	  ris_nan(BIGINT_DEN(x)) ||
	  (ris_zero(BIGINT_NUM(x)) && ris_zero(BIGINT_DEN(x))) ||
	  (ris_inf(BIGINT_NUM(x)) && ris_inf(BIGINT_DEN(x))));
}

int bigint_is_inf(bigint *x)
{
  return ((ris_inf(BIGINT_NUM(x)) && ris_number(BIGINT_DEN(x))) ||
	  (ris_number(BIGINT_NUM(x)) && ris_zero(BIGINT_DEN(x))));
}

int bigint_is_zero(bigint *x)
{
  return (ris_zero(BIGINT_NUM(x)) && ris_number(BIGINT_DEN(x)));
}

int bigint_is_one(bigint *x)
{
  return req(BIGINT_NUM(x),BIGINT_DEN(x));
}

int bigint_is_neg_one(bigint *x)
{
  return (rabs_eq(BIGINT_NUM(x),BIGINT_DEN(x)) && (rget_sgn(BIGINT_NUM(x))*rget_sgn(BIGINT_DEN(x)))<0);
}


/////////////////////////////////////////////


int bigint_sgn(bigint *x)
{
  return rget_sgn(BIGINT_NUM(x));
}

int bigint_get_si(bigint *x)
{
  if(bigint_is_integer(x)){
    return rget_si(BIGINT_NUM(x));
  }else{
    return 0;
  }
}

int bigint_get_rmulti(rmulti *z, bigint *x)
{
  return rdiv(z,BIGINT_NUM(x),BIGINT_DEN(x));
}

int bigint_get_cmulti(cmulti *z, bigint *x)
{
  int e=0;
  e+=ccopy_r(z,BIGINT_NUM(x));
  e+=cdiv_r2(z,z,BIGINT_DEN(x));
  return e;
}

/////////////////////////////////

void bigint_print(bigint *x)
{
  if(bigint_is_integer(x)){
    mpfr_printf("%.0Rf",BIGINT_NUM(x));
  }else{
    mpfr_printf("%.0Rf/%.0Rf",BIGINT_NUM(x),BIGINT_DEN(x));
  }
}

//EOF
