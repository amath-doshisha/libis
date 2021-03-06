#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"is_macros.h"
#include"is_rmulti.h"
#include"is_irmulti.h"
#include"is_icmulti.h"
#include"is_bigint.h"
#include"is_strings.h"

////////////////////////////////////////

// x=(bigint)value
void rset_bi(rmulti *x, bigint *value)
{
  rdiv_rr(x,BIGINT_NUM(value),BIGINT_DEN(value));
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

bigint *bigint_allocate_str(char *num, char *den)
{
  bigint *x=NULL;
  x=bigint_allocate();
  bigint_set_str(x,num,den);
  return x;
}

bigint *bigint_allocate_script(char *str)
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
  int s=1,e;
  rmulti *gcd=NULL;
  // special cases
  if     (ris_nan(BIGINT_NUM(x))  || ris_nan(BIGINT_DEN(x)))  { bigint_set_nan(x);  return; }
  else if(ris_zero(BIGINT_NUM(x)) && ris_zero(BIGINT_DEN(x))) { bigint_set_nan(x);  return; }
  else if(ris_zero(BIGINT_DEN(x)))                            { bigint_set_inf(x);  return; }
  else if(ris_zero(BIGINT_NUM(x)))                            { bigint_set_zero(x); return; }
  else if(eq_rr(BIGINT_NUM(x),BIGINT_DEN(x)))                   { bigint_set_one(x);  return; }
  // check sgn
  if(rget_sgn(BIGINT_NUM(x))<0)                               { s*=-1; rneg_r(BIGINT_NUM(x),BIGINT_NUM(x)); }
  if(rget_sgn(BIGINT_DEN(x))<0)                               { s*=-1; rneg_r(BIGINT_DEN(x),BIGINT_DEN(x)); }  
  // check common factors
  gcd=rallocate_prec(MAX2(rget_prec(BIGINT_NUM(x)),rget_prec(BIGINT_DEN(x))));
  bigint_gcd(gcd,BIGINT_NUM(x),BIGINT_DEN(x));
  if(!ris_zero(gcd) && !eq_rd(gcd,1)){
    if((e=rdiv_rr_rouding_check(BIGINT_NUM(x),BIGINT_NUM(x),gcd))!=0){ ERROR_AT; exit(-1); }
    if((e=rdiv_rr_rouding_check(BIGINT_DEN(x),BIGINT_DEN(x),gcd))!=0){ ERROR_AT; exit(-1); }
  }
  // set sgn
  if(s<0){ rneg_r(BIGINT_NUM(x),BIGINT_NUM(x)); }
  // done
  gcd=rmfree(gcd);
}

void bigint_gcd(rmulti *gcd, rmulti *m0, rmulti *n0)
{
  rmulti *q=NULL,*r=NULL,*m=NULL,*n=NULL;
  if(ris_zero(n0))    { rset_i(gcd,0); return; }
  if(eq_rd(n0,1))     { rset_i(gcd,1); return; }
  if     (lt_rr(m0,n0)) { bigint_gcd(gcd,n0,m0); return; }
  if(eq_rr(m0,n0))      { rclone_r(gcd,n0); return; }
  r=rallocate_prec(rget_prec(m0));
  q=rallocate_prec(rget_prec(m0));
  m=rallocate_prec(rget_prec(m0));
  n=rallocate_prec(rget_prec(n0));
  rclone_r(m,m0);                   // m:=m0
  rclone_r(n,n0);                   // n:=n0
  do{                             // repeat
    rdiv_rr(q,m,n); rfloor_r(q,q);     // q:=floor(m/n)
    rclone_r(r,m); rsub_mul_rr(r,q,n); // r:=m-q*n=m%n
    rclone_r(m,n);                  // m:=n
    rclone_r(n,r);                  // n:=r
  }while(!ris_zero(n));           // until n!=0
  rclone_r(gcd,m);                  // gcd:=m
  q=rmfree(q);
  r=rmfree(r);
  m=rmfree(m);
  n=rmfree(n);
}

bigint *bigint_free(bigint *x)
{
  if(x==NULL){ return NULL; }
  BIGINT_NUM(x)=rmfree(BIGINT_NUM(x));
  BIGINT_DEN(x)=rmfree(BIGINT_DEN(x));
  free(x);
  x=NULL;
  return x;
}

////////////////////////////////////////////

void bigint_set_nan(bigint *x)
{
  rset_nan(BIGINT_NUM(x));
  rset_i(BIGINT_DEN(x),1);
}

void bigint_set_inf(bigint *x)
{
  rset_inf(BIGINT_NUM(x),rget_sgn(BIGINT_NUM(x)));
  rset_i(BIGINT_DEN(x),1);
}

void bigint_set_zero(bigint *x)
{
  rset_i(BIGINT_NUM(x),0);
  rset_i(BIGINT_DEN(x),1);
}

void bigint_set_one(bigint *x)
{
  rset_i(BIGINT_NUM(x),1);
  rset_i(BIGINT_DEN(x),1);
}

void bigint_set_int(bigint *x, int num, int den)
{
  rset_i(BIGINT_NUM(x),num);
  rset_i(BIGINT_DEN(x),den);
  bigint_check(x);
}

void bigint_set_str(bigint *x, char *num, char *den)
{
  int p_num,p_den;
  p_num=(((int)(((int)strlen(num))/(log10(2.0))))/BIGINT_DEFAULT_PREC+1)*BIGINT_DEFAULT_PREC;
  p_den=(((int)(((int)strlen(den))/(log10(2.0))))/BIGINT_DEFAULT_PREC+1)*BIGINT_DEFAULT_PREC;
  rround(BIGINT_NUM(x),p_num);
  rround(BIGINT_DEN(x),p_den);
  rset_s(BIGINT_NUM(x),num);
  rset_s(BIGINT_DEN(x),den);
  bigint_check(x);
}

void bigint_set_script(bigint *x, char *str)
{
  strings *list=NULL;
  list=strings_split(str,"/",NULL,NULL," \t\n");
  if(list!=NULL && list->n==1){
    bigint_set_str(x,list->str[0],"1");
  }else if(list!=NULL && list->n==2){
    bigint_set_str(x,list->str[0],list->str[1]);
  }else{
    printf("Error in void bigint_set_script(bigint *x, char *str)\n");
    printf("str=%s",str);
    exit(0);
  }
  list=strings_del(list);
}


/**
 @brief bigint型から[z0,z1]へ型変換.
 */
void irset_bigint(rmulti *z0, rmulti *z1, bigint *x)
{
  mpfr_div(z0,BIGINT_NUM(x),BIGINT_DEN(x),MPFR_RNDD); // lower bound
  mpfr_div(z1,BIGINT_NUM(x),BIGINT_DEN(x),MPFR_RNDU); // upper bound
}

/**
 @brief bigint型から[z0,z1]へ型変換.
 */
void icset_bigint(cmulti *z0, cmulti *z1, bigint *x)
{
  irset_bigint(C_R(z0),C_R(z1),x);
  irset_d(C_I(z0),C_I(z1),0,0);
}


///////////////////////////////////////////////////////

// y=x
void bigint_clone(bigint *y, bigint *x)
{
  rclone_r(BIGINT_NUM(y),BIGINT_NUM(x));
  rclone_r(BIGINT_DEN(y),BIGINT_DEN(x));
}

// x<->y
void bigint_swap(bigint *x, bigint *y)
{
  rswap(BIGINT_NUM(x),BIGINT_NUM(y));
  rswap(BIGINT_DEN(x),BIGINT_DEN(y));
}

///////////////////////////////////////

// y=-x
void bigint_neg(bigint *y, bigint *x)
{
  bigint_clone(y,x);
  rneg_r(BIGINT_NUM(y),BIGINT_NUM(y));
  bigint_check(y);
}

// z=x*y
void bigint_mul(bigint *z, bigint *x, bigint *y)
{
  int e;
  bigint *a=NULL;
  a=bigint_allocate();
  if((e=rmul_rr_exact(BIGINT_NUM(a),BIGINT_NUM(x),BIGINT_NUM(y)))!=0){ ERROR_AT; exit(-1); }
  if((e=rmul_rr_exact(BIGINT_DEN(a),BIGINT_DEN(x),BIGINT_DEN(y)))!=0){ ERROR_AT; exit(-1); }
  bigint_check(a);
  bigint_clone(z,a);
  a=bigint_free(a);
}

// z=x/y
void bigint_div(bigint *z, bigint *x, bigint *y)
{
  int e;
  bigint *a=NULL;
  a=bigint_allocate();
  if((e=rmul_rr_exact(BIGINT_NUM(a),BIGINT_NUM(x),BIGINT_DEN(y))!=0)){ ERROR_AT; exit(-1); }
  if((e=rmul_rr_exact(BIGINT_DEN(a),BIGINT_DEN(x),BIGINT_NUM(y))!=0)){ ERROR_AT; exit(-1); }
  bigint_check(a);
  bigint_clone(z,a);
  a=bigint_free(a);
}

// z=1/x
void bigint_inv(bigint *z, bigint *x)
{
  bigint *a=NULL;
  a=bigint_allocate();
  rclone_r(BIGINT_NUM(a),BIGINT_DEN(x));
  rclone_r(BIGINT_DEN(a),BIGINT_NUM(x));
  bigint_check(a);  
  bigint_clone(z,a);
  a=bigint_free(a);
}

////////////////////////////////

// z=x+y
void bigint_add(bigint *z, bigint *x, bigint *y)
{
  int e;
  bigint *a=NULL;
  a=bigint_allocate();
  if((e=rmul_rr_exact    (BIGINT_NUM(a),BIGINT_NUM(x),BIGINT_DEN(y)))!=0){ ERROR_AT; exit(-1); }
  if((e=radd_mul_rr_exact(BIGINT_NUM(a),BIGINT_DEN(x),BIGINT_NUM(y)))!=0){ ERROR_AT; exit(-1); }
  if((e=rmul_rr_exact    (BIGINT_DEN(a),BIGINT_DEN(x),BIGINT_DEN(y)))!=0){ ERROR_AT; exit(-1); }
  bigint_check(a);
  bigint_clone(z,a);
  a=bigint_free(a);
}

// z=x-y
void bigint_sub(bigint *z, bigint *x, bigint *y)
{
  int e;
  bigint *a=NULL;
  a=bigint_allocate();
  if((e=rmul_rr_exact    (BIGINT_NUM(a),BIGINT_NUM(x),BIGINT_DEN(y)))!=0){ ERROR_AT; exit(-1); }
  if((e=rsub_mul_rr_exact(BIGINT_NUM(a),BIGINT_DEN(x),BIGINT_NUM(y)))!=0){ ERROR_AT; exit(-1); }
  if((e=rmul_rr_exact    (BIGINT_DEN(a),BIGINT_DEN(x),BIGINT_DEN(y)))!=0){ ERROR_AT; exit(-1); }
  bigint_check(a);
  bigint_clone(z,a);
  a=bigint_free(a);
}

// z=x^p
void bigint_pow_n(bigint *z, bigint *x, int p)
{
  int i,e;
  bigint *a=NULL;
  if(p==0){ bigint_set_one(z); return; }
  if(p==1){ bigint_clone(z,x); return; }
  a=bigint_allocate_int(1,1);
  if(p<0){
    for(i=0; i<(-p); i++){
      if((e=rmul_rr_exact(BIGINT_NUM(a),BIGINT_NUM(a),BIGINT_DEN(x)))!=0){ ERROR_AT; exit(-1); }
      if((e=rmul_rr_exact(BIGINT_DEN(a),BIGINT_DEN(a),BIGINT_NUM(x)))!=0){ ERROR_AT; exit(-1); }
    }
  }else{
    for(i=0; i<p; i++){
      if((e=rmul_rr_exact(BIGINT_NUM(a),BIGINT_NUM(a),BIGINT_NUM(x)))!=0){ ERROR_AT; exit(-1); }
      if((e=rmul_rr_exact(BIGINT_DEN(a),BIGINT_DEN(a),BIGINT_DEN(x)))!=0){ ERROR_AT; exit(-1); }
    }
  }
  bigint_check(a);
  bigint_clone(z,a);
  a=bigint_free(a);
}

int bigint_cmp(bigint *x, bigint *y)
{
  int value,e;
  rmulti *a=NULL,*b=NULL;
  a=rallocate();
  b=rallocate();
  if((e=rmul_rr_exact(a,BIGINT_NUM(x),BIGINT_DEN(y)))!=0){ ERROR_AT; exit(-1); }
  if((e=rmul_rr_exact(b,BIGINT_DEN(x),BIGINT_NUM(y)))!=0){ ERROR_AT; exit(-1); }
  value=cmp_rr(a,b);
  a=rmfree(a);
  b=rmfree(b);
  return value;
}

/////////////////////////////////////////////

int bigint_is_integer(bigint *x)
{
  return eq_rd(BIGINT_DEN(x),1);
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
  return eq_rr(BIGINT_NUM(x),BIGINT_DEN(x));
}

int bigint_is_neg_one(bigint *x)
{
  return (eq_abs_rr(BIGINT_NUM(x),BIGINT_DEN(x)) && (rget_sgn(BIGINT_NUM(x))*rget_sgn(BIGINT_DEN(x)))<0);
}


/////////////////////////////////////////////

int bigint_sgn(bigint *x)
{
  return rget_sgn(BIGINT_NUM(x));
}

int bigint_get_si(bigint *x)
{
  if(bigint_is_integer(x)){
    return (int)rget_i(BIGINT_NUM(x));
  }else{
    return 0;
  }
}

void bigint_get_rmulti(rmulti *z, bigint *x)
{
  rdiv_rr(z,BIGINT_NUM(x),BIGINT_DEN(x));
}

void bigint_get_cmulti(cmulti *z, bigint *x)
{
  cset_r(z,BIGINT_NUM(x));
  cdiv_cr(z,z,BIGINT_DEN(x));
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
