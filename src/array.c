#include "is_array.h"
#include "is_strings.h"
#include "is_svec.h"
#include "is_ivec.h"
#include "is_dvec.h"
#include "is_zvec.h"
#include "is_rvec.h"
#include "is_cvec.h"
#include "is_irvec.h"
#include "is_icvec.h"
#include "is_imat.h"
#include "is_dmat.h"
#include "is_zmat.h"
#include "is_rmat.h"
#include "is_cmat.h"
#include "is_irmat.h"
#include "is_icmat.h"

/**
 @file  array.c
 @brief 配列型arrayに関する関数の定義
 @detail
  type='s':文字列 ← (get_char/get_*) → 'i','d','z','r','c','R','C'

  type='i':整数
　 ⇅ (get_int/get_doubel)
  type='d':倍精度実数  ←(get_double/get_multi)→  type='r':倍精度実数 ←(get_multi/get_imulti)→ type='R':倍精度実数の機械区間
　 ⇅ (get_real/get_complex)                     ⇅ (get_real/get_complex)                   ⇅ (get_real/get_complex)
  type='z':倍精度実数  ←(get_double/get_multi)→  type='c':倍精度実数 ←(get_multi/get_imulti)→ type='C':倍精度実数の機械区間

  all←(get)→all
 */


/** @name array型のメモリ割当に関する関数 */
/** @{ */

/**
 @breif array型のメモリ割当
 */ 
array *array_allocate(int type, int ndim, int *dim)
{
  int i;
  array *x=NULL;
  // allocate struct
  x=malloc(sizeof(array));
  if(!array_type_check(type)){ ERROR_EXIT("Error at type=%c\n",type); }
  // set type
  ARRAY_TYPE(x)=type;
  // set number of dimensions
  if(ndim<0){ ERROR_EXIT("Error at ndim=%d\n",ndim); }
  ARRAY_NDIM(x)=ndim;
  // set dimensions
  if(ARRAY_NDIM(x)<=0){
    ARRAY_DIM_P(x)=NULL;
  }else{
    if(dim==NULL){ ARRAY_DIM_P(x)=ivec_allocate(ARRAY_NDIM(x)); ivec_set_zeros(ARRAY_NDIM(x),ARRAY_DIM_P(x)); }
    else         { ARRAY_DIM_P(x)=ivec_allocate_clone(ARRAY_NDIM(x),dim); }
    for(i=0; i<ARRAY_NDIM(x); i++){ if(ARRAY_DIM_P(x)[i]<0){ ARRAY_DIM_P(x)[i]=0; } }
  }  
  // set LD
  if(ARRAY_NDIM(x)<=0){
    ARRAY_LD_P(x)=NULL;
  }else{
    ARRAY_LD_P(x)=ivec_allocate(ARRAY_NDIM(x));
    ARRAY_LD_P(x)[0]=ARRAY_DIM(x,0);    
    for(i=1; i<ARRAY_NDIM(x); i++){ ARRAY_LD_P(x)[i]=ARRAY_LD(x,i-1)*ARRAY_DIM(x,i); }
  }
  // set size
  if(ARRAY_NDIM(x)<=0){ ARRAY_SIZE(x)=0; }else{ ARRAY_SIZE(x)=ARRAY_LD_END(x); }
  if(ARRAY_SIZE(x)<0){ ERROR_EXIT("size=%d\n",ARRAY_SIZE(x)); }
  // allocate #1
  ARRAY_P0(x)=NULL;
  if(ARRAY_NDIM(x)>0 && ARRAY_SIZE(x)>0){
    if(type=='a'){ ARRAY_P0(x)=(void*)avec_allocate(ARRAY_SIZE(x)); }
    if(type=='s'){ ARRAY_P0(x)=(void*)svec_allocate(ARRAY_SIZE(x)); }
    if(type=='i'){ ARRAY_P0(x)=(void*)ivec_allocate(ARRAY_SIZE(x)); }
    if(type=='d'){ ARRAY_P0(x)=(void*)dvec_allocate(ARRAY_SIZE(x)); }
    if(type=='z'){ ARRAY_P0(x)=(void*)zvec_allocate(ARRAY_SIZE(x)); }
    if(type=='r'){ ARRAY_P0(x)=(void*)rvec_allocate(ARRAY_SIZE(x)); }
    if(type=='c'){ ARRAY_P0(x)=(void*)cvec_allocate(ARRAY_SIZE(x)); }
    if(type=='R'){ ARRAY_P0(x)=(void*)rvec_allocate(ARRAY_SIZE(x)); }
    if(type=='C'){ ARRAY_P0(x)=(void*)cvec_allocate(ARRAY_SIZE(x)); }
  }
  // allocate #1
  ARRAY_P1(x)=NULL;
  if(ARRAY_NDIM(x)>0 && ARRAY_SIZE(x)>0){
    if(type=='R'){ ARRAY_P1(x)=(void*)rvec_allocate(ARRAY_SIZE(x)); }
    if(type=='C'){ ARRAY_P1(x)=(void*)cvec_allocate(ARRAY_SIZE(x)); }
  }
  return x;
}

/**
 @breif array型のメモリの解放
 */ 
array *array_free(array *x)
{
  if(x==NULL){ return NULL; }
  if(ARRAY_TYPE(x)=='a'){ ARRAY_P0(x)=avec_free(ARRAY_SIZE(x),ARRAY_P0_A(x)); }
  if(ARRAY_TYPE(x)=='s'){ ARRAY_P0(x)=svec_free(ARRAY_SIZE(x),ARRAY_P0_S(x)); }
  if(ARRAY_TYPE(x)=='i'){ ARRAY_P0(x)=ivec_free(ARRAY_P0_I(x)); }
  if(ARRAY_TYPE(x)=='d'){ ARRAY_P0(x)=dvec_free(ARRAY_P0_D(x)); }
  if(ARRAY_TYPE(x)=='z'){ ARRAY_P0(x)=zvec_free(ARRAY_P0_Z(x)); }
  if(ARRAY_TYPE(x)=='r'){ ARRAY_P0(x)=rvec_free(ARRAY_SIZE(x),ARRAY_P0_R(x)); }
  if(ARRAY_TYPE(x)=='c'){ ARRAY_P0(x)=cvec_free(ARRAY_SIZE(x),ARRAY_P0_C(x)); }
  if(ARRAY_TYPE(x)=='R'){ ARRAY_P0(x)=rvec_free(ARRAY_SIZE(x),ARRAY_P0_R(x)); ARRAY_P1(x)=rvec_free(ARRAY_SIZE(x),ARRAY_P1_R(x)); }
  if(ARRAY_TYPE(x)=='C'){ ARRAY_P0(x)=cvec_free(ARRAY_SIZE(x),ARRAY_P0_C(x)); ARRAY_P1(x)=cvec_free(ARRAY_SIZE(x),ARRAY_P1_C(x)); }
  ARRAY_LD_P(x)=ivec_free(ARRAY_LD_P(x));
  ARRAY_DIM_P(x)=ivec_free(ARRAY_DIM_P(x));
  free(x);
  x=NULL;
  return x;
}

/**
 @breif array型の配列のメモリ割当
 */
array **avec_allocate(int n)
{  
  int i;
  array **x=NULL;
  x=(array**)malloc(sizeof(array*)*n);
  for(i=0; i<n; i++){ x[i]=NULL; }
  return x;
}

/**
 @breif array型の配列のメモリ解放
 */
array **avec_free(int n, array **x)
{  
  int i;
  if(x==NULL){ return NULL; }
  for(i=0; i<n; i++){
    if(x[i]!=NULL){  x[i]=array_free(x[i]); }
  }
  free(x);
  x=NULL;
  return x;
}

/**
 @breif array型のコピー y=x
 */ 
array *array_copy(array *x)
{
  array *y=NULL;
  if(x==NULL){ return y; }
  y=array_allocate(ARRAY_TYPE(x),ARRAY_NDIM(x),ARRAY_DIM_P(x));
  if(ARRAY_TYPE(y)=='s'){ svec_copy(ARRAY_SIZE(y),ARRAY_P0_S(y),ARRAY_P0_S(x)); }
  if(ARRAY_TYPE(y)=='i'){ ivec_copy(ARRAY_SIZE(y),ARRAY_P0_I(y),ARRAY_P0_I(x)); }
  if(ARRAY_TYPE(y)=='d'){ dvec_copy(ARRAY_SIZE(y),ARRAY_P0_D(y),ARRAY_P0_D(x)); }
  if(ARRAY_TYPE(y)=='z'){ zvec_copy(ARRAY_SIZE(y),ARRAY_P0_Z(y),ARRAY_P0_Z(x)); }
  if(ARRAY_TYPE(y)=='r'){ rvec_copy(ARRAY_SIZE(y),ARRAY_P0_R(y),ARRAY_P0_R(x)); }
  if(ARRAY_TYPE(y)=='c'){ cvec_copy(ARRAY_SIZE(y),ARRAY_P0_C(y),ARRAY_P0_C(x)); }
  if(ARRAY_TYPE(y)=='R'){ rvec_copy(ARRAY_SIZE(y),ARRAY_P0_R(y),ARRAY_P0_R(x)); rvec_copy(ARRAY_SIZE(y),ARRAY_P1_R(y),ARRAY_P1_R(x)); }
  if(ARRAY_TYPE(y)=='C'){ cvec_copy(ARRAY_SIZE(y),ARRAY_P0_C(y),ARRAY_P0_C(x)); cvec_copy(ARRAY_SIZE(y),ARRAY_P1_C(y),ARRAY_P1_C(x)); }
  return y;
}

/**
 @breif array型のコピー y=x
 */ 
array *array_clone(array *x)
{
  array *y=NULL;
  if(x==NULL){ return y; }
  y=array_allocate(ARRAY_TYPE(x),ARRAY_NDIM(x),ARRAY_DIM_P(x));
  if(ARRAY_TYPE(y)=='s'){ svec_copy (ARRAY_SIZE(y),ARRAY_P0_S(y),ARRAY_P0_S(x)); }
  if(ARRAY_TYPE(y)=='i'){ ivec_copy (ARRAY_SIZE(y),ARRAY_P0_I(y),ARRAY_P0_I(x)); }
  if(ARRAY_TYPE(y)=='d'){ dvec_copy (ARRAY_SIZE(y),ARRAY_P0_D(y),ARRAY_P0_D(x)); }
  if(ARRAY_TYPE(y)=='z'){ zvec_copy (ARRAY_SIZE(y),ARRAY_P0_Z(y),ARRAY_P0_Z(x)); }
  if(ARRAY_TYPE(y)=='r'){ rvec_clone(ARRAY_SIZE(y),ARRAY_P0_R(y),ARRAY_P0_R(x)); }
  if(ARRAY_TYPE(y)=='c'){ cvec_clone(ARRAY_SIZE(y),ARRAY_P0_C(y),ARRAY_P0_C(x)); }
  if(ARRAY_TYPE(y)=='R'){ rvec_clone(ARRAY_SIZE(y),ARRAY_P0_R(y),ARRAY_P0_R(x)); rvec_clone(ARRAY_SIZE(y),ARRAY_P1_R(y),ARRAY_P1_R(x)); }
  if(ARRAY_TYPE(y)=='C'){ cvec_clone(ARRAY_SIZE(y),ARRAY_P0_C(y),ARRAY_P0_C(x)); cvec_clone(ARRAY_SIZE(y),ARRAY_P1_C(y),ARRAY_P1_C(x)); }
  return y;
}

/** @} */

///////////////////////////////////////////////

/** @name array型の値の設定に関する関数 */
/** @{ */

/**
 @breif array型の値をすべて同じ値に設定
 */ 
void array_set_all_d(array *x, double y)
{
  int i;
  if(x==NULL){ return; }
  if(ARRAY_TYPE(x)=='a'){ for(i=0; i<ARRAY_SIZE(x); i++){ array_set_all_d(ARRAY_AVEC(x,i),y); } }
  else if(ARRAY_TYPE(x)=='s'){ svec_set_all_d(ARRAY_SIZE(x),ARRAY_P0_S(x),y,"%.7g"); }
  else if(ARRAY_TYPE(x)=='i'){ ivec_set_all  (ARRAY_SIZE(x),ARRAY_P0_I(x),(int)y); }
  else if(ARRAY_TYPE(x)=='d'){ dvec_set_all  (ARRAY_SIZE(x),ARRAY_P0_D(x),y); }
  else if(ARRAY_TYPE(x)=='z'){ zvec_set_all_d(ARRAY_SIZE(x),ARRAY_P0_Z(x),y); }
  else if(ARRAY_TYPE(x)=='r'){ rvec_set_all_d(ARRAY_SIZE(x),ARRAY_P0_R(x),y); }
  else if(ARRAY_TYPE(x)=='c'){ cvec_set_all_d(ARRAY_SIZE(x),ARRAY_P0_C(x),y); }
  else if(ARRAY_TYPE(x)=='R'){ rvec_set_all_d(ARRAY_SIZE(x),ARRAY_P0_R(x),y); rvec_set_all_d(ARRAY_SIZE(x),ARRAY_P1_R(x),y); }
  else if(ARRAY_TYPE(x)=='C'){ cvec_set_all_d(ARRAY_SIZE(x),ARRAY_P0_C(x),y); cvec_set_all_d(ARRAY_SIZE(x),ARRAY_P1_C(x),y); }
}

/**
 @breif array型の値をすべて0に設定
 */ 
void array_set_zeros(array *x)
{
  int i;
  if(x==NULL){ return; }
  if(ARRAY_TYPE(x)=='a'){ for(i=0; i<ARRAY_SIZE(x); i++){ array_set_zeros(ARRAY_AVEC(x,i)); } }
  else if(ARRAY_TYPE(x)=='s'){ svec_set_all  (ARRAY_SIZE(x),ARRAY_P0_S(x),"0"); }
  else if(ARRAY_TYPE(x)=='i'){ ivec_set_zeros(ARRAY_SIZE(x),ARRAY_P0_I(x)); }
  else if(ARRAY_TYPE(x)=='d'){ dvec_set_zeros(ARRAY_SIZE(x),ARRAY_P0_D(x)); }
  else if(ARRAY_TYPE(x)=='z'){ zvec_set_zeros(ARRAY_SIZE(x),ARRAY_P0_Z(x)); }
  else if(ARRAY_TYPE(x)=='r'){ rvec_set_zeros(ARRAY_SIZE(x),ARRAY_P0_R(x)); }
  else if(ARRAY_TYPE(x)=='c'){ cvec_set_zeros(ARRAY_SIZE(x),ARRAY_P0_C(x)); }
  else if(ARRAY_TYPE(x)=='R'){ rvec_set_zeros(ARRAY_SIZE(x),ARRAY_P0_R(x)); rvec_set_zeros(ARRAY_SIZE(x),ARRAY_P1_R(x)); }
  else if(ARRAY_TYPE(x)=='C'){ cvec_set_zeros(ARRAY_SIZE(x),ARRAY_P0_C(x)); cvec_set_zeros(ARRAY_SIZE(x),ARRAY_P1_C(x)); }
}

/**
 @breif array型の値をすべて1に設定
 */ 
void array_set_ones(array *x)
{
  int i;
  if(x==NULL){ return; }
  if(ARRAY_TYPE(x)=='a'){ for(i=0; i<ARRAY_SIZE(x); i++){ array_set_ones(ARRAY_AVEC(x,i)); } }
  else if(ARRAY_TYPE(x)=='s'){ svec_set_all (ARRAY_SIZE(x),ARRAY_P0_S(x),"1"); }
  else if(ARRAY_TYPE(x)=='i'){ ivec_set_ones(ARRAY_SIZE(x),ARRAY_P0_I(x)); }
  else if(ARRAY_TYPE(x)=='d'){ dvec_set_ones(ARRAY_SIZE(x),ARRAY_P0_D(x)); }
  else if(ARRAY_TYPE(x)=='z'){ zvec_set_ones(ARRAY_SIZE(x),ARRAY_P0_Z(x)); }
  else if(ARRAY_TYPE(x)=='r'){ rvec_set_ones(ARRAY_SIZE(x),ARRAY_P0_R(x)); }
  else if(ARRAY_TYPE(x)=='c'){ cvec_set_ones(ARRAY_SIZE(x),ARRAY_P0_C(x)); }
  else if(ARRAY_TYPE(x)=='R'){ rvec_set_ones(ARRAY_SIZE(x),ARRAY_P0_R(x)); rvec_set_ones(ARRAY_SIZE(x),ARRAY_P1_R(x)); }
  else if(ARRAY_TYPE(x)=='C'){ cvec_set_ones(ARRAY_SIZE(x),ARRAY_P0_C(x)); cvec_set_ones(ARRAY_SIZE(x),ARRAY_P1_C(x)); }
}

/**
 @breif array型の値をすべてnanに設定
 */ 
void array_set_nan(array *x)
{
  int i;
  if(x==NULL){ return; }
  if(ARRAY_TYPE(x)=='a'){ for(i=0; i<ARRAY_SIZE(x); i++){ array_set_nan(ARRAY_AVEC(x,i)); } }
  else if(ARRAY_TYPE(x)=='s'){ svec_set_all  (ARRAY_SIZE(x),ARRAY_P0_S(x),"NaN"); }
  else if(ARRAY_TYPE(x)=='i'){ ivec_set_zeros(ARRAY_SIZE(x),ARRAY_P0_I(x)); }
  else if(ARRAY_TYPE(x)=='d'){ dvec_set_nan  (ARRAY_SIZE(x),ARRAY_P0_D(x)); }
  else if(ARRAY_TYPE(x)=='z'){ zvec_set_nan  (ARRAY_SIZE(x),ARRAY_P0_Z(x)); }
  else if(ARRAY_TYPE(x)=='r'){ rvec_set_nan  (ARRAY_SIZE(x),ARRAY_P0_R(x)); }
  else if(ARRAY_TYPE(x)=='c'){ cvec_set_nan  (ARRAY_SIZE(x),ARRAY_P0_C(x)); }
  else if(ARRAY_TYPE(x)=='R'){ rvec_set_nan  (ARRAY_SIZE(x),ARRAY_P0_R(x)); rvec_set_nan(ARRAY_SIZE(x),ARRAY_P1_R(x)); }
  else if(ARRAY_TYPE(x)=='C'){ cvec_set_nan  (ARRAY_SIZE(x),ARRAY_P0_C(x)); cvec_set_nan(ARRAY_SIZE(x),ARRAY_P1_C(x)); }
}

/**
 @breif array型の値をすべてinfに設定
 */ 
void array_set_inf(array *x)
{
  int i;
  if(x==NULL){ return; }
  if(ARRAY_TYPE(x)=='a'){ for(i=0; i<ARRAY_SIZE(x); i++){ array_set_inf(ARRAY_AVEC(x,i)); } }
  else if(ARRAY_TYPE(x)=='s'){ svec_set_all(ARRAY_SIZE(x),ARRAY_P0_S(x),"Inf"); }
  else if(ARRAY_TYPE(x)=='i'){ ivec_set_all(ARRAY_SIZE(x),ARRAY_P0_I(x),((int)1)<<(8*(sizeof(int)-1))); }
  else if(ARRAY_TYPE(x)=='d'){ dvec_set_inf(ARRAY_SIZE(x),ARRAY_P0_D(x),1); }
  else if(ARRAY_TYPE(x)=='z'){ zvec_set_inf(ARRAY_SIZE(x),ARRAY_P0_Z(x),1,1); }
  else if(ARRAY_TYPE(x)=='r'){ rvec_set_inf(ARRAY_SIZE(x),ARRAY_P0_R(x),1); }
  else if(ARRAY_TYPE(x)=='c'){ cvec_set_inf(ARRAY_SIZE(x),ARRAY_P0_C(x),1,1); }
  else if(ARRAY_TYPE(x)=='R'){ rvec_set_inf(ARRAY_SIZE(x),ARRAY_P0_R(x),1); rvec_set_inf(ARRAY_SIZE(x),ARRAY_P1_R(x),1); }
  else if(ARRAY_TYPE(x)=='C'){ cvec_set_inf(ARRAY_SIZE(x),ARRAY_P0_C(x),1,1); cvec_set_inf(ARRAY_SIZE(x),ARRAY_P1_C(x),1,1); }
}

/**
 @breif array型の値を乱数値[a0,a1]に設定
 */ 
void array_set_rand(array *x, double a0, double a1)
{
  int i,*y=NULL;
  if(x==NULL){ return; }
  if(ARRAY_TYPE(x)=='a'){ for(i=0; i<ARRAY_SIZE(x); i++){ array_set_rand(ARRAY_AVEC(x,i),a0,a1); } }
  else if(ARRAY_TYPE(x)=='s'){
    y=ivec_allocate(ARRAY_SIZE(x));
    ivec_set_rand(ARRAY_SIZE(x),y,a0,a1);   
    svec_set_si(ARRAY_SIZE(x),ARRAY_P0_S(x),y);
    y=ivec_free(y);
  }
  else if(ARRAY_TYPE(x)=='i'){ ivec_set_rand(ARRAY_SIZE(x),ARRAY_P0_I(x),a0,a1); }
  else if(ARRAY_TYPE(x)=='d'){ dvec_set_rand(ARRAY_SIZE(x),ARRAY_P0_D(x),a1-a0,a0); }
  else if(ARRAY_TYPE(x)=='z'){ zvec_set_rand(ARRAY_SIZE(x),ARRAY_P0_Z(x),a1-a0,a0); }
  else if(ARRAY_TYPE(x)=='r'){ rvec_set_rand(ARRAY_SIZE(x),ARRAY_P0_R(x),a1-a0,a0); }
  else if(ARRAY_TYPE(x)=='c'){ cvec_set_rand(ARRAY_SIZE(x),ARRAY_P0_C(x),a1-a0,a0); }
  else if(ARRAY_TYPE(x)=='R'){ rvec_set_rand(ARRAY_SIZE(x),ARRAY_P0_R(x),a1-a0,a0); rvec_copy(ARRAY_SIZE(x),ARRAY_P1_R(x),ARRAY_P0_R(x)); }
  else if(ARRAY_TYPE(x)=='C'){ cvec_set_rand(ARRAY_SIZE(x),ARRAY_P0_C(x),a1-a0,a0); cvec_copy(ARRAY_SIZE(x),ARRAY_P1_C(x),ARRAY_P0_C(x)); }
}

/**
 @breif array型の値を0,1,2,..に設定
 */ 
void array_set_grid(array *x)
{
  int i,*y=NULL;
  if(x==NULL){ return; }
  if(ARRAY_TYPE(x)=='a'){ for(i=0; i<ARRAY_SIZE(x); i++){ array_set_grid(ARRAY_AVEC(x,i)); } }
  if(ARRAY_TYPE(x)=='s'){
    y=ivec_allocate(ARRAY_SIZE(x));
    ivec_set_grid(ARRAY_SIZE(x),y);
    svec_set_si(ARRAY_SIZE(x),ARRAY_P0_S(x),y);
    y=ivec_free(y);
  }
  else if(ARRAY_TYPE(x)=='i'){ ivec_set_grid(ARRAY_SIZE(x),ARRAY_P0_I(x)); }
  else if(ARRAY_TYPE(x)=='d'){ dvec_set_grid(ARRAY_SIZE(x),ARRAY_P0_D(x)); }
  else if(ARRAY_TYPE(x)=='z'){ zvec_set_grid(ARRAY_SIZE(x),ARRAY_P0_Z(x)); }
  else if(ARRAY_TYPE(x)=='r'){ rvec_set_grid(ARRAY_SIZE(x),ARRAY_P0_R(x)); }
  else if(ARRAY_TYPE(x)=='c'){ cvec_set_grid(ARRAY_SIZE(x),ARRAY_P0_C(x)); }
  else if(ARRAY_TYPE(x)=='R'){ rvec_set_grid(ARRAY_SIZE(x),ARRAY_P0_R(x)); rvec_set_grid(ARRAY_SIZE(x),ARRAY_P1_R(x)); }
  else if(ARRAY_TYPE(x)=='C'){ cvec_set_grid(ARRAY_SIZE(x),ARRAY_P0_C(x)); cvec_set_grid(ARRAY_SIZE(x),ARRAY_P1_C(x)); }
}

/** @} */

///////////////////////////////////////////////

/** @name array型の型確認に関する関数 */
/** @{ */

/**
 @breif array型の内部型の確認
 */ 
int array_type_check(int type)
{
  return (type=='a' || type=='s' || type=='i' || type=='d' || type=='z' || type=='r' || type=='c' || type=='R' || type=='C');
}

/**
 @breif array型は空か？
 */ 
int array_is_empty(array *x)
{
  return (x==NULL || ARRAY_NDIM(x)<=0 || ARRAY_SIZE(x)<=0 || ARRAY_DIM_P(x)==NULL || ARRAY_LD_P(x)==NULL || ARRAY_P0(x)==NULL);
}

/*
 * array型は同じサイズ？
 */
int array_same_dim_check(array *x, array *y)
{
  int n,i;
  if(x==NULL || y==NULL){ return 0; }
  n=MAX2(ARRAY_NDIM(x),ARRAY_NDIM(y));
  for(i=0; i<n; i++){ if(ARRAY_DIM(x,i)!=ARRAY_DIM(y,i)){ return 0; } }
  return 1;
}

int array_get_ndim(array *x)
{
  int i;
  if(x==NULL){ return 0; }
  for(i=ARRAY_NDIM(x)-1; i>=0; i--){
    if(ARRAY_DIM(x,i)!=1){ return i+1; }
  }
  return 0;
}

int array_get_subdim(array *x, array *y)
{
  int k,i,n;
  if(x==NULL || y==NULL){ return 0; }
  n=array_get_ndim(y);
  for(i=0; i<n; i++){
    if(ARRAY_DIM(x,i)!=ARRAY_DIM(y,i)){ break; }
  }
  k=i;
  for(i=k; k>0 && i<ARRAY_NDIM(x); i++){
    if(ARRAY_DIM(x,i)!=1){ k=0; }
  }
  return k;
}

/*
 * array型はスカラー？
 */
int array_is_scalar(array *x)
{
  int i;
  if(x==NULL){ return 0; }
  for(i=0; i<ARRAY_NDIM(x); i++){ if(ARRAY_DIM(x,i)!=1){ return 0; } }
  return 1;
}

/** @} */

///////////////////////////////////////////////

/** @name array型の入出力に関する関数 */
/** @{ */

/**
 @breif array型のメンバの内容を出力
 */
void array_put(array *x)
{
  int i;
  if(x==NULL){ printf("array_struct: NULL\n"); return; }
  printf("array_struct {\n");
  printf("  type: '%c'\n",ARRAY_TYPE(x));
  printf("  size: %d\n",ARRAY_SIZE(x));
  printf("  ndim: %d\n",ARRAY_NDIM(x));
  if(ARRAY_DIM_P(x)==NULL){ printf("  dim: NULL\n"); }
  else{ printf("  dim: ["); ivec_put(ARRAY_NDIM(x),ARRAY_DIM_P(x),","); printf("]\n"); }
  if(ARRAY_LD_P(x)==NULL){ printf("  LD: NULL\n"); }
  else{ printf("  LD:  ["); ivec_put(ARRAY_NDIM(x),ARRAY_LD_P(x),","); printf("]\n"); }
  printf("  p0: %p\n",ARRAY_P0(x));
  printf("  p1: %p\n",ARRAY_P1(x));
  printf("}\n");
  for(i=0; i<ARRAY_SIZE(x); i++){
    if(ARRAY_TYPE(x)=='a'){
      printf("array[%d]=",i);
      if(ARRAY_AVEC(x,i)==NULL){ printf("NULL\n"); }
      else{
        printf("array('%c',%d,[",ARRAY_TYPE(ARRAY_AVEC(x,i)),ARRAY_NDIM(ARRAY_AVEC(x,i)));
        ivec_put(ARRAY_NDIM(ARRAY_AVEC(x,i)),ARRAY_DIM_P(ARRAY_AVEC(x,i)),",");
        printf("])\n");
      }
    }
    else if(ARRAY_TYPE(x)=='s'){ printf("array[%d]=%s\n",i,ARRAY_SVEC(x,i)); }
    else if(ARRAY_TYPE(x)=='i'){ printf("array[%d]=%d\n",i,ARRAY_IVEC(x,i)); }
    else if(ARRAY_TYPE(x)=='d'){ printf("array[%d]=%.7g\n",i,ARRAY_DVEC(x,i)); }
    else if(ARRAY_TYPE(x)=='z'){ printf("array[%d]=(%.7g, %.7g)\n",i,Z_R(ARRAY_ZVEC(x,i)),Z_I(ARRAY_ZVEC(x,i))); }
    else if(ARRAY_TYPE(x)=='r'){ mpfr_printf("array[%d]=%.7Rg\n",i,ARRAY_RVEC(x,i)); }
    else if(ARRAY_TYPE(x)=='c'){ mpfr_printf("array[%d]=(%.7Rg, %.7Rg)\n",i,C_R(ARRAY_CVEC(x,i)),C_I(ARRAY_CVEC(x,i))); }
    else if(ARRAY_TYPE(x)=='R'){ mpfr_printf("array[%d]=[%.7Rg, %.7Rg]\n",i,ARRAY_RVEC0(x,i),ARRAY_RVEC1(x,i)); }
    else if(ARRAY_TYPE(x)=='C'){ mpfr_printf("array[%d]=[(%.7Rg, %.7Rg), (%.7Rg, %.7Rg)]\n",i,C_R(ARRAY_CVEC0(x,i)),C_I(ARRAY_CVEC0(x,i)),C_R(ARRAY_CVEC1(x,i)),C_I(ARRAY_CVEC1(x,i))); }
  }
}


/**
 @breif array型の出力
 */
void array_print(array *x, char *name, char format, int digits)
{
  int i,k,*I=NULL;
  if(array_is_empty(x)){
    if(name!=NULL){ printf("%s=NULL\n",name); }else{ printf("NULL\n"); }
    return;
  }
  if(!(format=='f' || format=='e' || format=='g')){ format='g'; }
  if(digits<0){ digits=6; }
  if(ARRAY_TYPE(x)=='a'){
    I=ivec_allocate(ARRAY_NDIM(x));
    for(k=0; k<ARRAY_SIZE(x); k++){
      array_get_index(I,x,k);
      if(name!=NULL){ printf("%s",name); }
      for(i=0; i<ARRAY_NDIM(x); i++){ printf("[%d]",I[i]); }
      if(ARRAY_AVEC(x,k)==NULL){
        printf("=NULL\n");
      }else{
        printf("=array('%c',%d,[",ARRAY_TYPE(ARRAY_AVEC(x,k)),ARRAY_NDIM(ARRAY_AVEC(x,k)));
        ivec_put(ARRAY_NDIM(ARRAY_AVEC(x,k)),ARRAY_DIM_P(ARRAY_AVEC(x,k)),",");
        printf("])\n");
      }
    }
  }else{
    if(ARRAY_NDIM(x)==1){
      if(name!=NULL){ printf("%s=\n",name); }
      if(ARRAY_TYPE(x)=='s'){ svec_print(ARRAY_SIZE(x),&ARRAY_SVEC(x,0),NULL); }
      if(ARRAY_TYPE(x)=='i'){ ivec_print(ARRAY_SIZE(x),&ARRAY_IVEC(x,0),NULL); }
      if(ARRAY_TYPE(x)=='d'){ dvec_print(ARRAY_SIZE(x),&ARRAY_DVEC(x,0),NULL,format,digits); }
      if(ARRAY_TYPE(x)=='z'){ zvec_print(ARRAY_SIZE(x),&ARRAY_ZVEC(x,0),NULL,format,digits); }
      if(ARRAY_TYPE(x)=='r'){ rvec_print(ARRAY_SIZE(x),&ARRAY_RVEC(x,0),NULL,format,digits); }
      if(ARRAY_TYPE(x)=='c'){ cvec_print(ARRAY_SIZE(x),&ARRAY_CVEC(x,0),NULL,format,digits); }
      if(ARRAY_TYPE(x)=='R'){ irvec_print(ARRAY_SIZE(x),&ARRAY_RVEC0(x,0),&ARRAY_RVEC1(x,0),NULL,format,digits); }
      if(ARRAY_TYPE(x)=='C'){ icvec_print(ARRAY_SIZE(x),&ARRAY_CVEC0(x,0),&ARRAY_CVEC1(x,0),NULL,format,digits); }
    }else{
      I=ivec_allocate(ARRAY_NDIM(x));
      for(k=0; k<ARRAY_SIZE(x); k+=ARRAY_LD1(x)){
        array_get_index(I,x,k);
        printf("%s%s",(name==NULL?"":name),(ARRAY_NDIM(x)==2?"":"[][]"));
        for(i=2; i<ARRAY_NDIM(x); i++){ printf("[%d]",I[i]); }
        printf("=\n");
        if(ARRAY_TYPE(x)=='s'){ smat_print(ARRAY_M(x),ARRAY_N(x),ARRAY_P0_S(x)+k,ARRAY_LD0(x),NULL); }
        if(ARRAY_TYPE(x)=='i'){ imat_print(ARRAY_M(x),ARRAY_N(x),ARRAY_P0_I(x)+k,ARRAY_LD0(x),NULL); }
        if(ARRAY_TYPE(x)=='d'){ dmat_print(ARRAY_M(x),ARRAY_N(x),ARRAY_P0_D(x)+k,ARRAY_LD0(x),NULL,format,digits); }
        if(ARRAY_TYPE(x)=='z'){ zmat_print(ARRAY_M(x),ARRAY_N(x),ARRAY_P0_Z(x)+k,ARRAY_LD0(x),NULL,format,digits); }
        if(ARRAY_TYPE(x)=='r'){ rmat_print(ARRAY_M(x),ARRAY_N(x),ARRAY_P0_R(x)+k,ARRAY_LD0(x),NULL,format,digits); }
        if(ARRAY_TYPE(x)=='c'){ cmat_print(ARRAY_M(x),ARRAY_N(x),ARRAY_P0_C(x)+k,ARRAY_LD0(x),NULL,format,digits); }
        if(ARRAY_TYPE(x)=='R'){ irmat_print(ARRAY_M(x),ARRAY_N(x),ARRAY_P0_R(x)+k,ARRAY_LD0(x),ARRAY_P1_R(x)+k,ARRAY_LD0(x),NULL,format,digits); }
        if(ARRAY_TYPE(x)=='C'){ icmat_print(ARRAY_M(x),ARRAY_N(x),ARRAY_P0_C(x)+k,ARRAY_LD0(x),ARRAY_P1_C(x)+k,ARRAY_LD0(x),NULL,format,digits); }
      }
      I=ivec_free(I);
    }
  }
}


/** @} */

///////////////////////////////////////////////

/** @name array型の要素に関する関数 */
/** @{ */

/**
 @brief array型の配列の取得
 */
void array_get_index(int *index, array *x, int t)
{
  int i;
  NULL_EXC1(x);
  if(ARRAY_NDIM(x)==1){ index[0]=t; }
  else{
    for(i=ARRAY_NDIM(x)-1; i>=1; i--){
      index[i]=t/ARRAY_LD(x,i-1);
      index[i-1]=t%ARRAY_LD(x,i-1);
      t=index[i-1];
    }
  }
}
  

/** @} */

///////////////////////////////////////////////

/** @name array型のキャストに関する関数 */
/** @{ */

/**
 @breif array型の型変換 y=cast(x)
 */ 
array *array_get(int type, array *x)
{
  array *y=NULL;
  if(x==NULL){ return y; }
  if(type=='d'){  
    if(ARRAY_TYPE(x)=='d'){ y=array_allocate('d',ARRAY_NDIM(x),ARRAY_DIM_P(x)); dvec_set   (ARRAY_SIZE(y),ARRAY_P0_D(y),ARRAY_P0_D(x)); }
    if(ARRAY_TYPE(x)=='z'){ y=array_allocate('d',ARRAY_NDIM(x),ARRAY_DIM_P(x)); dvec_set_z (ARRAY_SIZE(y),ARRAY_P0_D(y),ARRAY_P0_Z(x)); }
    if(ARRAY_TYPE(x)=='r'){ y=array_allocate('d',ARRAY_NDIM(x),ARRAY_DIM_P(x)); dvec_set_r (ARRAY_SIZE(y),ARRAY_P0_D(y),ARRAY_P0_R(x)); }
    if(ARRAY_TYPE(x)=='c'){ y=array_allocate('d',ARRAY_NDIM(x),ARRAY_DIM_P(x)); dvec_set_c (ARRAY_SIZE(y),ARRAY_P0_D(y),ARRAY_P0_C(x)); }
    if(ARRAY_TYPE(x)=='R'){ y=array_allocate('d',ARRAY_NDIM(x),ARRAY_DIM_P(x)); dvec_set_ir(ARRAY_SIZE(y),ARRAY_P0_D(y),ARRAY_P0_R(x),ARRAY_P1_R(x)); }
    if(ARRAY_TYPE(x)=='C'){ y=array_allocate('d',ARRAY_NDIM(x),ARRAY_DIM_P(x)); dvec_set_ic(ARRAY_SIZE(y),ARRAY_P0_D(y),ARRAY_P0_C(x),ARRAY_P1_C(x)); }
  }
  if(type=='z'){  
    if(ARRAY_TYPE(x)=='d'){ y=array_allocate('z',ARRAY_NDIM(x),ARRAY_DIM_P(x)); zvec_set_d (ARRAY_SIZE(y),ARRAY_P0_Z(y),ARRAY_P0_D(x)); }
    if(ARRAY_TYPE(x)=='z'){ y=array_allocate('z',ARRAY_NDIM(x),ARRAY_DIM_P(x)); zvec_set   (ARRAY_SIZE(y),ARRAY_P0_Z(y),ARRAY_P0_Z(x)); }
    if(ARRAY_TYPE(x)=='r'){ y=array_allocate('z',ARRAY_NDIM(x),ARRAY_DIM_P(x)); zvec_set_r (ARRAY_SIZE(y),ARRAY_P0_Z(y),ARRAY_P0_R(x)); }
    if(ARRAY_TYPE(x)=='c'){ y=array_allocate('z',ARRAY_NDIM(x),ARRAY_DIM_P(x)); zvec_set_c (ARRAY_SIZE(y),ARRAY_P0_Z(y),ARRAY_P0_C(x)); }
    if(ARRAY_TYPE(x)=='R'){ y=array_allocate('z',ARRAY_NDIM(x),ARRAY_DIM_P(x)); zvec_set_ir(ARRAY_SIZE(y),ARRAY_P0_Z(y),ARRAY_P0_R(x),ARRAY_P1_R(x)); }
    if(ARRAY_TYPE(x)=='C'){ y=array_allocate('z',ARRAY_NDIM(x),ARRAY_DIM_P(x)); zvec_set_ic(ARRAY_SIZE(y),ARRAY_P0_Z(y),ARRAY_P0_C(x),ARRAY_P1_C(x)); }
  }
  if(type=='r'){  
    if(ARRAY_TYPE(x)=='d'){ y=array_allocate('r',ARRAY_NDIM(x),ARRAY_DIM_P(x)); rvec_set_d (ARRAY_SIZE(y),ARRAY_P0_R(y),ARRAY_P0_D(x)); }
    if(ARRAY_TYPE(x)=='z'){ y=array_allocate('r',ARRAY_NDIM(x),ARRAY_DIM_P(x)); rvec_set_z (ARRAY_SIZE(y),ARRAY_P0_R(y),ARRAY_P0_Z(x)); }
    if(ARRAY_TYPE(x)=='r'){ y=array_allocate('r',ARRAY_NDIM(x),ARRAY_DIM_P(x)); rvec_set   (ARRAY_SIZE(y),ARRAY_P0_R(y),ARRAY_P0_R(x)); }
    if(ARRAY_TYPE(x)=='c'){ y=array_allocate('r',ARRAY_NDIM(x),ARRAY_DIM_P(x)); rvec_set_c (ARRAY_SIZE(y),ARRAY_P0_R(y),ARRAY_P0_C(x)); }
    if(ARRAY_TYPE(x)=='R'){ y=array_allocate('r',ARRAY_NDIM(x),ARRAY_DIM_P(x)); rvec_set_ir(ARRAY_SIZE(y),ARRAY_P0_R(y),ARRAY_P0_R(x),ARRAY_P1_R(x)); }
    if(ARRAY_TYPE(x)=='C'){ y=array_allocate('r',ARRAY_NDIM(x),ARRAY_DIM_P(x)); rvec_set_ic(ARRAY_SIZE(y),ARRAY_P0_R(y),ARRAY_P0_C(x),ARRAY_P1_C(x)); }
  }
  if(type=='c'){  
    if(ARRAY_TYPE(x)=='d'){ y=array_allocate('c',ARRAY_NDIM(x),ARRAY_DIM_P(x)); cvec_set_d (ARRAY_SIZE(y),ARRAY_P0_C(y),ARRAY_P0_D(x)); }
    if(ARRAY_TYPE(x)=='z'){ y=array_allocate('c',ARRAY_NDIM(x),ARRAY_DIM_P(x)); cvec_set_z (ARRAY_SIZE(y),ARRAY_P0_C(y),ARRAY_P0_Z(x)); }
    if(ARRAY_TYPE(x)=='r'){ y=array_allocate('c',ARRAY_NDIM(x),ARRAY_DIM_P(x)); cvec_set_r (ARRAY_SIZE(y),ARRAY_P0_C(y),ARRAY_P0_R(x)); }
    if(ARRAY_TYPE(x)=='c'){ y=array_allocate('c',ARRAY_NDIM(x),ARRAY_DIM_P(x)); cvec_set   (ARRAY_SIZE(y),ARRAY_P0_C(y),ARRAY_P0_C(x)); }
    if(ARRAY_TYPE(x)=='R'){ y=array_allocate('c',ARRAY_NDIM(x),ARRAY_DIM_P(x)); cvec_set_ir(ARRAY_SIZE(y),ARRAY_P0_C(y),ARRAY_P0_R(x),ARRAY_P1_R(x)); }
    if(ARRAY_TYPE(x)=='C'){ y=array_allocate('c',ARRAY_NDIM(x),ARRAY_DIM_P(x)); cvec_set_ic(ARRAY_SIZE(y),ARRAY_P0_C(y),ARRAY_P0_C(x),ARRAY_P1_C(x)); }
  }
  if(type=='R'){  
    if(ARRAY_TYPE(x)=='d'){ y=array_allocate('R',ARRAY_NDIM(x),ARRAY_DIM_P(x)); irvec_set_d(ARRAY_SIZE(y),ARRAY_P0_R(y),ARRAY_P1_R(y),ARRAY_P0_D(x),ARRAY_P0_D(x)); }
    if(ARRAY_TYPE(x)=='z'){ y=array_allocate('R',ARRAY_NDIM(x),ARRAY_DIM_P(x)); irvec_set_z(ARRAY_SIZE(y),ARRAY_P0_R(y),ARRAY_P1_R(y),ARRAY_P0_Z(x),ARRAY_P0_Z(x)); }
    if(ARRAY_TYPE(x)=='r'){ y=array_allocate('R',ARRAY_NDIM(x),ARRAY_DIM_P(x)); irvec_set  (ARRAY_SIZE(y),ARRAY_P0_R(y),ARRAY_P1_R(y),ARRAY_P0_R(x),ARRAY_P0_R(x)); }
    if(ARRAY_TYPE(x)=='c'){ y=array_allocate('R',ARRAY_NDIM(x),ARRAY_DIM_P(x)); irvec_set_c(ARRAY_SIZE(y),ARRAY_P0_R(y),ARRAY_P1_R(y),ARRAY_P0_C(x),ARRAY_P0_C(x)); }
    if(ARRAY_TYPE(x)=='R'){ y=array_allocate('R',ARRAY_NDIM(x),ARRAY_DIM_P(x)); irvec_set  (ARRAY_SIZE(y),ARRAY_P0_R(y),ARRAY_P1_R(y),ARRAY_P0_R(x),ARRAY_P1_R(x)); }
    if(ARRAY_TYPE(x)=='C'){ y=array_allocate('R',ARRAY_NDIM(x),ARRAY_DIM_P(x)); irvec_set_c(ARRAY_SIZE(y),ARRAY_P0_R(y),ARRAY_P1_R(y),ARRAY_P0_C(x),ARRAY_P1_C(x)); }
  }
  if(type=='C'){  
    if(ARRAY_TYPE(x)=='d'){ y=array_allocate('C',ARRAY_NDIM(x),ARRAY_DIM_P(x)); icvec_set_d(ARRAY_SIZE(y),ARRAY_P0_C(y),ARRAY_P1_C(y),ARRAY_P0_D(x),ARRAY_P0_D(x)); }
    if(ARRAY_TYPE(x)=='z'){ y=array_allocate('C',ARRAY_NDIM(x),ARRAY_DIM_P(x)); icvec_set_z(ARRAY_SIZE(y),ARRAY_P0_C(y),ARRAY_P1_C(y),ARRAY_P0_Z(x),ARRAY_P0_Z(x)); }
    if(ARRAY_TYPE(x)=='r'){ y=array_allocate('C',ARRAY_NDIM(x),ARRAY_DIM_P(x)); icvec_set_r(ARRAY_SIZE(y),ARRAY_P0_C(y),ARRAY_P1_C(y),ARRAY_P0_R(x),ARRAY_P0_R(x)); }
    if(ARRAY_TYPE(x)=='c'){ y=array_allocate('C',ARRAY_NDIM(x),ARRAY_DIM_P(x)); icvec_set  (ARRAY_SIZE(y),ARRAY_P0_C(y),ARRAY_P1_C(y),ARRAY_P0_C(x),ARRAY_P0_C(x)); }
    if(ARRAY_TYPE(x)=='R'){ y=array_allocate('C',ARRAY_NDIM(x),ARRAY_DIM_P(x)); icvec_set_r(ARRAY_SIZE(y),ARRAY_P0_C(y),ARRAY_P1_C(y),ARRAY_P0_R(x),ARRAY_P1_R(x)); }
    if(ARRAY_TYPE(x)=='C'){ y=array_allocate('C',ARRAY_NDIM(x),ARRAY_DIM_P(x)); icvec_set  (ARRAY_SIZE(y),ARRAY_P0_C(y),ARRAY_P1_C(y),ARRAY_P0_C(x),ARRAY_P1_C(x)); }
  }
  return y;  
}

/**
 @breif array型の型変換 y=char(x)
 */ 
array *array_get_char(array *x, char format, int digits)
{
  array *y=NULL;
  if(!(format=='f' || format=='e' || format=='g')){ format='g'; }
  if(digits<0){ digits=6; }
  if(x==NULL){ return y; }
  if(ARRAY_TYPE(x)=='s'){ y=array_allocate('s',ARRAY_NDIM(x),ARRAY_DIM_P(x));  svec_copy (ARRAY_SIZE(y),ARRAY_P0_S(y),ARRAY_P0_S(x)); }
  if(ARRAY_TYPE(x)=='i'){ y=array_allocate('s',ARRAY_NDIM(x),ARRAY_DIM_P(x));  ivec_get_s(ARRAY_SIZE(y),ARRAY_P0_S(y),ARRAY_P0_I(x)); }
  if(ARRAY_TYPE(x)=='d'){ y=array_allocate('s',ARRAY_NDIM(x),ARRAY_DIM_P(x));  dvec_get_s(ARRAY_SIZE(y),ARRAY_P0_S(y),ARRAY_P0_D(x),format,digits); }
  if(ARRAY_TYPE(x)=='z'){ y=array_allocate('s',ARRAY_NDIM(x),ARRAY_DIM_P(x));  zvec_get_s(ARRAY_SIZE(y),ARRAY_P0_S(y),ARRAY_P0_Z(x),format,digits); }
  if(ARRAY_TYPE(x)=='r'){ y=array_allocate('s',ARRAY_NDIM(x),ARRAY_DIM_P(x));  rvec_get_s(ARRAY_SIZE(y),ARRAY_P0_S(y),ARRAY_P0_R(x),format,digits); }
  if(ARRAY_TYPE(x)=='c'){ y=array_allocate('s',ARRAY_NDIM(x),ARRAY_DIM_P(x));  cvec_get_s(ARRAY_SIZE(y),ARRAY_P0_S(y),ARRAY_P0_C(x),format,digits); }
  if(ARRAY_TYPE(x)=='R'){ y=array_allocate('s',ARRAY_NDIM(x),ARRAY_DIM_P(x)); irvec_get_s(ARRAY_SIZE(y),ARRAY_P0_S(y),ARRAY_P0_R(x),ARRAY_P1_R(x),format,digits); }
  if(ARRAY_TYPE(x)=='C'){ y=array_allocate('s',ARRAY_NDIM(x),ARRAY_DIM_P(x)); icvec_get_s(ARRAY_SIZE(y),ARRAY_P0_S(y),ARRAY_P0_C(x),ARRAY_P1_C(x),format,digits); }
  return y;
}

/**
 @breif array型の型変換 y=int(x)
 */ 
array *array_get_int(array *x)
{
  array *y=NULL;
  if(x==NULL){ return y; }
  if(ARRAY_TYPE(x)=='s'){ y=array_allocate('i',ARRAY_NDIM(x),ARRAY_DIM_P(x));  ivec_set_s (ARRAY_SIZE(y),ARRAY_P0_I(y),ARRAY_P0_S(x)); }
  if(ARRAY_TYPE(x)=='i'){ y=array_allocate('i',ARRAY_NDIM(x),ARRAY_DIM_P(x));  ivec_copy  (ARRAY_SIZE(y),ARRAY_P0_I(y),ARRAY_P0_I(x)); }
  if(ARRAY_TYPE(x)=='d'){ y=array_allocate('i',ARRAY_NDIM(x),ARRAY_DIM_P(x));  dvec_get_si(ARRAY_SIZE(y),ARRAY_P0_I(y),ARRAY_P0_D(x)); }
  if(ARRAY_TYPE(x)=='z'){ y=array_allocate('i',ARRAY_NDIM(x),ARRAY_DIM_P(x));  zvec_get_si(ARRAY_SIZE(y),ARRAY_P0_I(y),ARRAY_P0_Z(x)); }
  if(ARRAY_TYPE(x)=='r'){ y=array_allocate('i',ARRAY_NDIM(x),ARRAY_DIM_P(x));  rvec_get_si(ARRAY_SIZE(y),ARRAY_P0_I(y),ARRAY_P0_R(x)); }
  if(ARRAY_TYPE(x)=='c'){ y=array_allocate('i',ARRAY_NDIM(x),ARRAY_DIM_P(x));  cvec_get_si(ARRAY_SIZE(y),ARRAY_P0_I(y),ARRAY_P0_C(x)); }
  if(ARRAY_TYPE(x)=='R'){ y=array_allocate('i',ARRAY_NDIM(x),ARRAY_DIM_P(x)); irvec_get_si(ARRAY_SIZE(y),ARRAY_P0_I(y),ARRAY_P0_R(x),ARRAY_P1_R(x)); }
  if(ARRAY_TYPE(x)=='C'){ y=array_allocate('i',ARRAY_NDIM(x),ARRAY_DIM_P(x)); icvec_get_si(ARRAY_SIZE(y),ARRAY_P0_I(y),ARRAY_P0_C(x),ARRAY_P1_C(x)); }
  return y;
}

/**
 @breif array型の型変換 y=real(x)
 */ 
array *array_get_real(array *x)
{
  array *y=NULL;
  if(x==NULL){ return y; }
  if(ARRAY_TYPE(x)=='s'){ y=array_allocate('s',ARRAY_NDIM(x),ARRAY_DIM_P(x));  svec_copy(ARRAY_SIZE(y),ARRAY_P0_S(y),ARRAY_P0_S(x)); }
  if(ARRAY_TYPE(x)=='i'){ y=array_allocate('i',ARRAY_NDIM(x),ARRAY_DIM_P(x));  ivec_copy(ARRAY_SIZE(y),ARRAY_P0_I(y),ARRAY_P0_I(x)); }
  if(ARRAY_TYPE(x)=='d'){ y=array_allocate('d',ARRAY_NDIM(x),ARRAY_DIM_P(x));  dvec_copy(ARRAY_SIZE(y),ARRAY_P0_D(y),ARRAY_P0_D(x)); }
  if(ARRAY_TYPE(x)=='z'){ y=array_allocate('d',ARRAY_NDIM(x),ARRAY_DIM_P(x));  zvec_real(ARRAY_SIZE(y),ARRAY_P0_D(y),ARRAY_P0_Z(x)); }
  if(ARRAY_TYPE(x)=='r'){ y=array_allocate('r',ARRAY_NDIM(x),ARRAY_DIM_P(x));  rvec_copy(ARRAY_SIZE(y),ARRAY_P0_R(y),ARRAY_P0_R(x)); }
  if(ARRAY_TYPE(x)=='c'){ y=array_allocate('r',ARRAY_NDIM(x),ARRAY_DIM_P(x));  cvec_real(ARRAY_SIZE(y),ARRAY_P0_R(y),ARRAY_P0_C(x)); }
  if(ARRAY_TYPE(x)=='R'){ y=array_allocate('R',ARRAY_NDIM(x),ARRAY_DIM_P(x)); irvec_copy(ARRAY_SIZE(y),ARRAY_P0_R(y),ARRAY_P1_R(y),ARRAY_P0_R(x),ARRAY_P1_R(x)); }
  if(ARRAY_TYPE(x)=='C'){ y=array_allocate('R',ARRAY_NDIM(x),ARRAY_DIM_P(x)); icvec_real(ARRAY_SIZE(y),ARRAY_P0_R(y),ARRAY_P1_R(y),ARRAY_P0_C(x),ARRAY_P1_C(x)); }
  return y;
}

/**
 @breif array型の型変換 y=complex(x)
 */ 
array *array_get_complex(array *x)
{
  array *y=NULL;
  if(x==NULL){ return y; }
  if(ARRAY_TYPE(x)=='s'){ y=array_allocate('s',ARRAY_NDIM(x),ARRAY_DIM_P(x));  svec_copy  (ARRAY_SIZE(y),ARRAY_P0_S(y),ARRAY_P0_S(x)); }
  if(ARRAY_TYPE(x)=='i'){ y=array_allocate('i',ARRAY_NDIM(x),ARRAY_DIM_P(x));  ivec_copy  (ARRAY_SIZE(y),ARRAY_P0_I(y),ARRAY_P0_I(x)); }
  if(ARRAY_TYPE(x)=='d'){ y=array_allocate('z',ARRAY_NDIM(x),ARRAY_DIM_P(x));  zvec_copy_d(ARRAY_SIZE(y),ARRAY_P0_Z(y),ARRAY_P0_D(x)); }
  if(ARRAY_TYPE(x)=='z'){ y=array_allocate('z',ARRAY_NDIM(x),ARRAY_DIM_P(x));  zvec_copy  (ARRAY_SIZE(y),ARRAY_P0_Z(y),ARRAY_P0_Z(x)); }
  if(ARRAY_TYPE(x)=='r'){ y=array_allocate('c',ARRAY_NDIM(x),ARRAY_DIM_P(x));  cvec_copy_r(ARRAY_SIZE(y),ARRAY_P0_C(y),ARRAY_P0_R(x)); }
  if(ARRAY_TYPE(x)=='c'){ y=array_allocate('c',ARRAY_NDIM(x),ARRAY_DIM_P(x));  cvec_copy  (ARRAY_SIZE(y),ARRAY_P0_C(y),ARRAY_P0_C(x)); }
  if(ARRAY_TYPE(x)=='R'){ y=array_allocate('C',ARRAY_NDIM(x),ARRAY_DIM_P(x)); icvec_copy_r(ARRAY_SIZE(y),ARRAY_P0_C(y),ARRAY_P1_C(y),ARRAY_P0_R(x),ARRAY_P1_R(x)); }
  if(ARRAY_TYPE(x)=='C'){ y=array_allocate('C',ARRAY_NDIM(x),ARRAY_DIM_P(x)); icvec_copy  (ARRAY_SIZE(y),ARRAY_P0_C(y),ARRAY_P1_C(y),ARRAY_P0_C(x),ARRAY_P1_C(x)); }
  return y;
}
  
/**
 @breif array型の型変換 y=double(x)
 */ 
array *array_get_double(array *x)
{
  int flag=0;
  array *y=NULL;
  if(x==NULL){ return y; }
  if(ARRAY_TYPE(x)=='s'){ flag=svec_has(ARRAY_SIZE(x),ARRAY_P0_S(x),"iIjJ"); }
  if(ARRAY_TYPE(x)=='s' && !flag){ y=array_allocate('d',ARRAY_NDIM(x),ARRAY_DIM_P(x));  dvec_set_s (ARRAY_SIZE(y),ARRAY_P0_D(y),ARRAY_P0_S(x)); }
  if(ARRAY_TYPE(x)=='s' &&  flag){ y=array_allocate('z',ARRAY_NDIM(x),ARRAY_DIM_P(x));  zvec_set_s (ARRAY_SIZE(y),ARRAY_P0_Z(y),ARRAY_P0_S(x)); }
  if(ARRAY_TYPE(x)=='i')         { y=array_allocate('d',ARRAY_NDIM(x),ARRAY_DIM_P(x));  dvec_set_si(ARRAY_SIZE(y),ARRAY_P0_D(y),ARRAY_P0_I(x)); }
  if(ARRAY_TYPE(x)=='d')         { y=array_allocate('d',ARRAY_NDIM(x),ARRAY_DIM_P(x));  dvec_copy  (ARRAY_SIZE(y),ARRAY_P0_D(y),ARRAY_P0_D(x)); }
  if(ARRAY_TYPE(x)=='z')         { y=array_allocate('z',ARRAY_NDIM(x),ARRAY_DIM_P(x));  zvec_copy  (ARRAY_SIZE(y),ARRAY_P0_Z(y),ARRAY_P0_Z(x)); }
  if(ARRAY_TYPE(x)=='r')         { y=array_allocate('d',ARRAY_NDIM(x),ARRAY_DIM_P(x));  rvec_get_d (ARRAY_SIZE(y),ARRAY_P0_D(y),ARRAY_P0_R(x)); }
  if(ARRAY_TYPE(x)=='c')         { y=array_allocate('z',ARRAY_NDIM(x),ARRAY_DIM_P(x));  cvec_get_z (ARRAY_SIZE(y),ARRAY_P0_Z(y),ARRAY_P0_C(x)); }
  if(ARRAY_TYPE(x)=='R')         { y=array_allocate('d',ARRAY_NDIM(x),ARRAY_DIM_P(x)); irvec_get_d (ARRAY_SIZE(y),ARRAY_P0_D(y),ARRAY_P0_R(x),ARRAY_P1_R(x)); }
  if(ARRAY_TYPE(x)=='C')         { y=array_allocate('z',ARRAY_NDIM(x),ARRAY_DIM_P(x)); icvec_get_z (ARRAY_SIZE(y),ARRAY_P0_Z(y),ARRAY_P0_C(x),ARRAY_P1_C(x)); }
  return y;
}

/**
 @breif array型の型変換 y=multi(x)
 */ 
array *array_get_multi(array *x)
{
  int flag=0;
  array *y=NULL;
  if(x==NULL){ return y; }
  if(ARRAY_TYPE(x)=='s'){ flag=svec_has(ARRAY_SIZE(x),ARRAY_P0_S(x),"iIjJ"); }
  if(ARRAY_TYPE(x)=='s' &&  flag){ y=array_allocate('c',ARRAY_NDIM(x),ARRAY_DIM_P(x)); cvec_set_s (ARRAY_SIZE(y),ARRAY_P0_C(y),ARRAY_P0_S(x)); }
  if(ARRAY_TYPE(x)=='s' && !flag){ y=array_allocate('r',ARRAY_NDIM(x),ARRAY_DIM_P(x)); rvec_set_s (ARRAY_SIZE(y),ARRAY_P0_R(y),ARRAY_P0_S(x)); }
  if(ARRAY_TYPE(x)=='i')         { y=array_allocate('r',ARRAY_NDIM(x),ARRAY_DIM_P(x)); rvec_set_si(ARRAY_SIZE(y),ARRAY_P0_R(y),ARRAY_P0_I(x)); }
  if(ARRAY_TYPE(x)=='d')         { y=array_allocate('r',ARRAY_NDIM(x),ARRAY_DIM_P(x)); rvec_set_d (ARRAY_SIZE(y),ARRAY_P0_R(y),ARRAY_P0_D(x)); }
  if(ARRAY_TYPE(x)=='z')         { y=array_allocate('c',ARRAY_NDIM(x),ARRAY_DIM_P(x)); cvec_set_z (ARRAY_SIZE(y),ARRAY_P0_C(y),ARRAY_P0_Z(x)); }
  if(ARRAY_TYPE(x)=='r')         { y=array_allocate('r',ARRAY_NDIM(x),ARRAY_DIM_P(x)); rvec_set   (ARRAY_SIZE(y),ARRAY_P0_R(y),ARRAY_P0_R(x)); }
  if(ARRAY_TYPE(x)=='c')         { y=array_allocate('c',ARRAY_NDIM(x),ARRAY_DIM_P(x)); cvec_set   (ARRAY_SIZE(y),ARRAY_P0_C(y),ARRAY_P0_C(x)); }
  if(ARRAY_TYPE(x)=='R')         { y=array_allocate('r',ARRAY_NDIM(x),ARRAY_DIM_P(x)); rvec_set_ir(ARRAY_SIZE(y),ARRAY_P0_R(y),ARRAY_P0_R(x),ARRAY_P1_R(x)); }
  if(ARRAY_TYPE(x)=='C')         { y=array_allocate('c',ARRAY_NDIM(x),ARRAY_DIM_P(x)); cvec_set_ic(ARRAY_SIZE(y),ARRAY_P0_C(y),ARRAY_P0_C(x),ARRAY_P1_C(x)); }
  return y;
}

/**
 @breif array型の型変換 y=imulti(x)
 */ 
array *array_get_imulti(array *x)
{
  int flag=0;
  array *y=NULL;
  if(x==NULL){ return y; }
  if(ARRAY_TYPE(x)=='s'){ flag=svec_has(ARRAY_SIZE(x),ARRAY_P0_S(x),"iIjJ"); }
  if(ARRAY_TYPE(x)=='s' &&  flag){ y=array_allocate('C',ARRAY_NDIM(x),ARRAY_DIM_P(x)); icvec_set_s (ARRAY_SIZE(y),ARRAY_P0_C(y),ARRAY_P1_C(y),ARRAY_P0_S(x)); }
  if(ARRAY_TYPE(x)=='s' && !flag){ y=array_allocate('R',ARRAY_NDIM(x),ARRAY_DIM_P(x)); irvec_set_s (ARRAY_SIZE(y),ARRAY_P0_R(y),ARRAY_P1_R(y),ARRAY_P0_S(x)); }
  if(ARRAY_TYPE(x)=='i')         { y=array_allocate('R',ARRAY_NDIM(x),ARRAY_DIM_P(x)); irvec_set_si(ARRAY_SIZE(y),ARRAY_P0_R(y),ARRAY_P1_R(y),ARRAY_P0_I(x)); }
  if(ARRAY_TYPE(x)=='d')         { y=array_allocate('R',ARRAY_NDIM(x),ARRAY_DIM_P(x)); irvec_set_d (ARRAY_SIZE(y),ARRAY_P0_R(y),ARRAY_P1_R(y),ARRAY_P0_D(x),ARRAY_P0_D(x)); }
  if(ARRAY_TYPE(x)=='z')         { y=array_allocate('C',ARRAY_NDIM(x),ARRAY_DIM_P(x)); icvec_set_z (ARRAY_SIZE(y),ARRAY_P0_C(y),ARRAY_P1_C(y),ARRAY_P0_Z(x),ARRAY_P0_Z(x)); }
  if(ARRAY_TYPE(x)=='r')         { y=array_allocate('R',ARRAY_NDIM(x),ARRAY_DIM_P(x)); irvec_set   (ARRAY_SIZE(y),ARRAY_P0_R(y),ARRAY_P1_R(y),ARRAY_P0_R(x),ARRAY_P0_R(x)); }
  if(ARRAY_TYPE(x)=='c')         { y=array_allocate('C',ARRAY_NDIM(x),ARRAY_DIM_P(x)); icvec_set   (ARRAY_SIZE(y),ARRAY_P0_C(y),ARRAY_P1_C(y),ARRAY_P0_C(x),ARRAY_P0_C(x)); }
  if(ARRAY_TYPE(x)=='R')         { y=array_allocate('R',ARRAY_NDIM(x),ARRAY_DIM_P(x)); irvec_set   (ARRAY_SIZE(y),ARRAY_P0_R(y),ARRAY_P1_R(y),ARRAY_P0_R(x),ARRAY_P1_R(x)); }
  if(ARRAY_TYPE(x)=='C')         { y=array_allocate('C',ARRAY_NDIM(x),ARRAY_DIM_P(x)); icvec_set   (ARRAY_SIZE(y),ARRAY_P0_C(y),ARRAY_P1_C(y),ARRAY_P0_C(x),ARRAY_P1_C(x)); }
  return y;
}

/**
 @breif array型の型変換 y=imulti2(x0,x1)
 */ 
array *array_get_imulti2(array *x0, array *x1)
{
  array *y=NULL;
  if(x0==NULL || x1==NULL){ return y; }
  if(!array_same_dim_check(x0,x1)){ return y; }
  if(ARRAY_TYPE(x0)=='d' && ARRAY_TYPE(x1)=='d'){ y=array_allocate('R',ARRAY_NDIM(x0),ARRAY_DIM_P(x0)); irvec_set_d (ARRAY_SIZE(y),ARRAY_P0_R(y),ARRAY_P1_R(y),ARRAY_P0_D(x0),ARRAY_P0_D(x1)); }
  if(ARRAY_TYPE(x0)=='d' && ARRAY_TYPE(x1)=='z'){ y=array_allocate('C',ARRAY_NDIM(x0),ARRAY_DIM_P(x0)); icvec_set_dz(ARRAY_SIZE(y),ARRAY_P0_C(y),ARRAY_P1_C(y),ARRAY_P0_D(x0),ARRAY_P0_Z(x1)); }
  if(ARRAY_TYPE(x0)=='z' && ARRAY_TYPE(x1)=='d'){ y=array_allocate('C',ARRAY_NDIM(x0),ARRAY_DIM_P(x0)); icvec_set_zd(ARRAY_SIZE(y),ARRAY_P0_C(y),ARRAY_P1_C(y),ARRAY_P0_Z(x0),ARRAY_P0_D(x1)); }
  if(ARRAY_TYPE(x0)=='z' && ARRAY_TYPE(x1)=='z'){ y=array_allocate('C',ARRAY_NDIM(x0),ARRAY_DIM_P(x0)); icvec_set_z (ARRAY_SIZE(y),ARRAY_P0_C(y),ARRAY_P1_C(y),ARRAY_P0_Z(x0),ARRAY_P0_Z(x1)); }
  if(ARRAY_TYPE(x0)=='r' && ARRAY_TYPE(x1)=='r'){ y=array_allocate('R',ARRAY_NDIM(x0),ARRAY_DIM_P(x0)); irvec_set   (ARRAY_SIZE(y),ARRAY_P0_R(y),ARRAY_P1_R(y),ARRAY_P0_R(x0),ARRAY_P0_R(x1)); }
  if(ARRAY_TYPE(x0)=='r' && ARRAY_TYPE(x1)=='c'){ y=array_allocate('C',ARRAY_NDIM(x0),ARRAY_DIM_P(x0)); icvec_set_rc(ARRAY_SIZE(y),ARRAY_P0_C(y),ARRAY_P1_C(y),ARRAY_P0_R(x0),ARRAY_P0_C(x1)); }
  if(ARRAY_TYPE(x0)=='c' && ARRAY_TYPE(x1)=='r'){ y=array_allocate('C',ARRAY_NDIM(x0),ARRAY_DIM_P(x0)); icvec_set_cr(ARRAY_SIZE(y),ARRAY_P0_C(y),ARRAY_P1_C(y),ARRAY_P0_C(x0),ARRAY_P0_R(x1)); }
  if(ARRAY_TYPE(x0)=='c' && ARRAY_TYPE(x1)=='c'){ y=array_allocate('C',ARRAY_NDIM(x0),ARRAY_DIM_P(x0)); icvec_set   (ARRAY_SIZE(y),ARRAY_P0_C(y),ARRAY_P1_C(y),ARRAY_P0_C(x0),ARRAY_P0_C(x1)); }
  return y;
}

/**
 @breif array型の型変換 y=complex2(xr,xi)
 */ 
array *array_get_complex2(array *xr, array *xi)
{
  array *y=NULL;
  if(xr==NULL || xi==NULL){ return y; }
  if(!array_same_dim_check(xr,xi)){ return y; }
  if(ARRAY_TYPE(xr)=='d' && ARRAY_TYPE(xi)=='d'){ y=array_allocate('z',ARRAY_NDIM(xr),ARRAY_DIM_P(xr));  zvec_copy_Z(ARRAY_SIZE(y),ARRAY_P0_Z(y),ARRAY_P0_D(xr),ARRAY_P0_D(xi)); }
  if(ARRAY_TYPE(xr)=='r' && ARRAY_TYPE(xi)=='r'){ y=array_allocate('c',ARRAY_NDIM(xr),ARRAY_DIM_P(xr));  cvec_copy_C(ARRAY_SIZE(y),ARRAY_P0_C(y),ARRAY_P0_R(xr),ARRAY_P0_R(xi)); }
  if(ARRAY_TYPE(xr)=='R' && ARRAY_TYPE(xi)=='R'){ y=array_allocate('C',ARRAY_NDIM(xr),ARRAY_DIM_P(xr)); icvec_copy_C(ARRAY_SIZE(y),ARRAY_P0_C(y),ARRAY_P1_C(y),ARRAY_P0_R(xr),ARRAY_P0_R(xi),ARRAY_P1_R(xr),ARRAY_P1_R(xi)); }
  return y;
}

/** @} */

///////////////////////////////////////////////

/** @name array型の１入力演算に関する関数 */
/** @{ */


/**
 @breif y=conj(x)
 */ 
array *array_conj(array *x)
{
  array *y=NULL;
  if(x==NULL){ return y; }
  if(ARRAY_TYPE(x)=='d'){ y=array_allocate('d',ARRAY_NDIM(x),ARRAY_DIM_P(x));  dvec_copy(ARRAY_SIZE(y),ARRAY_P0_D(y),              ARRAY_P0_D(x)); }
  if(ARRAY_TYPE(x)=='z'){ y=array_allocate('z',ARRAY_NDIM(x),ARRAY_DIM_P(x));  zvec_conj(ARRAY_SIZE(y),ARRAY_P0_Z(y),              ARRAY_P0_Z(x)); }
  if(ARRAY_TYPE(x)=='r'){ y=array_allocate('r',ARRAY_NDIM(x),ARRAY_DIM_P(x));  rvec_copy(ARRAY_SIZE(y),ARRAY_P0_R(y),              ARRAY_P0_R(x)); }
  if(ARRAY_TYPE(x)=='c'){ y=array_allocate('c',ARRAY_NDIM(x),ARRAY_DIM_P(x));  cvec_conj(ARRAY_SIZE(y),ARRAY_P0_C(y),              ARRAY_P0_C(x)); }
  if(ARRAY_TYPE(x)=='R'){ y=array_allocate('R',ARRAY_NDIM(x),ARRAY_DIM_P(x)); irvec_copy(ARRAY_SIZE(y),ARRAY_P0_R(y),ARRAY_P1_R(y),ARRAY_P0_R(x),ARRAY_P1_R(x)); }
  if(ARRAY_TYPE(x)=='C'){ y=array_allocate('C',ARRAY_NDIM(x),ARRAY_DIM_P(x)); icvec_conj(ARRAY_SIZE(y),ARRAY_P0_C(y),ARRAY_P1_C(y),ARRAY_P0_C(x),ARRAY_P1_C(x)); }
  return y;
}

  
/** @} */

///////////////////////////////////////////////

/** @name array型の２入力演算に関する関数 */
/** @{ */

/**
 @breif z=x+y
 */ 
array *array_add(array *x, array *y)
{
  int k,i;
  array *z=NULL;
  if(x==NULL || y==NULL){ return z; }
  if(z==NULL && array_same_dim_check(x,y)){
    if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='d'){ z=array_allocate('d',ARRAY_NDIM(x),ARRAY_DIM_P(x));  dvec_add_dvec(ARRAY_SIZE(z),ARRAY_P0_D(z),ARRAY_P0_D(x),ARRAY_P0_D(y)); }
    if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='z'){ z=array_allocate('z',ARRAY_NDIM(x),ARRAY_DIM_P(x));  dvec_add_zvec(ARRAY_SIZE(z),ARRAY_P0_Z(z),ARRAY_P0_D(x),ARRAY_P0_Z(y)); }
    if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='d'){ z=array_allocate('z',ARRAY_NDIM(x),ARRAY_DIM_P(x));  zvec_add_dvec(ARRAY_SIZE(z),ARRAY_P0_Z(z),ARRAY_P0_Z(x),ARRAY_P0_D(y)); }
    if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='z'){ z=array_allocate('z',ARRAY_NDIM(x),ARRAY_DIM_P(x));  zvec_add_zvec(ARRAY_SIZE(z),ARRAY_P0_Z(z),ARRAY_P0_Z(x),ARRAY_P0_Z(y)); }
    if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='r'){ z=array_allocate('r',ARRAY_NDIM(x),ARRAY_DIM_P(x));  dvec_add_rvec(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P0_D(x),ARRAY_P0_R(y)); }
    if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='r'){ z=array_allocate('c',ARRAY_NDIM(x),ARRAY_DIM_P(x));  zvec_add_rvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_Z(x),ARRAY_P0_R(y)); }
    if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='d'){ z=array_allocate('r',ARRAY_NDIM(x),ARRAY_DIM_P(x));  rvec_add_dvec(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P0_R(x),ARRAY_P0_D(y)); }
    if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='z'){ z=array_allocate('c',ARRAY_NDIM(x),ARRAY_DIM_P(x));  rvec_add_zvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_R(x),ARRAY_P0_Z(y)); }
    if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='r'){ z=array_allocate('r',ARRAY_NDIM(x),ARRAY_DIM_P(x));  rvec_add_rvec(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P0_R(x),ARRAY_P0_R(y)); }
    if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='c'){ z=array_allocate('c',ARRAY_NDIM(x),ARRAY_DIM_P(x));  dvec_add_cvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_D(x),ARRAY_P0_C(y)); }
    if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='c'){ z=array_allocate('c',ARRAY_NDIM(x),ARRAY_DIM_P(x));  zvec_add_cvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_Z(x),ARRAY_P0_C(y)); }
    if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='c'){ z=array_allocate('c',ARRAY_NDIM(x),ARRAY_DIM_P(x));  rvec_add_cvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_R(x),ARRAY_P0_C(y)); }
    if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='d'){ z=array_allocate('c',ARRAY_NDIM(x),ARRAY_DIM_P(x));  cvec_add_dvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_C(x),ARRAY_P0_D(y)); }
    if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='z'){ z=array_allocate('c',ARRAY_NDIM(x),ARRAY_DIM_P(x));  cvec_add_zvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_C(x),ARRAY_P0_Z(y)); }
    if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='r'){ z=array_allocate('c',ARRAY_NDIM(x),ARRAY_DIM_P(x));  cvec_add_rvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_C(x),ARRAY_P0_R(y)); }
    if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='c'){ z=array_allocate('c',ARRAY_NDIM(x),ARRAY_DIM_P(x));  cvec_add_cvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_C(x),ARRAY_P0_C(y)); }
    if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='R'){ z=array_allocate('R',ARRAY_NDIM(x),ARRAY_DIM_P(x)); idvec_add_rvec(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P1_R(z),ARRAY_P0_D(x),ARRAY_P0_D(x),ARRAY_P0_R(y),ARRAY_P1_R(y)); }
    if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='R'){ z=array_allocate('C',ARRAY_NDIM(x),ARRAY_DIM_P(x)); izvec_add_rvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_Z(x),ARRAY_P0_Z(x),ARRAY_P0_R(y),ARRAY_P1_R(y)); }
    if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='R'){ z=array_allocate('R',ARRAY_NDIM(x),ARRAY_DIM_P(x)); irvec_add_rvec(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P1_R(z),ARRAY_P0_R(x),ARRAY_P0_R(x),ARRAY_P0_R(y),ARRAY_P1_R(y)); }
    if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='R'){ z=array_allocate('C',ARRAY_NDIM(x),ARRAY_DIM_P(x)); icvec_add_rvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(x),ARRAY_P0_C(x),ARRAY_P0_R(y),ARRAY_P1_R(y)); }
    if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='d'){ z=array_allocate('R',ARRAY_NDIM(x),ARRAY_DIM_P(x)); irvec_add_dvec(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P1_R(z),ARRAY_P0_R(x),ARRAY_P1_R(x),ARRAY_P0_D(y),ARRAY_P0_D(y)); }
    if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='z'){ z=array_allocate('C',ARRAY_NDIM(x),ARRAY_DIM_P(x)); irvec_add_zvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_R(x),ARRAY_P1_R(x),ARRAY_P0_Z(y),ARRAY_P0_Z(y)); }
    if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='r'){ z=array_allocate('R',ARRAY_NDIM(x),ARRAY_DIM_P(x)); irvec_add_rvec(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P1_R(z),ARRAY_P0_R(x),ARRAY_P1_R(x),ARRAY_P0_R(y),ARRAY_P0_R(y)); }
    if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='c'){ z=array_allocate('C',ARRAY_NDIM(x),ARRAY_DIM_P(x)); irvec_add_cvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_R(x),ARRAY_P1_R(x),ARRAY_P0_C(y),ARRAY_P0_C(y)); }
    if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='R'){ z=array_allocate('R',ARRAY_NDIM(x),ARRAY_DIM_P(x)); irvec_add_rvec(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P1_R(z),ARRAY_P0_R(x),ARRAY_P1_R(x),ARRAY_P0_R(y),ARRAY_P1_R(y)); }
    if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='C'){ z=array_allocate('C',ARRAY_NDIM(x),ARRAY_DIM_P(x)); idvec_add_cvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_D(x),ARRAY_P0_D(x),ARRAY_P0_C(y),ARRAY_P1_C(y)); }
    if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='C'){ z=array_allocate('C',ARRAY_NDIM(x),ARRAY_DIM_P(x)); izvec_add_cvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_Z(x),ARRAY_P0_Z(x),ARRAY_P0_C(y),ARRAY_P1_C(y)); }
    if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='C'){ z=array_allocate('C',ARRAY_NDIM(x),ARRAY_DIM_P(x)); irvec_add_cvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_R(x),ARRAY_P0_R(x),ARRAY_P0_C(y),ARRAY_P1_C(y)); }
    if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='C'){ z=array_allocate('C',ARRAY_NDIM(x),ARRAY_DIM_P(x)); icvec_add_cvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(x),ARRAY_P0_C(x),ARRAY_P0_C(y),ARRAY_P1_C(y)); }
    if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='C'){ z=array_allocate('C',ARRAY_NDIM(x),ARRAY_DIM_P(x)); irvec_add_cvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_R(x),ARRAY_P1_R(x),ARRAY_P0_C(y),ARRAY_P1_C(y)); }
    if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='d'){ z=array_allocate('C',ARRAY_NDIM(x),ARRAY_DIM_P(x)); icvec_add_dvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(x),ARRAY_P1_C(x),ARRAY_P0_D(y),ARRAY_P0_D(y)); }
    if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='z'){ z=array_allocate('C',ARRAY_NDIM(x),ARRAY_DIM_P(x)); icvec_add_zvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(x),ARRAY_P1_C(x),ARRAY_P0_Z(y),ARRAY_P0_Z(y)); }
    if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='r'){ z=array_allocate('C',ARRAY_NDIM(x),ARRAY_DIM_P(x)); icvec_add_rvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(x),ARRAY_P1_C(x),ARRAY_P0_R(y),ARRAY_P0_R(y)); }
    if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='c'){ z=array_allocate('C',ARRAY_NDIM(x),ARRAY_DIM_P(x)); icvec_add_cvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(x),ARRAY_P1_C(x),ARRAY_P0_C(y),ARRAY_P0_C(y)); }
    if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='R'){ z=array_allocate('C',ARRAY_NDIM(x),ARRAY_DIM_P(x)); icvec_add_rvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(x),ARRAY_P1_C(x),ARRAY_P0_R(y),ARRAY_P1_R(y)); }
    if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='C'){ z=array_allocate('C',ARRAY_NDIM(x),ARRAY_DIM_P(x)); icvec_add_cvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(x),ARRAY_P1_C(x),ARRAY_P0_C(y),ARRAY_P1_C(y)); }
  }
  if(z==NULL && array_is_scalar(y)){
    if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='d'){ z=array_allocate('d',ARRAY_NDIM(x),ARRAY_DIM_P(x));  dvec_add_dscalar(ARRAY_SIZE(z),ARRAY_P0_D(z),ARRAY_P0_D(x),ARRAY_P0_D(y)[0]); }
    if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='z'){ z=array_allocate('z',ARRAY_NDIM(x),ARRAY_DIM_P(x));  dvec_add_zscalar(ARRAY_SIZE(z),ARRAY_P0_Z(z),ARRAY_P0_D(x),ARRAY_P0_Z(y)[0]); }
    if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='d'){ z=array_allocate('z',ARRAY_NDIM(x),ARRAY_DIM_P(x));  zvec_add_dscalar(ARRAY_SIZE(z),ARRAY_P0_Z(z),ARRAY_P0_Z(x),ARRAY_P0_D(y)[0]); }
    if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='z'){ z=array_allocate('z',ARRAY_NDIM(x),ARRAY_DIM_P(x));  zvec_add_zscalar(ARRAY_SIZE(z),ARRAY_P0_Z(z),ARRAY_P0_Z(x),ARRAY_P0_Z(y)[0]); }
    if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='r'){ z=array_allocate('r',ARRAY_NDIM(x),ARRAY_DIM_P(x));  dvec_add_rscalar(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P0_D(x),ARRAY_P0_R(y)[0]); }
    if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='r'){ z=array_allocate('c',ARRAY_NDIM(x),ARRAY_DIM_P(x));  zvec_add_rscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_Z(x),ARRAY_P0_R(y)[0]); }
    if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='d'){ z=array_allocate('r',ARRAY_NDIM(x),ARRAY_DIM_P(x));  rvec_add_dscalar(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P0_R(x),ARRAY_P0_D(y)[0]); }
    if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='z'){ z=array_allocate('c',ARRAY_NDIM(x),ARRAY_DIM_P(x));  rvec_add_zscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_R(x),ARRAY_P0_Z(y)[0]); }
    if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='r'){ z=array_allocate('r',ARRAY_NDIM(x),ARRAY_DIM_P(x));  rvec_add_rscalar(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P0_R(x),ARRAY_P0_R(y)[0]); }
    if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='c'){ z=array_allocate('c',ARRAY_NDIM(x),ARRAY_DIM_P(x));  dvec_add_cscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_D(x),ARRAY_P0_C(y)[0]); }
    if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='c'){ z=array_allocate('c',ARRAY_NDIM(x),ARRAY_DIM_P(x));  zvec_add_cscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_Z(x),ARRAY_P0_C(y)[0]); }
    if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='c'){ z=array_allocate('c',ARRAY_NDIM(x),ARRAY_DIM_P(x));  rvec_add_cscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_R(x),ARRAY_P0_C(y)[0]); }
    if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='d'){ z=array_allocate('c',ARRAY_NDIM(x),ARRAY_DIM_P(x));  cvec_add_dscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_C(x),ARRAY_P0_D(y)[0]); }
    if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='z'){ z=array_allocate('c',ARRAY_NDIM(x),ARRAY_DIM_P(x));  cvec_add_zscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_C(x),ARRAY_P0_Z(y)[0]); }
    if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='r'){ z=array_allocate('c',ARRAY_NDIM(x),ARRAY_DIM_P(x));  cvec_add_rscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_C(x),ARRAY_P0_R(y)[0]); }
    if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='c'){ z=array_allocate('c',ARRAY_NDIM(x),ARRAY_DIM_P(x));  cvec_add_cscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_C(x),ARRAY_P0_C(y)[0]); }
    if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='R'){ z=array_allocate('R',ARRAY_NDIM(x),ARRAY_DIM_P(x)); idvec_add_rscalar(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P1_R(z),ARRAY_P0_D(x),ARRAY_P0_D(x),ARRAY_P0_R(y)[0],ARRAY_P1_R(y)[0]); }
    if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='R'){ z=array_allocate('C',ARRAY_NDIM(x),ARRAY_DIM_P(x)); izvec_add_rscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_Z(x),ARRAY_P0_Z(x),ARRAY_P0_R(y)[0],ARRAY_P1_R(y)[0]); }
    if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='R'){ z=array_allocate('R',ARRAY_NDIM(x),ARRAY_DIM_P(x)); irvec_add_rscalar(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P1_R(z),ARRAY_P0_R(x),ARRAY_P0_R(x),ARRAY_P0_R(y)[0],ARRAY_P1_R(y)[0]); }
    if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='R'){ z=array_allocate('C',ARRAY_NDIM(x),ARRAY_DIM_P(x)); icvec_add_rscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(x),ARRAY_P0_C(x),ARRAY_P0_R(y)[0],ARRAY_P1_R(y)[0]); }
    if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='d'){ z=array_allocate('R',ARRAY_NDIM(x),ARRAY_DIM_P(x)); irvec_add_dscalar(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P1_R(z),ARRAY_P0_R(x),ARRAY_P1_R(x),ARRAY_P0_D(y)[0],ARRAY_P0_D(y)[0]); }
    if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='z'){ z=array_allocate('C',ARRAY_NDIM(x),ARRAY_DIM_P(x)); irvec_add_zscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_R(x),ARRAY_P1_R(x),ARRAY_P0_Z(y)[0],ARRAY_P0_Z(y)[0]); }
    if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='r'){ z=array_allocate('R',ARRAY_NDIM(x),ARRAY_DIM_P(x)); irvec_add_rscalar(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P1_R(z),ARRAY_P0_R(x),ARRAY_P1_R(x),ARRAY_P0_R(y)[0],ARRAY_P0_R(y)[0]); }
    if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='c'){ z=array_allocate('C',ARRAY_NDIM(x),ARRAY_DIM_P(x)); irvec_add_cscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_R(x),ARRAY_P1_R(x),ARRAY_P0_C(y)[0],ARRAY_P0_C(y)[0]); }
    if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='R'){ z=array_allocate('R',ARRAY_NDIM(x),ARRAY_DIM_P(x)); irvec_add_rscalar(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P1_R(z),ARRAY_P0_R(x),ARRAY_P1_R(x),ARRAY_P0_R(y)[0],ARRAY_P1_R(y)[0]); }
    if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='C'){ z=array_allocate('C',ARRAY_NDIM(x),ARRAY_DIM_P(x)); idvec_add_cscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_D(x),ARRAY_P0_D(x),ARRAY_P0_C(y)[0],ARRAY_P1_C(y)[0]); }
    if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='C'){ z=array_allocate('C',ARRAY_NDIM(x),ARRAY_DIM_P(x)); izvec_add_cscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_Z(x),ARRAY_P0_Z(x),ARRAY_P0_C(y)[0],ARRAY_P1_C(y)[0]); }
    if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='C'){ z=array_allocate('C',ARRAY_NDIM(x),ARRAY_DIM_P(x)); irvec_add_cscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_R(x),ARRAY_P0_R(x),ARRAY_P0_C(y)[0],ARRAY_P1_C(y)[0]); }
    if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='C'){ z=array_allocate('C',ARRAY_NDIM(x),ARRAY_DIM_P(x)); icvec_add_cscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(x),ARRAY_P0_C(x),ARRAY_P0_C(y)[0],ARRAY_P1_C(y)[0]); }
    if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='C'){ z=array_allocate('C',ARRAY_NDIM(x),ARRAY_DIM_P(x)); irvec_add_cscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_R(x),ARRAY_P1_R(x),ARRAY_P0_C(y)[0],ARRAY_P1_C(y)[0]); }
    if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='d'){ z=array_allocate('C',ARRAY_NDIM(x),ARRAY_DIM_P(x)); icvec_add_dscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(x),ARRAY_P1_C(x),ARRAY_P0_D(y)[0],ARRAY_P0_D(y)[0]); }
    if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='z'){ z=array_allocate('C',ARRAY_NDIM(x),ARRAY_DIM_P(x)); icvec_add_zscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(x),ARRAY_P1_C(x),ARRAY_P0_Z(y)[0],ARRAY_P0_Z(y)[0]); }
    if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='r'){ z=array_allocate('C',ARRAY_NDIM(x),ARRAY_DIM_P(x)); icvec_add_rscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(x),ARRAY_P1_C(x),ARRAY_P0_R(y)[0],ARRAY_P0_R(y)[0]); }
    if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='c'){ z=array_allocate('C',ARRAY_NDIM(x),ARRAY_DIM_P(x)); icvec_add_cscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(x),ARRAY_P1_C(x),ARRAY_P0_C(y)[0],ARRAY_P0_C(y)[0]); }
    if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='R'){ z=array_allocate('C',ARRAY_NDIM(x),ARRAY_DIM_P(x)); icvec_add_rscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(x),ARRAY_P1_C(x),ARRAY_P0_R(y)[0],ARRAY_P1_R(y)[0]); }
    if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='C'){ z=array_allocate('C',ARRAY_NDIM(x),ARRAY_DIM_P(x)); icvec_add_cscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(x),ARRAY_P1_C(x),ARRAY_P0_C(y)[0],ARRAY_P1_C(y)[0]); }
  }
  if(z==NULL && array_is_scalar(x)){
    if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='d'){ z=array_allocate('d',ARRAY_NDIM(y),ARRAY_DIM_P(y));  dvec_add_dscalar(ARRAY_SIZE(z),ARRAY_P0_D(z),ARRAY_P0_D(y),ARRAY_P0_D(x)[0]); }
    if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='z'){ z=array_allocate('z',ARRAY_NDIM(y),ARRAY_DIM_P(y));  zvec_add_dscalar(ARRAY_SIZE(z),ARRAY_P0_Z(z),ARRAY_P0_Z(y),ARRAY_P0_D(x)[0]); }
    if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='d'){ z=array_allocate('z',ARRAY_NDIM(y),ARRAY_DIM_P(y));  dvec_add_zscalar(ARRAY_SIZE(z),ARRAY_P0_Z(z),ARRAY_P0_D(y),ARRAY_P0_Z(x)[0]); }
    if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='z'){ z=array_allocate('z',ARRAY_NDIM(y),ARRAY_DIM_P(y));  zvec_add_zscalar(ARRAY_SIZE(z),ARRAY_P0_Z(z),ARRAY_P0_Z(y),ARRAY_P0_Z(x)[0]); }
    if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='r'){ z=array_allocate('r',ARRAY_NDIM(y),ARRAY_DIM_P(y));  rvec_add_dscalar(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P0_R(y),ARRAY_P0_D(x)[0]); }
    if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='r'){ z=array_allocate('c',ARRAY_NDIM(y),ARRAY_DIM_P(y));  rvec_add_zscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_R(y),ARRAY_P0_Z(x)[0]); }
    if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='d'){ z=array_allocate('r',ARRAY_NDIM(y),ARRAY_DIM_P(y));  dvec_add_rscalar(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P0_D(y),ARRAY_P0_R(x)[0]); }
    if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='z'){ z=array_allocate('c',ARRAY_NDIM(y),ARRAY_DIM_P(y));  zvec_add_rscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_Z(y),ARRAY_P0_R(x)[0]); }
    if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='r'){ z=array_allocate('r',ARRAY_NDIM(y),ARRAY_DIM_P(y));  rvec_add_rscalar(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P0_R(y),ARRAY_P0_R(x)[0]); }
    if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='c'){ z=array_allocate('c',ARRAY_NDIM(y),ARRAY_DIM_P(y));  cvec_add_dscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_C(y),ARRAY_P0_D(x)[0]); }
    if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='c'){ z=array_allocate('c',ARRAY_NDIM(y),ARRAY_DIM_P(y));  cvec_add_zscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_C(y),ARRAY_P0_Z(x)[0]); }
    if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='c'){ z=array_allocate('c',ARRAY_NDIM(y),ARRAY_DIM_P(y));  cvec_add_rscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_C(y),ARRAY_P0_R(x)[0]); }
    if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='d'){ z=array_allocate('c',ARRAY_NDIM(y),ARRAY_DIM_P(y));  dvec_add_cscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_D(y),ARRAY_P0_C(x)[0]); }
    if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='z'){ z=array_allocate('c',ARRAY_NDIM(y),ARRAY_DIM_P(y));  zvec_add_cscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_Z(y),ARRAY_P0_C(x)[0]); }
    if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='r'){ z=array_allocate('c',ARRAY_NDIM(y),ARRAY_DIM_P(y));  rvec_add_cscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_R(y),ARRAY_P0_C(x)[0]); }
    if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='c'){ z=array_allocate('c',ARRAY_NDIM(y),ARRAY_DIM_P(y));  cvec_add_cscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_C(y),ARRAY_P0_C(x)[0]); }
    if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='R'){ z=array_allocate('R',ARRAY_NDIM(y),ARRAY_DIM_P(y)); irvec_add_dscalar(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P1_R(z),ARRAY_P0_R(y),ARRAY_P1_R(y),ARRAY_P0_D(x)[0],ARRAY_P0_D(x)[0]); }
    if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='R'){ z=array_allocate('C',ARRAY_NDIM(y),ARRAY_DIM_P(y)); irvec_add_zscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_R(y),ARRAY_P1_R(y),ARRAY_P0_Z(x)[0],ARRAY_P0_Z(x)[0]); }
    if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='R'){ z=array_allocate('R',ARRAY_NDIM(y),ARRAY_DIM_P(y)); irvec_add_rscalar(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P1_R(z),ARRAY_P0_R(y),ARRAY_P1_R(y),ARRAY_P0_R(x)[0],ARRAY_P0_R(x)[0]); }
    if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='R'){ z=array_allocate('C',ARRAY_NDIM(y),ARRAY_DIM_P(y)); irvec_add_cscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_R(y),ARRAY_P1_R(y),ARRAY_P0_C(x)[0],ARRAY_P0_C(x)[0]); }
    if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='d'){ z=array_allocate('R',ARRAY_NDIM(y),ARRAY_DIM_P(y)); idvec_add_rscalar(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P1_R(z),ARRAY_P0_D(y),ARRAY_P0_D(y),ARRAY_P0_R(x)[0],ARRAY_P1_R(x)[0]); }
    if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='z'){ z=array_allocate('C',ARRAY_NDIM(y),ARRAY_DIM_P(y)); izvec_add_rscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_Z(y),ARRAY_P0_Z(y),ARRAY_P0_R(x)[0],ARRAY_P1_R(x)[0]); }
    if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='r'){ z=array_allocate('R',ARRAY_NDIM(y),ARRAY_DIM_P(y)); irvec_add_rscalar(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P1_R(z),ARRAY_P0_R(y),ARRAY_P0_R(y),ARRAY_P0_R(x)[0],ARRAY_P1_R(x)[0]); }
    if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='c'){ z=array_allocate('C',ARRAY_NDIM(y),ARRAY_DIM_P(y)); icvec_add_rscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(y),ARRAY_P0_C(y),ARRAY_P0_R(x)[0],ARRAY_P1_R(x)[0]); }
    if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='R'){ z=array_allocate('R',ARRAY_NDIM(y),ARRAY_DIM_P(y)); irvec_add_rscalar(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P1_R(z),ARRAY_P0_R(y),ARRAY_P1_R(y),ARRAY_P0_R(x)[0],ARRAY_P1_R(x)[0]); }
    if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='C'){ z=array_allocate('C',ARRAY_NDIM(y),ARRAY_DIM_P(y)); icvec_add_dscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(y),ARRAY_P1_C(y),ARRAY_P0_D(x)[0],ARRAY_P0_D(x)[0]); }
    if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='C'){ z=array_allocate('C',ARRAY_NDIM(y),ARRAY_DIM_P(y)); icvec_add_zscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(y),ARRAY_P1_C(y),ARRAY_P0_Z(x)[0],ARRAY_P0_Z(x)[0]); }
    if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='C'){ z=array_allocate('C',ARRAY_NDIM(y),ARRAY_DIM_P(y)); icvec_add_rscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(y),ARRAY_P1_C(y),ARRAY_P0_R(x)[0],ARRAY_P0_R(x)[0]); }
    if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='C'){ z=array_allocate('C',ARRAY_NDIM(y),ARRAY_DIM_P(y)); icvec_add_cscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(y),ARRAY_P1_C(y),ARRAY_P0_C(x)[0],ARRAY_P0_C(x)[0]); }
    if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='C'){ z=array_allocate('C',ARRAY_NDIM(y),ARRAY_DIM_P(y)); icvec_add_rscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(y),ARRAY_P1_C(y),ARRAY_P0_R(x)[0],ARRAY_P1_R(x)[0]); }
    if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='d'){ z=array_allocate('C',ARRAY_NDIM(y),ARRAY_DIM_P(y)); idvec_add_cscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_D(y),ARRAY_P0_D(y),ARRAY_P0_C(x)[0],ARRAY_P1_C(x)[0]); }
    if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='z'){ z=array_allocate('C',ARRAY_NDIM(y),ARRAY_DIM_P(y)); izvec_add_cscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_Z(y),ARRAY_P0_Z(y),ARRAY_P0_C(x)[0],ARRAY_P1_C(x)[0]); }
    if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='r'){ z=array_allocate('C',ARRAY_NDIM(y),ARRAY_DIM_P(y)); irvec_add_cscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_R(y),ARRAY_P0_R(y),ARRAY_P0_C(x)[0],ARRAY_P1_C(x)[0]); }
    if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='c'){ z=array_allocate('C',ARRAY_NDIM(y),ARRAY_DIM_P(y)); icvec_add_cscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(y),ARRAY_P0_C(y),ARRAY_P0_C(x)[0],ARRAY_P1_C(x)[0]); }
    if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='R'){ z=array_allocate('C',ARRAY_NDIM(y),ARRAY_DIM_P(y)); irvec_add_cscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_R(y),ARRAY_P1_R(y),ARRAY_P0_C(x)[0],ARRAY_P1_C(x)[0]); }
    if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='C'){ z=array_allocate('C',ARRAY_NDIM(y),ARRAY_DIM_P(y)); icvec_add_cscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(y),ARRAY_P1_C(y),ARRAY_P0_C(x)[0],ARRAY_P1_C(x)[0]); }
  }
  if(z==NULL && (k=array_get_subdim(x,y))>0){
    if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='d'){ z=array_allocate('d',ARRAY_NDIM(y),ARRAY_DIM_P(y)); }
    if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='z'){ z=array_allocate('z',ARRAY_NDIM(y),ARRAY_DIM_P(y)); }
    if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='d'){ z=array_allocate('z',ARRAY_NDIM(y),ARRAY_DIM_P(y)); }
    if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='z'){ z=array_allocate('z',ARRAY_NDIM(y),ARRAY_DIM_P(y)); }
    if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='r'){ z=array_allocate('r',ARRAY_NDIM(y),ARRAY_DIM_P(y)); }
    if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='r'){ z=array_allocate('c',ARRAY_NDIM(y),ARRAY_DIM_P(y)); }
    if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='d'){ z=array_allocate('r',ARRAY_NDIM(y),ARRAY_DIM_P(y)); }
    if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='z'){ z=array_allocate('c',ARRAY_NDIM(y),ARRAY_DIM_P(y)); }
    if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='r'){ z=array_allocate('r',ARRAY_NDIM(y),ARRAY_DIM_P(y)); }
    if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='c'){ z=array_allocate('c',ARRAY_NDIM(y),ARRAY_DIM_P(y)); }
    if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='c'){ z=array_allocate('c',ARRAY_NDIM(y),ARRAY_DIM_P(y)); }
    if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='c'){ z=array_allocate('c',ARRAY_NDIM(y),ARRAY_DIM_P(y)); }
    if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='d'){ z=array_allocate('c',ARRAY_NDIM(y),ARRAY_DIM_P(y)); }
    if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='z'){ z=array_allocate('c',ARRAY_NDIM(y),ARRAY_DIM_P(y)); }
    if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='r'){ z=array_allocate('c',ARRAY_NDIM(y),ARRAY_DIM_P(y)); }
    if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='c'){ z=array_allocate('c',ARRAY_NDIM(y),ARRAY_DIM_P(y)); }
    if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='R'){ z=array_allocate('R',ARRAY_NDIM(y),ARRAY_DIM_P(y)); }
    if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='R'){ z=array_allocate('C',ARRAY_NDIM(y),ARRAY_DIM_P(y)); }
    if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='R'){ z=array_allocate('R',ARRAY_NDIM(y),ARRAY_DIM_P(y)); }
    if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='R'){ z=array_allocate('C',ARRAY_NDIM(y),ARRAY_DIM_P(y)); }
    if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='d'){ z=array_allocate('R',ARRAY_NDIM(y),ARRAY_DIM_P(y)); }
    if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='z'){ z=array_allocate('C',ARRAY_NDIM(y),ARRAY_DIM_P(y)); }
    if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='r'){ z=array_allocate('R',ARRAY_NDIM(y),ARRAY_DIM_P(y)); }
    if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='c'){ z=array_allocate('C',ARRAY_NDIM(y),ARRAY_DIM_P(y)); }
    if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='R'){ z=array_allocate('R',ARRAY_NDIM(y),ARRAY_DIM_P(y)); }
    if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='C'){ z=array_allocate('C',ARRAY_NDIM(y),ARRAY_DIM_P(y)); }
    if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='C'){ z=array_allocate('C',ARRAY_NDIM(y),ARRAY_DIM_P(y)); }
    if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='C'){ z=array_allocate('C',ARRAY_NDIM(y),ARRAY_DIM_P(y)); }
    if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='C'){ z=array_allocate('C',ARRAY_NDIM(y),ARRAY_DIM_P(y)); }
    if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='C'){ z=array_allocate('C',ARRAY_NDIM(y),ARRAY_DIM_P(y)); }
    if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='d'){ z=array_allocate('C',ARRAY_NDIM(y),ARRAY_DIM_P(y)); }
    if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='z'){ z=array_allocate('C',ARRAY_NDIM(y),ARRAY_DIM_P(y)); }
    if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='r'){ z=array_allocate('C',ARRAY_NDIM(y),ARRAY_DIM_P(y)); }
    if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='c'){ z=array_allocate('C',ARRAY_NDIM(y),ARRAY_DIM_P(y)); }
    if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='R'){ z=array_allocate('C',ARRAY_NDIM(y),ARRAY_DIM_P(y)); }
    if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='C'){ z=array_allocate('C',ARRAY_NDIM(y),ARRAY_DIM_P(y)); }
    for(i=0; i<ARRAY_LD_END(z)/ARRAY_LD(z,k-1); i++){
      if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='d'){  dvec_add_dvec(ARRAY_LD(z,k-1),ARRAY_P0_D(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_D(x),ARRAY_P0_D(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='z'){  dvec_add_zvec(ARRAY_LD(z,k-1),ARRAY_P0_Z(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_D(x),ARRAY_P0_Z(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='d'){  zvec_add_dvec(ARRAY_LD(z,k-1),ARRAY_P0_Z(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_Z(x),ARRAY_P0_D(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='z'){  zvec_add_zvec(ARRAY_LD(z,k-1),ARRAY_P0_Z(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_Z(x),ARRAY_P0_Z(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='r'){  dvec_add_rvec(ARRAY_LD(z,k-1),ARRAY_P0_R(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_D(x),ARRAY_P0_R(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='r'){  zvec_add_rvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_Z(x),ARRAY_P0_R(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='d'){  rvec_add_dvec(ARRAY_LD(z,k-1),ARRAY_P0_R(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_R(x),ARRAY_P0_D(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='z'){  rvec_add_zvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_R(x),ARRAY_P0_Z(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='r'){  rvec_add_rvec(ARRAY_LD(z,k-1),ARRAY_P0_R(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_R(x),ARRAY_P0_R(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='c'){  dvec_add_cvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_D(x),ARRAY_P0_C(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='c'){  zvec_add_cvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_Z(x),ARRAY_P0_C(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='c'){  rvec_add_cvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_R(x),ARRAY_P0_C(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='d'){  cvec_add_dvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_C(x),ARRAY_P0_D(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='z'){  cvec_add_zvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_C(x),ARRAY_P0_Z(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='r'){  cvec_add_rvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_C(x),ARRAY_P0_R(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='c'){  cvec_add_cvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_C(x),ARRAY_P0_C(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='R'){ idvec_add_rvec(ARRAY_LD(z,k-1),ARRAY_P0_R(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_R(z),ARRAY_P0_D(x),ARRAY_P0_D(x),ARRAY_P0_R(y)+i*ARRAY_LD(y,k-1),ARRAY_P1_R(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='R'){ izvec_add_rvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z),ARRAY_P0_Z(x),ARRAY_P0_Z(x),ARRAY_P0_R(y)+i*ARRAY_LD(y,k-1),ARRAY_P1_R(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='R'){ irvec_add_rvec(ARRAY_LD(z,k-1),ARRAY_P0_R(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_R(z),ARRAY_P0_R(x),ARRAY_P0_R(x),ARRAY_P0_R(y)+i*ARRAY_LD(y,k-1),ARRAY_P1_R(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='R'){ icvec_add_rvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z),ARRAY_P0_C(x),ARRAY_P0_C(x),ARRAY_P0_R(y)+i*ARRAY_LD(y,k-1),ARRAY_P1_R(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='d'){ irvec_add_dvec(ARRAY_LD(z,k-1),ARRAY_P0_R(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_R(z),ARRAY_P0_R(x),ARRAY_P1_R(x),ARRAY_P0_D(y)+i*ARRAY_LD(y,k-1),ARRAY_P0_D(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='z'){ irvec_add_zvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z),ARRAY_P0_R(x),ARRAY_P1_R(x),ARRAY_P0_Z(y)+i*ARRAY_LD(y,k-1),ARRAY_P0_Z(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='r'){ irvec_add_rvec(ARRAY_LD(z,k-1),ARRAY_P0_R(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_R(z),ARRAY_P0_R(x),ARRAY_P1_R(x),ARRAY_P0_R(y)+i*ARRAY_LD(y,k-1),ARRAY_P0_R(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='c'){ irvec_add_cvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z),ARRAY_P0_R(x),ARRAY_P1_R(x),ARRAY_P0_C(y)+i*ARRAY_LD(y,k-1),ARRAY_P0_C(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='R'){ irvec_add_rvec(ARRAY_LD(z,k-1),ARRAY_P0_R(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_R(z),ARRAY_P0_R(x),ARRAY_P1_R(x),ARRAY_P0_R(y)+i*ARRAY_LD(y,k-1),ARRAY_P1_R(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='C'){ idvec_add_cvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z),ARRAY_P0_D(x),ARRAY_P0_D(x),ARRAY_P0_C(y)+i*ARRAY_LD(y,k-1),ARRAY_P1_C(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='C'){ izvec_add_cvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z),ARRAY_P0_Z(x),ARRAY_P0_Z(x),ARRAY_P0_C(y)+i*ARRAY_LD(y,k-1),ARRAY_P1_C(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='C'){ irvec_add_cvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z),ARRAY_P0_R(x),ARRAY_P0_R(x),ARRAY_P0_C(y)+i*ARRAY_LD(y,k-1),ARRAY_P1_C(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='C'){ icvec_add_cvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z),ARRAY_P0_C(x),ARRAY_P0_C(x),ARRAY_P0_C(y)+i*ARRAY_LD(y,k-1),ARRAY_P1_C(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='C'){ irvec_add_cvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z),ARRAY_P0_R(x),ARRAY_P1_R(x),ARRAY_P0_C(y)+i*ARRAY_LD(y,k-1),ARRAY_P1_C(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='d'){ icvec_add_dvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z),ARRAY_P0_C(x),ARRAY_P1_C(x),ARRAY_P0_D(y)+i*ARRAY_LD(y,k-1),ARRAY_P0_D(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='z'){ icvec_add_zvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z),ARRAY_P0_C(x),ARRAY_P1_C(x),ARRAY_P0_Z(y)+i*ARRAY_LD(y,k-1),ARRAY_P0_Z(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='r'){ icvec_add_rvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z),ARRAY_P0_C(x),ARRAY_P1_C(x),ARRAY_P0_R(y)+i*ARRAY_LD(y,k-1),ARRAY_P0_R(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='c'){ icvec_add_cvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z),ARRAY_P0_C(x),ARRAY_P1_C(x),ARRAY_P0_C(y)+i*ARRAY_LD(y,k-1),ARRAY_P0_C(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='R'){ icvec_add_rvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z),ARRAY_P0_C(x),ARRAY_P1_C(x),ARRAY_P0_R(y)+i*ARRAY_LD(y,k-1),ARRAY_P1_R(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='C'){ icvec_add_cvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z),ARRAY_P0_C(x),ARRAY_P1_C(x),ARRAY_P0_C(y)+i*ARRAY_LD(y,k-1),ARRAY_P1_C(y)+i*ARRAY_LD(y,k-1)); }
    }
  }
  if(z==NULL && (k=array_get_subdim(y,x))>0){
    if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='d'){ z=array_allocate('d',ARRAY_NDIM(x),ARRAY_DIM_P(x)); }
    if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='z'){ z=array_allocate('z',ARRAY_NDIM(x),ARRAY_DIM_P(x)); }
    if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='d'){ z=array_allocate('z',ARRAY_NDIM(x),ARRAY_DIM_P(x)); }
    if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='z'){ z=array_allocate('z',ARRAY_NDIM(x),ARRAY_DIM_P(x)); }
    if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='r'){ z=array_allocate('r',ARRAY_NDIM(x),ARRAY_DIM_P(x)); }
    if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='r'){ z=array_allocate('c',ARRAY_NDIM(x),ARRAY_DIM_P(x)); }
    if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='d'){ z=array_allocate('r',ARRAY_NDIM(x),ARRAY_DIM_P(x)); }
    if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='z'){ z=array_allocate('c',ARRAY_NDIM(x),ARRAY_DIM_P(x)); }
    if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='r'){ z=array_allocate('r',ARRAY_NDIM(x),ARRAY_DIM_P(x)); }
    if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='c'){ z=array_allocate('c',ARRAY_NDIM(x),ARRAY_DIM_P(x)); }
    if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='c'){ z=array_allocate('c',ARRAY_NDIM(x),ARRAY_DIM_P(x)); }
    if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='c'){ z=array_allocate('c',ARRAY_NDIM(x),ARRAY_DIM_P(x)); }
    if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='d'){ z=array_allocate('c',ARRAY_NDIM(x),ARRAY_DIM_P(x)); }
    if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='z'){ z=array_allocate('c',ARRAY_NDIM(x),ARRAY_DIM_P(x)); }
    if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='r'){ z=array_allocate('c',ARRAY_NDIM(x),ARRAY_DIM_P(x)); }
    if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='c'){ z=array_allocate('c',ARRAY_NDIM(x),ARRAY_DIM_P(x)); }
    if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='R'){ z=array_allocate('R',ARRAY_NDIM(x),ARRAY_DIM_P(x)); }
    if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='R'){ z=array_allocate('C',ARRAY_NDIM(x),ARRAY_DIM_P(x)); }
    if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='R'){ z=array_allocate('R',ARRAY_NDIM(x),ARRAY_DIM_P(x)); }
    if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='R'){ z=array_allocate('C',ARRAY_NDIM(x),ARRAY_DIM_P(x)); }
    if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='d'){ z=array_allocate('R',ARRAY_NDIM(x),ARRAY_DIM_P(x)); }
    if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='z'){ z=array_allocate('C',ARRAY_NDIM(x),ARRAY_DIM_P(x)); }
    if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='r'){ z=array_allocate('R',ARRAY_NDIM(x),ARRAY_DIM_P(x)); }
    if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='c'){ z=array_allocate('C',ARRAY_NDIM(x),ARRAY_DIM_P(x)); }
    if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='R'){ z=array_allocate('R',ARRAY_NDIM(x),ARRAY_DIM_P(x)); }
    if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='C'){ z=array_allocate('C',ARRAY_NDIM(x),ARRAY_DIM_P(x)); }
    if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='C'){ z=array_allocate('C',ARRAY_NDIM(x),ARRAY_DIM_P(x)); }
    if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='C'){ z=array_allocate('C',ARRAY_NDIM(x),ARRAY_DIM_P(x)); }
    if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='C'){ z=array_allocate('C',ARRAY_NDIM(x),ARRAY_DIM_P(x)); }
    if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='C'){ z=array_allocate('C',ARRAY_NDIM(x),ARRAY_DIM_P(x)); }
    if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='d'){ z=array_allocate('C',ARRAY_NDIM(x),ARRAY_DIM_P(x)); }
    if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='z'){ z=array_allocate('C',ARRAY_NDIM(x),ARRAY_DIM_P(x)); }
    if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='r'){ z=array_allocate('C',ARRAY_NDIM(x),ARRAY_DIM_P(x)); }
    if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='c'){ z=array_allocate('C',ARRAY_NDIM(x),ARRAY_DIM_P(x)); }
    if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='R'){ z=array_allocate('C',ARRAY_NDIM(x),ARRAY_DIM_P(x)); }
    if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='C'){ z=array_allocate('C',ARRAY_NDIM(x),ARRAY_DIM_P(x)); }
    for(i=0; i<ARRAY_LD_END(z)/ARRAY_LD(z,k-1); i++){
      if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='d'){  dvec_add_dvec(ARRAY_LD(z,k-1),ARRAY_P0_D(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_D(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_D(y)); }
      if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='z'){  dvec_add_zvec(ARRAY_LD(z,k-1),ARRAY_P0_Z(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_D(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_Z(y)); }
      if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='d'){  zvec_add_dvec(ARRAY_LD(z,k-1),ARRAY_P0_Z(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_Z(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_D(y)); }
      if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='z'){  zvec_add_zvec(ARRAY_LD(z,k-1),ARRAY_P0_Z(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_Z(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_Z(y)); }
      if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='r'){  dvec_add_rvec(ARRAY_LD(z,k-1),ARRAY_P0_R(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_D(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_R(y)); }
      if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='r'){  zvec_add_rvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_Z(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_R(y)); }
      if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='d'){  rvec_add_dvec(ARRAY_LD(z,k-1),ARRAY_P0_R(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_R(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_D(y)); }
      if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='z'){  rvec_add_zvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_R(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_Z(y)); }
      if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='r'){  rvec_add_rvec(ARRAY_LD(z,k-1),ARRAY_P0_R(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_R(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_R(y)); }
      if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='c'){  dvec_add_cvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_D(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_C(y)); }
      if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='c'){  zvec_add_cvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_Z(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_C(y)); }
      if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='c'){  rvec_add_cvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_R(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_C(y)); }
      if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='d'){  cvec_add_dvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_C(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_D(y)); }
      if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='z'){  cvec_add_zvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_C(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_Z(y)); }
      if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='r'){  cvec_add_rvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_C(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_R(y)); }
      if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='c'){  cvec_add_cvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_C(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_C(y)); }
      if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='R'){ idvec_add_rvec(ARRAY_LD(z,k-1),ARRAY_P0_R(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_R(z),ARRAY_P0_D(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_D(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_R(y),ARRAY_P1_R(y)); }
      if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='R'){ izvec_add_rvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z),ARRAY_P0_Z(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_Z(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_R(y),ARRAY_P1_R(y)); }
      if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='R'){ irvec_add_rvec(ARRAY_LD(z,k-1),ARRAY_P0_R(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_R(z),ARRAY_P0_R(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_R(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_R(y),ARRAY_P1_R(y)); }
      if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='R'){ icvec_add_rvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z),ARRAY_P0_C(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_C(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_R(y),ARRAY_P1_R(y)); }
      if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='d'){ irvec_add_dvec(ARRAY_LD(z,k-1),ARRAY_P0_R(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_R(z),ARRAY_P0_R(x)+i*ARRAY_LD(x,k-1),ARRAY_P1_R(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_D(y),ARRAY_P0_D(y)); }
      if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='z'){ irvec_add_zvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z),ARRAY_P0_R(x)+i*ARRAY_LD(x,k-1),ARRAY_P1_R(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_Z(y),ARRAY_P0_Z(y)); }
      if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='r'){ irvec_add_rvec(ARRAY_LD(z,k-1),ARRAY_P0_R(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_R(z),ARRAY_P0_R(x)+i*ARRAY_LD(x,k-1),ARRAY_P1_R(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_R(y),ARRAY_P0_R(y)); }
      if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='c'){ irvec_add_cvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z),ARRAY_P0_R(x)+i*ARRAY_LD(x,k-1),ARRAY_P1_R(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_C(y),ARRAY_P0_C(y)); }
      if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='R'){ irvec_add_rvec(ARRAY_LD(z,k-1),ARRAY_P0_R(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_R(z),ARRAY_P0_R(x)+i*ARRAY_LD(x,k-1),ARRAY_P1_R(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_R(y),ARRAY_P1_R(y)); }
      if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='C'){ idvec_add_cvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z),ARRAY_P0_D(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_D(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_C(y),ARRAY_P1_C(y)); }
      if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='C'){ izvec_add_cvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z),ARRAY_P0_Z(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_Z(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_C(y),ARRAY_P1_C(y)); }
      if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='C'){ irvec_add_cvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z),ARRAY_P0_R(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_R(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_C(y),ARRAY_P1_C(y)); }
      if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='C'){ icvec_add_cvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z),ARRAY_P0_C(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_C(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_C(y),ARRAY_P1_C(y)); }
      if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='C'){ irvec_add_cvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z),ARRAY_P0_R(x)+i*ARRAY_LD(x,k-1),ARRAY_P1_R(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_C(y),ARRAY_P1_C(y)); }
      if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='d'){ icvec_add_dvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z),ARRAY_P0_C(x)+i*ARRAY_LD(x,k-1),ARRAY_P1_C(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_D(y),ARRAY_P0_D(y)); }
      if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='z'){ icvec_add_zvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z),ARRAY_P0_C(x)+i*ARRAY_LD(x,k-1),ARRAY_P1_C(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_Z(y),ARRAY_P0_Z(y)); }
      if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='r'){ icvec_add_rvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z),ARRAY_P0_C(x)+i*ARRAY_LD(x,k-1),ARRAY_P1_C(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_R(y),ARRAY_P0_R(y)); }
      if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='c'){ icvec_add_cvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z),ARRAY_P0_C(x)+i*ARRAY_LD(x,k-1),ARRAY_P1_C(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_C(y),ARRAY_P0_C(y)); }
      if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='R'){ icvec_add_rvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z),ARRAY_P0_C(x)+i*ARRAY_LD(x,k-1),ARRAY_P1_C(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_R(y),ARRAY_P1_R(y)); }
      if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='C'){ icvec_add_cvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z),ARRAY_P0_C(x)+i*ARRAY_LD(x,k-1),ARRAY_P1_C(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_C(y),ARRAY_P1_C(y)); }
    }
  }
  return z;
}

/**
 @breif z=x-y
 */
array *array_sub(array *x, array *y)
{
  array *z=NULL;
  if(x==NULL || y==NULL){ return z; }

  if(array_same_dim_check(x,y)){
    if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='d'){ z=array_allocate('d',ARRAY_NDIM(x),ARRAY_DIM_P(x));  dvec_sub_dvec(ARRAY_SIZE(z),ARRAY_P0_D(z),ARRAY_P0_D(x),ARRAY_P0_D(y)); }
    if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='z'){ z=array_allocate('z',ARRAY_NDIM(x),ARRAY_DIM_P(x));  dvec_sub_zvec(ARRAY_SIZE(z),ARRAY_P0_Z(z),ARRAY_P0_D(x),ARRAY_P0_Z(y)); }
    if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='d'){ z=array_allocate('z',ARRAY_NDIM(x),ARRAY_DIM_P(x));  zvec_sub_dvec(ARRAY_SIZE(z),ARRAY_P0_Z(z),ARRAY_P0_Z(x),ARRAY_P0_D(y)); }
    if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='z'){ z=array_allocate('z',ARRAY_NDIM(x),ARRAY_DIM_P(x));  zvec_sub_zvec(ARRAY_SIZE(z),ARRAY_P0_Z(z),ARRAY_P0_Z(x),ARRAY_P0_Z(y)); }
    if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='r'){ z=array_allocate('r',ARRAY_NDIM(x),ARRAY_DIM_P(x));  dvec_sub_rvec(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P0_D(x),ARRAY_P0_R(y)); }
    if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='r'){ z=array_allocate('c',ARRAY_NDIM(x),ARRAY_DIM_P(x));  zvec_sub_rvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_Z(x),ARRAY_P0_R(y)); }
    if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='d'){ z=array_allocate('r',ARRAY_NDIM(x),ARRAY_DIM_P(x));  rvec_sub_dvec(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P0_R(x),ARRAY_P0_D(y)); }
    if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='z'){ z=array_allocate('c',ARRAY_NDIM(x),ARRAY_DIM_P(x));  rvec_sub_zvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_R(x),ARRAY_P0_Z(y)); }
    if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='r'){ z=array_allocate('r',ARRAY_NDIM(x),ARRAY_DIM_P(x));  rvec_sub_rvec(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P0_R(x),ARRAY_P0_R(y)); }
    if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='c'){ z=array_allocate('c',ARRAY_NDIM(x),ARRAY_DIM_P(x));  dvec_sub_cvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_D(x),ARRAY_P0_C(y)); }
    if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='c'){ z=array_allocate('c',ARRAY_NDIM(x),ARRAY_DIM_P(x));  zvec_sub_cvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_Z(x),ARRAY_P0_C(y)); }
    if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='c'){ z=array_allocate('c',ARRAY_NDIM(x),ARRAY_DIM_P(x));  rvec_sub_cvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_R(x),ARRAY_P0_C(y)); }
    if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='d'){ z=array_allocate('c',ARRAY_NDIM(x),ARRAY_DIM_P(x));  cvec_sub_dvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_C(x),ARRAY_P0_D(y)); }
    if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='z'){ z=array_allocate('c',ARRAY_NDIM(x),ARRAY_DIM_P(x));  cvec_sub_zvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_C(x),ARRAY_P0_Z(y)); }
    if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='r'){ z=array_allocate('c',ARRAY_NDIM(x),ARRAY_DIM_P(x));  cvec_sub_rvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_C(x),ARRAY_P0_R(y)); }
    if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='c'){ z=array_allocate('c',ARRAY_NDIM(x),ARRAY_DIM_P(x));  cvec_sub_cvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_C(x),ARRAY_P0_C(y)); }
    if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='R'){ z=array_allocate('R',ARRAY_NDIM(x),ARRAY_DIM_P(x)); idvec_sub_rvec(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P1_R(z),ARRAY_P0_D(x),ARRAY_P0_D(x),ARRAY_P0_R(y),ARRAY_P1_R(y)); }
    if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='R'){ z=array_allocate('C',ARRAY_NDIM(x),ARRAY_DIM_P(x)); izvec_sub_rvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_Z(x),ARRAY_P0_Z(x),ARRAY_P0_R(y),ARRAY_P1_R(y)); }
    if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='R'){ z=array_allocate('R',ARRAY_NDIM(x),ARRAY_DIM_P(x)); irvec_sub_rvec(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P1_R(z),ARRAY_P0_R(x),ARRAY_P0_R(x),ARRAY_P0_R(y),ARRAY_P1_R(y)); }
    if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='R'){ z=array_allocate('C',ARRAY_NDIM(x),ARRAY_DIM_P(x)); icvec_sub_rvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(x),ARRAY_P0_C(x),ARRAY_P0_R(y),ARRAY_P1_R(y)); }
    if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='d'){ z=array_allocate('R',ARRAY_NDIM(x),ARRAY_DIM_P(x)); irvec_sub_dvec(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P1_R(z),ARRAY_P0_R(x),ARRAY_P1_R(x),ARRAY_P0_D(y),ARRAY_P0_D(y)); }
    if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='z'){ z=array_allocate('C',ARRAY_NDIM(x),ARRAY_DIM_P(x)); irvec_sub_zvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_R(x),ARRAY_P1_R(x),ARRAY_P0_Z(y),ARRAY_P0_Z(y)); }
    if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='r'){ z=array_allocate('R',ARRAY_NDIM(x),ARRAY_DIM_P(x)); irvec_sub_rvec(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P1_R(z),ARRAY_P0_R(x),ARRAY_P1_R(x),ARRAY_P0_R(y),ARRAY_P0_R(y)); }
    if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='c'){ z=array_allocate('C',ARRAY_NDIM(x),ARRAY_DIM_P(x)); irvec_sub_cvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_R(x),ARRAY_P1_R(x),ARRAY_P0_C(y),ARRAY_P0_C(y)); }
    if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='R'){ z=array_allocate('R',ARRAY_NDIM(x),ARRAY_DIM_P(x)); irvec_sub_rvec(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P1_R(z),ARRAY_P0_R(x),ARRAY_P1_R(x),ARRAY_P0_R(y),ARRAY_P1_R(y)); }
    if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='C'){ z=array_allocate('C',ARRAY_NDIM(x),ARRAY_DIM_P(x)); idvec_sub_cvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_D(x),ARRAY_P0_D(x),ARRAY_P0_C(y),ARRAY_P1_C(y)); }
    if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='C'){ z=array_allocate('C',ARRAY_NDIM(x),ARRAY_DIM_P(x)); izvec_sub_cvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_Z(x),ARRAY_P0_Z(x),ARRAY_P0_C(y),ARRAY_P1_C(y)); }
    if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='C'){ z=array_allocate('C',ARRAY_NDIM(x),ARRAY_DIM_P(x)); irvec_sub_cvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_R(x),ARRAY_P0_R(x),ARRAY_P0_C(y),ARRAY_P1_C(y)); }
    if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='C'){ z=array_allocate('C',ARRAY_NDIM(x),ARRAY_DIM_P(x)); icvec_sub_cvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(x),ARRAY_P0_C(x),ARRAY_P0_C(y),ARRAY_P1_C(y)); }
    if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='C'){ z=array_allocate('C',ARRAY_NDIM(x),ARRAY_DIM_P(x)); irvec_sub_cvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_R(x),ARRAY_P1_R(x),ARRAY_P0_C(y),ARRAY_P1_C(y)); }
    if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='d'){ z=array_allocate('C',ARRAY_NDIM(x),ARRAY_DIM_P(x)); icvec_sub_dvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(x),ARRAY_P1_C(x),ARRAY_P0_D(y),ARRAY_P0_D(y)); }
    if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='z'){ z=array_allocate('C',ARRAY_NDIM(x),ARRAY_DIM_P(x)); icvec_sub_zvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(x),ARRAY_P1_C(x),ARRAY_P0_Z(y),ARRAY_P0_Z(y)); }
    if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='r'){ z=array_allocate('C',ARRAY_NDIM(x),ARRAY_DIM_P(x)); icvec_sub_rvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(x),ARRAY_P1_C(x),ARRAY_P0_R(y),ARRAY_P0_R(y)); }
    if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='c'){ z=array_allocate('C',ARRAY_NDIM(x),ARRAY_DIM_P(x)); icvec_sub_cvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(x),ARRAY_P1_C(x),ARRAY_P0_C(y),ARRAY_P0_C(y)); }
    if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='R'){ z=array_allocate('C',ARRAY_NDIM(x),ARRAY_DIM_P(x)); icvec_sub_rvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(x),ARRAY_P1_C(x),ARRAY_P0_R(y),ARRAY_P1_R(y)); }
    if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='C'){ z=array_allocate('C',ARRAY_NDIM(x),ARRAY_DIM_P(x)); icvec_sub_cvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(x),ARRAY_P1_C(x),ARRAY_P0_C(y),ARRAY_P1_C(y)); }
  }else if(array_is_scalar(y)){
    if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='d'){ z=array_allocate('d',ARRAY_NDIM(x),ARRAY_DIM_P(x));  dvec_sub_dscalar(ARRAY_SIZE(z),ARRAY_P0_D(z),ARRAY_P0_D(x),ARRAY_P0_D(y)[0]); }
    if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='z'){ z=array_allocate('z',ARRAY_NDIM(x),ARRAY_DIM_P(x));  dvec_sub_zscalar(ARRAY_SIZE(z),ARRAY_P0_Z(z),ARRAY_P0_D(x),ARRAY_P0_Z(y)[0]); }
    if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='d'){ z=array_allocate('z',ARRAY_NDIM(x),ARRAY_DIM_P(x));  zvec_sub_dscalar(ARRAY_SIZE(z),ARRAY_P0_Z(z),ARRAY_P0_Z(x),ARRAY_P0_D(y)[0]); }
    if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='z'){ z=array_allocate('z',ARRAY_NDIM(x),ARRAY_DIM_P(x));  zvec_sub_zscalar(ARRAY_SIZE(z),ARRAY_P0_Z(z),ARRAY_P0_Z(x),ARRAY_P0_Z(y)[0]); }
    if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='r'){ z=array_allocate('r',ARRAY_NDIM(x),ARRAY_DIM_P(x));  dvec_sub_rscalar(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P0_D(x),ARRAY_P0_R(y)[0]); }
    if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='r'){ z=array_allocate('c',ARRAY_NDIM(x),ARRAY_DIM_P(x));  zvec_sub_rscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_Z(x),ARRAY_P0_R(y)[0]); }
    if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='d'){ z=array_allocate('r',ARRAY_NDIM(x),ARRAY_DIM_P(x));  rvec_sub_dscalar(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P0_R(x),ARRAY_P0_D(y)[0]); }
    if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='z'){ z=array_allocate('c',ARRAY_NDIM(x),ARRAY_DIM_P(x));  rvec_sub_zscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_R(x),ARRAY_P0_Z(y)[0]); }
    if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='r'){ z=array_allocate('r',ARRAY_NDIM(x),ARRAY_DIM_P(x));  rvec_sub_rscalar(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P0_R(x),ARRAY_P0_R(y)[0]); }
    if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='c'){ z=array_allocate('c',ARRAY_NDIM(x),ARRAY_DIM_P(x));  dvec_sub_cscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_D(x),ARRAY_P0_C(y)[0]); }
    if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='c'){ z=array_allocate('c',ARRAY_NDIM(x),ARRAY_DIM_P(x));  zvec_sub_cscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_Z(x),ARRAY_P0_C(y)[0]); }
    if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='c'){ z=array_allocate('c',ARRAY_NDIM(x),ARRAY_DIM_P(x));  rvec_sub_cscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_R(x),ARRAY_P0_C(y)[0]); }
    if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='d'){ z=array_allocate('c',ARRAY_NDIM(x),ARRAY_DIM_P(x));  cvec_sub_dscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_C(x),ARRAY_P0_D(y)[0]); }
    if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='z'){ z=array_allocate('c',ARRAY_NDIM(x),ARRAY_DIM_P(x));  cvec_sub_zscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_C(x),ARRAY_P0_Z(y)[0]); }
    if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='r'){ z=array_allocate('c',ARRAY_NDIM(x),ARRAY_DIM_P(x));  cvec_sub_rscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_C(x),ARRAY_P0_R(y)[0]); }
    if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='c'){ z=array_allocate('c',ARRAY_NDIM(x),ARRAY_DIM_P(x));  cvec_sub_cscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_C(x),ARRAY_P0_C(y)[0]); }
    if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='R'){ z=array_allocate('R',ARRAY_NDIM(x),ARRAY_DIM_P(x)); idvec_sub_rscalar(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P1_R(z),ARRAY_P0_D(x),ARRAY_P0_D(x),ARRAY_P0_R(y)[0],ARRAY_P1_R(y)[0]); }
    if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='R'){ z=array_allocate('C',ARRAY_NDIM(x),ARRAY_DIM_P(x)); izvec_sub_rscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_Z(x),ARRAY_P0_Z(x),ARRAY_P0_R(y)[0],ARRAY_P1_R(y)[0]); }
    if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='R'){ z=array_allocate('R',ARRAY_NDIM(x),ARRAY_DIM_P(x)); irvec_sub_rscalar(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P1_R(z),ARRAY_P0_R(x),ARRAY_P0_R(x),ARRAY_P0_R(y)[0],ARRAY_P1_R(y)[0]); }
    if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='R'){ z=array_allocate('C',ARRAY_NDIM(x),ARRAY_DIM_P(x)); icvec_sub_rscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(x),ARRAY_P0_C(x),ARRAY_P0_R(y)[0],ARRAY_P1_R(y)[0]); }
    if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='d'){ z=array_allocate('R',ARRAY_NDIM(x),ARRAY_DIM_P(x)); irvec_sub_dscalar(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P1_R(z),ARRAY_P0_R(x),ARRAY_P1_R(x),ARRAY_P0_D(y)[0],ARRAY_P0_D(y)[0]); }
    if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='z'){ z=array_allocate('C',ARRAY_NDIM(x),ARRAY_DIM_P(x)); irvec_sub_zscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_R(x),ARRAY_P1_R(x),ARRAY_P0_Z(y)[0],ARRAY_P0_Z(y)[0]); }
    if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='r'){ z=array_allocate('R',ARRAY_NDIM(x),ARRAY_DIM_P(x)); irvec_sub_rscalar(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P1_R(z),ARRAY_P0_R(x),ARRAY_P1_R(x),ARRAY_P0_R(y)[0],ARRAY_P0_R(y)[0]); }
    if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='c'){ z=array_allocate('C',ARRAY_NDIM(x),ARRAY_DIM_P(x)); irvec_sub_cscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_R(x),ARRAY_P1_R(x),ARRAY_P0_C(y)[0],ARRAY_P0_C(y)[0]); }
    if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='R'){ z=array_allocate('R',ARRAY_NDIM(x),ARRAY_DIM_P(x)); irvec_sub_rscalar(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P1_R(z),ARRAY_P0_R(x),ARRAY_P1_R(x),ARRAY_P0_R(y)[0],ARRAY_P1_R(y)[0]); }
    if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='C'){ z=array_allocate('C',ARRAY_NDIM(x),ARRAY_DIM_P(x)); idvec_sub_cscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_D(x),ARRAY_P0_D(x),ARRAY_P0_C(y)[0],ARRAY_P1_C(y)[0]); }
    if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='C'){ z=array_allocate('C',ARRAY_NDIM(x),ARRAY_DIM_P(x)); izvec_sub_cscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_Z(x),ARRAY_P0_Z(x),ARRAY_P0_C(y)[0],ARRAY_P1_C(y)[0]); }
    if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='C'){ z=array_allocate('C',ARRAY_NDIM(x),ARRAY_DIM_P(x)); irvec_sub_cscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_R(x),ARRAY_P0_R(x),ARRAY_P0_C(y)[0],ARRAY_P1_C(y)[0]); }
    if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='C'){ z=array_allocate('C',ARRAY_NDIM(x),ARRAY_DIM_P(x)); icvec_sub_cscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(x),ARRAY_P0_C(x),ARRAY_P0_C(y)[0],ARRAY_P1_C(y)[0]); }
    if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='C'){ z=array_allocate('C',ARRAY_NDIM(x),ARRAY_DIM_P(x)); irvec_sub_cscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_R(x),ARRAY_P1_R(x),ARRAY_P0_C(y)[0],ARRAY_P1_C(y)[0]); }
    if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='d'){ z=array_allocate('C',ARRAY_NDIM(x),ARRAY_DIM_P(x)); icvec_sub_dscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(x),ARRAY_P1_C(x),ARRAY_P0_D(y)[0],ARRAY_P0_D(y)[0]); }
    if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='z'){ z=array_allocate('C',ARRAY_NDIM(x),ARRAY_DIM_P(x)); icvec_sub_zscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(x),ARRAY_P1_C(x),ARRAY_P0_Z(y)[0],ARRAY_P0_Z(y)[0]); }
    if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='r'){ z=array_allocate('C',ARRAY_NDIM(x),ARRAY_DIM_P(x)); icvec_sub_rscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(x),ARRAY_P1_C(x),ARRAY_P0_R(y)[0],ARRAY_P0_R(y)[0]); }
    if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='c'){ z=array_allocate('C',ARRAY_NDIM(x),ARRAY_DIM_P(x)); icvec_sub_cscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(x),ARRAY_P1_C(x),ARRAY_P0_C(y)[0],ARRAY_P0_C(y)[0]); }
    if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='R'){ z=array_allocate('C',ARRAY_NDIM(x),ARRAY_DIM_P(x)); icvec_sub_rscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(x),ARRAY_P1_C(x),ARRAY_P0_R(y)[0],ARRAY_P1_R(y)[0]); }
    if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='C'){ z=array_allocate('C',ARRAY_NDIM(x),ARRAY_DIM_P(x)); icvec_sub_cscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(x),ARRAY_P1_C(x),ARRAY_P0_C(y)[0],ARRAY_P1_C(y)[0]); }
  }else if(array_is_scalar(x)){
    if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='d'){ z=array_allocate('d',ARRAY_NDIM(y),ARRAY_DIM_P(y));  dscalar_sub_dvec(ARRAY_SIZE(z),ARRAY_P0_D(z),ARRAY_P0_D(x)[0],ARRAY_P0_D(y)); }
    if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='z'){ z=array_allocate('z',ARRAY_NDIM(y),ARRAY_DIM_P(y));  dscalar_sub_zvec(ARRAY_SIZE(z),ARRAY_P0_Z(z),ARRAY_P0_D(x)[0],ARRAY_P0_Z(y)); }
    if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='d'){ z=array_allocate('z',ARRAY_NDIM(y),ARRAY_DIM_P(y));  zscalar_sub_dvec(ARRAY_SIZE(z),ARRAY_P0_Z(z),ARRAY_P0_Z(x)[0],ARRAY_P0_D(y)); }
    if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='z'){ z=array_allocate('z',ARRAY_NDIM(y),ARRAY_DIM_P(y));  zscalar_sub_zvec(ARRAY_SIZE(z),ARRAY_P0_Z(z),ARRAY_P0_Z(x)[0],ARRAY_P0_Z(y)); }
    if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='r'){ z=array_allocate('r',ARRAY_NDIM(y),ARRAY_DIM_P(y));  dscalar_sub_rvec(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P0_D(x)[0],ARRAY_P0_R(y)); }
    if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='r'){ z=array_allocate('c',ARRAY_NDIM(y),ARRAY_DIM_P(y));  zscalar_sub_rvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_Z(x)[0],ARRAY_P0_R(y)); }
    if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='d'){ z=array_allocate('r',ARRAY_NDIM(y),ARRAY_DIM_P(y));  rscalar_sub_dvec(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P0_R(x)[0],ARRAY_P0_D(y)); }
    if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='z'){ z=array_allocate('c',ARRAY_NDIM(y),ARRAY_DIM_P(y));  rscalar_sub_zvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_R(x)[0],ARRAY_P0_Z(y)); }
    if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='r'){ z=array_allocate('r',ARRAY_NDIM(y),ARRAY_DIM_P(y));  rscalar_sub_rvec(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P0_R(x)[0],ARRAY_P0_R(y)); }
    if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='c'){ z=array_allocate('c',ARRAY_NDIM(y),ARRAY_DIM_P(y));  dscalar_sub_cvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_D(x)[0],ARRAY_P0_C(y)); }
    if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='c'){ z=array_allocate('c',ARRAY_NDIM(y),ARRAY_DIM_P(y));  zscalar_sub_cvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_Z(x)[0],ARRAY_P0_C(y)); }
    if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='c'){ z=array_allocate('c',ARRAY_NDIM(y),ARRAY_DIM_P(y));  rscalar_sub_cvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_R(x)[0],ARRAY_P0_C(y)); }
    if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='d'){ z=array_allocate('c',ARRAY_NDIM(y),ARRAY_DIM_P(y));  cscalar_sub_dvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_C(x)[0],ARRAY_P0_D(y)); }
    if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='z'){ z=array_allocate('c',ARRAY_NDIM(y),ARRAY_DIM_P(y));  cscalar_sub_zvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_C(x)[0],ARRAY_P0_Z(y)); }
    if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='r'){ z=array_allocate('c',ARRAY_NDIM(y),ARRAY_DIM_P(y));  cscalar_sub_rvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_C(x)[0],ARRAY_P0_R(y)); }
    if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='c'){ z=array_allocate('c',ARRAY_NDIM(y),ARRAY_DIM_P(y));  cscalar_sub_cvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_C(x)[0],ARRAY_P0_C(y)); }
    if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='R'){ z=array_allocate('R',ARRAY_NDIM(y),ARRAY_DIM_P(y)); idscalar_sub_rvec(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P1_R(z),ARRAY_P0_D(x)[0],ARRAY_P0_D(x)[0],ARRAY_P0_R(y),ARRAY_P1_R(y)); }
    if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='R'){ z=array_allocate('C',ARRAY_NDIM(y),ARRAY_DIM_P(y)); izscalar_sub_rvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_Z(x)[0],ARRAY_P0_Z(x)[0],ARRAY_P0_R(y),ARRAY_P1_R(y)); }
    if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='R'){ z=array_allocate('R',ARRAY_NDIM(y),ARRAY_DIM_P(y)); irscalar_sub_rvec(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P1_R(z),ARRAY_P0_R(x)[0],ARRAY_P0_R(x)[0],ARRAY_P0_R(y),ARRAY_P1_R(y)); }
    if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='R'){ z=array_allocate('C',ARRAY_NDIM(y),ARRAY_DIM_P(y)); icscalar_sub_rvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(x)[0],ARRAY_P0_C(x)[0],ARRAY_P0_R(y),ARRAY_P1_R(y)); }
    if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='d'){ z=array_allocate('R',ARRAY_NDIM(y),ARRAY_DIM_P(y)); irscalar_sub_dvec(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P1_R(z),ARRAY_P0_R(x)[0],ARRAY_P1_R(x)[0],ARRAY_P0_D(y),ARRAY_P0_D(y)); }
    if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='z'){ z=array_allocate('C',ARRAY_NDIM(y),ARRAY_DIM_P(y)); irscalar_sub_zvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_R(x)[0],ARRAY_P1_R(x)[0],ARRAY_P0_Z(y),ARRAY_P0_Z(y)); }
    if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='r'){ z=array_allocate('R',ARRAY_NDIM(y),ARRAY_DIM_P(y)); irscalar_sub_rvec(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P1_R(z),ARRAY_P0_R(x)[0],ARRAY_P1_R(x)[0],ARRAY_P0_R(y),ARRAY_P0_R(y)); }
    if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='c'){ z=array_allocate('C',ARRAY_NDIM(y),ARRAY_DIM_P(y)); irscalar_sub_cvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_R(x)[0],ARRAY_P1_R(x)[0],ARRAY_P0_C(y),ARRAY_P0_C(y)); }
    if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='R'){ z=array_allocate('R',ARRAY_NDIM(y),ARRAY_DIM_P(y)); irscalar_sub_rvec(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P1_R(z),ARRAY_P0_R(x)[0],ARRAY_P1_R(x)[0],ARRAY_P0_R(y),ARRAY_P1_R(y)); }
    if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='C'){ z=array_allocate('C',ARRAY_NDIM(y),ARRAY_DIM_P(y)); idscalar_sub_cvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_D(x)[0],ARRAY_P0_D(x)[0],ARRAY_P0_C(y),ARRAY_P1_C(y)); }
    if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='C'){ z=array_allocate('C',ARRAY_NDIM(y),ARRAY_DIM_P(y)); izscalar_sub_cvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_Z(x)[0],ARRAY_P0_Z(x)[0],ARRAY_P0_C(y),ARRAY_P1_C(y)); }
    if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='C'){ z=array_allocate('C',ARRAY_NDIM(y),ARRAY_DIM_P(y)); irscalar_sub_cvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_R(x)[0],ARRAY_P0_R(x)[0],ARRAY_P0_C(y),ARRAY_P1_C(y)); }
    if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='C'){ z=array_allocate('C',ARRAY_NDIM(y),ARRAY_DIM_P(y)); icscalar_sub_cvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(x)[0],ARRAY_P0_C(x)[0],ARRAY_P0_C(y),ARRAY_P1_C(y)); }
    if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='C'){ z=array_allocate('C',ARRAY_NDIM(y),ARRAY_DIM_P(y)); irscalar_sub_cvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_R(x)[0],ARRAY_P1_R(x)[0],ARRAY_P0_C(y),ARRAY_P1_C(y)); }
    if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='d'){ z=array_allocate('C',ARRAY_NDIM(y),ARRAY_DIM_P(y)); icscalar_sub_dvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(x)[0],ARRAY_P1_C(x)[0],ARRAY_P0_D(y),ARRAY_P0_D(y)); }
    if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='z'){ z=array_allocate('C',ARRAY_NDIM(y),ARRAY_DIM_P(y)); icscalar_sub_zvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(x)[0],ARRAY_P1_C(x)[0],ARRAY_P0_Z(y),ARRAY_P0_Z(y)); }
    if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='r'){ z=array_allocate('C',ARRAY_NDIM(y),ARRAY_DIM_P(y)); icscalar_sub_rvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(x)[0],ARRAY_P1_C(x)[0],ARRAY_P0_R(y),ARRAY_P0_R(y)); }
    if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='c'){ z=array_allocate('C',ARRAY_NDIM(y),ARRAY_DIM_P(y)); icscalar_sub_cvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(x)[0],ARRAY_P1_C(x)[0],ARRAY_P0_C(y),ARRAY_P0_C(y)); }
    if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='R'){ z=array_allocate('C',ARRAY_NDIM(y),ARRAY_DIM_P(y)); icscalar_sub_rvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(x)[0],ARRAY_P1_C(x)[0],ARRAY_P0_R(y),ARRAY_P1_R(y)); }
    if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='C'){ z=array_allocate('C',ARRAY_NDIM(y),ARRAY_DIM_P(y)); icscalar_sub_cvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(x)[0],ARRAY_P1_C(x)[0],ARRAY_P0_C(y),ARRAY_P1_C(y)); }
  }
  return z;
}


/** @} */

///////////////////////////////////////////////

/** @name array型の３入力演算に関する関数 */
/** @{ */
/** @} */










//EOF
