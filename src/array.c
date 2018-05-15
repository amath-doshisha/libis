#include "is_array.h"
#include "is_strings.h"
#include "is_svec.h"
#include "is_ivec.h"
#include "is_dvec.h"
#include "is_zvec.h"
#include "is_rvec.h"
#include "is_cvec.h"
#include "is_irmulti.h"
#include "is_icmulti.h"
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
 @brief array型のメモリ割当
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
 @brief array型のメモリの解放
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
 @brief array型の配列のメモリ割当
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
 @brief array型の配列のメモリ解放
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
 @brief array型のコピー y=x
 */
array *array_copy(array *x)
{
  array *y=NULL;
  if(x==NULL){ return y; }
  y=array_allocate(ARRAY_TYPE(x),ARRAY_NDIM(x),ARRAY_DIM_P(x));
  if(ARRAY_TYPE(y)=='s'){ svec_copy(ARRAY_SIZE(y),ARRAY_P0_S(y),ARRAY_P0_S(x)); }
  if(ARRAY_TYPE(y)=='i'){ ivec_copy(ARRAY_SIZE(y),ARRAY_P0_I(y),ARRAY_P0_I(x)); }
  if(ARRAY_TYPE(y)=='d'){ dvec_copy_dvec(ARRAY_SIZE(y),ARRAY_P0_D(y),ARRAY_P0_D(x)); }
  if(ARRAY_TYPE(y)=='z'){ zvec_copy_zvec(ARRAY_SIZE(y),ARRAY_P0_Z(y),ARRAY_P0_Z(x)); }
  if(ARRAY_TYPE(y)=='r'){ rvec_copy_rvec(ARRAY_SIZE(y),ARRAY_P0_R(y),ARRAY_P0_R(x)); }
  if(ARRAY_TYPE(y)=='c'){ cvec_copy_cvec(ARRAY_SIZE(y),ARRAY_P0_C(y),ARRAY_P0_C(x)); }
  if(ARRAY_TYPE(y)=='R'){ rvec_copy_rvec(ARRAY_SIZE(y),ARRAY_P0_R(y),ARRAY_P0_R(x)); rvec_copy_rvec(ARRAY_SIZE(y),ARRAY_P1_R(y),ARRAY_P1_R(x)); }
  if(ARRAY_TYPE(y)=='C'){ cvec_copy_cvec(ARRAY_SIZE(y),ARRAY_P0_C(y),ARRAY_P0_C(x)); cvec_copy_cvec(ARRAY_SIZE(y),ARRAY_P1_C(y),ARRAY_P1_C(x)); }
  return y;
}

/**
 @brief array型のコピー y=x
 */
array *array_clone(array *x)
{
  array *y=NULL;
  if(x==NULL){ return y; }
  y=array_allocate(ARRAY_TYPE(x),ARRAY_NDIM(x),ARRAY_DIM_P(x));
  if(ARRAY_TYPE(y)=='s'){ svec_copy (ARRAY_SIZE(y),ARRAY_P0_S(y),ARRAY_P0_S(x)); }
  if(ARRAY_TYPE(y)=='i'){ ivec_copy (ARRAY_SIZE(y),ARRAY_P0_I(y),ARRAY_P0_I(x)); }
  if(ARRAY_TYPE(y)=='d'){ dvec_copy_dvec (ARRAY_SIZE(y),ARRAY_P0_D(y),ARRAY_P0_D(x)); }
  if(ARRAY_TYPE(y)=='z'){ zvec_copy_zvec (ARRAY_SIZE(y),ARRAY_P0_Z(y),ARRAY_P0_Z(x)); }
  if(ARRAY_TYPE(y)=='r'){ rvec_clone_rvec(ARRAY_SIZE(y),ARRAY_P0_R(y),ARRAY_P0_R(x)); }
  if(ARRAY_TYPE(y)=='c'){ cvec_clone_cvec(ARRAY_SIZE(y),ARRAY_P0_C(y),ARRAY_P0_C(x)); }
  if(ARRAY_TYPE(y)=='R'){ rvec_clone_rvec(ARRAY_SIZE(y),ARRAY_P0_R(y),ARRAY_P0_R(x)); rvec_clone_rvec(ARRAY_SIZE(y),ARRAY_P1_R(y),ARRAY_P1_R(x)); }
  if(ARRAY_TYPE(y)=='C'){ cvec_clone_cvec(ARRAY_SIZE(y),ARRAY_P0_C(y),ARRAY_P0_C(x)); cvec_clone_cvec(ARRAY_SIZE(y),ARRAY_P1_C(y),ARRAY_P1_C(x)); }
  return y;
}

/** @} */

///////////////////////////////////////////////

/** @name array型の値の設定に関する関数 */
/** @{ */

/**
 @brief array型の値をすべて同じ値に設定
 */
void array_set_all_d(array *x, double y)
{
  int i;
  if(x==NULL){ return; }
  if(ARRAY_TYPE(x)=='a'){ for(i=0; i<ARRAY_SIZE(x); i++){ array_set_all_d(ARRAY_AVEC(x,i),y); } }
  else if(ARRAY_TYPE(x)=='s'){ svec_set_all_d(ARRAY_SIZE(x),ARRAY_P0_S(x),y,"%.7g"); }
  else if(ARRAY_TYPE(x)=='i'){ ivec_set_all  (ARRAY_SIZE(x),ARRAY_P0_I(x),(int)y); }
  else if(ARRAY_TYPE(x)=='d'){ dvec_set_all_d  (ARRAY_SIZE(x),ARRAY_P0_D(x),y); }
  else if(ARRAY_TYPE(x)=='z'){ zvec_set_all_d(ARRAY_SIZE(x),ARRAY_P0_Z(x),y); }
  else if(ARRAY_TYPE(x)=='r'){ rvec_set_all_d(ARRAY_SIZE(x),ARRAY_P0_R(x),y); }
  else if(ARRAY_TYPE(x)=='c'){ cvec_set_all_d(ARRAY_SIZE(x),ARRAY_P0_C(x),y); }
  else if(ARRAY_TYPE(x)=='R'){ rvec_set_all_d(ARRAY_SIZE(x),ARRAY_P0_R(x),y); rvec_set_all_d(ARRAY_SIZE(x),ARRAY_P1_R(x),y); }
  else if(ARRAY_TYPE(x)=='C'){ cvec_set_all_d(ARRAY_SIZE(x),ARRAY_P0_C(x),y); cvec_set_all_d(ARRAY_SIZE(x),ARRAY_P1_C(x),y); }
}

/**
 @brief array型の値をすべて0に設定
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
 @brief array型の値をすべて1に設定
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
 @brief array型の値をすべてnanに設定
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
 @brief array型の値をすべてinfに設定
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
 @brief array型の値を乱数値[a0,a1]に設定
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
  else if(ARRAY_TYPE(x)=='R'){ rvec_set_rand(ARRAY_SIZE(x),ARRAY_P0_R(x),a1-a0,a0); rvec_copy_rvec(ARRAY_SIZE(x),ARRAY_P1_R(x),ARRAY_P0_R(x)); }
  else if(ARRAY_TYPE(x)=='C'){ cvec_set_rand(ARRAY_SIZE(x),ARRAY_P0_C(x),a1-a0,a0); cvec_copy_cvec(ARRAY_SIZE(x),ARRAY_P1_C(x),ARRAY_P0_C(x)); }
}

/**
 @brief array型の値を0,1,2,..に設定
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
 @brief array型の内部型の確認
 */
int array_type_check(int type)
{
  return (type=='a' || type=='s' || type=='i' || type=='d' || type=='z' || type=='r' || type=='c' || type=='R' || type=='C');
}

/**
 @brief array型は空か？
 */
int array_is_empty(array *x)
{
  return (x==NULL || ARRAY_NDIM(x)<=0 || ARRAY_SIZE(x)<=0 || ARRAY_DIM_P(x)==NULL || ARRAY_LD_P(x)==NULL || ARRAY_P0(x)==NULL);
}

/*
 @brief array型は同じサイズかの判定．
 */
int array_same_dim_check(array *x, array *y)
{
  int n,i;
  if(x==NULL || y==NULL){ return 0; }
  n=MAX2(ARRAY_NDIM(x),ARRAY_NDIM(y));
  for(i=0; i<n; i++){ if(ARRAY_DIM(x,i)!=ARRAY_DIM(y,i)){ return 0; } }
  return 1;
}

/*
 @brief array型は実質上のサイズの取得（サイズ1を無視する）．
 */
int array_get_ndim(array *x)
{
  int i;
  if(x==NULL){ return 0; }
  for(i=ARRAY_NDIM(x)-1; i>=0; i--){
    if(ARRAY_DIM(x,i)!=1){ return i+1; }
  }
  return 0;
}

/*
 @brief array型のxはyの部分集合であるかの判定．返値は部分集合の次元．
 */
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
 @brief array型のxとyは互換性のあるサイズであるかの判定（サイズ1は無視する）．
 */
int array_compatible_dim_check(array * x, array *y)
{
  int i,n;
  n=MAX2(ARRAY_NDIM(x),ARRAY_N(y));
  for(i=0; i<n; i++){
    if(ARRAY_DIM(x,i)!=ARRAY_DIM(y,i) && ARRAY_DIM(x,i)!=1 && ARRAY_DIM(y,i)!=1){ return 0; }
  }
  return 1;
}

/*
 @brief array型のxとyの両方が含まれる最小のサイズの取得．
 */
int *array_get_inclusion_dim(array * x, array *y)
{
  int i,n,*dim=NULL;
  n=MAX2(ARRAY_NDIM(x),ARRAY_N(y));
  dim=ivec_allocate(n);
  for(i=0; i<n; i++){
    dim[i]=MAX2(ARRAY_DIM(x,i),ARRAY_DIM(y,i));
  }
  return dim;
}

/*
 @brief array型のxとyを数学的な演算をしたとき，両方が含まれるタイプの取得．
 */
char array_get_inclusion_type(char x, char y)
{
  char z='-';
  if(x=='d' && y=='d'){ z='d'; }
  if(x=='d' && y=='z'){ z='z'; }
  if(x=='z' && y=='d'){ z='z'; }
  if(x=='z' && y=='z'){ z='z'; }
  if(x=='d' && y=='r'){ z='r'; }
  if(x=='z' && y=='r'){ z='c'; }
  if(x=='r' && y=='d'){ z='r'; }
  if(x=='r' && y=='z'){ z='c'; }
  if(x=='r' && y=='r'){ z='r'; }
  if(x=='d' && y=='c'){ z='c'; }
  if(x=='z' && y=='c'){ z='c'; }
  if(x=='r' && y=='c'){ z='c'; }
  if(x=='c' && y=='d'){ z='c'; }
  if(x=='c' && y=='z'){ z='c'; }
  if(x=='c' && y=='r'){ z='c'; }
  if(x=='c' && y=='c'){ z='c'; }
  if(x=='d' && y=='R'){ z='R'; }
  if(x=='z' && y=='R'){ z='C'; }
  if(x=='r' && y=='R'){ z='R'; }
  if(x=='c' && y=='R'){ z='C'; }
  if(x=='R' && y=='d'){ z='R'; }
  if(x=='R' && y=='z'){ z='C'; }
  if(x=='R' && y=='r'){ z='R'; }
  if(x=='R' && y=='c'){ z='C'; }
  if(x=='R' && y=='R'){ z='R'; }
  if(x=='d' && y=='C'){ z='C'; }
  if(x=='z' && y=='C'){ z='C'; }
  if(x=='r' && y=='C'){ z='C'; }
  if(x=='c' && y=='C'){ z='C'; }
  if(x=='R' && y=='C'){ z='C'; }
  if(x=='C' && y=='d'){ z='C'; }
  if(x=='C' && y=='z'){ z='C'; }
  if(x=='C' && y=='r'){ z='C'; }
  if(x=='C' && y=='c'){ z='C'; }
  if(x=='C' && y=='R'){ z='C'; }
  if(x=='C' && y=='C'){ z='C'; }
  return z;
}

/*
 @brief array型はスカラーであるかの判定．
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
 @brief array型のメンバの内容を出力
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
  else                    { printf("  dim: ["); ivec_put(ARRAY_NDIM(x),ARRAY_DIM_P(x),","); printf("]\n"); }
  if(ARRAY_LD_P(x)==NULL){ printf("  LD: NULL\n"); }
  else                   { printf("  LD:  ["); ivec_put(ARRAY_NDIM(x),ARRAY_LD_P(x),","); printf("]\n"); }
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
 @brief array型の出力
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
      if(ARRAY_TYPE(x)=='s'){  svec_print(ARRAY_SIZE(x),&ARRAY_SVEC(x,0),NULL); }
      if(ARRAY_TYPE(x)=='i'){  ivec_print(ARRAY_SIZE(x),&ARRAY_IVEC(x,0),NULL); }
      if(ARRAY_TYPE(x)=='d'){  dvec_print(ARRAY_SIZE(x),&ARRAY_DVEC(x,0),NULL,format,digits); }
      if(ARRAY_TYPE(x)=='z'){  zvec_print(ARRAY_SIZE(x),&ARRAY_ZVEC(x,0),NULL,format,digits); }
      if(ARRAY_TYPE(x)=='r'){  rvec_print(ARRAY_SIZE(x),&ARRAY_RVEC(x,0),NULL,format,digits); }
      if(ARRAY_TYPE(x)=='c'){  cvec_print(ARRAY_SIZE(x),&ARRAY_CVEC(x,0),NULL,format,digits); }
      if(ARRAY_TYPE(x)=='R'){ irvec_print(ARRAY_SIZE(x),&ARRAY_RVEC0(x,0),&ARRAY_RVEC1(x,0),NULL,format,digits); }
      if(ARRAY_TYPE(x)=='C'){ icvec_print(ARRAY_SIZE(x),&ARRAY_CVEC0(x,0),&ARRAY_CVEC1(x,0),NULL,format,digits); }
    }else{
      I=ivec_allocate(ARRAY_NDIM(x));
      for(k=0; k<ARRAY_SIZE(x); k+=ARRAY_LD1(x)){
        array_get_index(I,x,k);
        printf("%s%s",(name==NULL?"":name),(ARRAY_NDIM(x)==2?"":"[][]"));
        for(i=2; i<ARRAY_NDIM(x); i++){ printf("[%d]",I[i]); }
        printf("=\n");
        if(ARRAY_TYPE(x)=='s'){  smat_print(ARRAY_M(x),ARRAY_N(x),ARRAY_P0_S(x)+k,ARRAY_LD0(x),NULL); }
        if(ARRAY_TYPE(x)=='i'){  imat_print(ARRAY_M(x),ARRAY_N(x),ARRAY_P0_I(x)+k,ARRAY_LD0(x),NULL); }
        if(ARRAY_TYPE(x)=='d'){  dmat_print(ARRAY_M(x),ARRAY_N(x),ARRAY_P0_D(x)+k,ARRAY_LD0(x),NULL,format,digits); }
        if(ARRAY_TYPE(x)=='z'){  zmat_print(ARRAY_M(x),ARRAY_N(x),ARRAY_P0_Z(x)+k,ARRAY_LD0(x),NULL,format,digits); }
        if(ARRAY_TYPE(x)=='r'){  rmat_print(ARRAY_M(x),ARRAY_N(x),ARRAY_P0_R(x)+k,ARRAY_LD0(x),NULL,format,digits); }
        if(ARRAY_TYPE(x)=='c'){  cmat_print(ARRAY_M(x),ARRAY_N(x),ARRAY_P0_C(x)+k,ARRAY_LD0(x),NULL,format,digits); }
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

int array_get_position(int *index, array *x)
{
  int i,k=0;
  if(x==NULL){ return -1; }
  k=index[0];
  for(i=1; i<ARRAY_NDIM(x); i++){
    k=k+index[i]*ARRAY_LD(x,i-1);
  }
  return k;
}

/** @} */

///////////////////////////////////////////////

/** @name array型のキャストに関する関数 */
/** @{ */

/**
 @brief array型の型変換 y=cast(x)
 */
array *array_get(int type, array *x)
{
  array *y=NULL;
  if(x==NULL){ return y; }
  if(type=='d'){
    if(ARRAY_TYPE(x)=='d'){ y=array_allocate('d',ARRAY_NDIM(x),ARRAY_DIM_P(x)); dvec_set_dvec (ARRAY_SIZE(y),ARRAY_P0_D(y),ARRAY_P0_D(x)); }
    if(ARRAY_TYPE(x)=='z'){ y=array_allocate('d',ARRAY_NDIM(x),ARRAY_DIM_P(x)); dvec_set_zvec (ARRAY_SIZE(y),ARRAY_P0_D(y),ARRAY_P0_Z(x)); }
    if(ARRAY_TYPE(x)=='r'){ y=array_allocate('d',ARRAY_NDIM(x),ARRAY_DIM_P(x)); dvec_set_rvec (ARRAY_SIZE(y),ARRAY_P0_D(y),ARRAY_P0_R(x)); }
    if(ARRAY_TYPE(x)=='c'){ y=array_allocate('d',ARRAY_NDIM(x),ARRAY_DIM_P(x)); dvec_set_cvec (ARRAY_SIZE(y),ARRAY_P0_D(y),ARRAY_P0_C(x)); }
    if(ARRAY_TYPE(x)=='R'){ y=array_allocate('d',ARRAY_NDIM(x),ARRAY_DIM_P(x)); dvec_set_irvec(ARRAY_SIZE(y),ARRAY_P0_D(y),ARRAY_P0_R(x),ARRAY_P1_R(x)); }
    if(ARRAY_TYPE(x)=='C'){ y=array_allocate('d',ARRAY_NDIM(x),ARRAY_DIM_P(x)); dvec_set_icvec(ARRAY_SIZE(y),ARRAY_P0_D(y),ARRAY_P0_C(x),ARRAY_P1_C(x)); }
  }
  if(type=='z'){
    if(ARRAY_TYPE(x)=='d'){ y=array_allocate('z',ARRAY_NDIM(x),ARRAY_DIM_P(x)); zvec_set_dvec (ARRAY_SIZE(y),ARRAY_P0_Z(y),ARRAY_P0_D(x)); }
    if(ARRAY_TYPE(x)=='z'){ y=array_allocate('z',ARRAY_NDIM(x),ARRAY_DIM_P(x)); zvec_set_zvec (ARRAY_SIZE(y),ARRAY_P0_Z(y),ARRAY_P0_Z(x)); }
    if(ARRAY_TYPE(x)=='r'){ y=array_allocate('z',ARRAY_NDIM(x),ARRAY_DIM_P(x)); zvec_set_rvec (ARRAY_SIZE(y),ARRAY_P0_Z(y),ARRAY_P0_R(x)); }
    if(ARRAY_TYPE(x)=='c'){ y=array_allocate('z',ARRAY_NDIM(x),ARRAY_DIM_P(x)); zvec_set_cvec (ARRAY_SIZE(y),ARRAY_P0_Z(y),ARRAY_P0_C(x)); }
    if(ARRAY_TYPE(x)=='R'){ y=array_allocate('z',ARRAY_NDIM(x),ARRAY_DIM_P(x)); zvec_set_irvec(ARRAY_SIZE(y),ARRAY_P0_Z(y),ARRAY_P0_R(x),ARRAY_P1_R(x)); }
    if(ARRAY_TYPE(x)=='C'){ y=array_allocate('z',ARRAY_NDIM(x),ARRAY_DIM_P(x)); zvec_set_icvec(ARRAY_SIZE(y),ARRAY_P0_Z(y),ARRAY_P0_C(x),ARRAY_P1_C(x)); }
  }
  if(type=='r'){
    if(ARRAY_TYPE(x)=='d'){ y=array_allocate('r',ARRAY_NDIM(x),ARRAY_DIM_P(x)); rvec_set_dvec (ARRAY_SIZE(y),ARRAY_P0_R(y),ARRAY_P0_D(x)); }
    if(ARRAY_TYPE(x)=='z'){ y=array_allocate('r',ARRAY_NDIM(x),ARRAY_DIM_P(x)); rvec_set_zvec (ARRAY_SIZE(y),ARRAY_P0_R(y),ARRAY_P0_Z(x)); }
    if(ARRAY_TYPE(x)=='r'){ y=array_allocate('r',ARRAY_NDIM(x),ARRAY_DIM_P(x)); rvec_set_rvec (ARRAY_SIZE(y),ARRAY_P0_R(y),ARRAY_P0_R(x)); }
    if(ARRAY_TYPE(x)=='c'){ y=array_allocate('r',ARRAY_NDIM(x),ARRAY_DIM_P(x)); rvec_set_cvec (ARRAY_SIZE(y),ARRAY_P0_R(y),ARRAY_P0_C(x)); }
    if(ARRAY_TYPE(x)=='R'){ y=array_allocate('r',ARRAY_NDIM(x),ARRAY_DIM_P(x)); rvec_set_irvec(ARRAY_SIZE(y),ARRAY_P0_R(y),ARRAY_P0_R(x),ARRAY_P1_R(x)); }
    if(ARRAY_TYPE(x)=='C'){ y=array_allocate('r',ARRAY_NDIM(x),ARRAY_DIM_P(x)); rvec_set_icvec(ARRAY_SIZE(y),ARRAY_P0_R(y),ARRAY_P0_C(x),ARRAY_P1_C(x)); }
  }
  if(type=='c'){
    if(ARRAY_TYPE(x)=='d'){ y=array_allocate('c',ARRAY_NDIM(x),ARRAY_DIM_P(x)); cvec_set_dvec (ARRAY_SIZE(y),ARRAY_P0_C(y),ARRAY_P0_D(x)); }
    if(ARRAY_TYPE(x)=='z'){ y=array_allocate('c',ARRAY_NDIM(x),ARRAY_DIM_P(x)); cvec_set_zvec (ARRAY_SIZE(y),ARRAY_P0_C(y),ARRAY_P0_Z(x)); }
    if(ARRAY_TYPE(x)=='r'){ y=array_allocate('c',ARRAY_NDIM(x),ARRAY_DIM_P(x)); cvec_set_rvec (ARRAY_SIZE(y),ARRAY_P0_C(y),ARRAY_P0_R(x)); }
    if(ARRAY_TYPE(x)=='c'){ y=array_allocate('c',ARRAY_NDIM(x),ARRAY_DIM_P(x)); cvec_set_cvec (ARRAY_SIZE(y),ARRAY_P0_C(y),ARRAY_P0_C(x)); }
    if(ARRAY_TYPE(x)=='R'){ y=array_allocate('c',ARRAY_NDIM(x),ARRAY_DIM_P(x)); cvec_set_irvec(ARRAY_SIZE(y),ARRAY_P0_C(y),ARRAY_P0_R(x),ARRAY_P1_R(x)); }
    if(ARRAY_TYPE(x)=='C'){ y=array_allocate('c',ARRAY_NDIM(x),ARRAY_DIM_P(x)); cvec_set_icvec(ARRAY_SIZE(y),ARRAY_P0_C(y),ARRAY_P0_C(x),ARRAY_P1_C(x)); }
  }
  if(type=='R'){
    if(ARRAY_TYPE(x)=='d'){ y=array_allocate('R',ARRAY_NDIM(x),ARRAY_DIM_P(x)); irvec_set_dvec(ARRAY_SIZE(y),ARRAY_P0_R(y),ARRAY_P1_R(y),ARRAY_P0_D(x),ARRAY_P0_D(x)); }
    if(ARRAY_TYPE(x)=='z'){ y=array_allocate('R',ARRAY_NDIM(x),ARRAY_DIM_P(x)); irvec_set_zvec(ARRAY_SIZE(y),ARRAY_P0_R(y),ARRAY_P1_R(y),ARRAY_P0_Z(x),ARRAY_P0_Z(x)); }
    if(ARRAY_TYPE(x)=='r'){ y=array_allocate('R',ARRAY_NDIM(x),ARRAY_DIM_P(x)); irvec_set_rvec(ARRAY_SIZE(y),ARRAY_P0_R(y),ARRAY_P1_R(y),ARRAY_P0_R(x),ARRAY_P0_R(x)); }
    if(ARRAY_TYPE(x)=='c'){ y=array_allocate('R',ARRAY_NDIM(x),ARRAY_DIM_P(x)); irvec_set_cvec(ARRAY_SIZE(y),ARRAY_P0_R(y),ARRAY_P1_R(y),ARRAY_P0_C(x),ARRAY_P0_C(x)); }
    if(ARRAY_TYPE(x)=='R'){ y=array_allocate('R',ARRAY_NDIM(x),ARRAY_DIM_P(x)); irvec_set_rvec(ARRAY_SIZE(y),ARRAY_P0_R(y),ARRAY_P1_R(y),ARRAY_P0_R(x),ARRAY_P1_R(x)); }
    if(ARRAY_TYPE(x)=='C'){ y=array_allocate('R',ARRAY_NDIM(x),ARRAY_DIM_P(x)); irvec_set_cvec(ARRAY_SIZE(y),ARRAY_P0_R(y),ARRAY_P1_R(y),ARRAY_P0_C(x),ARRAY_P1_C(x)); }
  }
  if(type=='C'){
    if(ARRAY_TYPE(x)=='d'){ y=array_allocate('C',ARRAY_NDIM(x),ARRAY_DIM_P(x)); icvec_set_dvec(ARRAY_SIZE(y),ARRAY_P0_C(y),ARRAY_P1_C(y),ARRAY_P0_D(x),ARRAY_P0_D(x)); }
    if(ARRAY_TYPE(x)=='z'){ y=array_allocate('C',ARRAY_NDIM(x),ARRAY_DIM_P(x)); icvec_set_zvec(ARRAY_SIZE(y),ARRAY_P0_C(y),ARRAY_P1_C(y),ARRAY_P0_Z(x),ARRAY_P0_Z(x)); }
    if(ARRAY_TYPE(x)=='r'){ y=array_allocate('C',ARRAY_NDIM(x),ARRAY_DIM_P(x)); icvec_set_rvec(ARRAY_SIZE(y),ARRAY_P0_C(y),ARRAY_P1_C(y),ARRAY_P0_R(x),ARRAY_P0_R(x)); }
    if(ARRAY_TYPE(x)=='c'){ y=array_allocate('C',ARRAY_NDIM(x),ARRAY_DIM_P(x)); icvec_set_cvec(ARRAY_SIZE(y),ARRAY_P0_C(y),ARRAY_P1_C(y),ARRAY_P0_C(x),ARRAY_P0_C(x)); }
    if(ARRAY_TYPE(x)=='R'){ y=array_allocate('C',ARRAY_NDIM(x),ARRAY_DIM_P(x)); icvec_set_rvec(ARRAY_SIZE(y),ARRAY_P0_C(y),ARRAY_P1_C(y),ARRAY_P0_R(x),ARRAY_P1_R(x)); }
    if(ARRAY_TYPE(x)=='C'){ y=array_allocate('C',ARRAY_NDIM(x),ARRAY_DIM_P(x)); icvec_set_cvec(ARRAY_SIZE(y),ARRAY_P0_C(y),ARRAY_P1_C(y),ARRAY_P0_C(x),ARRAY_P1_C(x)); }
  }
  return y;
}

/**
 @brief array型の型変換 y=char(x)
 */
array *array_get_char(array *x, char format, int digits)
{
  array *y=NULL;
  if(!(format=='f' || format=='e' || format=='g')){ format='g'; }
  if(digits<0){ digits=6; }
  if(x==NULL){ return y; }
  if(ARRAY_TYPE(x)=='s'){ y=array_allocate('s',ARRAY_NDIM(x),ARRAY_DIM_P(x));  svec_copy (ARRAY_SIZE(y),ARRAY_P0_S(y),ARRAY_P0_S(x)); }
  if(ARRAY_TYPE(x)=='i'){ y=array_allocate('s',ARRAY_NDIM(x),ARRAY_DIM_P(x));  ivec_get_s(ARRAY_SIZE(y),ARRAY_P0_S(y),ARRAY_P0_I(x)); }
  if(ARRAY_TYPE(x)=='d'){ y=array_allocate('s',ARRAY_NDIM(x),ARRAY_DIM_P(x));  dvec_get_svec(ARRAY_SIZE(y),ARRAY_P0_S(y),ARRAY_P0_D(x),format,digits); }
  if(ARRAY_TYPE(x)=='z'){ y=array_allocate('s',ARRAY_NDIM(x),ARRAY_DIM_P(x));  zvec_get_svec(ARRAY_SIZE(y),ARRAY_P0_S(y),ARRAY_P0_Z(x),format,digits); }
  if(ARRAY_TYPE(x)=='r'){ y=array_allocate('s',ARRAY_NDIM(x),ARRAY_DIM_P(x));  rvec_get_svec(ARRAY_SIZE(y),ARRAY_P0_S(y),ARRAY_P0_R(x),format,digits); }
  if(ARRAY_TYPE(x)=='c'){ y=array_allocate('s',ARRAY_NDIM(x),ARRAY_DIM_P(x));  cvec_get_svec(ARRAY_SIZE(y),ARRAY_P0_S(y),ARRAY_P0_C(x),format,digits); }
  if(ARRAY_TYPE(x)=='R'){ y=array_allocate('s',ARRAY_NDIM(x),ARRAY_DIM_P(x)); irvec_get_svec(ARRAY_SIZE(y),ARRAY_P0_S(y),ARRAY_P0_R(x),ARRAY_P1_R(x),format,digits); }
  if(ARRAY_TYPE(x)=='C'){ y=array_allocate('s',ARRAY_NDIM(x),ARRAY_DIM_P(x)); icvec_get_svec(ARRAY_SIZE(y),ARRAY_P0_S(y),ARRAY_P0_C(x),ARRAY_P1_C(x),format,digits); }
  return y;
}

/**
 @brief array型の型変換 y=int(x)
 */
array *array_get_int(array *x)
{
  array *y=NULL;
  if(x==NULL){ return y; }
  if(ARRAY_TYPE(x)=='s'){ y=array_allocate('i',ARRAY_NDIM(x),ARRAY_DIM_P(x));  ivec_set_s (ARRAY_SIZE(y),ARRAY_P0_I(y),ARRAY_P0_S(x)); }
  if(ARRAY_TYPE(x)=='i'){ y=array_allocate('i',ARRAY_NDIM(x),ARRAY_DIM_P(x));  ivec_copy  (ARRAY_SIZE(y),ARRAY_P0_I(y),ARRAY_P0_I(x)); }
  if(ARRAY_TYPE(x)=='d'){ y=array_allocate('i',ARRAY_NDIM(x),ARRAY_DIM_P(x));  dvec_get_ivec(ARRAY_SIZE(y),ARRAY_P0_I(y),ARRAY_P0_D(x)); }
  if(ARRAY_TYPE(x)=='z'){ y=array_allocate('i',ARRAY_NDIM(x),ARRAY_DIM_P(x));  zvec_get_ivec(ARRAY_SIZE(y),ARRAY_P0_I(y),ARRAY_P0_Z(x)); }
  if(ARRAY_TYPE(x)=='r'){ y=array_allocate('i',ARRAY_NDIM(x),ARRAY_DIM_P(x));  rvec_get_ivec(ARRAY_SIZE(y),ARRAY_P0_I(y),ARRAY_P0_R(x)); }
  if(ARRAY_TYPE(x)=='c'){ y=array_allocate('i',ARRAY_NDIM(x),ARRAY_DIM_P(x));  cvec_get_ivec(ARRAY_SIZE(y),ARRAY_P0_I(y),ARRAY_P0_C(x)); }
  if(ARRAY_TYPE(x)=='R'){ y=array_allocate('i',ARRAY_NDIM(x),ARRAY_DIM_P(x)); irvec_get_ivec(ARRAY_SIZE(y),ARRAY_P0_I(y),ARRAY_P0_R(x),ARRAY_P1_R(x)); }
  if(ARRAY_TYPE(x)=='C'){ y=array_allocate('i',ARRAY_NDIM(x),ARRAY_DIM_P(x)); icvec_get_ivec(ARRAY_SIZE(y),ARRAY_P0_I(y),ARRAY_P0_C(x),ARRAY_P1_C(x)); }
  return y;
}

/**
 @brief array型の型変換 y=real(x)
 */
array *array_get_real(array *x)
{
  array *y=NULL;
  if(x==NULL){ return y; }
  if(ARRAY_TYPE(x)=='s'){ y=array_allocate('s',ARRAY_NDIM(x),ARRAY_DIM_P(x));  svec_copy(ARRAY_SIZE(y),ARRAY_P0_S(y),ARRAY_P0_S(x)); }
  if(ARRAY_TYPE(x)=='i'){ y=array_allocate('i',ARRAY_NDIM(x),ARRAY_DIM_P(x));  ivec_copy(ARRAY_SIZE(y),ARRAY_P0_I(y),ARRAY_P0_I(x)); }
  if(ARRAY_TYPE(x)=='d'){ y=array_allocate('d',ARRAY_NDIM(x),ARRAY_DIM_P(x));  dvec_copy_dvec(ARRAY_SIZE(y),ARRAY_P0_D(y),ARRAY_P0_D(x)); }
  if(ARRAY_TYPE(x)=='z'){ y=array_allocate('d',ARRAY_NDIM(x),ARRAY_DIM_P(x));  dvec_real_zvec(ARRAY_SIZE(y),ARRAY_P0_D(y),ARRAY_P0_Z(x)); }
  if(ARRAY_TYPE(x)=='r'){ y=array_allocate('r',ARRAY_NDIM(x),ARRAY_DIM_P(x));  rvec_copy_rvec(ARRAY_SIZE(y),ARRAY_P0_R(y),ARRAY_P0_R(x)); }
  if(ARRAY_TYPE(x)=='c'){ y=array_allocate('r',ARRAY_NDIM(x),ARRAY_DIM_P(x));  rvec_real_cvec(ARRAY_SIZE(y),ARRAY_P0_R(y),ARRAY_P0_C(x)); }
  if(ARRAY_TYPE(x)=='R'){ y=array_allocate('R',ARRAY_NDIM(x),ARRAY_DIM_P(x)); irvec_copy_rvec(ARRAY_SIZE(y),ARRAY_P0_R(y),ARRAY_P1_R(y),ARRAY_P0_R(x),ARRAY_P1_R(x)); }
  if(ARRAY_TYPE(x)=='C'){ y=array_allocate('R',ARRAY_NDIM(x),ARRAY_DIM_P(x)); irvec_real_cvec(ARRAY_SIZE(y),ARRAY_P0_R(y),ARRAY_P1_R(y),ARRAY_P0_C(x),ARRAY_P1_C(x)); }
  return y;
}

/**
 @brief array型の型変換 y=complex(x)
 */
array *array_get_complex(array *x)
{
  array *y=NULL;
  if(x==NULL){ return y; }
  if(ARRAY_TYPE(x)=='s'){ y=array_allocate('s',ARRAY_NDIM(x),ARRAY_DIM_P(x));  svec_copy  (ARRAY_SIZE(y),ARRAY_P0_S(y),ARRAY_P0_S(x)); }
  if(ARRAY_TYPE(x)=='i'){ y=array_allocate('i',ARRAY_NDIM(x),ARRAY_DIM_P(x));  ivec_copy  (ARRAY_SIZE(y),ARRAY_P0_I(y),ARRAY_P0_I(x)); }
  if(ARRAY_TYPE(x)=='d'){ y=array_allocate('z',ARRAY_NDIM(x),ARRAY_DIM_P(x));  zvec_set_dvec (ARRAY_SIZE(y),ARRAY_P0_Z(y),ARRAY_P0_D(x)); }
  if(ARRAY_TYPE(x)=='z'){ y=array_allocate('z',ARRAY_NDIM(x),ARRAY_DIM_P(x));  zvec_copy_zvec(ARRAY_SIZE(y),ARRAY_P0_Z(y),ARRAY_P0_Z(x)); }
  if(ARRAY_TYPE(x)=='r'){ y=array_allocate('c',ARRAY_NDIM(x),ARRAY_DIM_P(x));  cvec_set_rvec(ARRAY_SIZE(y),ARRAY_P0_C(y),ARRAY_P0_R(x)); }
  if(ARRAY_TYPE(x)=='c'){ y=array_allocate('c',ARRAY_NDIM(x),ARRAY_DIM_P(x));  cvec_copy_cvec  (ARRAY_SIZE(y),ARRAY_P0_C(y),ARRAY_P0_C(x)); }
  if(ARRAY_TYPE(x)=='R'){ y=array_allocate('C',ARRAY_NDIM(x),ARRAY_DIM_P(x)); icvec_copy_rvec(ARRAY_SIZE(y),ARRAY_P0_C(y),ARRAY_P1_C(y),ARRAY_P0_R(x),ARRAY_P1_R(x)); }
  if(ARRAY_TYPE(x)=='C'){ y=array_allocate('C',ARRAY_NDIM(x),ARRAY_DIM_P(x)); icvec_copy_cvec  (ARRAY_SIZE(y),ARRAY_P0_C(y),ARRAY_P1_C(y),ARRAY_P0_C(x),ARRAY_P1_C(x)); }
  return y;
}

/**
 @brief array型の型変換 y=double(x)
 */
array *array_get_double(array *x)
{
  int flag=0;
  array *y=NULL;
  if(x==NULL){ return y; }
  if(ARRAY_TYPE(x)=='s'){ flag=svec_has(ARRAY_SIZE(x),ARRAY_P0_S(x),"iIjJ"); }
  if(ARRAY_TYPE(x)=='s' && !flag){ y=array_allocate('d',ARRAY_NDIM(x),ARRAY_DIM_P(x));  dvec_set_svec (ARRAY_SIZE(y),ARRAY_P0_D(y),ARRAY_P0_S(x)); }
  if(ARRAY_TYPE(x)=='s' &&  flag){ y=array_allocate('z',ARRAY_NDIM(x),ARRAY_DIM_P(x));  zvec_set_svec (ARRAY_SIZE(y),ARRAY_P0_Z(y),ARRAY_P0_S(x)); }
  if(ARRAY_TYPE(x)=='i')         { y=array_allocate('d',ARRAY_NDIM(x),ARRAY_DIM_P(x));  dvec_set_ivec(ARRAY_SIZE(y),ARRAY_P0_D(y),ARRAY_P0_I(x)); }
  if(ARRAY_TYPE(x)=='d')         { y=array_allocate('d',ARRAY_NDIM(x),ARRAY_DIM_P(x));  dvec_set_dvec (ARRAY_SIZE(y),ARRAY_P0_D(y),ARRAY_P0_D(x)); }
  if(ARRAY_TYPE(x)=='z')         { y=array_allocate('z',ARRAY_NDIM(x),ARRAY_DIM_P(x));  zvec_set_zvec (ARRAY_SIZE(y),ARRAY_P0_Z(y),ARRAY_P0_Z(x)); }
  if(ARRAY_TYPE(x)=='r')         { y=array_allocate('d',ARRAY_NDIM(x),ARRAY_DIM_P(x));  rvec_get_dvec (ARRAY_SIZE(y),ARRAY_P0_D(y),ARRAY_P0_R(x)); }
  if(ARRAY_TYPE(x)=='c')         { y=array_allocate('z',ARRAY_NDIM(x),ARRAY_DIM_P(x));  cvec_get_zvec (ARRAY_SIZE(y),ARRAY_P0_Z(y),ARRAY_P0_C(x)); }
  if(ARRAY_TYPE(x)=='R')         { y=array_allocate('d',ARRAY_NDIM(x),ARRAY_DIM_P(x)); irvec_get_dvec (ARRAY_SIZE(y),ARRAY_P0_D(y),ARRAY_P0_R(x),ARRAY_P1_R(x)); }
  if(ARRAY_TYPE(x)=='C')         { y=array_allocate('z',ARRAY_NDIM(x),ARRAY_DIM_P(x)); icvec_get_zvec (ARRAY_SIZE(y),ARRAY_P0_Z(y),ARRAY_P0_C(x),ARRAY_P1_C(x)); }
  return y;
}

/**
 @brief array型の型変換 y=multi(x)
 */
array *array_get_multi(array *x)
{
  int flag=0;
  array *y=NULL;
  if(x==NULL){ return y; }
  if(ARRAY_TYPE(x)=='s'){ flag=svec_has(ARRAY_SIZE(x),ARRAY_P0_S(x),"iIjJ"); }
  if(ARRAY_TYPE(x)=='s' &&  flag){ y=array_allocate('c',ARRAY_NDIM(x),ARRAY_DIM_P(x)); cvec_set_svec (ARRAY_SIZE(y),ARRAY_P0_C(y),ARRAY_P0_S(x)); }
  if(ARRAY_TYPE(x)=='s' && !flag){ y=array_allocate('r',ARRAY_NDIM(x),ARRAY_DIM_P(x)); rvec_set_svec (ARRAY_SIZE(y),ARRAY_P0_R(y),ARRAY_P0_S(x)); }
  if(ARRAY_TYPE(x)=='i')         { y=array_allocate('r',ARRAY_NDIM(x),ARRAY_DIM_P(x)); rvec_set_ivec(ARRAY_SIZE(y),ARRAY_P0_R(y),ARRAY_P0_I(x)); }
  if(ARRAY_TYPE(x)=='d')         { y=array_allocate('r',ARRAY_NDIM(x),ARRAY_DIM_P(x)); rvec_set_dvec (ARRAY_SIZE(y),ARRAY_P0_R(y),ARRAY_P0_D(x)); }
  if(ARRAY_TYPE(x)=='z')         { y=array_allocate('c',ARRAY_NDIM(x),ARRAY_DIM_P(x)); cvec_set_zvec (ARRAY_SIZE(y),ARRAY_P0_C(y),ARRAY_P0_Z(x)); }
  if(ARRAY_TYPE(x)=='r')         { y=array_allocate('r',ARRAY_NDIM(x),ARRAY_DIM_P(x)); rvec_set_rvec (ARRAY_SIZE(y),ARRAY_P0_R(y),ARRAY_P0_R(x)); }
  if(ARRAY_TYPE(x)=='c')         { y=array_allocate('c',ARRAY_NDIM(x),ARRAY_DIM_P(x)); cvec_set_cvec (ARRAY_SIZE(y),ARRAY_P0_C(y),ARRAY_P0_C(x)); }
  if(ARRAY_TYPE(x)=='R')         { y=array_allocate('r',ARRAY_NDIM(x),ARRAY_DIM_P(x)); rvec_set_irvec(ARRAY_SIZE(y),ARRAY_P0_R(y),ARRAY_P0_R(x),ARRAY_P1_R(x)); }
  if(ARRAY_TYPE(x)=='C')         { y=array_allocate('c',ARRAY_NDIM(x),ARRAY_DIM_P(x)); cvec_set_icvec(ARRAY_SIZE(y),ARRAY_P0_C(y),ARRAY_P0_C(x),ARRAY_P1_C(x)); }
  return y;
}

/**
 @brief array型の型変換 y=imulti(x)
 */
array *array_get_imulti(array *x)
{
  int flag=0;
  array *y=NULL;
  if(x==NULL){ return y; }
  if(ARRAY_TYPE(x)=='s'){ flag=svec_has(ARRAY_SIZE(x),ARRAY_P0_S(x),"iIjJ"); }
  if(ARRAY_TYPE(x)=='s' &&  flag){ y=array_allocate('C',ARRAY_NDIM(x),ARRAY_DIM_P(x)); icvec_set_svec (ARRAY_SIZE(y),ARRAY_P0_C(y),ARRAY_P1_C(y),ARRAY_P0_S(x)); }
  if(ARRAY_TYPE(x)=='s' && !flag){ y=array_allocate('R',ARRAY_NDIM(x),ARRAY_DIM_P(x)); irvec_set_svec (ARRAY_SIZE(y),ARRAY_P0_R(y),ARRAY_P1_R(y),ARRAY_P0_S(x)); }
  if(ARRAY_TYPE(x)=='i')         { y=array_allocate('R',ARRAY_NDIM(x),ARRAY_DIM_P(x)); irvec_set_ivec(ARRAY_SIZE(y),ARRAY_P0_R(y),ARRAY_P1_R(y),ARRAY_P0_I(x),ARRAY_P0_I(x)); }
  if(ARRAY_TYPE(x)=='d')         { y=array_allocate('R',ARRAY_NDIM(x),ARRAY_DIM_P(x)); irvec_set_dvec (ARRAY_SIZE(y),ARRAY_P0_R(y),ARRAY_P1_R(y),ARRAY_P0_D(x),ARRAY_P0_D(x)); }
  if(ARRAY_TYPE(x)=='z')         { y=array_allocate('C',ARRAY_NDIM(x),ARRAY_DIM_P(x)); icvec_set_zvec (ARRAY_SIZE(y),ARRAY_P0_C(y),ARRAY_P1_C(y),ARRAY_P0_Z(x),ARRAY_P0_Z(x)); }
  if(ARRAY_TYPE(x)=='r')         { y=array_allocate('R',ARRAY_NDIM(x),ARRAY_DIM_P(x)); irvec_set_rvec (ARRAY_SIZE(y),ARRAY_P0_R(y),ARRAY_P1_R(y),ARRAY_P0_R(x),ARRAY_P0_R(x)); }
  if(ARRAY_TYPE(x)=='c')         { y=array_allocate('C',ARRAY_NDIM(x),ARRAY_DIM_P(x)); icvec_set_cvec (ARRAY_SIZE(y),ARRAY_P0_C(y),ARRAY_P1_C(y),ARRAY_P0_C(x),ARRAY_P0_C(x)); }
  if(ARRAY_TYPE(x)=='R')         { y=array_allocate('R',ARRAY_NDIM(x),ARRAY_DIM_P(x)); irvec_set_rvec (ARRAY_SIZE(y),ARRAY_P0_R(y),ARRAY_P1_R(y),ARRAY_P0_R(x),ARRAY_P1_R(x)); }
  if(ARRAY_TYPE(x)=='C')         { y=array_allocate('C',ARRAY_NDIM(x),ARRAY_DIM_P(x)); icvec_set_cvec (ARRAY_SIZE(y),ARRAY_P0_C(y),ARRAY_P1_C(y),ARRAY_P0_C(x),ARRAY_P1_C(x)); }
  return y;
}

/**
 @brief array型の型変換 y=imulti2(x0,x1)
 */
array *array_get_imulti2(array *x0, array *x1)
{
  array *y=NULL;
  if(x0==NULL || x1==NULL){ return y; }
  if(!array_same_dim_check(x0,x1)){ return y; }
  if(ARRAY_TYPE(x0)=='d' && ARRAY_TYPE(x1)=='d'){ y=array_allocate('R',ARRAY_NDIM(x0),ARRAY_DIM_P(x0)); irvec_set_dvec (ARRAY_SIZE(y),ARRAY_P0_R(y),ARRAY_P1_R(y),ARRAY_P0_D(x0),ARRAY_P0_D(x1)); }
  if(ARRAY_TYPE(x0)=='d' && ARRAY_TYPE(x1)=='z'){ y=array_allocate('C',ARRAY_NDIM(x0),ARRAY_DIM_P(x0)); icvec_set_dvec_zvec(ARRAY_SIZE(y),ARRAY_P0_C(y),ARRAY_P1_C(y),ARRAY_P0_D(x0),ARRAY_P0_Z(x1)); }
  if(ARRAY_TYPE(x0)=='z' && ARRAY_TYPE(x1)=='d'){ y=array_allocate('C',ARRAY_NDIM(x0),ARRAY_DIM_P(x0)); icvec_set_zvec_dvec(ARRAY_SIZE(y),ARRAY_P0_C(y),ARRAY_P1_C(y),ARRAY_P0_Z(x0),ARRAY_P0_D(x1)); }
  if(ARRAY_TYPE(x0)=='z' && ARRAY_TYPE(x1)=='z'){ y=array_allocate('C',ARRAY_NDIM(x0),ARRAY_DIM_P(x0)); icvec_set_zvec (ARRAY_SIZE(y),ARRAY_P0_C(y),ARRAY_P1_C(y),ARRAY_P0_Z(x0),ARRAY_P0_Z(x1)); }
  if(ARRAY_TYPE(x0)=='r' && ARRAY_TYPE(x1)=='r'){ y=array_allocate('R',ARRAY_NDIM(x0),ARRAY_DIM_P(x0)); irvec_set_rvec (ARRAY_SIZE(y),ARRAY_P0_R(y),ARRAY_P1_R(y),ARRAY_P0_R(x0),ARRAY_P0_R(x1)); }
  if(ARRAY_TYPE(x0)=='r' && ARRAY_TYPE(x1)=='c'){ y=array_allocate('C',ARRAY_NDIM(x0),ARRAY_DIM_P(x0)); icvec_set_rvec_cvec(ARRAY_SIZE(y),ARRAY_P0_C(y),ARRAY_P1_C(y),ARRAY_P0_R(x0),ARRAY_P0_C(x1)); }
  if(ARRAY_TYPE(x0)=='c' && ARRAY_TYPE(x1)=='r'){ y=array_allocate('C',ARRAY_NDIM(x0),ARRAY_DIM_P(x0)); icvec_set_cvec_rvec(ARRAY_SIZE(y),ARRAY_P0_C(y),ARRAY_P1_C(y),ARRAY_P0_C(x0),ARRAY_P0_R(x1)); }
  if(ARRAY_TYPE(x0)=='c' && ARRAY_TYPE(x1)=='c'){ y=array_allocate('C',ARRAY_NDIM(x0),ARRAY_DIM_P(x0)); icvec_set_cvec (ARRAY_SIZE(y),ARRAY_P0_C(y),ARRAY_P1_C(y),ARRAY_P0_C(x0),ARRAY_P0_C(x1)); }
  return y;
}

/**
 @brief array型の型変換 y=complex2(xr,xi)
 */
array *array_get_complex2(array *xr, array *xi)
{
  array *y=NULL;
  if(xr==NULL || xi==NULL){ return y; }
  if(!array_same_dim_check(xr,xi)){ return y; }
  if(ARRAY_TYPE(xr)=='d' && ARRAY_TYPE(xi)=='d'){ y=array_allocate('z',ARRAY_NDIM(xr),ARRAY_DIM_P(xr));  zvec_set_dvec_dvec(ARRAY_SIZE(y),ARRAY_P0_Z(y),ARRAY_P0_D(xr),ARRAY_P0_D(xi)); }
  if(ARRAY_TYPE(xr)=='r' && ARRAY_TYPE(xi)=='r'){ y=array_allocate('c',ARRAY_NDIM(xr),ARRAY_DIM_P(xr));  cvec_set_rvec_rvec(ARRAY_SIZE(y),ARRAY_P0_C(y),ARRAY_P0_R(xr),ARRAY_P0_R(xi)); }
  if(ARRAY_TYPE(xr)=='R' && ARRAY_TYPE(xi)=='R'){ y=array_allocate('C',ARRAY_NDIM(xr),ARRAY_DIM_P(xr)); icvec_set_rvec_rvec(ARRAY_SIZE(y),ARRAY_P0_C(y),ARRAY_P1_C(y),ARRAY_P0_R(xr),ARRAY_P0_R(xi),ARRAY_P1_R(xr),ARRAY_P1_R(xi)); }
  return y;
}

/** @} */

///////////////////////////////////////////////

/** @name array型の１入力演算に関する関数 */
/** @{ */


/**
 @brief y=conj(x)
 */
array *array_conj(array *x)
{
  array *y=NULL;
  if(x==NULL){ return y; }
  if(ARRAY_TYPE(x)=='d'){ y=array_allocate('d',ARRAY_NDIM(x),ARRAY_DIM_P(x));  dvec_copy_dvec(ARRAY_SIZE(y),ARRAY_P0_D(y),              ARRAY_P0_D(x)); }
  if(ARRAY_TYPE(x)=='z'){ y=array_allocate('z',ARRAY_NDIM(x),ARRAY_DIM_P(x));  zvec_conj_zvec(ARRAY_SIZE(y),ARRAY_P0_Z(y),              ARRAY_P0_Z(x)); }
  if(ARRAY_TYPE(x)=='r'){ y=array_allocate('r',ARRAY_NDIM(x),ARRAY_DIM_P(x));  rvec_copy_rvec(ARRAY_SIZE(y),ARRAY_P0_R(y),              ARRAY_P0_R(x)); }
  if(ARRAY_TYPE(x)=='c'){ y=array_allocate('c',ARRAY_NDIM(x),ARRAY_DIM_P(x));  cvec_conj_cvec(ARRAY_SIZE(y),ARRAY_P0_C(y),              ARRAY_P0_C(x)); }
  if(ARRAY_TYPE(x)=='R'){ y=array_allocate('R',ARRAY_NDIM(x),ARRAY_DIM_P(x)); irvec_copy_rvec(ARRAY_SIZE(y),ARRAY_P0_R(y),ARRAY_P1_R(y),ARRAY_P0_R(x),ARRAY_P1_R(x)); }
  if(ARRAY_TYPE(x)=='C'){ y=array_allocate('C',ARRAY_NDIM(x),ARRAY_DIM_P(x)); icvec_conj_cvec(ARRAY_SIZE(y),ARRAY_P0_C(y),ARRAY_P1_C(y),ARRAY_P0_C(x),ARRAY_P1_C(x)); }
  return y;
}


/** @} */

///////////////////////////////////////////////

/** @name array型の２入力演算に関する関数 */
/** @{ */

/**
 @brief z=x+y
 */
array *array_add(array *x, array *y)
{
  int i,j,k,l,ndim,*dim=NULL,*I=NULL,*K=NULL,*L=NULL;
  array *z=NULL;
  if(x==NULL || y==NULL){ return z; }
  // サイズが同じ場合
  if(z==NULL && array_same_dim_check(x,y)){
    z=array_allocate(array_get_inclusion_type(ARRAY_TYPE(x),ARRAY_TYPE(y)),ARRAY_NDIM(x),ARRAY_DIM_P(x));
    if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='d'){  dvec_add_dvec(ARRAY_SIZE(z),ARRAY_P0_D(z),ARRAY_P0_D(x),ARRAY_P0_D(y)); }
    if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='z'){  dvec_add_zvec(ARRAY_SIZE(z),ARRAY_P0_Z(z),ARRAY_P0_D(x),ARRAY_P0_Z(y)); }
    if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='d'){  zvec_add_dvec(ARRAY_SIZE(z),ARRAY_P0_Z(z),ARRAY_P0_Z(x),ARRAY_P0_D(y)); }
    if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='z'){  zvec_add_zvec(ARRAY_SIZE(z),ARRAY_P0_Z(z),ARRAY_P0_Z(x),ARRAY_P0_Z(y)); }
    if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='r'){  dvec_add_rvec(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P0_D(x),ARRAY_P0_R(y)); }
    if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='r'){  zvec_add_rvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_Z(x),ARRAY_P0_R(y)); }
    if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='d'){  rvec_add_dvec(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P0_R(x),ARRAY_P0_D(y)); }
    if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='z'){  rvec_add_zvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_R(x),ARRAY_P0_Z(y)); }
    if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='r'){  rvec_add_rvec(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P0_R(x),ARRAY_P0_R(y)); }
    if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='c'){  dvec_add_cvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_D(x),ARRAY_P0_C(y)); }
    if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='c'){  zvec_add_cvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_Z(x),ARRAY_P0_C(y)); }
    if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='c'){  rvec_add_cvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_R(x),ARRAY_P0_C(y)); }
    if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='d'){  cvec_add_dvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_C(x),ARRAY_P0_D(y)); }
    if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='z'){  cvec_add_zvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_C(x),ARRAY_P0_Z(y)); }
    if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='r'){  cvec_add_rvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_C(x),ARRAY_P0_R(y)); }
    if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='c'){  cvec_add_cvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_C(x),ARRAY_P0_C(y)); }
    if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='R'){ idvec_add_rvec(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P1_R(z),ARRAY_P0_D(x),ARRAY_P0_D(x),ARRAY_P0_R(y),ARRAY_P1_R(y)); }
    if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='R'){ izvec_add_rvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_Z(x),ARRAY_P0_Z(x),ARRAY_P0_R(y),ARRAY_P1_R(y)); }
    if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='R'){ irvec_add_rvec(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P1_R(z),ARRAY_P0_R(x),ARRAY_P0_R(x),ARRAY_P0_R(y),ARRAY_P1_R(y)); }
    if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='R'){ icvec_add_rvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(x),ARRAY_P0_C(x),ARRAY_P0_R(y),ARRAY_P1_R(y)); }
    if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='d'){ irvec_add_dvec(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P1_R(z),ARRAY_P0_R(x),ARRAY_P1_R(x),ARRAY_P0_D(y),ARRAY_P0_D(y)); }
    if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='z'){ irvec_add_zvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_R(x),ARRAY_P1_R(x),ARRAY_P0_Z(y),ARRAY_P0_Z(y)); }
    if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='r'){ irvec_add_rvec(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P1_R(z),ARRAY_P0_R(x),ARRAY_P1_R(x),ARRAY_P0_R(y),ARRAY_P0_R(y)); }
    if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='c'){ irvec_add_cvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_R(x),ARRAY_P1_R(x),ARRAY_P0_C(y),ARRAY_P0_C(y)); }
    if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='R'){ irvec_add_rvec(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P1_R(z),ARRAY_P0_R(x),ARRAY_P1_R(x),ARRAY_P0_R(y),ARRAY_P1_R(y)); }
    if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='C'){ idvec_add_cvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_D(x),ARRAY_P0_D(x),ARRAY_P0_C(y),ARRAY_P1_C(y)); }
    if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='C'){ izvec_add_cvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_Z(x),ARRAY_P0_Z(x),ARRAY_P0_C(y),ARRAY_P1_C(y)); }
    if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='C'){ irvec_add_cvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_R(x),ARRAY_P0_R(x),ARRAY_P0_C(y),ARRAY_P1_C(y)); }
    if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='C'){ icvec_add_cvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(x),ARRAY_P0_C(x),ARRAY_P0_C(y),ARRAY_P1_C(y)); }
    if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='C'){ irvec_add_cvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_R(x),ARRAY_P1_R(x),ARRAY_P0_C(y),ARRAY_P1_C(y)); }
    if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='d'){ icvec_add_dvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(x),ARRAY_P1_C(x),ARRAY_P0_D(y),ARRAY_P0_D(y)); }
    if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='z'){ icvec_add_zvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(x),ARRAY_P1_C(x),ARRAY_P0_Z(y),ARRAY_P0_Z(y)); }
    if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='r'){ icvec_add_rvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(x),ARRAY_P1_C(x),ARRAY_P0_R(y),ARRAY_P0_R(y)); }
    if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='c'){ icvec_add_cvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(x),ARRAY_P1_C(x),ARRAY_P0_C(y),ARRAY_P0_C(y)); }
    if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='R'){ icvec_add_rvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(x),ARRAY_P1_C(x),ARRAY_P0_R(y),ARRAY_P1_R(y)); }
    if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='C'){ icvec_add_cvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(x),ARRAY_P1_C(x),ARRAY_P0_C(y),ARRAY_P1_C(y)); }
  }
  // yがスカラーの場合
  if(z==NULL && array_is_scalar(y)){
    z=array_allocate(array_get_inclusion_type(ARRAY_TYPE(x),ARRAY_TYPE(y)),ARRAY_NDIM(x),ARRAY_DIM_P(x));
    if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='d'){  dvec_add_dscalar(ARRAY_SIZE(z),ARRAY_P0_D(z),ARRAY_P0_D(x),ARRAY_P0_D(y)[0]); }
    if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='z'){  dvec_add_zscalar(ARRAY_SIZE(z),ARRAY_P0_Z(z),ARRAY_P0_D(x),ARRAY_P0_Z(y)[0]); }
    if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='d'){  zvec_add_dscalar(ARRAY_SIZE(z),ARRAY_P0_Z(z),ARRAY_P0_Z(x),ARRAY_P0_D(y)[0]); }
    if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='z'){  zvec_add_zscalar(ARRAY_SIZE(z),ARRAY_P0_Z(z),ARRAY_P0_Z(x),ARRAY_P0_Z(y)[0]); }
    if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='r'){  dvec_add_rscalar(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P0_D(x),ARRAY_P0_R(y)[0]); }
    if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='r'){  zvec_add_rscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_Z(x),ARRAY_P0_R(y)[0]); }
    if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='d'){  rvec_add_dscalar(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P0_R(x),ARRAY_P0_D(y)[0]); }
    if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='z'){  rvec_add_zscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_R(x),ARRAY_P0_Z(y)[0]); }
    if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='r'){  rvec_add_rscalar(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P0_R(x),ARRAY_P0_R(y)[0]); }
    if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='c'){  dvec_add_cscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_D(x),ARRAY_P0_C(y)[0]); }
    if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='c'){  zvec_add_cscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_Z(x),ARRAY_P0_C(y)[0]); }
    if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='c'){  rvec_add_cscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_R(x),ARRAY_P0_C(y)[0]); }
    if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='d'){  cvec_add_dscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_C(x),ARRAY_P0_D(y)[0]); }
    if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='z'){  cvec_add_zscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_C(x),ARRAY_P0_Z(y)[0]); }
    if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='r'){  cvec_add_rscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_C(x),ARRAY_P0_R(y)[0]); }
    if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='c'){  cvec_add_cscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_C(x),ARRAY_P0_C(y)[0]); }
    if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='R'){ idvec_add_rscalar(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P1_R(z),ARRAY_P0_D(x),ARRAY_P0_D(x),ARRAY_P0_R(y)[0],ARRAY_P1_R(y)[0]); }
    if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='R'){ izvec_add_rscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_Z(x),ARRAY_P0_Z(x),ARRAY_P0_R(y)[0],ARRAY_P1_R(y)[0]); }
    if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='R'){ irvec_add_rscalar(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P1_R(z),ARRAY_P0_R(x),ARRAY_P0_R(x),ARRAY_P0_R(y)[0],ARRAY_P1_R(y)[0]); }
    if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='R'){ icvec_add_rscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(x),ARRAY_P0_C(x),ARRAY_P0_R(y)[0],ARRAY_P1_R(y)[0]); }
    if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='d'){ irvec_add_dscalar(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P1_R(z),ARRAY_P0_R(x),ARRAY_P1_R(x),ARRAY_P0_D(y)[0],ARRAY_P0_D(y)[0]); }
    if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='z'){ irvec_add_zscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_R(x),ARRAY_P1_R(x),ARRAY_P0_Z(y)[0],ARRAY_P0_Z(y)[0]); }
    if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='r'){ irvec_add_rscalar(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P1_R(z),ARRAY_P0_R(x),ARRAY_P1_R(x),ARRAY_P0_R(y)[0],ARRAY_P0_R(y)[0]); }
    if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='c'){ irvec_add_cscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_R(x),ARRAY_P1_R(x),ARRAY_P0_C(y)[0],ARRAY_P0_C(y)[0]); }
    if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='R'){ irvec_add_rscalar(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P1_R(z),ARRAY_P0_R(x),ARRAY_P1_R(x),ARRAY_P0_R(y)[0],ARRAY_P1_R(y)[0]); }
    if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='C'){ idvec_add_cscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_D(x),ARRAY_P0_D(x),ARRAY_P0_C(y)[0],ARRAY_P1_C(y)[0]); }
    if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='C'){ izvec_add_cscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_Z(x),ARRAY_P0_Z(x),ARRAY_P0_C(y)[0],ARRAY_P1_C(y)[0]); }
    if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='C'){ irvec_add_cscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_R(x),ARRAY_P0_R(x),ARRAY_P0_C(y)[0],ARRAY_P1_C(y)[0]); }
    if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='C'){ icvec_add_cscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(x),ARRAY_P0_C(x),ARRAY_P0_C(y)[0],ARRAY_P1_C(y)[0]); }
    if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='C'){ irvec_add_cscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_R(x),ARRAY_P1_R(x),ARRAY_P0_C(y)[0],ARRAY_P1_C(y)[0]); }
    if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='d'){ icvec_add_dscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(x),ARRAY_P1_C(x),ARRAY_P0_D(y)[0],ARRAY_P0_D(y)[0]); }
    if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='z'){ icvec_add_zscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(x),ARRAY_P1_C(x),ARRAY_P0_Z(y)[0],ARRAY_P0_Z(y)[0]); }
    if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='r'){ icvec_add_rscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(x),ARRAY_P1_C(x),ARRAY_P0_R(y)[0],ARRAY_P0_R(y)[0]); }
    if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='c'){ icvec_add_cscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(x),ARRAY_P1_C(x),ARRAY_P0_C(y)[0],ARRAY_P0_C(y)[0]); }
    if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='R'){ icvec_add_rscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(x),ARRAY_P1_C(x),ARRAY_P0_R(y)[0],ARRAY_P1_R(y)[0]); }
    if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='C'){ icvec_add_cscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(x),ARRAY_P1_C(x),ARRAY_P0_C(y)[0],ARRAY_P1_C(y)[0]); }
  }
  // xがスカラーの場合
  if(z==NULL && array_is_scalar(x)){
    z=array_allocate(array_get_inclusion_type(ARRAY_TYPE(x),ARRAY_TYPE(y)),ARRAY_NDIM(y),ARRAY_DIM_P(y));
    if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='d'){  dvec_add_dscalar(ARRAY_SIZE(z),ARRAY_P0_D(z),ARRAY_P0_D(y),ARRAY_P0_D(x)[0]); }
    if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='z'){  zvec_add_dscalar(ARRAY_SIZE(z),ARRAY_P0_Z(z),ARRAY_P0_Z(y),ARRAY_P0_D(x)[0]); }
    if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='d'){  dvec_add_zscalar(ARRAY_SIZE(z),ARRAY_P0_Z(z),ARRAY_P0_D(y),ARRAY_P0_Z(x)[0]); }
    if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='z'){  zvec_add_zscalar(ARRAY_SIZE(z),ARRAY_P0_Z(z),ARRAY_P0_Z(y),ARRAY_P0_Z(x)[0]); }
    if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='r'){  rvec_add_dscalar(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P0_R(y),ARRAY_P0_D(x)[0]); }
    if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='r'){  rvec_add_zscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_R(y),ARRAY_P0_Z(x)[0]); }
    if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='d'){  dvec_add_rscalar(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P0_D(y),ARRAY_P0_R(x)[0]); }
    if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='z'){  zvec_add_rscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_Z(y),ARRAY_P0_R(x)[0]); }
    if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='r'){  rvec_add_rscalar(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P0_R(y),ARRAY_P0_R(x)[0]); }
    if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='c'){  cvec_add_dscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_C(y),ARRAY_P0_D(x)[0]); }
    if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='c'){  cvec_add_zscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_C(y),ARRAY_P0_Z(x)[0]); }
    if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='c'){  cvec_add_rscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_C(y),ARRAY_P0_R(x)[0]); }
    if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='d'){  dvec_add_cscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_D(y),ARRAY_P0_C(x)[0]); }
    if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='z'){  zvec_add_cscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_Z(y),ARRAY_P0_C(x)[0]); }
    if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='r'){  rvec_add_cscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_R(y),ARRAY_P0_C(x)[0]); }
    if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='c'){  cvec_add_cscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_C(y),ARRAY_P0_C(x)[0]); }
    if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='R'){ irvec_add_dscalar(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P1_R(z),ARRAY_P0_R(y),ARRAY_P1_R(y),ARRAY_P0_D(x)[0],ARRAY_P0_D(x)[0]); }
    if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='R'){ irvec_add_zscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_R(y),ARRAY_P1_R(y),ARRAY_P0_Z(x)[0],ARRAY_P0_Z(x)[0]); }
    if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='R'){ irvec_add_rscalar(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P1_R(z),ARRAY_P0_R(y),ARRAY_P1_R(y),ARRAY_P0_R(x)[0],ARRAY_P0_R(x)[0]); }
    if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='R'){ irvec_add_cscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_R(y),ARRAY_P1_R(y),ARRAY_P0_C(x)[0],ARRAY_P0_C(x)[0]); }
    if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='d'){ idvec_add_rscalar(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P1_R(z),ARRAY_P0_D(y),ARRAY_P0_D(y),ARRAY_P0_R(x)[0],ARRAY_P1_R(x)[0]); }
    if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='z'){ izvec_add_rscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_Z(y),ARRAY_P0_Z(y),ARRAY_P0_R(x)[0],ARRAY_P1_R(x)[0]); }
    if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='r'){ irvec_add_rscalar(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P1_R(z),ARRAY_P0_R(y),ARRAY_P0_R(y),ARRAY_P0_R(x)[0],ARRAY_P1_R(x)[0]); }
    if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='c'){ icvec_add_rscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(y),ARRAY_P0_C(y),ARRAY_P0_R(x)[0],ARRAY_P1_R(x)[0]); }
    if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='R'){ irvec_add_rscalar(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P1_R(z),ARRAY_P0_R(y),ARRAY_P1_R(y),ARRAY_P0_R(x)[0],ARRAY_P1_R(x)[0]); }
    if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='C'){ icvec_add_dscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(y),ARRAY_P1_C(y),ARRAY_P0_D(x)[0],ARRAY_P0_D(x)[0]); }
    if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='C'){ icvec_add_zscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(y),ARRAY_P1_C(y),ARRAY_P0_Z(x)[0],ARRAY_P0_Z(x)[0]); }
    if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='C'){ icvec_add_rscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(y),ARRAY_P1_C(y),ARRAY_P0_R(x)[0],ARRAY_P0_R(x)[0]); }
    if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='C'){ icvec_add_cscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(y),ARRAY_P1_C(y),ARRAY_P0_C(x)[0],ARRAY_P0_C(x)[0]); }
    if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='C'){ icvec_add_rscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(y),ARRAY_P1_C(y),ARRAY_P0_R(x)[0],ARRAY_P1_R(x)[0]); }
    if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='d'){ idvec_add_cscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_D(y),ARRAY_P0_D(y),ARRAY_P0_C(x)[0],ARRAY_P1_C(x)[0]); }
    if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='z'){ izvec_add_cscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_Z(y),ARRAY_P0_Z(y),ARRAY_P0_C(x)[0],ARRAY_P1_C(x)[0]); }
    if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='r'){ irvec_add_cscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_R(y),ARRAY_P0_R(y),ARRAY_P0_C(x)[0],ARRAY_P1_C(x)[0]); }
    if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='c'){ icvec_add_cscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(y),ARRAY_P0_C(y),ARRAY_P0_C(x)[0],ARRAY_P1_C(x)[0]); }
    if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='R'){ irvec_add_cscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_R(y),ARRAY_P1_R(y),ARRAY_P0_C(x)[0],ARRAY_P1_C(x)[0]); }
    if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='C'){ icvec_add_cscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(y),ARRAY_P1_C(y),ARRAY_P0_C(x)[0],ARRAY_P1_C(x)[0]); }
  }
  // 互換性のあるサイズでありxが小さい場合
  if(z==NULL && (k=array_get_subdim(x,y))>0){
    z=array_allocate(array_get_inclusion_type(ARRAY_TYPE(x),ARRAY_TYPE(y)),ARRAY_NDIM(y),ARRAY_DIM_P(y));
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
      if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='R'){ idvec_add_rvec(ARRAY_LD(z,k-1),ARRAY_P0_R(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_R(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_D(x),ARRAY_P0_D(x),ARRAY_P0_R(y)+i*ARRAY_LD(y,k-1),ARRAY_P1_R(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='R'){ izvec_add_rvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_Z(x),ARRAY_P0_Z(x),ARRAY_P0_R(y)+i*ARRAY_LD(y,k-1),ARRAY_P1_R(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='R'){ irvec_add_rvec(ARRAY_LD(z,k-1),ARRAY_P0_R(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_R(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_R(x),ARRAY_P0_R(x),ARRAY_P0_R(y)+i*ARRAY_LD(y,k-1),ARRAY_P1_R(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='R'){ icvec_add_rvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_C(x),ARRAY_P0_C(x),ARRAY_P0_R(y)+i*ARRAY_LD(y,k-1),ARRAY_P1_R(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='d'){ irvec_add_dvec(ARRAY_LD(z,k-1),ARRAY_P0_R(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_R(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_R(x),ARRAY_P1_R(x),ARRAY_P0_D(y)+i*ARRAY_LD(y,k-1),ARRAY_P0_D(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='z'){ irvec_add_zvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_R(x),ARRAY_P1_R(x),ARRAY_P0_Z(y)+i*ARRAY_LD(y,k-1),ARRAY_P0_Z(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='r'){ irvec_add_rvec(ARRAY_LD(z,k-1),ARRAY_P0_R(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_R(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_R(x),ARRAY_P1_R(x),ARRAY_P0_R(y)+i*ARRAY_LD(y,k-1),ARRAY_P0_R(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='c'){ irvec_add_cvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_R(x),ARRAY_P1_R(x),ARRAY_P0_C(y)+i*ARRAY_LD(y,k-1),ARRAY_P0_C(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='R'){ irvec_add_rvec(ARRAY_LD(z,k-1),ARRAY_P0_R(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_R(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_R(x),ARRAY_P1_R(x),ARRAY_P0_R(y)+i*ARRAY_LD(y,k-1),ARRAY_P1_R(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='C'){ idvec_add_cvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_D(x),ARRAY_P0_D(x),ARRAY_P0_C(y)+i*ARRAY_LD(y,k-1),ARRAY_P1_C(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='C'){ izvec_add_cvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_Z(x),ARRAY_P0_Z(x),ARRAY_P0_C(y)+i*ARRAY_LD(y,k-1),ARRAY_P1_C(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='C'){ irvec_add_cvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_R(x),ARRAY_P0_R(x),ARRAY_P0_C(y)+i*ARRAY_LD(y,k-1),ARRAY_P1_C(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='C'){ icvec_add_cvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_C(x),ARRAY_P0_C(x),ARRAY_P0_C(y)+i*ARRAY_LD(y,k-1),ARRAY_P1_C(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='C'){ irvec_add_cvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_R(x),ARRAY_P1_R(x),ARRAY_P0_C(y)+i*ARRAY_LD(y,k-1),ARRAY_P1_C(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='d'){ icvec_add_dvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_C(x),ARRAY_P1_C(x),ARRAY_P0_D(y)+i*ARRAY_LD(y,k-1),ARRAY_P0_D(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='z'){ icvec_add_zvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_C(x),ARRAY_P1_C(x),ARRAY_P0_Z(y)+i*ARRAY_LD(y,k-1),ARRAY_P0_Z(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='r'){ icvec_add_rvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_C(x),ARRAY_P1_C(x),ARRAY_P0_R(y)+i*ARRAY_LD(y,k-1),ARRAY_P0_R(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='c'){ icvec_add_cvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_C(x),ARRAY_P1_C(x),ARRAY_P0_C(y)+i*ARRAY_LD(y,k-1),ARRAY_P0_C(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='R'){ icvec_add_rvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_C(x),ARRAY_P1_C(x),ARRAY_P0_R(y)+i*ARRAY_LD(y,k-1),ARRAY_P1_R(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='C'){ icvec_add_cvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_C(x),ARRAY_P1_C(x),ARRAY_P0_C(y)+i*ARRAY_LD(y,k-1),ARRAY_P1_C(y)+i*ARRAY_LD(y,k-1)); }
    }
  }
  // 互換性のあるサイズでありyが小さい場合
  if(z==NULL && (k=array_get_subdim(y,x))>0){
    z=array_allocate(array_get_inclusion_type(ARRAY_TYPE(x),ARRAY_TYPE(y)),ARRAY_NDIM(x),ARRAY_DIM_P(x));
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
      if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='R'){ idvec_add_rvec(ARRAY_LD(z,k-1),ARRAY_P0_R(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_R(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_D(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_D(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_R(y),ARRAY_P1_R(y)); }
      if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='R'){ izvec_add_rvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_Z(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_Z(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_R(y),ARRAY_P1_R(y)); }
      if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='R'){ irvec_add_rvec(ARRAY_LD(z,k-1),ARRAY_P0_R(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_R(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_R(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_R(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_R(y),ARRAY_P1_R(y)); }
      if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='R'){ icvec_add_rvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_C(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_C(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_R(y),ARRAY_P1_R(y)); }
      if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='d'){ irvec_add_dvec(ARRAY_LD(z,k-1),ARRAY_P0_R(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_R(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_R(x)+i*ARRAY_LD(x,k-1),ARRAY_P1_R(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_D(y),ARRAY_P0_D(y)); }
      if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='z'){ irvec_add_zvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_R(x)+i*ARRAY_LD(x,k-1),ARRAY_P1_R(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_Z(y),ARRAY_P0_Z(y)); }
      if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='r'){ irvec_add_rvec(ARRAY_LD(z,k-1),ARRAY_P0_R(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_R(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_R(x)+i*ARRAY_LD(x,k-1),ARRAY_P1_R(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_R(y),ARRAY_P0_R(y)); }
      if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='c'){ irvec_add_cvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_R(x)+i*ARRAY_LD(x,k-1),ARRAY_P1_R(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_C(y),ARRAY_P0_C(y)); }
      if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='R'){ irvec_add_rvec(ARRAY_LD(z,k-1),ARRAY_P0_R(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_R(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_R(x)+i*ARRAY_LD(x,k-1),ARRAY_P1_R(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_R(y),ARRAY_P1_R(y)); }
      if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='C'){ idvec_add_cvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_D(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_D(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_C(y),ARRAY_P1_C(y)); }
      if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='C'){ izvec_add_cvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_Z(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_Z(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_C(y),ARRAY_P1_C(y)); }
      if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='C'){ irvec_add_cvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_R(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_R(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_C(y),ARRAY_P1_C(y)); }
      if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='C'){ icvec_add_cvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_C(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_C(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_C(y),ARRAY_P1_C(y)); }
      if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='C'){ irvec_add_cvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_R(x)+i*ARRAY_LD(x,k-1),ARRAY_P1_R(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_C(y),ARRAY_P1_C(y)); }
      if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='d'){ icvec_add_dvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_C(x)+i*ARRAY_LD(x,k-1),ARRAY_P1_C(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_D(y),ARRAY_P0_D(y)); }
      if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='z'){ icvec_add_zvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_C(x)+i*ARRAY_LD(x,k-1),ARRAY_P1_C(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_Z(y),ARRAY_P0_Z(y)); }
      if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='r'){ icvec_add_rvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_C(x)+i*ARRAY_LD(x,k-1),ARRAY_P1_C(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_R(y),ARRAY_P0_R(y)); }
      if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='c'){ icvec_add_cvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_C(x)+i*ARRAY_LD(x,k-1),ARRAY_P1_C(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_C(y),ARRAY_P0_C(y)); }
      if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='R'){ icvec_add_rvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_C(x)+i*ARRAY_LD(x,k-1),ARRAY_P1_C(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_R(y),ARRAY_P1_R(y)); }
      if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='C'){ icvec_add_cvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_C(x)+i*ARRAY_LD(x,k-1),ARRAY_P1_C(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_C(y),ARRAY_P1_C(y)); }
    }
  }
  // 互換性のあるサイズの場合
  if(z==NULL && array_compatible_dim_check(x,y)){
    ndim=MAX2(ARRAY_NDIM(x),ARRAY_NDIM(y));
    dim=array_get_inclusion_dim(x,y);
    z=array_allocate(array_get_inclusion_type(ARRAY_TYPE(x),ARRAY_TYPE(y)),ndim,dim);
    // loop
    I=ivec_allocate(ARRAY_NDIM(z));
    K=ivec_allocate(ARRAY_NDIM(z));
    L=ivec_allocate(ARRAY_NDIM(z));
    for(i=0; i<ARRAY_SIZE(z); i++){
      array_get_index(I,z,i);
      ivec_copy(ndim,K,I);
      ivec_copy(ndim,L,I);
      for(j=0; j<ndim; j++){
        if(K[j]>=ARRAY_DIM(x,j)){ K[j]=ARRAY_DIM(x,j)-1; }
        if(L[j]>=ARRAY_DIM(y,j)){ L[j]=ARRAY_DIM(y,j)-1; }
      }
      k=array_get_position(K,x);
      l=array_get_position(L,y);
      if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='d'){ ARRAY_DVEC0(z,i)=ARRAY_DVEC0(x,k)+ARRAY_DVEC0(y,l); }
      if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='z'){ ARRAY_ZVEC0(z,i)=zadd_dz(ARRAY_DVEC0(x,k),ARRAY_ZVEC0(y,l)); }
      if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='d'){ ARRAY_ZVEC0(z,i)=zadd_zd(ARRAY_ZVEC0(x,k),ARRAY_DVEC0(y,l)); }
      if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='z'){ ARRAY_ZVEC0(z,i)=zadd_zz(ARRAY_ZVEC0(x,k),ARRAY_ZVEC0(y,l)); }
      if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='r'){  radd_dr(ARRAY_RVEC0(z,i),ARRAY_DVEC0(x,k),ARRAY_RVEC0(y,l)); }
      if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='r'){  cadd_zr(ARRAY_CVEC0(z,i),ARRAY_ZVEC0(x,k),ARRAY_RVEC0(y,l)); }
      if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='d'){  radd_rd(ARRAY_RVEC0(z,i),ARRAY_RVEC0(x,k),ARRAY_DVEC0(y,l)); }
      if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='z'){  cadd_rz(ARRAY_CVEC0(z,i),ARRAY_RVEC0(x,k),ARRAY_ZVEC0(y,l)); }
      if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='r'){  radd_rr(ARRAY_RVEC0(z,i),ARRAY_RVEC0(x,k),ARRAY_RVEC0(y,l)); }
      if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='c'){  cadd_dc(ARRAY_CVEC0(z,i),ARRAY_DVEC0(x,k),ARRAY_CVEC0(y,l)); }
      if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='c'){  cadd_zc(ARRAY_CVEC0(z,i),ARRAY_ZVEC0(x,k),ARRAY_CVEC0(y,l)); }
      if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='c'){  cadd_rc(ARRAY_CVEC0(z,i),ARRAY_RVEC0(x,k),ARRAY_CVEC0(y,l)); }
      if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='d'){  cadd_cd(ARRAY_CVEC0(z,i),ARRAY_CVEC0(x,k),ARRAY_DVEC0(y,l)); }
      if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='z'){  cadd_cz(ARRAY_CVEC0(z,i),ARRAY_CVEC0(x,k),ARRAY_ZVEC0(y,l)); }
      if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='r'){  cadd_cr(ARRAY_CVEC0(z,i),ARRAY_CVEC0(x,k),ARRAY_RVEC0(y,l)); }
      if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='c'){  cadd_cc(ARRAY_CVEC0(z,i),ARRAY_CVEC0(x,k),ARRAY_CVEC0(y,l)); }
      if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='R'){ iradd_dr(ARRAY_RVEC0(z,i),ARRAY_RVEC1(z,i),ARRAY_DVEC0(x,k),ARRAY_DVEC0(x,k),ARRAY_RVEC0(y,l),ARRAY_RVEC1(y,l)); }
      if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='R'){ icadd_zr(ARRAY_CVEC0(z,i),ARRAY_CVEC1(z,i),ARRAY_ZVEC0(x,k),ARRAY_ZVEC0(x,k),ARRAY_RVEC0(y,l),ARRAY_RVEC1(y,l)); }
      if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='R'){ iradd_rr(ARRAY_RVEC0(z,i),ARRAY_RVEC1(z,i),ARRAY_RVEC0(x,k),ARRAY_RVEC0(x,k),ARRAY_RVEC0(y,l),ARRAY_RVEC1(y,l)); }
      if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='R'){ icadd_cr(ARRAY_CVEC0(z,i),ARRAY_CVEC1(z,i),ARRAY_CVEC0(x,k),ARRAY_CVEC0(x,k),ARRAY_RVEC0(y,l),ARRAY_RVEC1(y,l)); }
      if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='d'){ iradd_rd(ARRAY_RVEC0(z,i),ARRAY_RVEC1(z,i),ARRAY_RVEC0(x,k),ARRAY_RVEC1(x,k),ARRAY_DVEC0(y,l),ARRAY_DVEC0(y,l)); }
      if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='z'){ icadd_rz(ARRAY_CVEC0(z,i),ARRAY_CVEC1(z,i),ARRAY_RVEC0(x,k),ARRAY_RVEC1(x,k),ARRAY_ZVEC0(y,l),ARRAY_ZVEC0(y,l)); }
      if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='r'){ iradd_rr(ARRAY_RVEC0(z,i),ARRAY_RVEC1(z,i),ARRAY_RVEC0(x,k),ARRAY_RVEC1(x,k),ARRAY_RVEC0(y,l),ARRAY_RVEC0(y,l)); }
      if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='c'){ icadd_rc(ARRAY_CVEC0(z,i),ARRAY_CVEC1(z,i),ARRAY_RVEC0(x,k),ARRAY_RVEC1(x,k),ARRAY_CVEC0(y,l),ARRAY_CVEC0(y,l)); }
      if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='R'){ iradd_rr(ARRAY_RVEC0(z,i),ARRAY_RVEC1(z,i),ARRAY_RVEC0(x,k),ARRAY_RVEC1(x,k),ARRAY_RVEC0(y,l),ARRAY_RVEC1(y,l)); }
      if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='C'){ icadd_dc(ARRAY_CVEC0(z,i),ARRAY_CVEC1(z,i),ARRAY_DVEC0(x,k),ARRAY_DVEC0(x,k),ARRAY_CVEC0(y,l),ARRAY_CVEC1(y,l)); }
      if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='C'){ icadd_zc(ARRAY_CVEC0(z,i),ARRAY_CVEC1(z,i),ARRAY_ZVEC0(x,k),ARRAY_ZVEC0(x,k),ARRAY_CVEC0(y,l),ARRAY_CVEC1(y,l)); }
      if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='C'){ icadd_rc(ARRAY_CVEC0(z,i),ARRAY_CVEC1(z,i),ARRAY_RVEC0(x,k),ARRAY_RVEC0(x,k),ARRAY_CVEC0(y,l),ARRAY_CVEC1(y,l)); }
      if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='C'){ icadd_cc(ARRAY_CVEC0(z,i),ARRAY_CVEC1(z,i),ARRAY_CVEC0(x,k),ARRAY_CVEC0(x,k),ARRAY_CVEC0(y,l),ARRAY_CVEC1(y,l)); }
      if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='C'){ icadd_rc(ARRAY_CVEC0(z,i),ARRAY_CVEC1(z,i),ARRAY_RVEC0(x,k),ARRAY_RVEC1(x,k),ARRAY_CVEC0(y,l),ARRAY_CVEC1(y,l)); }
      if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='d'){ icadd_cd(ARRAY_CVEC0(z,i),ARRAY_CVEC1(z,i),ARRAY_CVEC0(x,k),ARRAY_CVEC1(x,k),ARRAY_DVEC0(y,l),ARRAY_DVEC0(y,l)); }
      if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='z'){ icadd_cz(ARRAY_CVEC0(z,i),ARRAY_CVEC1(z,i),ARRAY_CVEC0(x,k),ARRAY_CVEC1(x,k),ARRAY_ZVEC0(y,l),ARRAY_ZVEC0(y,l)); }
      if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='r'){ icadd_cr(ARRAY_CVEC0(z,i),ARRAY_CVEC1(z,i),ARRAY_CVEC0(x,k),ARRAY_CVEC1(x,k),ARRAY_RVEC0(y,l),ARRAY_RVEC0(y,l)); }
      if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='c'){ icadd_cc(ARRAY_CVEC0(z,i),ARRAY_CVEC1(z,i),ARRAY_CVEC0(x,k),ARRAY_CVEC1(x,k),ARRAY_CVEC0(y,l),ARRAY_CVEC0(y,l)); }
      if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='R'){ icadd_cr(ARRAY_CVEC0(z,i),ARRAY_CVEC1(z,i),ARRAY_CVEC0(x,k),ARRAY_CVEC1(x,k),ARRAY_RVEC0(y,l),ARRAY_RVEC1(y,l)); }
      if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='C'){ icadd_cc(ARRAY_CVEC0(z,i),ARRAY_CVEC1(z,i),ARRAY_CVEC0(x,k),ARRAY_CVEC1(x,k),ARRAY_CVEC0(y,l),ARRAY_CVEC1(y,l)); }
    }
    dim=ivec_free(dim);
    I=ivec_free(I);
    K=ivec_free(K);
    L=ivec_free(L);
  }
  return z;
}

/**
 @brief z=x-y
 */
array *array_sub(array *x, array *y)
{
  int i,j,k,l,ndim,*dim=NULL,*I=NULL,*K=NULL,*L=NULL;
  array *z=NULL;
  if(x==NULL || y==NULL){ return z; }
  // サイズが同じ場合
  if(z==NULL && array_same_dim_check(x,y))
  {
    z=array_allocate(array_get_inclusion_type(ARRAY_TYPE(x),ARRAY_TYPE(y)),ARRAY_NDIM(x),ARRAY_DIM_P(x));
    if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='d'){  dvec_sub_dvec(ARRAY_SIZE(z),ARRAY_P0_D(z),ARRAY_P0_D(x),ARRAY_P0_D(y)); }
    if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='z'){  dvec_sub_zvec(ARRAY_SIZE(z),ARRAY_P0_Z(z),ARRAY_P0_D(x),ARRAY_P0_Z(y)); }
    if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='d'){  zvec_sub_dvec(ARRAY_SIZE(z),ARRAY_P0_Z(z),ARRAY_P0_Z(x),ARRAY_P0_D(y)); }
    if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='z'){  zvec_sub_zvec(ARRAY_SIZE(z),ARRAY_P0_Z(z),ARRAY_P0_Z(x),ARRAY_P0_Z(y)); }
    if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='r'){  dvec_sub_rvec(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P0_D(x),ARRAY_P0_R(y)); }
    if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='r'){  zvec_sub_rvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_Z(x),ARRAY_P0_R(y)); }
    if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='d'){  rvec_sub_dvec(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P0_R(x),ARRAY_P0_D(y)); }
    if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='z'){  rvec_sub_zvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_R(x),ARRAY_P0_Z(y)); }
    if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='r'){  rvec_sub_rvec(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P0_R(x),ARRAY_P0_R(y)); }
    if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='c'){  dvec_sub_cvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_D(x),ARRAY_P0_C(y)); }
    if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='c'){  zvec_sub_cvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_Z(x),ARRAY_P0_C(y)); }
    if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='c'){  rvec_sub_cvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_R(x),ARRAY_P0_C(y)); }
    if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='d'){  cvec_sub_dvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_C(x),ARRAY_P0_D(y)); }
    if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='z'){  cvec_sub_zvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_C(x),ARRAY_P0_Z(y)); }
    if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='r'){  cvec_sub_rvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_C(x),ARRAY_P0_R(y)); }
    if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='c'){  cvec_sub_cvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_C(x),ARRAY_P0_C(y)); }
    if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='R'){ idvec_sub_rvec(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P1_R(z),ARRAY_P0_D(x),ARRAY_P0_D(x),ARRAY_P0_R(y),ARRAY_P1_R(y)); }
    if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='R'){ izvec_sub_rvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_Z(x),ARRAY_P0_Z(x),ARRAY_P0_R(y),ARRAY_P1_R(y)); }
    if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='R'){ irvec_sub_rvec(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P1_R(z),ARRAY_P0_R(x),ARRAY_P0_R(x),ARRAY_P0_R(y),ARRAY_P1_R(y)); }
    if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='R'){ icvec_sub_rvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(x),ARRAY_P0_C(x),ARRAY_P0_R(y),ARRAY_P1_R(y)); }
    if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='d'){ irvec_sub_dvec(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P1_R(z),ARRAY_P0_R(x),ARRAY_P1_R(x),ARRAY_P0_D(y),ARRAY_P0_D(y)); }
    if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='z'){ irvec_sub_zvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_R(x),ARRAY_P1_R(x),ARRAY_P0_Z(y),ARRAY_P0_Z(y)); }
    if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='r'){ irvec_sub_rvec(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P1_R(z),ARRAY_P0_R(x),ARRAY_P1_R(x),ARRAY_P0_R(y),ARRAY_P0_R(y)); }
    if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='c'){ irvec_sub_cvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_R(x),ARRAY_P1_R(x),ARRAY_P0_C(y),ARRAY_P0_C(y)); }
    if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='R'){ irvec_sub_rvec(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P1_R(z),ARRAY_P0_R(x),ARRAY_P1_R(x),ARRAY_P0_R(y),ARRAY_P1_R(y)); }
    if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='C'){ idvec_sub_cvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_D(x),ARRAY_P0_D(x),ARRAY_P0_C(y),ARRAY_P1_C(y)); }
    if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='C'){ izvec_sub_cvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_Z(x),ARRAY_P0_Z(x),ARRAY_P0_C(y),ARRAY_P1_C(y)); }
    if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='C'){ irvec_sub_cvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_R(x),ARRAY_P0_R(x),ARRAY_P0_C(y),ARRAY_P1_C(y)); }
    if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='C'){ icvec_sub_cvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(x),ARRAY_P0_C(x),ARRAY_P0_C(y),ARRAY_P1_C(y)); }
    if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='C'){ irvec_sub_cvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_R(x),ARRAY_P1_R(x),ARRAY_P0_C(y),ARRAY_P1_C(y)); }
    if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='d'){ icvec_sub_dvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(x),ARRAY_P1_C(x),ARRAY_P0_D(y),ARRAY_P0_D(y)); }
    if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='z'){ icvec_sub_zvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(x),ARRAY_P1_C(x),ARRAY_P0_Z(y),ARRAY_P0_Z(y)); }
    if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='r'){ icvec_sub_rvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(x),ARRAY_P1_C(x),ARRAY_P0_R(y),ARRAY_P0_R(y)); }
    if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='c'){ icvec_sub_cvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(x),ARRAY_P1_C(x),ARRAY_P0_C(y),ARRAY_P0_C(y)); }
    if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='R'){ icvec_sub_rvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(x),ARRAY_P1_C(x),ARRAY_P0_R(y),ARRAY_P1_R(y)); }
    if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='C'){ icvec_sub_cvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(x),ARRAY_P1_C(x),ARRAY_P0_C(y),ARRAY_P1_C(y)); }
  }
  // yがスカラーの場合
  if(z==NULL && array_is_scalar(y))
  {
    z=array_allocate(array_get_inclusion_type(ARRAY_TYPE(x),ARRAY_TYPE(y)),ARRAY_NDIM(x),ARRAY_DIM_P(x));
    if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='d'){  dvec_sub_dscalar(ARRAY_SIZE(z),ARRAY_P0_D(z),ARRAY_P0_D(x),ARRAY_P0_D(y)[0]); }
    if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='z'){  dvec_sub_zscalar(ARRAY_SIZE(z),ARRAY_P0_Z(z),ARRAY_P0_D(x),ARRAY_P0_Z(y)[0]); }
    if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='d'){  zvec_sub_dscalar(ARRAY_SIZE(z),ARRAY_P0_Z(z),ARRAY_P0_Z(x),ARRAY_P0_D(y)[0]); }
    if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='z'){  zvec_sub_zscalar(ARRAY_SIZE(z),ARRAY_P0_Z(z),ARRAY_P0_Z(x),ARRAY_P0_Z(y)[0]); }
    if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='r'){  dvec_sub_rscalar(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P0_D(x),ARRAY_P0_R(y)[0]); }
    if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='r'){  zvec_sub_rscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_Z(x),ARRAY_P0_R(y)[0]); }
    if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='d'){  rvec_sub_dscalar(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P0_R(x),ARRAY_P0_D(y)[0]); }
    if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='z'){  rvec_sub_zscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_R(x),ARRAY_P0_Z(y)[0]); }
    if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='r'){  rvec_sub_rscalar(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P0_R(x),ARRAY_P0_R(y)[0]); }
    if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='c'){  dvec_sub_cscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_D(x),ARRAY_P0_C(y)[0]); }
    if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='c'){  zvec_sub_cscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_Z(x),ARRAY_P0_C(y)[0]); }
    if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='c'){  rvec_sub_cscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_R(x),ARRAY_P0_C(y)[0]); }
    if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='d'){  cvec_sub_dscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_C(x),ARRAY_P0_D(y)[0]); }
    if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='z'){  cvec_sub_zscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_C(x),ARRAY_P0_Z(y)[0]); }
    if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='r'){  cvec_sub_rscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_C(x),ARRAY_P0_R(y)[0]); }
    if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='c'){  cvec_sub_cscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_C(x),ARRAY_P0_C(y)[0]); }
    if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='R'){ idvec_sub_rscalar(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P1_R(z),ARRAY_P0_D(x),ARRAY_P0_D(x),ARRAY_P0_R(y)[0],ARRAY_P1_R(y)[0]); }
    if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='R'){ izvec_sub_rscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_Z(x),ARRAY_P0_Z(x),ARRAY_P0_R(y)[0],ARRAY_P1_R(y)[0]); }
    if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='R'){ irvec_sub_rscalar(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P1_R(z),ARRAY_P0_R(x),ARRAY_P0_R(x),ARRAY_P0_R(y)[0],ARRAY_P1_R(y)[0]); }
    if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='R'){ icvec_sub_rscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(x),ARRAY_P0_C(x),ARRAY_P0_R(y)[0],ARRAY_P1_R(y)[0]); }
    if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='d'){ irvec_sub_dscalar(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P1_R(z),ARRAY_P0_R(x),ARRAY_P1_R(x),ARRAY_P0_D(y)[0],ARRAY_P0_D(y)[0]); }
    if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='z'){ irvec_sub_zscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_R(x),ARRAY_P1_R(x),ARRAY_P0_Z(y)[0],ARRAY_P0_Z(y)[0]); }
    if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='r'){ irvec_sub_rscalar(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P1_R(z),ARRAY_P0_R(x),ARRAY_P1_R(x),ARRAY_P0_R(y)[0],ARRAY_P0_R(y)[0]); }
    if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='c'){ irvec_sub_cscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_R(x),ARRAY_P1_R(x),ARRAY_P0_C(y)[0],ARRAY_P0_C(y)[0]); }
    if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='R'){ irvec_sub_rscalar(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P1_R(z),ARRAY_P0_R(x),ARRAY_P1_R(x),ARRAY_P0_R(y)[0],ARRAY_P1_R(y)[0]); }
    if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='C'){ idvec_sub_cscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_D(x),ARRAY_P0_D(x),ARRAY_P0_C(y)[0],ARRAY_P1_C(y)[0]); }
    if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='C'){ izvec_sub_cscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_Z(x),ARRAY_P0_Z(x),ARRAY_P0_C(y)[0],ARRAY_P1_C(y)[0]); }
    if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='C'){ irvec_sub_cscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_R(x),ARRAY_P0_R(x),ARRAY_P0_C(y)[0],ARRAY_P1_C(y)[0]); }
    if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='C'){ icvec_sub_cscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(x),ARRAY_P0_C(x),ARRAY_P0_C(y)[0],ARRAY_P1_C(y)[0]); }
    if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='C'){ irvec_sub_cscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_R(x),ARRAY_P1_R(x),ARRAY_P0_C(y)[0],ARRAY_P1_C(y)[0]); }
    if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='d'){ icvec_sub_dscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(x),ARRAY_P1_C(x),ARRAY_P0_D(y)[0],ARRAY_P0_D(y)[0]); }
    if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='z'){ icvec_sub_zscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(x),ARRAY_P1_C(x),ARRAY_P0_Z(y)[0],ARRAY_P0_Z(y)[0]); }
    if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='r'){ icvec_sub_rscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(x),ARRAY_P1_C(x),ARRAY_P0_R(y)[0],ARRAY_P0_R(y)[0]); }
    if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='c'){ icvec_sub_cscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(x),ARRAY_P1_C(x),ARRAY_P0_C(y)[0],ARRAY_P0_C(y)[0]); }
    if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='R'){ icvec_sub_rscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(x),ARRAY_P1_C(x),ARRAY_P0_R(y)[0],ARRAY_P1_R(y)[0]); }
    if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='C'){ icvec_sub_cscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(x),ARRAY_P1_C(x),ARRAY_P0_C(y)[0],ARRAY_P1_C(y)[0]); }
  }
  // xがスカラーの場合
  if(z==NULL && array_is_scalar(x))
  {
    z=array_allocate(array_get_inclusion_type(ARRAY_TYPE(x),ARRAY_TYPE(y)),ARRAY_NDIM(y),ARRAY_DIM_P(y));
    if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='d'){  dscalar_sub_dvec(ARRAY_SIZE(z),ARRAY_P0_D(z),ARRAY_P0_D(x)[0],ARRAY_P0_D(y)); }
    if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='z'){  dscalar_sub_zvec(ARRAY_SIZE(z),ARRAY_P0_Z(z),ARRAY_P0_D(x)[0],ARRAY_P0_Z(y)); }
    if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='d'){  zscalar_sub_dvec(ARRAY_SIZE(z),ARRAY_P0_Z(z),ARRAY_P0_Z(x)[0],ARRAY_P0_D(y)); }
    if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='z'){  zscalar_sub_zvec(ARRAY_SIZE(z),ARRAY_P0_Z(z),ARRAY_P0_Z(x)[0],ARRAY_P0_Z(y)); }
    if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='r'){  dscalar_sub_rvec(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P0_D(x)[0],ARRAY_P0_R(y)); }
    if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='r'){  zscalar_sub_rvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_Z(x)[0],ARRAY_P0_R(y)); }
    if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='d'){  rscalar_sub_dvec(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P0_R(x)[0],ARRAY_P0_D(y)); }
    if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='z'){  rscalar_sub_zvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_R(x)[0],ARRAY_P0_Z(y)); }
    if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='r'){  rscalar_sub_rvec(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P0_R(x)[0],ARRAY_P0_R(y)); }
    if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='c'){  dscalar_sub_cvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_D(x)[0],ARRAY_P0_C(y)); }
    if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='c'){  zscalar_sub_cvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_Z(x)[0],ARRAY_P0_C(y)); }
    if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='c'){  rscalar_sub_cvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_R(x)[0],ARRAY_P0_C(y)); }
    if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='d'){  cscalar_sub_dvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_C(x)[0],ARRAY_P0_D(y)); }
    if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='z'){  cscalar_sub_zvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_C(x)[0],ARRAY_P0_Z(y)); }
    if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='r'){  cscalar_sub_rvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_C(x)[0],ARRAY_P0_R(y)); }
    if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='c'){  cscalar_sub_cvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_C(x)[0],ARRAY_P0_C(y)); }
    if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='R'){ idscalar_sub_rvec(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P1_R(z),ARRAY_P0_D(x)[0],ARRAY_P0_D(x)[0],ARRAY_P0_R(y),ARRAY_P1_R(y)); }
    if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='R'){ izscalar_sub_rvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_Z(x)[0],ARRAY_P0_Z(x)[0],ARRAY_P0_R(y),ARRAY_P1_R(y)); }
    if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='R'){ irscalar_sub_rvec(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P1_R(z),ARRAY_P0_R(x)[0],ARRAY_P0_R(x)[0],ARRAY_P0_R(y),ARRAY_P1_R(y)); }
    if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='R'){ icscalar_sub_rvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(x)[0],ARRAY_P0_C(x)[0],ARRAY_P0_R(y),ARRAY_P1_R(y)); }
    if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='d'){ irscalar_sub_dvec(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P1_R(z),ARRAY_P0_R(x)[0],ARRAY_P1_R(x)[0],ARRAY_P0_D(y),ARRAY_P0_D(y)); }
    if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='z'){ irscalar_sub_zvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_R(x)[0],ARRAY_P1_R(x)[0],ARRAY_P0_Z(y),ARRAY_P0_Z(y)); }
    if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='r'){ irscalar_sub_rvec(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P1_R(z),ARRAY_P0_R(x)[0],ARRAY_P1_R(x)[0],ARRAY_P0_R(y),ARRAY_P0_R(y)); }
    if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='c'){ irscalar_sub_cvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_R(x)[0],ARRAY_P1_R(x)[0],ARRAY_P0_C(y),ARRAY_P0_C(y)); }
    if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='R'){ irscalar_sub_rvec(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P1_R(z),ARRAY_P0_R(x)[0],ARRAY_P1_R(x)[0],ARRAY_P0_R(y),ARRAY_P1_R(y)); }
    if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='C'){ idscalar_sub_cvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_D(x)[0],ARRAY_P0_D(x)[0],ARRAY_P0_C(y),ARRAY_P1_C(y)); }
    if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='C'){ izscalar_sub_cvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_Z(x)[0],ARRAY_P0_Z(x)[0],ARRAY_P0_C(y),ARRAY_P1_C(y)); }
    if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='C'){ irscalar_sub_cvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_R(x)[0],ARRAY_P0_R(x)[0],ARRAY_P0_C(y),ARRAY_P1_C(y)); }
    if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='C'){ icscalar_sub_cvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(x)[0],ARRAY_P0_C(x)[0],ARRAY_P0_C(y),ARRAY_P1_C(y)); }
    if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='C'){ irscalar_sub_cvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_R(x)[0],ARRAY_P1_R(x)[0],ARRAY_P0_C(y),ARRAY_P1_C(y)); }
    if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='d'){ icscalar_sub_dvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(x)[0],ARRAY_P1_C(x)[0],ARRAY_P0_D(y),ARRAY_P0_D(y)); }
    if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='z'){ icscalar_sub_zvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(x)[0],ARRAY_P1_C(x)[0],ARRAY_P0_Z(y),ARRAY_P0_Z(y)); }
    if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='r'){ icscalar_sub_rvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(x)[0],ARRAY_P1_C(x)[0],ARRAY_P0_R(y),ARRAY_P0_R(y)); }
    if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='c'){ icscalar_sub_cvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(x)[0],ARRAY_P1_C(x)[0],ARRAY_P0_C(y),ARRAY_P0_C(y)); }
    if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='R'){ icscalar_sub_rvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(x)[0],ARRAY_P1_C(x)[0],ARRAY_P0_R(y),ARRAY_P1_R(y)); }
    if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='C'){ icscalar_sub_cvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(x)[0],ARRAY_P1_C(x)[0],ARRAY_P0_C(y),ARRAY_P1_C(y)); }
  }
  // 互換性のあるサイズでありxが小さい場合
  if(z==NULL && (k=array_get_subdim(x,y))>0){
    z=array_allocate(array_get_inclusion_type(ARRAY_TYPE(x),ARRAY_TYPE(y)),ARRAY_NDIM(y),ARRAY_DIM_P(y));
    for(i=0; i<ARRAY_LD_END(z)/ARRAY_LD(z,k-1); i++){
      if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='d'){  dvec_sub_dvec(ARRAY_LD(z,k-1),ARRAY_P0_D(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_D(x),ARRAY_P0_D(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='z'){  dvec_sub_zvec(ARRAY_LD(z,k-1),ARRAY_P0_Z(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_D(x),ARRAY_P0_Z(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='d'){  zvec_sub_dvec(ARRAY_LD(z,k-1),ARRAY_P0_Z(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_Z(x),ARRAY_P0_D(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='z'){  zvec_sub_zvec(ARRAY_LD(z,k-1),ARRAY_P0_Z(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_Z(x),ARRAY_P0_Z(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='r'){  dvec_sub_rvec(ARRAY_LD(z,k-1),ARRAY_P0_R(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_D(x),ARRAY_P0_R(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='r'){  zvec_sub_rvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_Z(x),ARRAY_P0_R(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='d'){  rvec_sub_dvec(ARRAY_LD(z,k-1),ARRAY_P0_R(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_R(x),ARRAY_P0_D(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='z'){  rvec_sub_zvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_R(x),ARRAY_P0_Z(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='r'){  rvec_sub_rvec(ARRAY_LD(z,k-1),ARRAY_P0_R(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_R(x),ARRAY_P0_R(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='c'){  dvec_sub_cvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_D(x),ARRAY_P0_C(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='c'){  zvec_sub_cvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_Z(x),ARRAY_P0_C(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='c'){  rvec_sub_cvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_R(x),ARRAY_P0_C(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='d'){  cvec_sub_dvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_C(x),ARRAY_P0_D(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='z'){  cvec_sub_zvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_C(x),ARRAY_P0_Z(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='r'){  cvec_sub_rvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_C(x),ARRAY_P0_R(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='c'){  cvec_sub_cvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_C(x),ARRAY_P0_C(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='R'){ idvec_sub_rvec(ARRAY_LD(z,k-1),ARRAY_P0_R(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_R(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_D(x),ARRAY_P0_D(x),ARRAY_P0_R(y)+i*ARRAY_LD(y,k-1),ARRAY_P1_R(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='R'){ izvec_sub_rvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_Z(x),ARRAY_P0_Z(x),ARRAY_P0_R(y)+i*ARRAY_LD(y,k-1),ARRAY_P1_R(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='R'){ irvec_sub_rvec(ARRAY_LD(z,k-1),ARRAY_P0_R(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_R(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_R(x),ARRAY_P0_R(x),ARRAY_P0_R(y)+i*ARRAY_LD(y,k-1),ARRAY_P1_R(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='R'){ icvec_sub_rvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_C(x),ARRAY_P0_C(x),ARRAY_P0_R(y)+i*ARRAY_LD(y,k-1),ARRAY_P1_R(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='d'){ irvec_sub_dvec(ARRAY_LD(z,k-1),ARRAY_P0_R(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_R(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_R(x),ARRAY_P1_R(x),ARRAY_P0_D(y)+i*ARRAY_LD(y,k-1),ARRAY_P0_D(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='z'){ irvec_sub_zvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_R(x),ARRAY_P1_R(x),ARRAY_P0_Z(y)+i*ARRAY_LD(y,k-1),ARRAY_P0_Z(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='r'){ irvec_sub_rvec(ARRAY_LD(z,k-1),ARRAY_P0_R(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_R(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_R(x),ARRAY_P1_R(x),ARRAY_P0_R(y)+i*ARRAY_LD(y,k-1),ARRAY_P0_R(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='c'){ irvec_sub_cvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_R(x),ARRAY_P1_R(x),ARRAY_P0_C(y)+i*ARRAY_LD(y,k-1),ARRAY_P0_C(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='R'){ irvec_sub_rvec(ARRAY_LD(z,k-1),ARRAY_P0_R(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_R(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_R(x),ARRAY_P1_R(x),ARRAY_P0_R(y)+i*ARRAY_LD(y,k-1),ARRAY_P1_R(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='C'){ idvec_sub_cvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_D(x),ARRAY_P0_D(x),ARRAY_P0_C(y)+i*ARRAY_LD(y,k-1),ARRAY_P1_C(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='C'){ izvec_sub_cvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_Z(x),ARRAY_P0_Z(x),ARRAY_P0_C(y)+i*ARRAY_LD(y,k-1),ARRAY_P1_C(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='C'){ irvec_sub_cvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_R(x),ARRAY_P0_R(x),ARRAY_P0_C(y)+i*ARRAY_LD(y,k-1),ARRAY_P1_C(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='C'){ icvec_sub_cvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_C(x),ARRAY_P0_C(x),ARRAY_P0_C(y)+i*ARRAY_LD(y,k-1),ARRAY_P1_C(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='C'){ irvec_sub_cvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_R(x),ARRAY_P1_R(x),ARRAY_P0_C(y)+i*ARRAY_LD(y,k-1),ARRAY_P1_C(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='d'){ icvec_sub_dvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_C(x),ARRAY_P1_C(x),ARRAY_P0_D(y)+i*ARRAY_LD(y,k-1),ARRAY_P0_D(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='z'){ icvec_sub_zvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_C(x),ARRAY_P1_C(x),ARRAY_P0_Z(y)+i*ARRAY_LD(y,k-1),ARRAY_P0_Z(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='r'){ icvec_sub_rvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_C(x),ARRAY_P1_C(x),ARRAY_P0_R(y)+i*ARRAY_LD(y,k-1),ARRAY_P0_R(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='c'){ icvec_sub_cvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_C(x),ARRAY_P1_C(x),ARRAY_P0_C(y)+i*ARRAY_LD(y,k-1),ARRAY_P0_C(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='R'){ icvec_sub_rvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_C(x),ARRAY_P1_C(x),ARRAY_P0_R(y)+i*ARRAY_LD(y,k-1),ARRAY_P1_R(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='C'){ icvec_sub_cvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_C(x),ARRAY_P1_C(x),ARRAY_P0_C(y)+i*ARRAY_LD(y,k-1),ARRAY_P1_C(y)+i*ARRAY_LD(y,k-1)); }
    }
  }
  // 互換性のあるサイズでありyが小さい場合
  if(z==NULL && (k=array_get_subdim(y,x))>0){
    z=array_allocate(array_get_inclusion_type(ARRAY_TYPE(x),ARRAY_TYPE(y)),ARRAY_NDIM(x),ARRAY_DIM_P(x));
    for(i=0; i<ARRAY_LD_END(z)/ARRAY_LD(z,k-1); i++){
      if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='d'){  dvec_sub_dvec(ARRAY_LD(z,k-1),ARRAY_P0_D(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_D(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_D(y)); }
      if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='z'){  dvec_sub_zvec(ARRAY_LD(z,k-1),ARRAY_P0_Z(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_D(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_Z(y)); }
      if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='d'){  zvec_sub_dvec(ARRAY_LD(z,k-1),ARRAY_P0_Z(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_Z(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_D(y)); }
      if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='z'){  zvec_sub_zvec(ARRAY_LD(z,k-1),ARRAY_P0_Z(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_Z(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_Z(y)); }
      if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='r'){  dvec_sub_rvec(ARRAY_LD(z,k-1),ARRAY_P0_R(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_D(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_R(y)); }
      if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='r'){  zvec_sub_rvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_Z(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_R(y)); }
      if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='d'){  rvec_sub_dvec(ARRAY_LD(z,k-1),ARRAY_P0_R(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_R(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_D(y)); }
      if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='z'){  rvec_sub_zvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_R(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_Z(y)); }
      if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='r'){  rvec_sub_rvec(ARRAY_LD(z,k-1),ARRAY_P0_R(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_R(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_R(y)); }
      if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='c'){  dvec_sub_cvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_D(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_C(y)); }
      if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='c'){  zvec_sub_cvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_Z(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_C(y)); }
      if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='c'){  rvec_sub_cvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_R(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_C(y)); }
      if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='d'){  cvec_sub_dvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_C(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_D(y)); }
      if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='z'){  cvec_sub_zvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_C(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_Z(y)); }
      if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='r'){  cvec_sub_rvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_C(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_R(y)); }
      if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='c'){  cvec_sub_cvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_C(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_C(y)); }
      if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='R'){ idvec_sub_rvec(ARRAY_LD(z,k-1),ARRAY_P0_R(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_R(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_D(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_D(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_R(y),ARRAY_P1_R(y)); }
      if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='R'){ izvec_sub_rvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_Z(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_Z(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_R(y),ARRAY_P1_R(y)); }
      if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='R'){ irvec_sub_rvec(ARRAY_LD(z,k-1),ARRAY_P0_R(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_R(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_R(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_R(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_R(y),ARRAY_P1_R(y)); }
      if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='R'){ icvec_sub_rvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_C(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_C(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_R(y),ARRAY_P1_R(y)); }
      if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='d'){ irvec_sub_dvec(ARRAY_LD(z,k-1),ARRAY_P0_R(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_R(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_R(x)+i*ARRAY_LD(x,k-1),ARRAY_P1_R(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_D(y),ARRAY_P0_D(y)); }
      if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='z'){ irvec_sub_zvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_R(x)+i*ARRAY_LD(x,k-1),ARRAY_P1_R(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_Z(y),ARRAY_P0_Z(y)); }
      if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='r'){ irvec_sub_rvec(ARRAY_LD(z,k-1),ARRAY_P0_R(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_R(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_R(x)+i*ARRAY_LD(x,k-1),ARRAY_P1_R(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_R(y),ARRAY_P0_R(y)); }
      if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='c'){ irvec_sub_cvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_R(x)+i*ARRAY_LD(x,k-1),ARRAY_P1_R(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_C(y),ARRAY_P0_C(y)); }
      if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='R'){ irvec_sub_rvec(ARRAY_LD(z,k-1),ARRAY_P0_R(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_R(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_R(x)+i*ARRAY_LD(x,k-1),ARRAY_P1_R(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_R(y),ARRAY_P1_R(y)); }
      if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='C'){ idvec_sub_cvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_D(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_D(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_C(y),ARRAY_P1_C(y)); }
      if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='C'){ izvec_sub_cvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_Z(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_Z(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_C(y),ARRAY_P1_C(y)); }
      if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='C'){ irvec_sub_cvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_R(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_R(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_C(y),ARRAY_P1_C(y)); }
      if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='C'){ icvec_sub_cvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_C(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_C(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_C(y),ARRAY_P1_C(y)); }
      if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='C'){ irvec_sub_cvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_R(x)+i*ARRAY_LD(x,k-1),ARRAY_P1_R(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_C(y),ARRAY_P1_C(y)); }
      if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='d'){ icvec_sub_dvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_C(x)+i*ARRAY_LD(x,k-1),ARRAY_P1_C(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_D(y),ARRAY_P0_D(y)); }
      if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='z'){ icvec_sub_zvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_C(x)+i*ARRAY_LD(x,k-1),ARRAY_P1_C(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_Z(y),ARRAY_P0_Z(y)); }
      if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='r'){ icvec_sub_rvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_C(x)+i*ARRAY_LD(x,k-1),ARRAY_P1_C(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_R(y),ARRAY_P0_R(y)); }
      if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='c'){ icvec_sub_cvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_C(x)+i*ARRAY_LD(x,k-1),ARRAY_P1_C(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_C(y),ARRAY_P0_C(y)); }
      if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='R'){ icvec_sub_rvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_C(x)+i*ARRAY_LD(x,k-1),ARRAY_P1_C(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_R(y),ARRAY_P1_R(y)); }
      if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='C'){ icvec_sub_cvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_C(x)+i*ARRAY_LD(x,k-1),ARRAY_P1_C(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_C(y),ARRAY_P1_C(y)); }
    }
  }
  // 互換性のあるサイズの場合
  if(z==NULL && array_compatible_dim_check(x,y)){
    ndim=MAX2(ARRAY_NDIM(x),ARRAY_NDIM(y));
    dim=array_get_inclusion_dim(x,y);
    z=array_allocate(array_get_inclusion_type(ARRAY_TYPE(x),ARRAY_TYPE(y)),ndim,dim);
    // loop
    I=ivec_allocate(ARRAY_NDIM(z));
    K=ivec_allocate(ARRAY_NDIM(z));
    L=ivec_allocate(ARRAY_NDIM(z));
    for(i=0; i<ARRAY_SIZE(z); i++){
      array_get_index(I,z,i);
      ivec_copy(ndim,K,I);
      ivec_copy(ndim,L,I);
      for(j=0; j<ndim; j++){
        if(K[j]>=ARRAY_DIM(x,j)){ K[j]=ARRAY_DIM(x,j)-1; }
        if(L[j]>=ARRAY_DIM(y,j)){ L[j]=ARRAY_DIM(y,j)-1; }
      }
      k=array_get_position(K,x);
      l=array_get_position(L,y);
      if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='d'){ ARRAY_DVEC0(z,i)=ARRAY_DVEC0(x,k)-ARRAY_DVEC0(y,l); }
      if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='z'){ ARRAY_ZVEC0(z,i)=zsub_dz(ARRAY_DVEC0(x,k),ARRAY_ZVEC0(y,l)); }
      if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='d'){ ARRAY_ZVEC0(z,i)=zsub_zd(ARRAY_ZVEC0(x,k),ARRAY_DVEC0(y,l)); }
      if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='z'){ ARRAY_ZVEC0(z,i)=zsub_zz(ARRAY_ZVEC0(x,k),ARRAY_ZVEC0(y,l)); }
      if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='r'){  rsub_dr(ARRAY_RVEC0(z,i),ARRAY_DVEC0(x,k),ARRAY_RVEC0(y,l)); }
      if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='r'){  csub_zr(ARRAY_CVEC0(z,i),ARRAY_ZVEC0(x,k),ARRAY_RVEC0(y,l)); }
      if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='d'){  rsub_rd(ARRAY_RVEC0(z,i),ARRAY_RVEC0(x,k),ARRAY_DVEC0(y,l)); }
      if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='z'){  csub_rz(ARRAY_CVEC0(z,i),ARRAY_RVEC0(x,k),ARRAY_ZVEC0(y,l)); }
      if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='r'){  rsub_rr(ARRAY_RVEC0(z,i),ARRAY_RVEC0(x,k),ARRAY_RVEC0(y,l)); }
      if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='c'){  csub_dc(ARRAY_CVEC0(z,i),ARRAY_DVEC0(x,k),ARRAY_CVEC0(y,l)); }
      if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='c'){  csub_zc(ARRAY_CVEC0(z,i),ARRAY_ZVEC0(x,k),ARRAY_CVEC0(y,l)); }
      if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='c'){  csub_rc(ARRAY_CVEC0(z,i),ARRAY_RVEC0(x,k),ARRAY_CVEC0(y,l)); }
      if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='d'){  csub_cd(ARRAY_CVEC0(z,i),ARRAY_CVEC0(x,k),ARRAY_DVEC0(y,l)); }
      if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='z'){  csub_cz(ARRAY_CVEC0(z,i),ARRAY_CVEC0(x,k),ARRAY_ZVEC0(y,l)); }
      if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='r'){  csub_cr(ARRAY_CVEC0(z,i),ARRAY_CVEC0(x,k),ARRAY_RVEC0(y,l)); }
      if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='c'){  csub_cc(ARRAY_CVEC0(z,i),ARRAY_CVEC0(x,k),ARRAY_CVEC0(y,l)); }
      if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='R'){ irsub_dr(ARRAY_RVEC0(z,i),ARRAY_RVEC1(z,i),ARRAY_DVEC0(x,k),ARRAY_DVEC0(x,k),ARRAY_RVEC0(y,l),ARRAY_RVEC1(y,l)); }
      if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='R'){ icsub_zr(ARRAY_CVEC0(z,i),ARRAY_CVEC1(z,i),ARRAY_ZVEC0(x,k),ARRAY_ZVEC0(x,k),ARRAY_RVEC0(y,l),ARRAY_RVEC1(y,l)); }
      if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='R'){ irsub_rr(ARRAY_RVEC0(z,i),ARRAY_RVEC1(z,i),ARRAY_RVEC0(x,k),ARRAY_RVEC0(x,k),ARRAY_RVEC0(y,l),ARRAY_RVEC1(y,l)); }
      if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='R'){ icsub_cr(ARRAY_CVEC0(z,i),ARRAY_CVEC1(z,i),ARRAY_CVEC0(x,k),ARRAY_CVEC0(x,k),ARRAY_RVEC0(y,l),ARRAY_RVEC1(y,l)); }
      if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='d'){ irsub_rd(ARRAY_RVEC0(z,i),ARRAY_RVEC1(z,i),ARRAY_RVEC0(x,k),ARRAY_RVEC1(x,k),ARRAY_DVEC0(y,l),ARRAY_DVEC0(y,l)); }
      if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='z'){ icsub_rz(ARRAY_CVEC0(z,i),ARRAY_CVEC1(z,i),ARRAY_RVEC0(x,k),ARRAY_RVEC1(x,k),ARRAY_ZVEC0(y,l),ARRAY_ZVEC0(y,l)); }
      if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='r'){ irsub_rr(ARRAY_RVEC0(z,i),ARRAY_RVEC1(z,i),ARRAY_RVEC0(x,k),ARRAY_RVEC1(x,k),ARRAY_RVEC0(y,l),ARRAY_RVEC0(y,l)); }
      if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='c'){ icsub_rc(ARRAY_CVEC0(z,i),ARRAY_CVEC1(z,i),ARRAY_RVEC0(x,k),ARRAY_RVEC1(x,k),ARRAY_CVEC0(y,l),ARRAY_CVEC0(y,l)); }
      if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='R'){ irsub_rr(ARRAY_RVEC0(z,i),ARRAY_RVEC1(z,i),ARRAY_RVEC0(x,k),ARRAY_RVEC1(x,k),ARRAY_RVEC0(y,l),ARRAY_RVEC1(y,l)); }
      if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='C'){ icsub_dc(ARRAY_CVEC0(z,i),ARRAY_CVEC1(z,i),ARRAY_DVEC0(x,k),ARRAY_DVEC0(x,k),ARRAY_CVEC0(y,l),ARRAY_CVEC1(y,l)); }
      if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='C'){ icsub_zc(ARRAY_CVEC0(z,i),ARRAY_CVEC1(z,i),ARRAY_ZVEC0(x,k),ARRAY_ZVEC0(x,k),ARRAY_CVEC0(y,l),ARRAY_CVEC1(y,l)); }
      if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='C'){ icsub_rc(ARRAY_CVEC0(z,i),ARRAY_CVEC1(z,i),ARRAY_RVEC0(x,k),ARRAY_RVEC0(x,k),ARRAY_CVEC0(y,l),ARRAY_CVEC1(y,l)); }
      if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='C'){ icsub_cc(ARRAY_CVEC0(z,i),ARRAY_CVEC1(z,i),ARRAY_CVEC0(x,k),ARRAY_CVEC0(x,k),ARRAY_CVEC0(y,l),ARRAY_CVEC1(y,l)); }
      if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='C'){ icsub_rc(ARRAY_CVEC0(z,i),ARRAY_CVEC1(z,i),ARRAY_RVEC0(x,k),ARRAY_RVEC1(x,k),ARRAY_CVEC0(y,l),ARRAY_CVEC1(y,l)); }
      if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='d'){ icsub_cd(ARRAY_CVEC0(z,i),ARRAY_CVEC1(z,i),ARRAY_CVEC0(x,k),ARRAY_CVEC1(x,k),ARRAY_DVEC0(y,l),ARRAY_DVEC0(y,l)); }
      if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='z'){ icsub_cz(ARRAY_CVEC0(z,i),ARRAY_CVEC1(z,i),ARRAY_CVEC0(x,k),ARRAY_CVEC1(x,k),ARRAY_ZVEC0(y,l),ARRAY_ZVEC0(y,l)); }
      if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='r'){ icsub_cr(ARRAY_CVEC0(z,i),ARRAY_CVEC1(z,i),ARRAY_CVEC0(x,k),ARRAY_CVEC1(x,k),ARRAY_RVEC0(y,l),ARRAY_RVEC0(y,l)); }
      if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='c'){ icsub_cc(ARRAY_CVEC0(z,i),ARRAY_CVEC1(z,i),ARRAY_CVEC0(x,k),ARRAY_CVEC1(x,k),ARRAY_CVEC0(y,l),ARRAY_CVEC0(y,l)); }
      if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='R'){ icsub_cr(ARRAY_CVEC0(z,i),ARRAY_CVEC1(z,i),ARRAY_CVEC0(x,k),ARRAY_CVEC1(x,k),ARRAY_RVEC0(y,l),ARRAY_RVEC1(y,l)); }
      if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='C'){ icsub_cc(ARRAY_CVEC0(z,i),ARRAY_CVEC1(z,i),ARRAY_CVEC0(x,k),ARRAY_CVEC1(x,k),ARRAY_CVEC0(y,l),ARRAY_CVEC1(y,l)); }
    }
    dim=ivec_free(dim);
    I=ivec_free(I);
    K=ivec_free(K);
    L=ivec_free(L);
  }
  return z;
}

/**
 @brief z=x.*y
 */
array *array_mul(array *x, array *y)
{
  int i,j,k,l,ndim,*dim=NULL,*I=NULL,*K=NULL,*L=NULL;
  array *z=NULL;
  if(x==NULL || y==NULL){ return z; }
  // サイズが同じ場合
  if(z==NULL && array_same_dim_check(x,y))
  {
    z=array_allocate(array_get_inclusion_type(ARRAY_TYPE(x),ARRAY_TYPE(y)),ARRAY_NDIM(x),ARRAY_DIM_P(x));
    if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='d'){  dvec_mul_dvec(ARRAY_SIZE(z),ARRAY_P0_D(z),ARRAY_P0_D(x),ARRAY_P0_D(y)); }
    if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='z'){  dvec_mul_zvec(ARRAY_SIZE(z),ARRAY_P0_Z(z),ARRAY_P0_D(x),ARRAY_P0_Z(y)); }
    if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='d'){  zvec_mul_dvec(ARRAY_SIZE(z),ARRAY_P0_Z(z),ARRAY_P0_Z(x),ARRAY_P0_D(y)); }
    if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='z'){  zvec_mul_zvec(ARRAY_SIZE(z),ARRAY_P0_Z(z),ARRAY_P0_Z(x),ARRAY_P0_Z(y)); }
    if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='r'){  dvec_mul_rvec(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P0_D(x),ARRAY_P0_R(y)); }
    if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='r'){  zvec_mul_rvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_Z(x),ARRAY_P0_R(y)); }
    if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='d'){  rvec_mul_dvec(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P0_R(x),ARRAY_P0_D(y)); }
    if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='z'){  rvec_mul_zvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_R(x),ARRAY_P0_Z(y)); }
    if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='r'){  rvec_mul_rvec(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P0_R(x),ARRAY_P0_R(y)); }
    if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='c'){  dvec_mul_cvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_D(x),ARRAY_P0_C(y)); }
    if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='c'){  zvec_mul_cvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_Z(x),ARRAY_P0_C(y)); }
    if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='c'){  rvec_mul_cvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_R(x),ARRAY_P0_C(y)); }
    if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='d'){  cvec_mul_dvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_C(x),ARRAY_P0_D(y)); }
    if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='z'){  cvec_mul_zvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_C(x),ARRAY_P0_Z(y)); }
    if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='r'){  cvec_mul_rvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_C(x),ARRAY_P0_R(y)); }
    if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='c'){  cvec_mul_cvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_C(x),ARRAY_P0_C(y)); }
    if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='R'){ idvec_mul_rvec(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P1_R(z),ARRAY_P0_D(x),ARRAY_P0_D(x),ARRAY_P0_R(y),ARRAY_P1_R(y)); }
    if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='R'){ izvec_mul_rvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_Z(x),ARRAY_P0_Z(x),ARRAY_P0_R(y),ARRAY_P1_R(y)); }
    if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='R'){ irvec_mul_rvec(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P1_R(z),ARRAY_P0_R(x),ARRAY_P0_R(x),ARRAY_P0_R(y),ARRAY_P1_R(y)); }
    if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='R'){ icvec_mul_rvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(x),ARRAY_P0_C(x),ARRAY_P0_R(y),ARRAY_P1_R(y)); }
    if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='d'){ irvec_mul_dvec(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P1_R(z),ARRAY_P0_R(x),ARRAY_P1_R(x),ARRAY_P0_D(y),ARRAY_P0_D(y)); }
    if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='z'){ irvec_mul_zvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_R(x),ARRAY_P1_R(x),ARRAY_P0_Z(y),ARRAY_P0_Z(y)); }
    if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='r'){ irvec_mul_rvec(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P1_R(z),ARRAY_P0_R(x),ARRAY_P1_R(x),ARRAY_P0_R(y),ARRAY_P0_R(y)); }
    if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='c'){ irvec_mul_cvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_R(x),ARRAY_P1_R(x),ARRAY_P0_C(y),ARRAY_P0_C(y)); }
    if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='R'){ irvec_mul_rvec(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P1_R(z),ARRAY_P0_R(x),ARRAY_P1_R(x),ARRAY_P0_R(y),ARRAY_P1_R(y)); }
    if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='C'){ idvec_mul_cvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_D(x),ARRAY_P0_D(x),ARRAY_P0_C(y),ARRAY_P1_C(y)); }
    if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='C'){ izvec_mul_cvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_Z(x),ARRAY_P0_Z(x),ARRAY_P0_C(y),ARRAY_P1_C(y)); }
    if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='C'){ irvec_mul_cvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_R(x),ARRAY_P0_R(x),ARRAY_P0_C(y),ARRAY_P1_C(y)); }
    if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='C'){ icvec_mul_cvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(x),ARRAY_P0_C(x),ARRAY_P0_C(y),ARRAY_P1_C(y)); }
    if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='C'){ irvec_mul_cvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_R(x),ARRAY_P1_R(x),ARRAY_P0_C(y),ARRAY_P1_C(y)); }
    if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='d'){ icvec_mul_dvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(x),ARRAY_P1_C(x),ARRAY_P0_D(y),ARRAY_P0_D(y)); }
    if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='z'){ icvec_mul_zvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(x),ARRAY_P1_C(x),ARRAY_P0_Z(y),ARRAY_P0_Z(y)); }
    if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='r'){ icvec_mul_rvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(x),ARRAY_P1_C(x),ARRAY_P0_R(y),ARRAY_P0_R(y)); }
    if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='c'){ icvec_mul_cvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(x),ARRAY_P1_C(x),ARRAY_P0_C(y),ARRAY_P0_C(y)); }
    if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='R'){ icvec_mul_rvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(x),ARRAY_P1_C(x),ARRAY_P0_R(y),ARRAY_P1_R(y)); }
    if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='C'){ icvec_mul_cvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(x),ARRAY_P1_C(x),ARRAY_P0_C(y),ARRAY_P1_C(y)); }
  }
  // yがスカラーの場合
  if(z==NULL && array_is_scalar(y))
  {
    z=array_allocate(array_get_inclusion_type(ARRAY_TYPE(x),ARRAY_TYPE(y)),ARRAY_NDIM(x),ARRAY_DIM_P(x));
    if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='d'){  dvec_mul_dscalar(ARRAY_SIZE(z),ARRAY_P0_D(z),ARRAY_P0_D(x),ARRAY_P0_D(y)[0]); }
    if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='z'){  dvec_mul_zscalar(ARRAY_SIZE(z),ARRAY_P0_Z(z),ARRAY_P0_D(x),ARRAY_P0_Z(y)[0]); }
    if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='d'){  zvec_mul_dscalar(ARRAY_SIZE(z),ARRAY_P0_Z(z),ARRAY_P0_Z(x),ARRAY_P0_D(y)[0]); }
    if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='z'){  zvec_mul_zscalar(ARRAY_SIZE(z),ARRAY_P0_Z(z),ARRAY_P0_Z(x),ARRAY_P0_Z(y)[0]); }
    if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='r'){  dvec_mul_rscalar(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P0_D(x),ARRAY_P0_R(y)[0]); }
    if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='r'){  zvec_mul_rscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_Z(x),ARRAY_P0_R(y)[0]); }
    if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='d'){  rvec_mul_dscalar(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P0_R(x),ARRAY_P0_D(y)[0]); }
    if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='z'){  rvec_mul_zscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_R(x),ARRAY_P0_Z(y)[0]); }
    if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='r'){  rvec_mul_rscalar(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P0_R(x),ARRAY_P0_R(y)[0]); }
    if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='c'){  dvec_mul_cscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_D(x),ARRAY_P0_C(y)[0]); }
    if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='c'){  zvec_mul_cscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_Z(x),ARRAY_P0_C(y)[0]); }
    if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='c'){  rvec_mul_cscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_R(x),ARRAY_P0_C(y)[0]); }
    if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='d'){  cvec_mul_dscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_C(x),ARRAY_P0_D(y)[0]); }
    if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='z'){  cvec_mul_zscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_C(x),ARRAY_P0_Z(y)[0]); }
    if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='r'){  cvec_mul_rscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_C(x),ARRAY_P0_R(y)[0]); }
    if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='c'){  cvec_mul_cscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_C(x),ARRAY_P0_C(y)[0]); }
    if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='R'){ idvec_mul_rscalar(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P1_R(z),ARRAY_P0_D(x),ARRAY_P0_D(x),ARRAY_P0_R(y)[0],ARRAY_P1_R(y)[0]); }
    if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='R'){ izvec_mul_rscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_Z(x),ARRAY_P0_Z(x),ARRAY_P0_R(y)[0],ARRAY_P1_R(y)[0]); }
    if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='R'){ irvec_mul_rscalar(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P1_R(z),ARRAY_P0_R(x),ARRAY_P0_R(x),ARRAY_P0_R(y)[0],ARRAY_P1_R(y)[0]); }
    if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='R'){ icvec_mul_rscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(x),ARRAY_P0_C(x),ARRAY_P0_R(y)[0],ARRAY_P1_R(y)[0]); }
    if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='d'){ irvec_mul_dscalar(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P1_R(z),ARRAY_P0_R(x),ARRAY_P1_R(x),ARRAY_P0_D(y)[0],ARRAY_P0_D(y)[0]); }
    if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='z'){ irvec_mul_zscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_R(x),ARRAY_P1_R(x),ARRAY_P0_Z(y)[0],ARRAY_P0_Z(y)[0]); }
    if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='r'){ irvec_mul_rscalar(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P1_R(z),ARRAY_P0_R(x),ARRAY_P1_R(x),ARRAY_P0_R(y)[0],ARRAY_P0_R(y)[0]); }
    if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='c'){ irvec_mul_cscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_R(x),ARRAY_P1_R(x),ARRAY_P0_C(y)[0],ARRAY_P0_C(y)[0]); }
    if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='R'){ irvec_mul_rscalar(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P1_R(z),ARRAY_P0_R(x),ARRAY_P1_R(x),ARRAY_P0_R(y)[0],ARRAY_P1_R(y)[0]); }
    if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='C'){ idvec_mul_cscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_D(x),ARRAY_P0_D(x),ARRAY_P0_C(y)[0],ARRAY_P1_C(y)[0]); }
    if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='C'){ izvec_mul_cscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_Z(x),ARRAY_P0_Z(x),ARRAY_P0_C(y)[0],ARRAY_P1_C(y)[0]); }
    if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='C'){ irvec_mul_cscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_R(x),ARRAY_P0_R(x),ARRAY_P0_C(y)[0],ARRAY_P1_C(y)[0]); }
    if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='C'){ icvec_mul_cscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(x),ARRAY_P0_C(x),ARRAY_P0_C(y)[0],ARRAY_P1_C(y)[0]); }
    if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='C'){ irvec_mul_cscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_R(x),ARRAY_P1_R(x),ARRAY_P0_C(y)[0],ARRAY_P1_C(y)[0]); }
    if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='d'){ icvec_mul_dscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(x),ARRAY_P1_C(x),ARRAY_P0_D(y)[0],ARRAY_P0_D(y)[0]); }
    if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='z'){ icvec_mul_zscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(x),ARRAY_P1_C(x),ARRAY_P0_Z(y)[0],ARRAY_P0_Z(y)[0]); }
    if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='r'){ icvec_mul_rscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(x),ARRAY_P1_C(x),ARRAY_P0_R(y)[0],ARRAY_P0_R(y)[0]); }
    if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='c'){ icvec_mul_cscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(x),ARRAY_P1_C(x),ARRAY_P0_C(y)[0],ARRAY_P0_C(y)[0]); }
    if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='R'){ icvec_mul_rscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(x),ARRAY_P1_C(x),ARRAY_P0_R(y)[0],ARRAY_P1_R(y)[0]); }
    if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='C'){ icvec_mul_cscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(x),ARRAY_P1_C(x),ARRAY_P0_C(y)[0],ARRAY_P1_C(y)[0]); }
  }
  // xがスカラーの場合
  if(z==NULL && array_is_scalar(x))
  {
    z=array_allocate(array_get_inclusion_type(ARRAY_TYPE(x),ARRAY_TYPE(y)),ARRAY_NDIM(y),ARRAY_DIM_P(y));
    if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='d'){  dscalar_mul_dvec(ARRAY_SIZE(z),ARRAY_P0_D(z),ARRAY_P0_D(x)[0],ARRAY_P0_D(y)); }
    if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='z'){  dscalar_mul_zvec(ARRAY_SIZE(z),ARRAY_P0_Z(z),ARRAY_P0_D(x)[0],ARRAY_P0_Z(y)); }
    if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='d'){  zscalar_mul_dvec(ARRAY_SIZE(z),ARRAY_P0_Z(z),ARRAY_P0_Z(x)[0],ARRAY_P0_D(y)); }
    if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='z'){  zscalar_mul_zvec(ARRAY_SIZE(z),ARRAY_P0_Z(z),ARRAY_P0_Z(x)[0],ARRAY_P0_Z(y)); }
    if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='r'){  dscalar_mul_rvec(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P0_D(x)[0],ARRAY_P0_R(y)); }
    if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='r'){  zscalar_mul_rvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_Z(x)[0],ARRAY_P0_R(y)); }
    if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='d'){  rscalar_mul_dvec(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P0_R(x)[0],ARRAY_P0_D(y)); }
    if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='z'){  rscalar_mul_zvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_R(x)[0],ARRAY_P0_Z(y)); }
    if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='r'){  rscalar_mul_rvec(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P0_R(x)[0],ARRAY_P0_R(y)); }
    if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='c'){  dscalar_mul_cvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_D(x)[0],ARRAY_P0_C(y)); }
    if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='c'){  zscalar_mul_cvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_Z(x)[0],ARRAY_P0_C(y)); }
    if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='c'){  rscalar_mul_cvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_R(x)[0],ARRAY_P0_C(y)); }
    if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='d'){  cscalar_mul_dvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_C(x)[0],ARRAY_P0_D(y)); }
    if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='z'){  cscalar_mul_zvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_C(x)[0],ARRAY_P0_Z(y)); }
    if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='r'){  cscalar_mul_rvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_C(x)[0],ARRAY_P0_R(y)); }
    if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='c'){  cscalar_mul_cvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_C(x)[0],ARRAY_P0_C(y)); }
    if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='R'){ idscalar_mul_rvec(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P1_R(z),ARRAY_P0_D(x)[0],ARRAY_P0_D(x)[0],ARRAY_P0_R(y),ARRAY_P1_R(y)); }
    if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='R'){ izscalar_mul_rvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_Z(x)[0],ARRAY_P0_Z(x)[0],ARRAY_P0_R(y),ARRAY_P1_R(y)); }
    if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='R'){ irscalar_mul_rvec(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P1_R(z),ARRAY_P0_R(x)[0],ARRAY_P0_R(x)[0],ARRAY_P0_R(y),ARRAY_P1_R(y)); }
    if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='R'){ icscalar_mul_rvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(x)[0],ARRAY_P0_C(x)[0],ARRAY_P0_R(y),ARRAY_P1_R(y)); }
    if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='d'){ irscalar_mul_dvec(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P1_R(z),ARRAY_P0_R(x)[0],ARRAY_P1_R(x)[0],ARRAY_P0_D(y),ARRAY_P0_D(y)); }
    if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='z'){ irscalar_mul_zvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_R(x)[0],ARRAY_P1_R(x)[0],ARRAY_P0_Z(y),ARRAY_P0_Z(y)); }
    if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='r'){ irscalar_mul_rvec(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P1_R(z),ARRAY_P0_R(x)[0],ARRAY_P1_R(x)[0],ARRAY_P0_R(y),ARRAY_P0_R(y)); }
    if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='c'){ irscalar_mul_cvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_R(x)[0],ARRAY_P1_R(x)[0],ARRAY_P0_C(y),ARRAY_P0_C(y)); }
    if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='R'){ irscalar_mul_rvec(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P1_R(z),ARRAY_P0_R(x)[0],ARRAY_P1_R(x)[0],ARRAY_P0_R(y),ARRAY_P1_R(y)); }
    if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='C'){ idscalar_mul_cvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_D(x)[0],ARRAY_P0_D(x)[0],ARRAY_P0_C(y),ARRAY_P1_C(y)); }
    if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='C'){ izscalar_mul_cvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_Z(x)[0],ARRAY_P0_Z(x)[0],ARRAY_P0_C(y),ARRAY_P1_C(y)); }
    if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='C'){ irscalar_mul_cvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_R(x)[0],ARRAY_P0_R(x)[0],ARRAY_P0_C(y),ARRAY_P1_C(y)); }
    if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='C'){ icscalar_mul_cvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(x)[0],ARRAY_P0_C(x)[0],ARRAY_P0_C(y),ARRAY_P1_C(y)); }
    if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='C'){ irscalar_mul_cvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_R(x)[0],ARRAY_P1_R(x)[0],ARRAY_P0_C(y),ARRAY_P1_C(y)); }
    if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='d'){ icscalar_mul_dvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(x)[0],ARRAY_P1_C(x)[0],ARRAY_P0_D(y),ARRAY_P0_D(y)); }
    if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='z'){ icscalar_mul_zvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(x)[0],ARRAY_P1_C(x)[0],ARRAY_P0_Z(y),ARRAY_P0_Z(y)); }
    if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='r'){ icscalar_mul_rvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(x)[0],ARRAY_P1_C(x)[0],ARRAY_P0_R(y),ARRAY_P0_R(y)); }
    if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='c'){ icscalar_mul_cvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(x)[0],ARRAY_P1_C(x)[0],ARRAY_P0_C(y),ARRAY_P0_C(y)); }
    if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='R'){ icscalar_mul_rvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(x)[0],ARRAY_P1_C(x)[0],ARRAY_P0_R(y),ARRAY_P1_R(y)); }
    if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='C'){ icscalar_mul_cvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(x)[0],ARRAY_P1_C(x)[0],ARRAY_P0_C(y),ARRAY_P1_C(y)); }
  }
  // 互換性のあるサイズでありxが小さい場合
  if(z==NULL && (k=array_get_subdim(x,y))>0){
    z=array_allocate(array_get_inclusion_type(ARRAY_TYPE(x),ARRAY_TYPE(y)),ARRAY_NDIM(y),ARRAY_DIM_P(y));
    for(i=0; i<ARRAY_LD_END(z)/ARRAY_LD(z,k-1); i++){
      if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='d'){  dvec_mul_dvec(ARRAY_LD(z,k-1),ARRAY_P0_D(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_D(x),ARRAY_P0_D(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='z'){  dvec_mul_zvec(ARRAY_LD(z,k-1),ARRAY_P0_Z(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_D(x),ARRAY_P0_Z(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='d'){  zvec_mul_dvec(ARRAY_LD(z,k-1),ARRAY_P0_Z(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_Z(x),ARRAY_P0_D(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='z'){  zvec_mul_zvec(ARRAY_LD(z,k-1),ARRAY_P0_Z(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_Z(x),ARRAY_P0_Z(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='r'){  dvec_mul_rvec(ARRAY_LD(z,k-1),ARRAY_P0_R(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_D(x),ARRAY_P0_R(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='r'){  zvec_mul_rvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_Z(x),ARRAY_P0_R(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='d'){  rvec_mul_dvec(ARRAY_LD(z,k-1),ARRAY_P0_R(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_R(x),ARRAY_P0_D(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='z'){  rvec_mul_zvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_R(x),ARRAY_P0_Z(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='r'){  rvec_mul_rvec(ARRAY_LD(z,k-1),ARRAY_P0_R(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_R(x),ARRAY_P0_R(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='c'){  dvec_mul_cvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_D(x),ARRAY_P0_C(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='c'){  zvec_mul_cvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_Z(x),ARRAY_P0_C(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='c'){  rvec_mul_cvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_R(x),ARRAY_P0_C(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='d'){  cvec_mul_dvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_C(x),ARRAY_P0_D(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='z'){  cvec_mul_zvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_C(x),ARRAY_P0_Z(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='r'){  cvec_mul_rvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_C(x),ARRAY_P0_R(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='c'){  cvec_mul_cvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_C(x),ARRAY_P0_C(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='R'){ idvec_mul_rvec(ARRAY_LD(z,k-1),ARRAY_P0_R(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_R(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_D(x),ARRAY_P0_D(x),ARRAY_P0_R(y)+i*ARRAY_LD(y,k-1),ARRAY_P1_R(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='R'){ izvec_mul_rvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_Z(x),ARRAY_P0_Z(x),ARRAY_P0_R(y)+i*ARRAY_LD(y,k-1),ARRAY_P1_R(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='R'){ irvec_mul_rvec(ARRAY_LD(z,k-1),ARRAY_P0_R(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_R(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_R(x),ARRAY_P0_R(x),ARRAY_P0_R(y)+i*ARRAY_LD(y,k-1),ARRAY_P1_R(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='R'){ icvec_mul_rvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_C(x),ARRAY_P0_C(x),ARRAY_P0_R(y)+i*ARRAY_LD(y,k-1),ARRAY_P1_R(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='d'){ irvec_mul_dvec(ARRAY_LD(z,k-1),ARRAY_P0_R(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_R(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_R(x),ARRAY_P1_R(x),ARRAY_P0_D(y)+i*ARRAY_LD(y,k-1),ARRAY_P0_D(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='z'){ irvec_mul_zvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_R(x),ARRAY_P1_R(x),ARRAY_P0_Z(y)+i*ARRAY_LD(y,k-1),ARRAY_P0_Z(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='r'){ irvec_mul_rvec(ARRAY_LD(z,k-1),ARRAY_P0_R(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_R(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_R(x),ARRAY_P1_R(x),ARRAY_P0_R(y)+i*ARRAY_LD(y,k-1),ARRAY_P0_R(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='c'){ irvec_mul_cvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_R(x),ARRAY_P1_R(x),ARRAY_P0_C(y)+i*ARRAY_LD(y,k-1),ARRAY_P0_C(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='R'){ irvec_mul_rvec(ARRAY_LD(z,k-1),ARRAY_P0_R(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_R(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_R(x),ARRAY_P1_R(x),ARRAY_P0_R(y)+i*ARRAY_LD(y,k-1),ARRAY_P1_R(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='C'){ idvec_mul_cvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_D(x),ARRAY_P0_D(x),ARRAY_P0_C(y)+i*ARRAY_LD(y,k-1),ARRAY_P1_C(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='C'){ izvec_mul_cvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_Z(x),ARRAY_P0_Z(x),ARRAY_P0_C(y)+i*ARRAY_LD(y,k-1),ARRAY_P1_C(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='C'){ irvec_mul_cvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_R(x),ARRAY_P0_R(x),ARRAY_P0_C(y)+i*ARRAY_LD(y,k-1),ARRAY_P1_C(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='C'){ icvec_mul_cvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_C(x),ARRAY_P0_C(x),ARRAY_P0_C(y)+i*ARRAY_LD(y,k-1),ARRAY_P1_C(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='C'){ irvec_mul_cvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_R(x),ARRAY_P1_R(x),ARRAY_P0_C(y)+i*ARRAY_LD(y,k-1),ARRAY_P1_C(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='d'){ icvec_mul_dvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_C(x),ARRAY_P1_C(x),ARRAY_P0_D(y)+i*ARRAY_LD(y,k-1),ARRAY_P0_D(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='z'){ icvec_mul_zvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_C(x),ARRAY_P1_C(x),ARRAY_P0_Z(y)+i*ARRAY_LD(y,k-1),ARRAY_P0_Z(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='r'){ icvec_mul_rvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_C(x),ARRAY_P1_C(x),ARRAY_P0_R(y)+i*ARRAY_LD(y,k-1),ARRAY_P0_R(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='c'){ icvec_mul_cvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_C(x),ARRAY_P1_C(x),ARRAY_P0_C(y)+i*ARRAY_LD(y,k-1),ARRAY_P0_C(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='R'){ icvec_mul_rvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_C(x),ARRAY_P1_C(x),ARRAY_P0_R(y)+i*ARRAY_LD(y,k-1),ARRAY_P1_R(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='C'){ icvec_mul_cvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_C(x),ARRAY_P1_C(x),ARRAY_P0_C(y)+i*ARRAY_LD(y,k-1),ARRAY_P1_C(y)+i*ARRAY_LD(y,k-1)); }
    }
  }
  // 互換性のあるサイズでありyが小さい場合
  if(z==NULL && (k=array_get_subdim(y,x))>0){
    z=array_allocate(array_get_inclusion_type(ARRAY_TYPE(x),ARRAY_TYPE(y)),ARRAY_NDIM(x),ARRAY_DIM_P(x));
    for(i=0; i<ARRAY_LD_END(z)/ARRAY_LD(z,k-1); i++){
      if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='d'){  dvec_mul_dvec(ARRAY_LD(z,k-1),ARRAY_P0_D(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_D(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_D(y)); }
      if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='z'){  dvec_mul_zvec(ARRAY_LD(z,k-1),ARRAY_P0_Z(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_D(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_Z(y)); }
      if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='d'){  zvec_mul_dvec(ARRAY_LD(z,k-1),ARRAY_P0_Z(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_Z(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_D(y)); }
      if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='z'){  zvec_mul_zvec(ARRAY_LD(z,k-1),ARRAY_P0_Z(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_Z(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_Z(y)); }
      if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='r'){  dvec_mul_rvec(ARRAY_LD(z,k-1),ARRAY_P0_R(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_D(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_R(y)); }
      if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='r'){  zvec_mul_rvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_Z(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_R(y)); }
      if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='d'){  rvec_mul_dvec(ARRAY_LD(z,k-1),ARRAY_P0_R(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_R(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_D(y)); }
      if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='z'){  rvec_mul_zvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_R(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_Z(y)); }
      if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='r'){  rvec_mul_rvec(ARRAY_LD(z,k-1),ARRAY_P0_R(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_R(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_R(y)); }
      if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='c'){  dvec_mul_cvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_D(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_C(y)); }
      if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='c'){  zvec_mul_cvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_Z(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_C(y)); }
      if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='c'){  rvec_mul_cvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_R(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_C(y)); }
      if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='d'){  cvec_mul_dvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_C(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_D(y)); }
      if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='z'){  cvec_mul_zvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_C(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_Z(y)); }
      if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='r'){  cvec_mul_rvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_C(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_R(y)); }
      if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='c'){  cvec_mul_cvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_C(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_C(y)); }
      if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='R'){ idvec_mul_rvec(ARRAY_LD(z,k-1),ARRAY_P0_R(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_R(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_D(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_D(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_R(y),ARRAY_P1_R(y)); }
      if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='R'){ izvec_mul_rvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_Z(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_Z(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_R(y),ARRAY_P1_R(y)); }
      if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='R'){ irvec_mul_rvec(ARRAY_LD(z,k-1),ARRAY_P0_R(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_R(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_R(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_R(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_R(y),ARRAY_P1_R(y)); }
      if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='R'){ icvec_mul_rvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_C(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_C(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_R(y),ARRAY_P1_R(y)); }
      if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='d'){ irvec_mul_dvec(ARRAY_LD(z,k-1),ARRAY_P0_R(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_R(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_R(x)+i*ARRAY_LD(x,k-1),ARRAY_P1_R(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_D(y),ARRAY_P0_D(y)); }
      if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='z'){ irvec_mul_zvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_R(x)+i*ARRAY_LD(x,k-1),ARRAY_P1_R(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_Z(y),ARRAY_P0_Z(y)); }
      if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='r'){ irvec_mul_rvec(ARRAY_LD(z,k-1),ARRAY_P0_R(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_R(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_R(x)+i*ARRAY_LD(x,k-1),ARRAY_P1_R(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_R(y),ARRAY_P0_R(y)); }
      if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='c'){ irvec_mul_cvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_R(x)+i*ARRAY_LD(x,k-1),ARRAY_P1_R(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_C(y),ARRAY_P0_C(y)); }
      if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='R'){ irvec_mul_rvec(ARRAY_LD(z,k-1),ARRAY_P0_R(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_R(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_R(x)+i*ARRAY_LD(x,k-1),ARRAY_P1_R(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_R(y),ARRAY_P1_R(y)); }
      if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='C'){ idvec_mul_cvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_D(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_D(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_C(y),ARRAY_P1_C(y)); }
      if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='C'){ izvec_mul_cvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_Z(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_Z(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_C(y),ARRAY_P1_C(y)); }
      if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='C'){ irvec_mul_cvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_R(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_R(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_C(y),ARRAY_P1_C(y)); }
      if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='C'){ icvec_mul_cvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_C(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_C(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_C(y),ARRAY_P1_C(y)); }
      if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='C'){ irvec_mul_cvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_R(x)+i*ARRAY_LD(x,k-1),ARRAY_P1_R(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_C(y),ARRAY_P1_C(y)); }
      if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='d'){ icvec_mul_dvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_C(x)+i*ARRAY_LD(x,k-1),ARRAY_P1_C(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_D(y),ARRAY_P0_D(y)); }
      if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='z'){ icvec_mul_zvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_C(x)+i*ARRAY_LD(x,k-1),ARRAY_P1_C(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_Z(y),ARRAY_P0_Z(y)); }
      if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='r'){ icvec_mul_rvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_C(x)+i*ARRAY_LD(x,k-1),ARRAY_P1_C(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_R(y),ARRAY_P0_R(y)); }
      if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='c'){ icvec_mul_cvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_C(x)+i*ARRAY_LD(x,k-1),ARRAY_P1_C(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_C(y),ARRAY_P0_C(y)); }
      if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='R'){ icvec_mul_rvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_C(x)+i*ARRAY_LD(x,k-1),ARRAY_P1_C(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_R(y),ARRAY_P1_R(y)); }
      if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='C'){ icvec_mul_cvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_C(x)+i*ARRAY_LD(x,k-1),ARRAY_P1_C(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_C(y),ARRAY_P1_C(y)); }
    }
  }
  // 互換性のあるサイズの場合
  if(z==NULL && array_compatible_dim_check(x,y)){
    ndim=MAX2(ARRAY_NDIM(x),ARRAY_NDIM(y));
    dim=array_get_inclusion_dim(x,y);
    z=array_allocate(array_get_inclusion_type(ARRAY_TYPE(x),ARRAY_TYPE(y)),ndim,dim);
    // loop
    I=ivec_allocate(ARRAY_NDIM(z));
    K=ivec_allocate(ARRAY_NDIM(z));
    L=ivec_allocate(ARRAY_NDIM(z));
    for(i=0; i<ARRAY_SIZE(z); i++){
      array_get_index(I,z,i);
      ivec_copy(ndim,K,I);
      ivec_copy(ndim,L,I);
      for(j=0; j<ndim; j++){
        if(K[j]>=ARRAY_DIM(x,j)){ K[j]=ARRAY_DIM(x,j)-1; }
        if(L[j]>=ARRAY_DIM(y,j)){ L[j]=ARRAY_DIM(y,j)-1; }
      }
      k=array_get_position(K,x);
      l=array_get_position(L,y);
      if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='d'){ ARRAY_DVEC0(z,i)=ARRAY_DVEC0(x,k)*ARRAY_DVEC0(y,l); }
      if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='z'){ ARRAY_ZVEC0(z,i)=zmul_dz(ARRAY_DVEC0(x,k),ARRAY_ZVEC0(y,l)); }
      if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='d'){ ARRAY_ZVEC0(z,i)=zmul_zd(ARRAY_ZVEC0(x,k),ARRAY_DVEC0(y,l)); }
      if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='z'){ ARRAY_ZVEC0(z,i)=zmul_zz(ARRAY_ZVEC0(x,k),ARRAY_ZVEC0(y,l)); }
      if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='r'){  rmul_dr(ARRAY_RVEC0(z,i),ARRAY_DVEC0(x,k),ARRAY_RVEC0(y,l)); }
      if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='r'){  cmul_zr(ARRAY_CVEC0(z,i),ARRAY_ZVEC0(x,k),ARRAY_RVEC0(y,l)); }
      if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='d'){  rmul_rd(ARRAY_RVEC0(z,i),ARRAY_RVEC0(x,k),ARRAY_DVEC0(y,l)); }
      if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='z'){  cmul_rz(ARRAY_CVEC0(z,i),ARRAY_RVEC0(x,k),ARRAY_ZVEC0(y,l)); }
      if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='r'){  rmul_rr(ARRAY_RVEC0(z,i),ARRAY_RVEC0(x,k),ARRAY_RVEC0(y,l)); }
      if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='c'){  cmul_dc(ARRAY_CVEC0(z,i),ARRAY_DVEC0(x,k),ARRAY_CVEC0(y,l)); }
      if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='c'){  cmul_zc(ARRAY_CVEC0(z,i),ARRAY_ZVEC0(x,k),ARRAY_CVEC0(y,l)); }
      if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='c'){  cmul_rc(ARRAY_CVEC0(z,i),ARRAY_RVEC0(x,k),ARRAY_CVEC0(y,l)); }
      if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='d'){  cmul_cd(ARRAY_CVEC0(z,i),ARRAY_CVEC0(x,k),ARRAY_DVEC0(y,l)); }
      if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='z'){  cmul_cz(ARRAY_CVEC0(z,i),ARRAY_CVEC0(x,k),ARRAY_ZVEC0(y,l)); }
      if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='r'){  cmul_cr(ARRAY_CVEC0(z,i),ARRAY_CVEC0(x,k),ARRAY_RVEC0(y,l)); }
      if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='c'){  cmul_cc(ARRAY_CVEC0(z,i),ARRAY_CVEC0(x,k),ARRAY_CVEC0(y,l)); }
      if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='R'){ irmul_dr(ARRAY_RVEC0(z,i),ARRAY_RVEC1(z,i),ARRAY_DVEC0(x,k),ARRAY_DVEC0(x,k),ARRAY_RVEC0(y,l),ARRAY_RVEC1(y,l)); }
      if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='R'){ icmul_zr(ARRAY_CVEC0(z,i),ARRAY_CVEC1(z,i),ARRAY_ZVEC0(x,k),ARRAY_ZVEC0(x,k),ARRAY_RVEC0(y,l),ARRAY_RVEC1(y,l)); }
      if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='R'){ irmul_rr(ARRAY_RVEC0(z,i),ARRAY_RVEC1(z,i),ARRAY_RVEC0(x,k),ARRAY_RVEC0(x,k),ARRAY_RVEC0(y,l),ARRAY_RVEC1(y,l)); }
      if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='R'){ icmul_cr(ARRAY_CVEC0(z,i),ARRAY_CVEC1(z,i),ARRAY_CVEC0(x,k),ARRAY_CVEC0(x,k),ARRAY_RVEC0(y,l),ARRAY_RVEC1(y,l)); }
      if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='d'){ irmul_rd(ARRAY_RVEC0(z,i),ARRAY_RVEC1(z,i),ARRAY_RVEC0(x,k),ARRAY_RVEC1(x,k),ARRAY_DVEC0(y,l),ARRAY_DVEC0(y,l)); }
      if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='z'){ icmul_rz(ARRAY_CVEC0(z,i),ARRAY_CVEC1(z,i),ARRAY_RVEC0(x,k),ARRAY_RVEC1(x,k),ARRAY_ZVEC0(y,l),ARRAY_ZVEC0(y,l)); }
      if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='r'){ irmul_rr(ARRAY_RVEC0(z,i),ARRAY_RVEC1(z,i),ARRAY_RVEC0(x,k),ARRAY_RVEC1(x,k),ARRAY_RVEC0(y,l),ARRAY_RVEC0(y,l)); }
      if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='c'){ icmul_rc(ARRAY_CVEC0(z,i),ARRAY_CVEC1(z,i),ARRAY_RVEC0(x,k),ARRAY_RVEC1(x,k),ARRAY_CVEC0(y,l),ARRAY_CVEC0(y,l)); }
      if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='R'){ irmul_rr(ARRAY_RVEC0(z,i),ARRAY_RVEC1(z,i),ARRAY_RVEC0(x,k),ARRAY_RVEC1(x,k),ARRAY_RVEC0(y,l),ARRAY_RVEC1(y,l)); }
      if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='C'){ icmul_dc(ARRAY_CVEC0(z,i),ARRAY_CVEC1(z,i),ARRAY_DVEC0(x,k),ARRAY_DVEC0(x,k),ARRAY_CVEC0(y,l),ARRAY_CVEC1(y,l)); }
      if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='C'){ icmul_zc(ARRAY_CVEC0(z,i),ARRAY_CVEC1(z,i),ARRAY_ZVEC0(x,k),ARRAY_ZVEC0(x,k),ARRAY_CVEC0(y,l),ARRAY_CVEC1(y,l)); }
      if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='C'){ icmul_rc(ARRAY_CVEC0(z,i),ARRAY_CVEC1(z,i),ARRAY_RVEC0(x,k),ARRAY_RVEC0(x,k),ARRAY_CVEC0(y,l),ARRAY_CVEC1(y,l)); }
      if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='C'){ icmul_cc(ARRAY_CVEC0(z,i),ARRAY_CVEC1(z,i),ARRAY_CVEC0(x,k),ARRAY_CVEC0(x,k),ARRAY_CVEC0(y,l),ARRAY_CVEC1(y,l)); }
      if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='C'){ icmul_rc(ARRAY_CVEC0(z,i),ARRAY_CVEC1(z,i),ARRAY_RVEC0(x,k),ARRAY_RVEC1(x,k),ARRAY_CVEC0(y,l),ARRAY_CVEC1(y,l)); }
      if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='d'){ icmul_cd(ARRAY_CVEC0(z,i),ARRAY_CVEC1(z,i),ARRAY_CVEC0(x,k),ARRAY_CVEC1(x,k),ARRAY_DVEC0(y,l),ARRAY_DVEC0(y,l)); }
      if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='z'){ icmul_cz(ARRAY_CVEC0(z,i),ARRAY_CVEC1(z,i),ARRAY_CVEC0(x,k),ARRAY_CVEC1(x,k),ARRAY_ZVEC0(y,l),ARRAY_ZVEC0(y,l)); }
      if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='r'){ icmul_cr(ARRAY_CVEC0(z,i),ARRAY_CVEC1(z,i),ARRAY_CVEC0(x,k),ARRAY_CVEC1(x,k),ARRAY_RVEC0(y,l),ARRAY_RVEC0(y,l)); }
      if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='c'){ icmul_cc(ARRAY_CVEC0(z,i),ARRAY_CVEC1(z,i),ARRAY_CVEC0(x,k),ARRAY_CVEC1(x,k),ARRAY_CVEC0(y,l),ARRAY_CVEC0(y,l)); }
      if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='R'){ icmul_cr(ARRAY_CVEC0(z,i),ARRAY_CVEC1(z,i),ARRAY_CVEC0(x,k),ARRAY_CVEC1(x,k),ARRAY_RVEC0(y,l),ARRAY_RVEC1(y,l)); }
      if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='C'){ icmul_cc(ARRAY_CVEC0(z,i),ARRAY_CVEC1(z,i),ARRAY_CVEC0(x,k),ARRAY_CVEC1(x,k),ARRAY_CVEC0(y,l),ARRAY_CVEC1(y,l)); }
    }
    dim=ivec_free(dim);
    I=ivec_free(I);
    K=ivec_free(K);
    L=ivec_free(L);
  }
  return z;
}

/**
 @brief z=x./y
 */
array *array_div(array *x, array *y)
{
  int i,j,k,l,ndim,*dim=NULL,*I=NULL,*K=NULL,*L=NULL;
  array *z=NULL;
  if(x==NULL || y==NULL){ return z; }
  // サイズが同じ場合
  if(z==NULL && array_same_dim_check(x,y))
  {
    z=array_allocate(array_get_inclusion_type(ARRAY_TYPE(x),ARRAY_TYPE(y)),ARRAY_NDIM(x),ARRAY_DIM_P(x));
    if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='d'){  dvec_div_dvec(ARRAY_SIZE(z),ARRAY_P0_D(z),ARRAY_P0_D(x),ARRAY_P0_D(y)); }
    if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='z'){  dvec_div_zvec(ARRAY_SIZE(z),ARRAY_P0_Z(z),ARRAY_P0_D(x),ARRAY_P0_Z(y)); }
    if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='d'){  zvec_div_dvec(ARRAY_SIZE(z),ARRAY_P0_Z(z),ARRAY_P0_Z(x),ARRAY_P0_D(y)); }
    if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='z'){  zvec_div_zvec(ARRAY_SIZE(z),ARRAY_P0_Z(z),ARRAY_P0_Z(x),ARRAY_P0_Z(y)); }
    if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='r'){  dvec_div_rvec(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P0_D(x),ARRAY_P0_R(y)); }
    if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='r'){  zvec_div_rvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_Z(x),ARRAY_P0_R(y)); }
    if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='d'){  rvec_div_dvec(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P0_R(x),ARRAY_P0_D(y)); }
    if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='z'){  rvec_div_zvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_R(x),ARRAY_P0_Z(y)); }
    if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='r'){  rvec_div_rvec(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P0_R(x),ARRAY_P0_R(y)); }
    if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='c'){  dvec_div_cvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_D(x),ARRAY_P0_C(y)); }
    if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='c'){  zvec_div_cvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_Z(x),ARRAY_P0_C(y)); }
    if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='c'){  rvec_div_cvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_R(x),ARRAY_P0_C(y)); }
    if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='d'){  cvec_div_dvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_C(x),ARRAY_P0_D(y)); }
    if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='z'){  cvec_div_zvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_C(x),ARRAY_P0_Z(y)); }
    if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='r'){  cvec_div_rvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_C(x),ARRAY_P0_R(y)); }
    if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='c'){  cvec_div_cvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_C(x),ARRAY_P0_C(y)); }
    if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='R'){ idvec_div_rvec(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P1_R(z),ARRAY_P0_D(x),ARRAY_P0_D(x),ARRAY_P0_R(y),ARRAY_P1_R(y)); }
    if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='R'){ izvec_div_rvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_Z(x),ARRAY_P0_Z(x),ARRAY_P0_R(y),ARRAY_P1_R(y)); }
    if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='R'){ irvec_div_rvec(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P1_R(z),ARRAY_P0_R(x),ARRAY_P0_R(x),ARRAY_P0_R(y),ARRAY_P1_R(y)); }
    if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='R'){ icvec_div_rvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(x),ARRAY_P0_C(x),ARRAY_P0_R(y),ARRAY_P1_R(y)); }
    if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='d'){ irvec_div_dvec(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P1_R(z),ARRAY_P0_R(x),ARRAY_P1_R(x),ARRAY_P0_D(y),ARRAY_P0_D(y)); }
    if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='z'){ irvec_div_zvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_R(x),ARRAY_P1_R(x),ARRAY_P0_Z(y),ARRAY_P0_Z(y)); }
    if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='r'){ irvec_div_rvec(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P1_R(z),ARRAY_P0_R(x),ARRAY_P1_R(x),ARRAY_P0_R(y),ARRAY_P0_R(y)); }
    if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='c'){ irvec_div_cvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_R(x),ARRAY_P1_R(x),ARRAY_P0_C(y),ARRAY_P0_C(y)); }
    if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='R'){ irvec_div_rvec(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P1_R(z),ARRAY_P0_R(x),ARRAY_P1_R(x),ARRAY_P0_R(y),ARRAY_P1_R(y)); }
    if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='C'){ idvec_div_cvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_D(x),ARRAY_P0_D(x),ARRAY_P0_C(y),ARRAY_P1_C(y)); }
    if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='C'){ izvec_div_cvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_Z(x),ARRAY_P0_Z(x),ARRAY_P0_C(y),ARRAY_P1_C(y)); }
    if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='C'){ irvec_div_cvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_R(x),ARRAY_P0_R(x),ARRAY_P0_C(y),ARRAY_P1_C(y)); }
    if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='C'){ icvec_div_cvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(x),ARRAY_P0_C(x),ARRAY_P0_C(y),ARRAY_P1_C(y)); }
    if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='C'){ irvec_div_cvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_R(x),ARRAY_P1_R(x),ARRAY_P0_C(y),ARRAY_P1_C(y)); }
    if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='d'){ icvec_div_dvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(x),ARRAY_P1_C(x),ARRAY_P0_D(y),ARRAY_P0_D(y)); }
    if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='z'){ icvec_div_zvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(x),ARRAY_P1_C(x),ARRAY_P0_Z(y),ARRAY_P0_Z(y)); }
    if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='r'){ icvec_div_rvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(x),ARRAY_P1_C(x),ARRAY_P0_R(y),ARRAY_P0_R(y)); }
    if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='c'){ icvec_div_cvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(x),ARRAY_P1_C(x),ARRAY_P0_C(y),ARRAY_P0_C(y)); }
    if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='R'){ icvec_div_rvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(x),ARRAY_P1_C(x),ARRAY_P0_R(y),ARRAY_P1_R(y)); }
    if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='C'){ icvec_div_cvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(x),ARRAY_P1_C(x),ARRAY_P0_C(y),ARRAY_P1_C(y)); }
  }
  // yがスカラーの場合
  if(z==NULL && array_is_scalar(y))
  {
    z=array_allocate(array_get_inclusion_type(ARRAY_TYPE(x),ARRAY_TYPE(y)),ARRAY_NDIM(x),ARRAY_DIM_P(x));
    if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='d'){  dvec_div_dscalar(ARRAY_SIZE(z),ARRAY_P0_D(z),ARRAY_P0_D(x),ARRAY_P0_D(y)[0]); }
    if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='z'){  dvec_div_zscalar(ARRAY_SIZE(z),ARRAY_P0_Z(z),ARRAY_P0_D(x),ARRAY_P0_Z(y)[0]); }
    if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='d'){  zvec_div_dscalar(ARRAY_SIZE(z),ARRAY_P0_Z(z),ARRAY_P0_Z(x),ARRAY_P0_D(y)[0]); }
    if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='z'){  zvec_div_zscalar(ARRAY_SIZE(z),ARRAY_P0_Z(z),ARRAY_P0_Z(x),ARRAY_P0_Z(y)[0]); }
    if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='r'){  dvec_div_rscalar(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P0_D(x),ARRAY_P0_R(y)[0]); }
    if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='r'){  zvec_div_rscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_Z(x),ARRAY_P0_R(y)[0]); }
    if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='d'){  rvec_div_dscalar(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P0_R(x),ARRAY_P0_D(y)[0]); }
    if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='z'){  rvec_div_zscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_R(x),ARRAY_P0_Z(y)[0]); }
    if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='r'){  rvec_div_rscalar(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P0_R(x),ARRAY_P0_R(y)[0]); }
    if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='c'){  dvec_div_cscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_D(x),ARRAY_P0_C(y)[0]); }
    if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='c'){  zvec_div_cscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_Z(x),ARRAY_P0_C(y)[0]); }
    if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='c'){  rvec_div_cscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_R(x),ARRAY_P0_C(y)[0]); }
    if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='d'){  cvec_div_dscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_C(x),ARRAY_P0_D(y)[0]); }
    if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='z'){  cvec_div_zscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_C(x),ARRAY_P0_Z(y)[0]); }
    if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='r'){  cvec_div_rscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_C(x),ARRAY_P0_R(y)[0]); }
    if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='c'){  cvec_div_cscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_C(x),ARRAY_P0_C(y)[0]); }
    if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='R'){ idvec_div_rscalar(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P1_R(z),ARRAY_P0_D(x),ARRAY_P0_D(x),ARRAY_P0_R(y)[0],ARRAY_P1_R(y)[0]); }
    if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='R'){ izvec_div_rscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_Z(x),ARRAY_P0_Z(x),ARRAY_P0_R(y)[0],ARRAY_P1_R(y)[0]); }
    if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='R'){ irvec_div_rscalar(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P1_R(z),ARRAY_P0_R(x),ARRAY_P0_R(x),ARRAY_P0_R(y)[0],ARRAY_P1_R(y)[0]); }
    if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='R'){ icvec_div_rscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(x),ARRAY_P0_C(x),ARRAY_P0_R(y)[0],ARRAY_P1_R(y)[0]); }
    if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='d'){ irvec_div_dscalar(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P1_R(z),ARRAY_P0_R(x),ARRAY_P1_R(x),ARRAY_P0_D(y)[0],ARRAY_P0_D(y)[0]); }
    if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='z'){ irvec_div_zscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_R(x),ARRAY_P1_R(x),ARRAY_P0_Z(y)[0],ARRAY_P0_Z(y)[0]); }
    if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='r'){ irvec_div_rscalar(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P1_R(z),ARRAY_P0_R(x),ARRAY_P1_R(x),ARRAY_P0_R(y)[0],ARRAY_P0_R(y)[0]); }
    if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='c'){ irvec_div_cscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_R(x),ARRAY_P1_R(x),ARRAY_P0_C(y)[0],ARRAY_P0_C(y)[0]); }
    if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='R'){ irvec_div_rscalar(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P1_R(z),ARRAY_P0_R(x),ARRAY_P1_R(x),ARRAY_P0_R(y)[0],ARRAY_P1_R(y)[0]); }
    if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='C'){ idvec_div_cscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_D(x),ARRAY_P0_D(x),ARRAY_P0_C(y)[0],ARRAY_P1_C(y)[0]); }
    if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='C'){ izvec_div_cscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_Z(x),ARRAY_P0_Z(x),ARRAY_P0_C(y)[0],ARRAY_P1_C(y)[0]); }
    if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='C'){ irvec_div_cscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_R(x),ARRAY_P0_R(x),ARRAY_P0_C(y)[0],ARRAY_P1_C(y)[0]); }
    if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='C'){ icvec_div_cscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(x),ARRAY_P0_C(x),ARRAY_P0_C(y)[0],ARRAY_P1_C(y)[0]); }
    if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='C'){ irvec_div_cscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_R(x),ARRAY_P1_R(x),ARRAY_P0_C(y)[0],ARRAY_P1_C(y)[0]); }
    if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='d'){ icvec_div_dscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(x),ARRAY_P1_C(x),ARRAY_P0_D(y)[0],ARRAY_P0_D(y)[0]); }
    if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='z'){ icvec_div_zscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(x),ARRAY_P1_C(x),ARRAY_P0_Z(y)[0],ARRAY_P0_Z(y)[0]); }
    if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='r'){ icvec_div_rscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(x),ARRAY_P1_C(x),ARRAY_P0_R(y)[0],ARRAY_P0_R(y)[0]); }
    if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='c'){ icvec_div_cscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(x),ARRAY_P1_C(x),ARRAY_P0_C(y)[0],ARRAY_P0_C(y)[0]); }
    if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='R'){ icvec_div_rscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(x),ARRAY_P1_C(x),ARRAY_P0_R(y)[0],ARRAY_P1_R(y)[0]); }
    if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='C'){ icvec_div_cscalar(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(x),ARRAY_P1_C(x),ARRAY_P0_C(y)[0],ARRAY_P1_C(y)[0]); }
  }
  // xがスカラーの場合
  if(z==NULL && array_is_scalar(x))
  {
    z=array_allocate(array_get_inclusion_type(ARRAY_TYPE(x),ARRAY_TYPE(y)),ARRAY_NDIM(y),ARRAY_DIM_P(y));
    if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='d'){  dscalar_div_dvec(ARRAY_SIZE(z),ARRAY_P0_D(z),ARRAY_P0_D(x)[0],ARRAY_P0_D(y)); }
    if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='z'){  dscalar_div_zvec(ARRAY_SIZE(z),ARRAY_P0_Z(z),ARRAY_P0_D(x)[0],ARRAY_P0_Z(y)); }
    if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='d'){  zscalar_div_dvec(ARRAY_SIZE(z),ARRAY_P0_Z(z),ARRAY_P0_Z(x)[0],ARRAY_P0_D(y)); }
    if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='z'){  zscalar_div_zvec(ARRAY_SIZE(z),ARRAY_P0_Z(z),ARRAY_P0_Z(x)[0],ARRAY_P0_Z(y)); }
    if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='r'){  dscalar_div_rvec(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P0_D(x)[0],ARRAY_P0_R(y)); }
    if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='r'){  zscalar_div_rvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_Z(x)[0],ARRAY_P0_R(y)); }
    if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='d'){  rscalar_div_dvec(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P0_R(x)[0],ARRAY_P0_D(y)); }
    if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='z'){  rscalar_div_zvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_R(x)[0],ARRAY_P0_Z(y)); }
    if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='r'){  rscalar_div_rvec(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P0_R(x)[0],ARRAY_P0_R(y)); }
    if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='c'){  dscalar_div_cvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_D(x)[0],ARRAY_P0_C(y)); }
    if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='c'){  zscalar_div_cvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_Z(x)[0],ARRAY_P0_C(y)); }
    if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='c'){  rscalar_div_cvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_R(x)[0],ARRAY_P0_C(y)); }
    if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='d'){  cscalar_div_dvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_C(x)[0],ARRAY_P0_D(y)); }
    if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='z'){  cscalar_div_zvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_C(x)[0],ARRAY_P0_Z(y)); }
    if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='r'){  cscalar_div_rvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_C(x)[0],ARRAY_P0_R(y)); }
    if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='c'){  cscalar_div_cvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P0_C(x)[0],ARRAY_P0_C(y)); }
    if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='R'){ idscalar_div_rvec(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P1_R(z),ARRAY_P0_D(x)[0],ARRAY_P0_D(x)[0],ARRAY_P0_R(y),ARRAY_P1_R(y)); }
    if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='R'){ izscalar_div_rvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_Z(x)[0],ARRAY_P0_Z(x)[0],ARRAY_P0_R(y),ARRAY_P1_R(y)); }
    if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='R'){ irscalar_div_rvec(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P1_R(z),ARRAY_P0_R(x)[0],ARRAY_P0_R(x)[0],ARRAY_P0_R(y),ARRAY_P1_R(y)); }
    if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='R'){ icscalar_div_rvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(x)[0],ARRAY_P0_C(x)[0],ARRAY_P0_R(y),ARRAY_P1_R(y)); }
    if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='d'){ irscalar_div_dvec(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P1_R(z),ARRAY_P0_R(x)[0],ARRAY_P1_R(x)[0],ARRAY_P0_D(y),ARRAY_P0_D(y)); }
    if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='z'){ irscalar_div_zvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_R(x)[0],ARRAY_P1_R(x)[0],ARRAY_P0_Z(y),ARRAY_P0_Z(y)); }
    if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='r'){ irscalar_div_rvec(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P1_R(z),ARRAY_P0_R(x)[0],ARRAY_P1_R(x)[0],ARRAY_P0_R(y),ARRAY_P0_R(y)); }
    if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='c'){ irscalar_div_cvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_R(x)[0],ARRAY_P1_R(x)[0],ARRAY_P0_C(y),ARRAY_P0_C(y)); }
    if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='R'){ irscalar_div_rvec(ARRAY_SIZE(z),ARRAY_P0_R(z),ARRAY_P1_R(z),ARRAY_P0_R(x)[0],ARRAY_P1_R(x)[0],ARRAY_P0_R(y),ARRAY_P1_R(y)); }
    if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='C'){ idscalar_div_cvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_D(x)[0],ARRAY_P0_D(x)[0],ARRAY_P0_C(y),ARRAY_P1_C(y)); }
    if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='C'){ izscalar_div_cvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_Z(x)[0],ARRAY_P0_Z(x)[0],ARRAY_P0_C(y),ARRAY_P1_C(y)); }
    if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='C'){ irscalar_div_cvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_R(x)[0],ARRAY_P0_R(x)[0],ARRAY_P0_C(y),ARRAY_P1_C(y)); }
    if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='C'){ icscalar_div_cvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(x)[0],ARRAY_P0_C(x)[0],ARRAY_P0_C(y),ARRAY_P1_C(y)); }
    if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='C'){ irscalar_div_cvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_R(x)[0],ARRAY_P1_R(x)[0],ARRAY_P0_C(y),ARRAY_P1_C(y)); }
    if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='d'){ icscalar_div_dvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(x)[0],ARRAY_P1_C(x)[0],ARRAY_P0_D(y),ARRAY_P0_D(y)); }
    if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='z'){ icscalar_div_zvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(x)[0],ARRAY_P1_C(x)[0],ARRAY_P0_Z(y),ARRAY_P0_Z(y)); }
    if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='r'){ icscalar_div_rvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(x)[0],ARRAY_P1_C(x)[0],ARRAY_P0_R(y),ARRAY_P0_R(y)); }
    if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='c'){ icscalar_div_cvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(x)[0],ARRAY_P1_C(x)[0],ARRAY_P0_C(y),ARRAY_P0_C(y)); }
    if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='R'){ icscalar_div_rvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(x)[0],ARRAY_P1_C(x)[0],ARRAY_P0_R(y),ARRAY_P1_R(y)); }
    if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='C'){ icscalar_div_cvec(ARRAY_SIZE(z),ARRAY_P0_C(z),ARRAY_P1_C(z),ARRAY_P0_C(x)[0],ARRAY_P1_C(x)[0],ARRAY_P0_C(y),ARRAY_P1_C(y)); }
  }
  // 互換性のあるサイズでありxが小さい場合
  if(z==NULL && (k=array_get_subdim(x,y))>0){
    z=array_allocate(array_get_inclusion_type(ARRAY_TYPE(x),ARRAY_TYPE(y)),ARRAY_NDIM(y),ARRAY_DIM_P(y));
    for(i=0; i<ARRAY_LD_END(z)/ARRAY_LD(z,k-1); i++){
      if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='d'){  dvec_div_dvec(ARRAY_LD(z,k-1),ARRAY_P0_D(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_D(x),ARRAY_P0_D(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='z'){  dvec_div_zvec(ARRAY_LD(z,k-1),ARRAY_P0_Z(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_D(x),ARRAY_P0_Z(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='d'){  zvec_div_dvec(ARRAY_LD(z,k-1),ARRAY_P0_Z(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_Z(x),ARRAY_P0_D(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='z'){  zvec_div_zvec(ARRAY_LD(z,k-1),ARRAY_P0_Z(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_Z(x),ARRAY_P0_Z(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='r'){  dvec_div_rvec(ARRAY_LD(z,k-1),ARRAY_P0_R(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_D(x),ARRAY_P0_R(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='r'){  zvec_div_rvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_Z(x),ARRAY_P0_R(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='d'){  rvec_div_dvec(ARRAY_LD(z,k-1),ARRAY_P0_R(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_R(x),ARRAY_P0_D(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='z'){  rvec_div_zvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_R(x),ARRAY_P0_Z(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='r'){  rvec_div_rvec(ARRAY_LD(z,k-1),ARRAY_P0_R(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_R(x),ARRAY_P0_R(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='c'){  dvec_div_cvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_D(x),ARRAY_P0_C(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='c'){  zvec_div_cvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_Z(x),ARRAY_P0_C(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='c'){  rvec_div_cvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_R(x),ARRAY_P0_C(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='d'){  cvec_div_dvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_C(x),ARRAY_P0_D(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='z'){  cvec_div_zvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_C(x),ARRAY_P0_Z(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='r'){  cvec_div_rvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_C(x),ARRAY_P0_R(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='c'){  cvec_div_cvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_C(x),ARRAY_P0_C(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='R'){ idvec_div_rvec(ARRAY_LD(z,k-1),ARRAY_P0_R(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_R(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_D(x),ARRAY_P0_D(x),ARRAY_P0_R(y)+i*ARRAY_LD(y,k-1),ARRAY_P1_R(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='R'){ izvec_div_rvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_Z(x),ARRAY_P0_Z(x),ARRAY_P0_R(y)+i*ARRAY_LD(y,k-1),ARRAY_P1_R(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='R'){ irvec_div_rvec(ARRAY_LD(z,k-1),ARRAY_P0_R(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_R(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_R(x),ARRAY_P0_R(x),ARRAY_P0_R(y)+i*ARRAY_LD(y,k-1),ARRAY_P1_R(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='R'){ icvec_div_rvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_C(x),ARRAY_P0_C(x),ARRAY_P0_R(y)+i*ARRAY_LD(y,k-1),ARRAY_P1_R(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='d'){ irvec_div_dvec(ARRAY_LD(z,k-1),ARRAY_P0_R(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_R(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_R(x),ARRAY_P1_R(x),ARRAY_P0_D(y)+i*ARRAY_LD(y,k-1),ARRAY_P0_D(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='z'){ irvec_div_zvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_R(x),ARRAY_P1_R(x),ARRAY_P0_Z(y)+i*ARRAY_LD(y,k-1),ARRAY_P0_Z(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='r'){ irvec_div_rvec(ARRAY_LD(z,k-1),ARRAY_P0_R(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_R(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_R(x),ARRAY_P1_R(x),ARRAY_P0_R(y)+i*ARRAY_LD(y,k-1),ARRAY_P0_R(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='c'){ irvec_div_cvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_R(x),ARRAY_P1_R(x),ARRAY_P0_C(y)+i*ARRAY_LD(y,k-1),ARRAY_P0_C(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='R'){ irvec_div_rvec(ARRAY_LD(z,k-1),ARRAY_P0_R(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_R(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_R(x),ARRAY_P1_R(x),ARRAY_P0_R(y)+i*ARRAY_LD(y,k-1),ARRAY_P1_R(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='C'){ idvec_div_cvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_D(x),ARRAY_P0_D(x),ARRAY_P0_C(y)+i*ARRAY_LD(y,k-1),ARRAY_P1_C(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='C'){ izvec_div_cvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_Z(x),ARRAY_P0_Z(x),ARRAY_P0_C(y)+i*ARRAY_LD(y,k-1),ARRAY_P1_C(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='C'){ irvec_div_cvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_R(x),ARRAY_P0_R(x),ARRAY_P0_C(y)+i*ARRAY_LD(y,k-1),ARRAY_P1_C(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='C'){ icvec_div_cvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_C(x),ARRAY_P0_C(x),ARRAY_P0_C(y)+i*ARRAY_LD(y,k-1),ARRAY_P1_C(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='C'){ irvec_div_cvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_R(x),ARRAY_P1_R(x),ARRAY_P0_C(y)+i*ARRAY_LD(y,k-1),ARRAY_P1_C(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='d'){ icvec_div_dvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_C(x),ARRAY_P1_C(x),ARRAY_P0_D(y)+i*ARRAY_LD(y,k-1),ARRAY_P0_D(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='z'){ icvec_div_zvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_C(x),ARRAY_P1_C(x),ARRAY_P0_Z(y)+i*ARRAY_LD(y,k-1),ARRAY_P0_Z(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='r'){ icvec_div_rvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_C(x),ARRAY_P1_C(x),ARRAY_P0_R(y)+i*ARRAY_LD(y,k-1),ARRAY_P0_R(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='c'){ icvec_div_cvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_C(x),ARRAY_P1_C(x),ARRAY_P0_C(y)+i*ARRAY_LD(y,k-1),ARRAY_P0_C(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='R'){ icvec_div_rvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_C(x),ARRAY_P1_C(x),ARRAY_P0_R(y)+i*ARRAY_LD(y,k-1),ARRAY_P1_R(y)+i*ARRAY_LD(y,k-1)); }
      if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='C'){ icvec_div_cvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_C(x),ARRAY_P1_C(x),ARRAY_P0_C(y)+i*ARRAY_LD(y,k-1),ARRAY_P1_C(y)+i*ARRAY_LD(y,k-1)); }
    }
  }
  // 互換性のあるサイズでありyが小さい場合
  if(z==NULL && (k=array_get_subdim(y,x))>0){
    z=array_allocate(array_get_inclusion_type(ARRAY_TYPE(x),ARRAY_TYPE(y)),ARRAY_NDIM(x),ARRAY_DIM_P(x));
    for(i=0; i<ARRAY_LD_END(z)/ARRAY_LD(z,k-1); i++){
      if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='d'){  dvec_div_dvec(ARRAY_LD(z,k-1),ARRAY_P0_D(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_D(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_D(y)); }
      if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='z'){  dvec_div_zvec(ARRAY_LD(z,k-1),ARRAY_P0_Z(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_D(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_Z(y)); }
      if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='d'){  zvec_div_dvec(ARRAY_LD(z,k-1),ARRAY_P0_Z(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_Z(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_D(y)); }
      if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='z'){  zvec_div_zvec(ARRAY_LD(z,k-1),ARRAY_P0_Z(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_Z(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_Z(y)); }
      if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='r'){  dvec_div_rvec(ARRAY_LD(z,k-1),ARRAY_P0_R(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_D(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_R(y)); }
      if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='r'){  zvec_div_rvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_Z(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_R(y)); }
      if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='d'){  rvec_div_dvec(ARRAY_LD(z,k-1),ARRAY_P0_R(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_R(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_D(y)); }
      if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='z'){  rvec_div_zvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_R(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_Z(y)); }
      if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='r'){  rvec_div_rvec(ARRAY_LD(z,k-1),ARRAY_P0_R(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_R(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_R(y)); }
      if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='c'){  dvec_div_cvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_D(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_C(y)); }
      if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='c'){  zvec_div_cvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_Z(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_C(y)); }
      if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='c'){  rvec_div_cvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_R(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_C(y)); }
      if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='d'){  cvec_div_dvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_C(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_D(y)); }
      if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='z'){  cvec_div_zvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_C(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_Z(y)); }
      if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='r'){  cvec_div_rvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_C(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_R(y)); }
      if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='c'){  cvec_div_cvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_C(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_C(y)); }
      if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='R'){ idvec_div_rvec(ARRAY_LD(z,k-1),ARRAY_P0_R(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_R(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_D(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_D(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_R(y),ARRAY_P1_R(y)); }
      if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='R'){ izvec_div_rvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_Z(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_Z(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_R(y),ARRAY_P1_R(y)); }
      if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='R'){ irvec_div_rvec(ARRAY_LD(z,k-1),ARRAY_P0_R(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_R(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_R(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_R(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_R(y),ARRAY_P1_R(y)); }
      if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='R'){ icvec_div_rvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_C(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_C(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_R(y),ARRAY_P1_R(y)); }
      if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='d'){ irvec_div_dvec(ARRAY_LD(z,k-1),ARRAY_P0_R(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_R(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_R(x)+i*ARRAY_LD(x,k-1),ARRAY_P1_R(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_D(y),ARRAY_P0_D(y)); }
      if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='z'){ irvec_div_zvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_R(x)+i*ARRAY_LD(x,k-1),ARRAY_P1_R(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_Z(y),ARRAY_P0_Z(y)); }
      if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='r'){ irvec_div_rvec(ARRAY_LD(z,k-1),ARRAY_P0_R(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_R(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_R(x)+i*ARRAY_LD(x,k-1),ARRAY_P1_R(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_R(y),ARRAY_P0_R(y)); }
      if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='c'){ irvec_div_cvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_R(x)+i*ARRAY_LD(x,k-1),ARRAY_P1_R(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_C(y),ARRAY_P0_C(y)); }
      if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='R'){ irvec_div_rvec(ARRAY_LD(z,k-1),ARRAY_P0_R(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_R(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_R(x)+i*ARRAY_LD(x,k-1),ARRAY_P1_R(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_R(y),ARRAY_P1_R(y)); }
      if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='C'){ idvec_div_cvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_D(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_D(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_C(y),ARRAY_P1_C(y)); }
      if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='C'){ izvec_div_cvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_Z(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_Z(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_C(y),ARRAY_P1_C(y)); }
      if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='C'){ irvec_div_cvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_R(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_R(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_C(y),ARRAY_P1_C(y)); }
      if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='C'){ icvec_div_cvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_C(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_C(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_C(y),ARRAY_P1_C(y)); }
      if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='C'){ irvec_div_cvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_R(x)+i*ARRAY_LD(x,k-1),ARRAY_P1_R(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_C(y),ARRAY_P1_C(y)); }
      if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='d'){ icvec_div_dvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_C(x)+i*ARRAY_LD(x,k-1),ARRAY_P1_C(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_D(y),ARRAY_P0_D(y)); }
      if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='z'){ icvec_div_zvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_C(x)+i*ARRAY_LD(x,k-1),ARRAY_P1_C(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_Z(y),ARRAY_P0_Z(y)); }
      if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='r'){ icvec_div_rvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_C(x)+i*ARRAY_LD(x,k-1),ARRAY_P1_C(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_R(y),ARRAY_P0_R(y)); }
      if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='c'){ icvec_div_cvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_C(x)+i*ARRAY_LD(x,k-1),ARRAY_P1_C(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_C(y),ARRAY_P0_C(y)); }
      if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='R'){ icvec_div_rvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_C(x)+i*ARRAY_LD(x,k-1),ARRAY_P1_C(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_R(y),ARRAY_P1_R(y)); }
      if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='C'){ icvec_div_cvec(ARRAY_LD(z,k-1),ARRAY_P0_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P1_C(z)+i*ARRAY_LD(z,k-1),ARRAY_P0_C(x)+i*ARRAY_LD(x,k-1),ARRAY_P1_C(x)+i*ARRAY_LD(x,k-1),ARRAY_P0_C(y),ARRAY_P1_C(y)); }
    }
  }
  // 互換性のあるサイズの場合
  if(z==NULL && array_compatible_dim_check(x,y)){
    ndim=MAX2(ARRAY_NDIM(x),ARRAY_NDIM(y));
    dim=array_get_inclusion_dim(x,y);
    z=array_allocate(array_get_inclusion_type(ARRAY_TYPE(x),ARRAY_TYPE(y)),ndim,dim);
    // loop
    I=ivec_allocate(ARRAY_NDIM(z));
    K=ivec_allocate(ARRAY_NDIM(z));
    L=ivec_allocate(ARRAY_NDIM(z));
    for(i=0; i<ARRAY_SIZE(z); i++){
      array_get_index(I,z,i);
      ivec_copy(ndim,K,I);
      ivec_copy(ndim,L,I);
      for(j=0; j<ndim; j++){
        if(K[j]>=ARRAY_DIM(x,j)){ K[j]=ARRAY_DIM(x,j)-1; }
        if(L[j]>=ARRAY_DIM(y,j)){ L[j]=ARRAY_DIM(y,j)-1; }
      }
      k=array_get_position(K,x);
      l=array_get_position(L,y);
      if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='d'){ ARRAY_DVEC0(z,i)=ARRAY_DVEC0(x,k)/ARRAY_DVEC0(y,l); }
      if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='z'){ ARRAY_ZVEC0(z,i)=zdiv_dz(ARRAY_DVEC0(x,k),ARRAY_ZVEC0(y,l)); }
      if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='d'){ ARRAY_ZVEC0(z,i)=zdiv_zd(ARRAY_ZVEC0(x,k),ARRAY_DVEC0(y,l)); }
      if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='z'){ ARRAY_ZVEC0(z,i)=zdiv_zz(ARRAY_ZVEC0(x,k),ARRAY_ZVEC0(y,l)); }
      if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='r'){  rdiv_dr(ARRAY_RVEC0(z,i),ARRAY_DVEC0(x,k),ARRAY_RVEC0(y,l)); }
      if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='r'){  cdiv_zr(ARRAY_CVEC0(z,i),ARRAY_ZVEC0(x,k),ARRAY_RVEC0(y,l)); }
      if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='d'){  rdiv_rd(ARRAY_RVEC0(z,i),ARRAY_RVEC0(x,k),ARRAY_DVEC0(y,l)); }
      if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='z'){  cdiv_rz(ARRAY_CVEC0(z,i),ARRAY_RVEC0(x,k),ARRAY_ZVEC0(y,l)); }
      if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='r'){  rdiv_rr(ARRAY_RVEC0(z,i),ARRAY_RVEC0(x,k),ARRAY_RVEC0(y,l)); }
      if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='c'){  cdiv_dc(ARRAY_CVEC0(z,i),ARRAY_DVEC0(x,k),ARRAY_CVEC0(y,l)); }
      if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='c'){  cdiv_zc(ARRAY_CVEC0(z,i),ARRAY_ZVEC0(x,k),ARRAY_CVEC0(y,l)); }
      if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='c'){  cdiv_rc(ARRAY_CVEC0(z,i),ARRAY_RVEC0(x,k),ARRAY_CVEC0(y,l)); }
      if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='d'){  cdiv_cd(ARRAY_CVEC0(z,i),ARRAY_CVEC0(x,k),ARRAY_DVEC0(y,l)); }
      if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='z'){  cdiv_cz(ARRAY_CVEC0(z,i),ARRAY_CVEC0(x,k),ARRAY_ZVEC0(y,l)); }
      if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='r'){  cdiv_cr(ARRAY_CVEC0(z,i),ARRAY_CVEC0(x,k),ARRAY_RVEC0(y,l)); }
      if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='c'){  cdiv_cc(ARRAY_CVEC0(z,i),ARRAY_CVEC0(x,k),ARRAY_CVEC0(y,l)); }
      if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='R'){ irdiv_dr(ARRAY_RVEC0(z,i),ARRAY_RVEC1(z,i),ARRAY_DVEC0(x,k),ARRAY_DVEC0(x,k),ARRAY_RVEC0(y,l),ARRAY_RVEC1(y,l)); }
      if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='R'){ icdiv_zr(ARRAY_CVEC0(z,i),ARRAY_CVEC1(z,i),ARRAY_ZVEC0(x,k),ARRAY_ZVEC0(x,k),ARRAY_RVEC0(y,l),ARRAY_RVEC1(y,l)); }
      if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='R'){ irdiv_rr(ARRAY_RVEC0(z,i),ARRAY_RVEC1(z,i),ARRAY_RVEC0(x,k),ARRAY_RVEC0(x,k),ARRAY_RVEC0(y,l),ARRAY_RVEC1(y,l)); }
      if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='R'){ icdiv_cr(ARRAY_CVEC0(z,i),ARRAY_CVEC1(z,i),ARRAY_CVEC0(x,k),ARRAY_CVEC0(x,k),ARRAY_RVEC0(y,l),ARRAY_RVEC1(y,l)); }
      if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='d'){ irdiv_rd(ARRAY_RVEC0(z,i),ARRAY_RVEC1(z,i),ARRAY_RVEC0(x,k),ARRAY_RVEC1(x,k),ARRAY_DVEC0(y,l),ARRAY_DVEC0(y,l)); }
      if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='z'){ icdiv_rz(ARRAY_CVEC0(z,i),ARRAY_CVEC1(z,i),ARRAY_RVEC0(x,k),ARRAY_RVEC1(x,k),ARRAY_ZVEC0(y,l),ARRAY_ZVEC0(y,l)); }
      if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='r'){ irdiv_rr(ARRAY_RVEC0(z,i),ARRAY_RVEC1(z,i),ARRAY_RVEC0(x,k),ARRAY_RVEC1(x,k),ARRAY_RVEC0(y,l),ARRAY_RVEC0(y,l)); }
      if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='c'){ icdiv_rc(ARRAY_CVEC0(z,i),ARRAY_CVEC1(z,i),ARRAY_RVEC0(x,k),ARRAY_RVEC1(x,k),ARRAY_CVEC0(y,l),ARRAY_CVEC0(y,l)); }
      if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='R'){ irdiv_rr(ARRAY_RVEC0(z,i),ARRAY_RVEC1(z,i),ARRAY_RVEC0(x,k),ARRAY_RVEC1(x,k),ARRAY_RVEC0(y,l),ARRAY_RVEC1(y,l)); }
      if(ARRAY_TYPE(x)=='d' && ARRAY_TYPE(y)=='C'){ icdiv_dc(ARRAY_CVEC0(z,i),ARRAY_CVEC1(z,i),ARRAY_DVEC0(x,k),ARRAY_DVEC0(x,k),ARRAY_CVEC0(y,l),ARRAY_CVEC1(y,l)); }
      if(ARRAY_TYPE(x)=='z' && ARRAY_TYPE(y)=='C'){ icdiv_zc(ARRAY_CVEC0(z,i),ARRAY_CVEC1(z,i),ARRAY_ZVEC0(x,k),ARRAY_ZVEC0(x,k),ARRAY_CVEC0(y,l),ARRAY_CVEC1(y,l)); }
      if(ARRAY_TYPE(x)=='r' && ARRAY_TYPE(y)=='C'){ icdiv_rc(ARRAY_CVEC0(z,i),ARRAY_CVEC1(z,i),ARRAY_RVEC0(x,k),ARRAY_RVEC0(x,k),ARRAY_CVEC0(y,l),ARRAY_CVEC1(y,l)); }
      if(ARRAY_TYPE(x)=='c' && ARRAY_TYPE(y)=='C'){ icdiv_cc(ARRAY_CVEC0(z,i),ARRAY_CVEC1(z,i),ARRAY_CVEC0(x,k),ARRAY_CVEC0(x,k),ARRAY_CVEC0(y,l),ARRAY_CVEC1(y,l)); }
      if(ARRAY_TYPE(x)=='R' && ARRAY_TYPE(y)=='C'){ icdiv_rc(ARRAY_CVEC0(z,i),ARRAY_CVEC1(z,i),ARRAY_RVEC0(x,k),ARRAY_RVEC1(x,k),ARRAY_CVEC0(y,l),ARRAY_CVEC1(y,l)); }
      if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='d'){ icdiv_cd(ARRAY_CVEC0(z,i),ARRAY_CVEC1(z,i),ARRAY_CVEC0(x,k),ARRAY_CVEC1(x,k),ARRAY_DVEC0(y,l),ARRAY_DVEC0(y,l)); }
      if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='z'){ icdiv_cz(ARRAY_CVEC0(z,i),ARRAY_CVEC1(z,i),ARRAY_CVEC0(x,k),ARRAY_CVEC1(x,k),ARRAY_ZVEC0(y,l),ARRAY_ZVEC0(y,l)); }
      if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='r'){ icdiv_cr(ARRAY_CVEC0(z,i),ARRAY_CVEC1(z,i),ARRAY_CVEC0(x,k),ARRAY_CVEC1(x,k),ARRAY_RVEC0(y,l),ARRAY_RVEC0(y,l)); }
      if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='c'){ icdiv_cc(ARRAY_CVEC0(z,i),ARRAY_CVEC1(z,i),ARRAY_CVEC0(x,k),ARRAY_CVEC1(x,k),ARRAY_CVEC0(y,l),ARRAY_CVEC0(y,l)); }
      if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='R'){ icdiv_cr(ARRAY_CVEC0(z,i),ARRAY_CVEC1(z,i),ARRAY_CVEC0(x,k),ARRAY_CVEC1(x,k),ARRAY_RVEC0(y,l),ARRAY_RVEC1(y,l)); }
      if(ARRAY_TYPE(x)=='C' && ARRAY_TYPE(y)=='C'){ icdiv_cc(ARRAY_CVEC0(z,i),ARRAY_CVEC1(z,i),ARRAY_CVEC0(x,k),ARRAY_CVEC1(x,k),ARRAY_CVEC0(y,l),ARRAY_CVEC1(y,l)); }
    }
    dim=ivec_free(dim);
    I=ivec_free(I);
    K=ivec_free(K);
    L=ivec_free(L);
  }
  return z;
}



/** @} */

///////////////////////////////////////////////

/** @name array型の３入力演算に関する関数 */
/** @{ */
/** @} */










//EOF
