#include"is_cmulti.h"
#include"is_cmat.h"
#include"is_irmulti.h"
#include"is_icmulti.h"
#include"is_icvec.h"
#include"is_icmat.h"
#include"is_func.h"

/**
 @file  icmat.c
 @brief 多倍長精度実数型cmultiの機械区間演算の行列に関する関数の定義.
 @details スカラーに関しては@link icmulti.c@endlinkを参照のこと.
          ベクトルに関しては@link icvec.c@endlinkを参照のこと.
          行列に関しては@link icmat.c@endlinkを参照のこと.
 */

/** @name 基本操作 */
/** @{ */


/**
 @brief [B0,B1]=[A0,A1]
*/
int icmat_copy(int m, int n, cmulti **B0, int LDB0, cmulti **B1, int LDB1, cmulti **A0, int LDA0, cmulti **A1, int LDA1)
{
  int i,j,e=0;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      e+=iccopy(MAT(B0,i,j,LDB0),MAT(B1,i,j,LDB1),MAT(A0,i,j,LDA0),MAT(A1,i,j,LDA1));
    }
  }
  return e;
}

/**
 @brief 転置のコピー [B0,B1]=[A0,A1]'
*/
int icmat_copy_t(int m, int n, cmulti **B0, int LDB0, cmulti **B1, int LDB1, cmulti **A0, int LDA0, cmulti **A1, int LDA1)
{
  int i,j,e=0;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      e+=iccopy(MAT(B0,j,i,LDB0),MAT(B1,j,i,LDB1),MAT(A0,i,j,LDA0),MAT(A1,i,j,LDA1));
    }
  }
  return e;
}

/**
 @brief 共役転置のコピー [B0,B1]=[A0,A1]'
*/
int icmat_copy_ct(int m, int n, cmulti **B0, int LDB0, cmulti **B1, int LDB1, cmulti **A0, int LDA0, cmulti **A1, int LDA1)
{
  int i,j,e=0;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      e+=icconj(MAT(B0,j,i,LDB0),MAT(B1,j,i,LDB1),MAT(A0,i,j,LDA0),MAT(A1,i,j,LDA1));
    }
  }
  return e;
}

/**
 @brief 表示
 */
void icmat_print(int m, int n, cmulti **A0, int LDA0, cmulti **A1, int LDA1, const char *name, const char *f, int digits)
{
  int i,j,k;
  char format[128];
  if(STR_EQ(f,"f") || STR_EQ(f,"F")){
    k=3;
    sprintf(format,"[%%%d.%dR%s, %%%d.%dR%s] [%%%d.%dR%s, %%%d.%dR%s]  ",digits+k,digits,f,digits+k,digits,f,digits+k,digits,f,digits+k,digits,f);
  }else if((strcmp(f,"e")==0) || (strcmp(f,"E")==0)){
    k=7;
    sprintf(format,"[%%%d.%dR%s, %%%d.%dR%s] [%%%d.%dR%s, %%%d.%dR%s]  ",digits+k,digits,f,digits+k,digits,f,digits+k,digits,f,digits+k,digits,f);
  }
  if(A0==NULL || A1==NULL){
    if(name!=NULL){ printf("%s=NULL\n",name); }
    else          { printf("NULL\n"); }
    return;
  }
  if(name!=NULL){ printf("%s\n",name); }
  if(A0==NULL || A1==NULL || m<=0 || n<=0) return;
  for(i=0; i<m; i++){
    for(j=0; j<n; j++){
      mpfr_printf(format,C_R(MAT(A0,i,j,LDA0)),C_R(MAT(A1,i,j,LDA1)),C_I(MAT(A0,i,j,LDA0)),C_I(MAT(A1,i,j,LDA1)));
    }
    printf("\n");
  }
}

/**
 @brief 区間の中心 Ac=(A1+A0)/2 と半径 Ar=A1-A0
*/
int icmat_center_radius(int m, int n, cmulti **Ac, int LDAc, cmulti **Ar, int LDAr, cmulti **A0, int LDA0, cmulti **A1, int LDA1)
{
  int e=0;
  mpfr_rnd_t mode;
  mode=get_round_mode();
  set_round_mode(MPFR_RNDU);             // up
  e+=cmat_sub(m,n,Ac,LDAc,A1,LDA1,A0,LDA0); // Ac=A1-A0
  e+=cmat_mul_d(m,n,Ac,LDAc,Ac,LDAc,0.5);   // Ac=(A1-A0)/2
  e+=cmat_add(m,n,Ac,LDAc,Ac,LDAc,A0,LDA0); // Ac=(A1-A0)/2+A0
  e+=cmat_sub(m,n,Ar,LDAr,Ac,LDAc,A0,LDA0); // Ar=Ac-A0
  set_round_mode(mode);                  // back
  return e;
}


/**
 @brief [A0,A1]=ones(m,n)*a
*/
int icmat_set_all_d(int m, int n, cmulti **A0, int LDA0, cmulti **A1, int LDA1, double a)
{
  int i,j,e=0;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      e+=icset_d(MAT(A0,i,j,LDA0),MAT(A1,i,j,LDA1),a);
    }
  }
  return e;
}


/**
 @brief [A0,A1]=[zeros(m,n),zeros(m,n)]
*/
int icmat_set_zeros(int m, int n, cmulti **A0, int LDA0, cmulti **A1, int LDA1)
{
  int e=0;
  e+=icmat_set_all_d(m,n,A0,LDA0,A1,LDA1,0);
  return e;
}

/** @} */

//////////////////////////////////////////////////////////////////////////

/** @name 四則演算 */
/** @{ */

/**
 @brief 積 [C0,C1]=[A0,A1]*[B0,B1]
*/
int icmat_prod(int l, int m, int n, cmulti **C0, int LDC0, cmulti **C1, int LDC1, cmulti **A0, int LDA0, cmulti **A1, int LDA1, cmulti **B0, int LDB0, cmulti **B1, int LDB1)
{
  int p0,p1,prec,i,j,k,e=0;
  cmulti **Z0=NULL,**Z1=NULL;
  p0=cmat_get_prec_max(l,n,C0,LDC0); p1=cmat_get_prec_max(l,n,C1,LDC1); prec=MAX2(p0,p1);
  Z0=cmat_allocate_prec(l,n,prec); Z1=cmat_allocate_prec(l,n,prec);
  for(i=0; i<l; i++){
    for(j=0; j<n; j++){
      e+=icset_d(MAT(Z0,i,j,l),MAT(Z1,i,j,l),0);
      for(k=0; k<m; k++){
	e+=icadd_mul(MAT(Z0,i,j,l),MAT(Z1,i,j,l),MAT(A0,i,k,LDA0),MAT(A1,i,k,LDA1),MAT(B0,k,j,LDB0),MAT(B1,k,j,LDB1));
      }
    }
  }
  e+=icmat_copy(l,n,C0,LDC0,C1,LDC1,Z0,l,Z1,l);
  Z0=cmat_free(l,n,Z0);
  Z1=cmat_free(l,n,Z1);
  return e;
}

/**
 @brief 積の加算 [C0,C1]=[C0,C1]+[A0,A1]*[B0,B1]
*/
int icmat_add_prod(int l, int m, int n, cmulti **C0, int LDC0, cmulti **C1, int LDC1, cmulti **A0, int LDA0, cmulti **A1, int LDA1, cmulti **B0, int LDB0, cmulti **B1, int LDB1)
{
  int p0,p1,prec,i,j,k,e=0;
  cmulti **Z0=NULL,**Z1=NULL;
  p0=cmat_get_prec_max(l,n,C0,LDC0); p1=cmat_get_prec_max(l,n,C1,LDC1); prec=MAX2(p0,p1);
  Z0=cmat_allocate_prec(l,n,prec); Z1=cmat_allocate_prec(l,n,prec);
  for(i=0; i<l; i++){
    for(j=0; j<n; j++){
      e+=iccopy(MAT(Z0,i,j,l),MAT(Z1,i,j,l),MAT(C0,i,j,LDC0),MAT(C1,i,j,LDC1));
      for(k=0; k<m; k++){
	e+=icadd_mul(MAT(Z0,i,j,l),MAT(Z1,i,j,l),MAT(A0,i,k,LDA0),MAT(A1,i,k,LDA1),MAT(B0,k,j,LDB0),MAT(B1,k,j,LDB1));
      }
    }
  }
  e+=icmat_copy(l,n,C0,LDC0,C1,LDC1,Z0,l,Z1,l);
  Z0=cmat_free(l,n,Z0); Z1=cmat_free(l,n,Z1);
  return e;
}

/**
 @brief 積の減算 [C0,C1]=[C0,C1]-[A0,A1]*[B0,B1]
*/
int icmat_sub_prod(int l, int m, int n, cmulti **C0, int LDC0, cmulti **C1, int LDC1, cmulti **A0, int LDA0, cmulti **A1, int LDA1, cmulti **B0, int LDB0, cmulti **B1, int LDB1)
{
  int p0,p1,prec,i,j,k,e=0;
  cmulti **Z0=NULL,**Z1=NULL;
  p0=cmat_get_prec_max(l,n,C0,LDC0); p1=cmat_get_prec_max(l,n,C1,LDC1); prec=MAX2(p0,p1);
  Z0=cmat_allocate_prec(l,n,prec); Z1=cmat_allocate_prec(l,n,prec);
  for(i=0; i<l; i++){
    for(j=0; j<n; j++){
      e+=iccopy(MAT(Z0,i,j,l),MAT(Z1,i,j,l),MAT(C0,i,j,LDC0),MAT(C1,i,j,LDC1));
      for(k=0; k<m; k++){
	e+=icsub_mul(MAT(Z0,i,j,l),MAT(Z1,i,j,l),MAT(A0,i,k,LDA0),MAT(A1,i,k,LDA1),MAT(B0,k,j,LDB0),MAT(B1,k,j,LDB1));
      }
    }
  }
  e+=icmat_copy(l,n,C0,LDC0,C1,LDC1,Z0,l,Z1,l);
  Z0=cmat_free(l,n,Z0); Z1=cmat_free(l,n,Z1);
  return e;
}

/** @} */

//////////////////////////////////////////////////////////////////////////

/** @name 写像 */
/** @{ */

/**
 @brief 行列写像 [A0,A1]=F([x0,x1])
*/
int icmat_func_list2(int m, int n, cmulti **A0, cmulti **A1, int LDA, func_t *f, int l, cmulti **x0, cmulti **x1)
{
  int i,j,e=0;
  for(i=0; i<m; i++){
    for(j=0; j<n; j++){
      if(func_is_list(f) && i<func_asize(f) && func_is_list(f->a[i]) && j<func_asize(func_aget(f,i))){
	e+=icvec_func(MAT(A0,i,j,LDA),MAT(A1,i,j,LDA),func_aget(func_aget(f,i),j),l,x0,x1);
      }else{ cset_nan(MAT(A0,i,j,LDA)); cset_nan(MAT(A1,i,j,LDA)); }
    }
  }
  return e;
}

/** @} */

//EOF
