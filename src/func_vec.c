#include<stdio.h>
#include<stdlib.h>
#include"is_macros.h"
#include"is_ivec.h"
#include"is_strings.h"
#include"is_print.h"
#include"is_func.h"
#include"is_rmulti.h"
#include"is_cmulti.h"
#include"is_rvec.h"
#include"is_cvec.h"
#include"is_rmat.h"
#include"is_cmat.h"
#include"is_rsolve.h"
#include"is_csolve.h"

#define FR(F) func_retain(F)


////////////////////////////////////////////

int func_rows(func_t *f)
{
  if(func_is(f,"rmat")){ return func_rmat_rows(f); }
  else if(func_is(f,"cmat")){ return func_cmat_rows(f); }
  else{ FUNC_ERROR_ARG1("func_rows",f); }
}

int func_cols(func_t *f)
{
  if(func_is(f,"rmat")){ return func_rmat_cols(f); }
  else if(func_is(f,"cmat")){ return func_cmat_cols(f); }
  else{ FUNC_ERROR_ARG1("func_cols",f); }
}

////////////////////////////////////////////

int func_is_rvec(func_t *f){ return func_is(f,"rvec"); }
int func_is_cvec(func_t *f){ return func_is(f,"cvec"); }
int func_in_vec(func_t *f){ return (func_is_rvec(f) || func_is_cvec(f)); }
int func_is_rmat(func_t *f){ return func_is(f,"rmat"); }
int func_is_cmat(func_t *f){ return func_is(f,"cmat"); }
int func_in_mat(func_t *f){ return (func_is_rmat(f) || func_is_cmat(f)); }

////////////////////////////////////////////


#define FRV(f) (f)->p.rvec->x
#define FCV(f) (f)->p.cvec->x
#define FRS(f) func_rvec_size(f)
#define FCS(f) func_cvec_size(f)

func_t *func_vec_add(func_t *f1, func_t *f2)
{
  func_t *f=NULL;
  if     (func_is(f1,"rvec") && func_is(f2,"rvec") && FRS(f1)==FRS(f2)){ f=func_rvec(FRS(f1)); rvec_add  (FRS(f),FRV(f),FRV(f1),FRV(f2)); }
  else if(func_is(f1,"rvec") && func_is(f2,"cvec") && FRS(f1)==FCS(f2)){ f=func_cvec(FRS(f1)); cvec_add_rvec(FCS(f),FCV(f),FCV(f2),FRV(f1)); }
  else if(func_is(f1,"cvec") && func_is(f2,"rvec") && FCS(f1)==FRS(f2)){ f=func_cvec(FCS(f1)); cvec_add_rvec(FCS(f),FCV(f),FCV(f1),FRV(f2)); }
  else if(func_is(f1,"cvec") && func_is(f2,"cvec") && FCS(f1)==FCS(f2)){ f=func_cvec(FCS(f1)); cvec_add  (FCS(f),FCV(f),FCV(f1),FCV(f2)); }
  else { FUNC_ERROR_ARG2("func_vec_add",f1,f2); }
  f1=func_del(f1);
  f2=func_del(f2);
  return f;
}

func_t *func_vec_sub(func_t *f1, func_t *f2)
{
  func_t *f=NULL;
  if     (func_is(f1,"rvec") && func_is(f2,"rvec") && FRS(f1)==FRS(f2)){ f=func_rvec(FRS(f1)); rvec_sub   (FRS(f),FRV(f),FRV(f1),FRV(f2)); }
  else if(func_is(f1,"rvec") && func_is(f2,"cvec") && FRS(f1)==FCS(f2)){ f=func_cvec(FRS(f1)); cvec_sub_rvec1(FCS(f),FCV(f),FRV(f1),FCV(f2)); }
  else if(func_is(f1,"cvec") && func_is(f2,"rvec") && FCS(f1)==FRS(f2)){ f=func_cvec(FCS(f1)); cvec_sub_rvec2(FCS(f),FCV(f),FCV(f1),FRV(f2)); }
  else if(func_is(f1,"cvec") && func_is(f2,"cvec") && FCS(f1)==FCS(f2)){ f=func_cvec(FCS(f1)); cvec_sub   (FCS(f),FCV(f),FCV(f1),FCV(f2)); }
  else { FUNC_ERROR_ARG2("func_vec_sub",f1,f2); }
  f1=func_del(f1);
  f2=func_del(f2);
  return f;
}

////////////////////////////////////////////////////////////////////////////////////////////////

#define FRR(f) f->p.rmat->m
#define FCR(f) func_cmat_rows(f)
#define FRC(f) f->p.rmat->n
#define FCC(f) func_cmat_cols(f)
#define FRD(f) f->p.rmat->LD
#define FCD(f) func_cmat_ld(f)
#define FRM(f) f->p.rmat->A
#define FCM(f) func_cmat_p(f)

func_t *func_mat_add(func_t *f1, func_t *f2)
{
  func_t *f=NULL;
  if     (func_is(f1,"rmat") && func_is(f2,"rmat") && FRR(f1)==FRR(f2) && FRC(f1)==FRC(f2)){ f=func_rmat(FRR(f1),FRC(f1)); rmat_add  (FRR(f),FRC(f),FRM(f),FRD(f),FRM(f1),FRD(f1),FRM(f2),FRD(f2)); }
  else if(func_is(f1,"rmat") && func_is(f2,"cmat") && FRR(f1)==FCR(f2) && FRC(f1)==FCC(f2)){ f=func_cmat(FRR(f1),FRC(f1)); cmat_add_r(FCR(f),FCC(f),FCM(f),FCD(f),FCM(f2),FCD(f2),FRM(f1),FRD(f1)); }
  else if(func_is(f1,"cmat") && func_is(f2,"rmat") && FCR(f1)==FRR(f2) && FCC(f1)==FRC(f2)){ f=func_cmat(FCR(f1),FCC(f1)); cmat_add_r(FCR(f),FCC(f),FCM(f),FCD(f),FCM(f1),FCD(f1),FRM(f2),FRD(f2)); }
  else if(func_is(f1,"cmat") && func_is(f2,"cmat") && FCR(f1)==FCR(f2) && FCC(f1)==FCC(f2)){ f=func_cmat(FCR(f1),FCC(f1)); cmat_add  (FCR(f),FCC(f),FCM(f),FCD(f),FCM(f1),FCD(f1),FCM(f2),FCD(f2)); }
  else { FUNC_ERROR_ARG2("func_mat_add",f1,f2); }
  f1=func_del(f1);
  f2=func_del(f2);
  return f;
}

func_t *func_mat_sub(func_t *f1, func_t *f2)
{
  func_t *f=NULL;
  if     (func_is(f1,"rmat") && func_is(f2,"rmat") && FRR(f1)==FRR(f2) && FRC(f1)==FRC(f2)){ f=func_rmat(FRR(f1),FRC(f1)); rmat_sub   (FRR(f),FRC(f),FRM(f),FRD(f),FRM(f1),FRD(f1),FRM(f2),FRD(f2)); }
  else if(func_is(f1,"rmat") && func_is(f2,"cmat") && FRR(f1)==FCR(f2) && FRC(f1)==FCC(f2)){ f=func_cmat(FRR(f1),FRC(f1)); cmat_sub_r1(FCR(f),FCC(f),FCM(f),FCD(f),FRM(f1),FRD(f1),FCM(f2),FCD(f2)); }
  else if(func_is(f1,"cmat") && func_is(f2,"rmat") && FCR(f1)==FRR(f2) && FCC(f1)==FRC(f2)){ f=func_cmat(FCR(f1),FCC(f1)); cmat_sub_r2(FCR(f),FCC(f),FCM(f),FCD(f),FCM(f1),FCD(f1),FRM(f2),FRD(f2)); }
  else if(func_is(f1,"cmat") && func_is(f2,"cmat") && FCR(f1)==FCR(f2) && FCC(f1)==FCC(f2)){ f=func_cmat(FCR(f1),FCC(f1)); cmat_sub   (FCR(f),FCC(f),FCM(f),FCD(f),FCM(f1),FCD(f1),FCM(f2),FCD(f2)); }
  else { FUNC_ERROR_ARG2("func_mat_sub",f1,f2); }
  f1=func_del(f1);
  f2=func_del(f2);
  return f;
}

//EOF
