#include "mex.h"
#include <string.h>
#include <isys.h>

#define IS_DUBL(N,P,I)  ((N>I && mxIsDouble(P[I])))
#define IS_NUMR(N,P,I)  ((N>I && mxIsDouble(P[I]) && mxGetM(P[I])==1 && mxGetN(P[I])==1))
#define IS_ROW(N,P,I)   ((N>I && mxIsDouble(P[I]) && mxGetM(P[I])==1 && mxGetN(P[I])>1 && (mxGetNumberOfDimensions(P[I])<=2 || (mxGetNumberOfDimensions(P[I])==3 && mxGetDimensions(P[I])[2]==1))))
#define IS_CHAR(N,P,I)  ((N>I && mxIsChar(P[I])))
#define IS_STRT(N,P,I)  ((N>I && mxIsStruct(P[I])))
#define IS_CELL(N,P,I)  ((N>I && mxIsCell(P[I])))
#define GET_DOUBLE(P)   ((double*)mxGetData(P))
#define MATLAB_ERROR(S) (mexErrMsgIdAndTxt("MATLAB:multi_mex",S));
#define N0 3

#define _T(A)   (A->type)
#define _M(A)   (A->M)
#define _N(A)   (A->N)
#define _L(A)   (A->L)
#define _LD1(A) (A->LD1)
#define _LD2(A) (A->LD2)
#define _R(A)   ((rmulti**)(A->p0))
#define _R0(A)  ((rmulti**)(A->p0))
#define _R1(A)  ((rmulti**)(A->p1))
#define _C(A)   ((cmulti**)(A->p0))
#define _C0(A)  ((cmulti**)(A->p0))
#define _C1(A)  ((cmulti**)(A->p1))
#define _D(A)   ((double*)(A->p0))
#define _Z(A)   ((dcomplex*)(A->p0))
#define _I(A)   ((int*)(A->p0))
#define _S(A)   ((char**)(A->p0))

typedef struct {
  char type;     // type='r','c','R','C','i'
  int LD1,LD2;   // memory=(LD1,LD2,l)
  int M,N,L;     // size=(m,n,l), 
  void *p0,*p1;  // p0 or [p0,p1]
} multi_struct;
typedef multi_struct multi;

typedef struct {
  int ndim;
  int *dims;
  int **index;
} subs_index_t;

/**
 * @breif The field names of struct for MATLAB
 */
const char *r_field_names[]={"r_prec","r_sign","r_exp","r_digits"};  
const char *c_field_names[]={"cr_prec","cr_sign","cr_exp","cr_digits",
			     "ci_prec","ci_sign","ci_exp","ci_digits"};
const char *cr_field_names[]={"cr_prec","cr_sign","cr_exp","cr_digits"};
const char *ci_field_names[]={"ci_prec","ci_sign","ci_exp","ci_digits"};  
const char *R_field_names[]={"R0_prec","R0_sign","R0_exp","R0_digits",
			     "R1_prec","R1_sign","R1_exp","R1_digits"};  
const char *R0_field_names[]={"R0_prec","R0_sign","R0_exp","R0_digits"};
const char *R1_field_names[]={"R1_prec","R1_sign","R1_exp","R1_digits"};
const char *C_field_names[]={"C0r_prec","C0r_sign","C0r_exp","C0r_digits",
			     "C0i_prec","C0i_sign","C0i_exp","C0i_digits",
			     "C1r_prec","C1r_sign","C1r_exp","C1r_digits",
			     "C1i_prec","C1i_sign","C1i_exp","C1i_digits"};
const char *C0r_field_names[]={"C0r_prec","C0r_sign","C0r_exp","C0r_digits"};
const char *C0i_field_names[]={"C0i_prec","C0i_sign","C0i_exp","C0i_digits"};
const char *C1r_field_names[]={"C1r_prec","C1r_sign","C1r_exp","C1r_digits"};
const char *C1i_field_names[]={"C1i_prec","C1i_sign","C1i_exp","C1i_digits"};

/**
 * functions
 */
#include"./c/mxCreateStructMultiValue.c"
#include"./c/mxCreateStructMulti.c"
#include"./c/mxArrayToMulti.c"
#include"./c/subs_allocate.c"
#include"./c/multi_free.c"
#include"./c/multi_allocate.c"
#include"./c/multi_allocate_mxArray.c"
#include"./c/multi_get_prec.c"
#include"./c/multi_get_sign.c"
#include"./c/multi_get_exp.c"
#include"./c/multi_set_zeros.c"
#include"./c/multi_set_ones.c"
#include"./c/multi_set_eye.c"
#include"./c/multi_set_d.c"
#include"./c/multi_get_d.c"
#include"./c/multi_set_s.c"
#include"./c/multi_get_s.c"
#include"./c/multi_copy.c"
#include"./c/multi_uminus.c"
#include"./c/multi_ctranspose.c"
#include"./c/multi_transpose.c"
#include"./c/multi_plus.c"
#include"./c/multi_minus.c"
#include"./c/multi_times.c"
#include"./c/multi_mtimes.c"
#include"./c/multi_rdivide.c"
#include"./c/multi_mrdivide.c"
#include"./c/multi_mldivide.c"
#include"./c/multi_power.c"
#include"./c/multi_mpower.c"
#include"./c/multi_eq.c"
#include"./c/multi_ne.c"
#include"./c/multi_ge.c"
#include"./c/multi_gt.c"
#include"./c/multi_le.c"
#include"./c/multi_lt.c"
#include"./c/multi_subsref.c"
#include"./c/multi_subsasgn.c"
#include"./c/multi_horzcat.c"
#include"./c/multi_vertcat.c"
#include"./c/multi_inv.c"
#include"./c/multi_real.c"
#include"./c/multi_imag.c"
#include"./c/multi_conj.c"
#include"./c/multi_abs.c"
#include"./c/multi_angle.c"
#include"./c/multi_max.c"
#include"./c/multi_min.c"
#include"./c/multi_eig.c"
#include"./c/multi_eig_verify.c"
#include"./c/multi_matgen_dhToda.c"

/**
 * @breif mexFunction() for multi_mex
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  int prec=64,ap_mode=0;
  char *cmd=NULL;
  if(!(IS_CHAR(nrhs,prhs,0))){ MATLAB_ERROR("mexFunction: The arg0 should be Char."); }
  if(!(IS_NUMR(nrhs,prhs,1))){ MATLAB_ERROR("mexFunction: The arg1 should be Double."); }
  if(!(IS_NUMR(nrhs,prhs,2))){ MATLAB_ERROR("mexFunction: The arg2 should be Double."); }
  if(IS_CHAR(nrhs,prhs,0)){ cmd=mxArrayToString(prhs[0]); }
  if(IS_NUMR(nrhs,prhs,1)){ prec=GET_DOUBLE(prhs[1])[0]; }
  if(IS_NUMR(nrhs,prhs,2)){ ap_mode=GET_DOUBLE(prhs[2])[0]; }
  set_default_prec(prec);
  set_auto_prec_mode(ap_mode);
  //  mexPrintf("cmd=%s prec=%d ap_mode=%d\n",cmd,prec,ap_mode);  
       if(STR_EQ(cmd,"get_prec"))  { multi_get_prec  (nlhs,plhs,nrhs,prhs); }  // y=get_prec(x)
  else if(STR_EQ(cmd,"get_sign"))  { multi_get_sign  (nlhs,plhs,nrhs,prhs); }  // y=get_sign(x)
  else if(STR_EQ(cmd,"get_exp"))   { multi_get_exp   (nlhs,plhs,nrhs,prhs); }  // y=get_exp(x)
  else if(STR_EQ(cmd,"set_zeros")) { multi_set_zeros (nlhs,plhs,nrhs,prhs); }  // x=zeros(M,N,L)
  else if(STR_EQ(cmd,"set_ones"))  { multi_set_ones  (nlhs,plhs,nrhs,prhs); }  // x=ones(M,N,L)
  else if(STR_EQ(cmd,"set_eye"))   { multi_set_eye   (nlhs,plhs,nrhs,prhs); }  // x=eye(M,N)
  else if(STR_EQ(cmd,"set_d"))     { multi_set_d     (nlhs,plhs,nrhs,prhs); }  // y=multi(x), where x is double
  else if(STR_EQ(cmd,"get_d"))     { multi_get_d     (nlhs,plhs,nrhs,prhs); }  // y=double(x)
  else if(STR_EQ(cmd,"set_s"))     { multi_set_s     (nlhs,plhs,nrhs,prhs); }  // y=multi(x), where x is cell of char
  else if(STR_EQ(cmd,"get_s"))     { multi_get_s     (nlhs,plhs,nrhs,prhs); }  // y=char(x)
  else if(STR_EQ(cmd,"copy"))      { multi_copy      (nlhs,plhs,nrhs,prhs); }  // y=x
  else if(STR_EQ(cmd,"uminus"))    { multi_uminus    (nlhs,plhs,nrhs,prhs); }  // y=-x
  else if(STR_EQ(cmd,"ctranspose")){ multi_ctranspose(nlhs,plhs,nrhs,prhs); }  // y=x'
  else if(STR_EQ(cmd,"transpose")) { multi_transpose (nlhs,plhs,nrhs,prhs); }  // y=x.'
  else if(STR_EQ(cmd,"uminus"))    { multi_uminus    (nlhs,plhs,nrhs,prhs); }  // y=-x
  else if(STR_EQ(cmd,"plus"))      { multi_plus      (nlhs,plhs,nrhs,prhs); }  // z=x+y
  else if(STR_EQ(cmd,"minus"))     { multi_minus     (nlhs,plhs,nrhs,prhs); }  // z=x-y
  else if(STR_EQ(cmd,"times"))     { multi_times     (nlhs,plhs,nrhs,prhs); }  // z=x.*y
  else if(STR_EQ(cmd,"mtimes"))    { multi_mtimes    (nlhs,plhs,nrhs,prhs); }  // z=x*y
  else if(STR_EQ(cmd,"rdivide"))   { multi_rdivide   (nlhs,plhs,nrhs,prhs); }  // z=x./y
  else if(STR_EQ(cmd,"mrdivide"))  { multi_mrdivide  (nlhs,plhs,nrhs,prhs); }  // z=x/y
  else if(STR_EQ(cmd,"mldivide"))  { multi_mldivide  (nlhs,plhs,nrhs,prhs); }  // z=x\y
  else if(STR_EQ(cmd,"power"))     { multi_power     (nlhs,plhs,nrhs,prhs); }  // z=x.^y
  else if(STR_EQ(cmd,"mpower"))    { multi_mpower    (nlhs,plhs,nrhs,prhs); }  // z=x^y
  else if(STR_EQ(cmd,"eq"))        { multi_eq        (nlhs,plhs,nrhs,prhs); }  // z=(x==y)
  else if(STR_EQ(cmd,"ne"))        { multi_ne        (nlhs,plhs,nrhs,prhs); }  // z=(x~=y)
  else if(STR_EQ(cmd,"ge"))        { multi_ge        (nlhs,plhs,nrhs,prhs); }  // z=(x>=y)
  else if(STR_EQ(cmd,"gt"))        { multi_gt        (nlhs,plhs,nrhs,prhs); }  // z=(x>y)
  else if(STR_EQ(cmd,"le"))        { multi_le        (nlhs,plhs,nrhs,prhs); }  // z=(x<=y)
  else if(STR_EQ(cmd,"lt"))        { multi_lt        (nlhs,plhs,nrhs,prhs); }  // z=(x<y)
  else if(STR_EQ(cmd,"subsref"))   { multi_subsref   (nlhs,plhs,nrhs,prhs); }  // y=x(s)
  else if(STR_EQ(cmd,"subsasgn"))  { multi_subsasgn  (nlhs,plhs,nrhs,prhs); }  // y(s)=x
  else if(STR_EQ(cmd,"horzcat"))   { multi_horzcat   (nlhs,plhs,nrhs,prhs); }  // y=[x1 x2 ...]
  else if(STR_EQ(cmd,"vertcat"))   { multi_vertcat   (nlhs,plhs,nrhs,prhs); }  // y=[x1; x2; ...]
  else if(STR_EQ(cmd,"inv"))       { multi_inv       (nlhs,plhs,nrhs,prhs); }  // y=inv(x)
  else if(STR_EQ(cmd,"real"))      { multi_real      (nlhs,plhs,nrhs,prhs); }  // y=real(x)
  else if(STR_EQ(cmd,"imag"))      { multi_imag      (nlhs,plhs,nrhs,prhs); }  // y=imag(x)
  else if(STR_EQ(cmd,"conj"))      { multi_conj      (nlhs,plhs,nrhs,prhs); }  // y=conj(x)
  else if(STR_EQ(cmd,"abs"))       { multi_abs       (nlhs,plhs,nrhs,prhs); }  // y=abs(x)
  else if(STR_EQ(cmd,"angle"))     { multi_angle     (nlhs,plhs,nrhs,prhs); }  // y=angle(x)      
  else if(STR_EQ(cmd,"max"))       { multi_max       (nlhs,plhs,nrhs,prhs); }  // y=max(x)
  else if(STR_EQ(cmd,"min"))       { multi_min       (nlhs,plhs,nrhs,prhs); }  // y=min(x)
  else if(STR_EQ(cmd,"eig"))       { multi_eig       (nlhs,plhs,nrhs,prhs); }  // lambda=eig(A), [V,D]=eig(A)
  else if(STR_EQ(cmd,"eig_verify")){ multi_eig_verify(nlhs,plhs,nrhs,prhs); }  // [V,D,EV,ED]=eig_verify(A)
  else if(STR_EQ(cmd,"matgen_dhToda")){ multi_matgen_dhToda(nlhs,plhs,nrhs,prhs); } // A=matgen_dhToda(...)       
  else{
    mexPrintf("\n\n\nError!\nmulti_mex(cmd='%s',....)\n",cmd);
    mexErrMsgIdAndTxt("MATLAB:multi_mex","Unknown command.");
  }
  if(cmd!=NULL){ mxFree(cmd); }
  return;
}

// EOF
