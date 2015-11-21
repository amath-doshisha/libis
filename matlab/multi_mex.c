#include "mex.h"
#include <string.h>
#include <isys.h>

#define IS_DUBL(N,P,I)  ((N>I && mxIsDouble(P[I])))
#define IS_NUMR(N,P,I)  ((N>I && mxIsDouble(P[I]) && mxGetM(P[I])==1 && mxGetN(P[I])==1))
#define IS_ROW(N,P,I)   ((N>I && mxIsDouble(P[I]) && mxGetM(P[I])==1 && mxGetN(P[I])>1 && (mxGetNumberOfDimensions(P[I])<=2 || (mxGetNumberOfDimensions(P[I])==3 && mxGetDimensions(P[I])[2]==1))))
#define IS_CHAR(N,P,I)  ((N>I && mxIsChar(P[I])))
#define IS_STRT(N,P,I)  ((N>I && mxIsStruct(P[I])))
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

typedef struct {
  char type;     // type='r','c','R','C','i'
  int LD1,LD2;   // memory=(LD1,LD2,l)
  int M,N,L;     // size=(m,n,l), 
  void *p0,*p1;  // p0 or [p0,p1]
} multi_struct;
typedef multi_struct multi;

/**
 * @breif MATLABの構造体のフィールド名
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
#include"./c/multi_free.c"
#include"./c/multi_allocate.c"
#include"./c/multi_allocate_mxArray.c"
#include"./c/multi_get_prec.c"
#include"./c/multi_get_sign.c"
#include"./c/multi_get_exp.c"
#include"./c/multi_double.c"
#include"./c/multi_set_zeros.c"
#include"./c/multi_set_ones.c"
#include"./c/multi_set_d.c"
#include"./c/multi_copy.c"
#include"./c/multi_uminus.c"
#include"./c/multi_plus.c"
#include"./c/multi_eq.c"
#include"./c/multi_subsref.c"

/**
 * @breif メイン
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
       if(STR_EQ(cmd,"get_prec")) { multi_get_prec (nlhs,plhs,nrhs,prhs); } // y=get_prec(x)
  else if(STR_EQ(cmd,"get_sign")) { multi_get_sign (nlhs,plhs,nrhs,prhs); } // y=get_sign(x)
  else if(STR_EQ(cmd,"get_exp"))  { multi_get_exp  (nlhs,plhs,nrhs,prhs); } // y=get_exp(x)
  else if(STR_EQ(cmd,"set_zeros")){ multi_set_zeros(nlhs,plhs,nrhs,prhs); } // x=zeros(M,N,L)
  else if(STR_EQ(cmd,"set_ones")) { multi_set_ones (nlhs,plhs,nrhs,prhs); } // x=ones(M,N,L)
  else if(STR_EQ(cmd,"set_d"))    { multi_set_d    (nlhs,plhs,nrhs,prhs); } // y=multi(x), where x is double
  else if(STR_EQ(cmd,"double"))   { multi_double   (nlhs,plhs,nrhs,prhs); } // y=double(x), where x is multi
  else if(STR_EQ(cmd,"copy"))     { multi_copy     (nlhs,plhs,nrhs,prhs); } // y=x
  else if(STR_EQ(cmd,"uminus"))   { multi_uminus   (nlhs,plhs,nrhs,prhs); } // y=-x
  else if(STR_EQ(cmd,"plus"))     { multi_plus     (nlhs,plhs,nrhs,prhs); } // z=x+y
  else if(STR_EQ(cmd,"eq"))       { multi_eq       (nlhs,plhs,nrhs,prhs); } // z=(x==y)
  else if(STR_EQ(cmd,"subsref"))  { multi_subsref  (nlhs,plhs,nrhs,prhs); } // y=x(s)
  else{
    mexPrintf("\n\n\nError!\nmulti_mex(cmd='%s',....)\n",cmd);
    mexErrMsgIdAndTxt("MATLAB:multi_mex","Unknown command.");
  }
  if(cmd!=NULL){ mxFree(cmd); }
  return;
}

// EOF
