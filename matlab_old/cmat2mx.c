#include "mex.h"
#include <isys.h>

mxArray *cmat2mx(int m, int n, cmulti **A, int LDA)
{
  const char *field_names[]={"r_prec","r_sign","r_exp","r_digits","i_prec","i_sign","i_exp","i_digits"};
  mwSize dims[2]={m,n},scalar[2]={1,1},size[2]={1,1};
  mxArray *ret=NULL,*value=NULL;
  int i,j,k;
  ret=mxCreateStructArray(2,dims,8,field_names);
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      // real part
      // prec
      value=mxCreateNumericArray(2,scalar,mxINT64_CLASS,mxREAL);
      (*(int64_t*)mxGetData(value))=C_R(MAT(A,i,j,LDA))->_mpfr_prec;
      mxSetField(ret,j*m+i,"r_prec",value);
      // sign
      value=mxCreateNumericArray(2,scalar,mxINT32_CLASS,mxREAL);
      (*(int32_t*)mxGetData(value))=C_R(MAT(A,i,j,LDA))->_mpfr_sign;
      mxSetField(ret,j*m+i,"r_sign",value);
      // exp
      value=mxCreateNumericArray(2,scalar,mxINT64_CLASS,mxREAL);
      (*(int64_t*)mxGetData(value))=C_R(MAT(A,i,j,LDA))->_mpfr_exp;
      mxSetField(ret,j*m+i,"r_exp",value);
      // digits
      size[1]=rget_size(C_R(MAT(A,i,j,LDA)));
      value=mxCreateNumericArray(2,size,mxUINT64_CLASS,mxREAL);
      for(k=0; k<size[1]; k++){ ((uint64_t*)mxGetData(value))[k]=C_R(MAT(A,i,j,LDA))->_mpfr_d[k]; }
      mxSetField(ret,j*m+i,"r_digits",value);
      // imaginary part
      // prec
      value=mxCreateNumericArray(2,scalar,mxINT64_CLASS,mxREAL);
      (*(int64_t*)mxGetData(value))=C_I(MAT(A,i,j,LDA))->_mpfr_prec;
      mxSetField(ret,j*m+i,"i_prec",value);
      // sign
      value=mxCreateNumericArray(2,scalar,mxINT32_CLASS,mxREAL);
      (*(int32_t*)mxGetData(value))=C_I(MAT(A,i,j,LDA))->_mpfr_sign;
      mxSetField(ret,j*m+i,"i_sign",value);
      // exp
      value=mxCreateNumericArray(2,scalar,mxINT64_CLASS,mxREAL);
      (*(int64_t*)mxGetData(value))=C_I(MAT(A,i,j,LDA))->_mpfr_exp;
      mxSetField(ret,j*m+i,"i_exp",value);
      // digits
      size[1]=rget_size(C_I(MAT(A,i,j,LDA)));
      value=mxCreateNumericArray(2,size,mxUINT64_CLASS,mxREAL);
      for(k=0; k<size[1]; k++){ ((uint64_t*)mxGetData(value))[k]=C_I(MAT(A,i,j,LDA))->_mpfr_d[k]; }
      mxSetField(ret,j*m+i,"i_digits",value);
    }
  }
  return ret;
}
