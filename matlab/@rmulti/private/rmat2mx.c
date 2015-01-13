#include "mex.h"
#include <isys.h>

mxArray *rmat2mx(int m, int n, rmulti **A, int LDA)
{
  const char *field_names[]={"prec","sign","exp","digits"};
  mwSize dims[2]={m,n},scalar[2]={1,1},size[2]={1,1};
  mxArray *ret=NULL,*value=NULL;
  int i,j,k;
  ret=mxCreateStructArray(2,dims,4,field_names);
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      // prec
      value=mxCreateNumericArray(2,scalar,mxINT64_CLASS,mxREAL);
      (*(int64_t*)mxGetData(value))=MAT(A,i,j,LDA)->_mpfr_prec;
      mxSetField(ret,j*m+i,"prec",value);
      // sign
      value=mxCreateNumericArray(2,scalar,mxINT32_CLASS,mxREAL);
      (*(int32_t*)mxGetData(value))=MAT(A,i,j,LDA)->_mpfr_sign;
      mxSetField(ret,j*m+i,"sign",value);
      // exp
      value=mxCreateNumericArray(2,scalar,mxINT64_CLASS,mxREAL);
      (*(int64_t*)mxGetData(value))=MAT(A,i,j,LDA)->_mpfr_exp;
      mxSetField(ret,j*m+i,"exp",value);
      // digits
      size[1]=rget_size(MAT(A,i,j,LDA));
      value=mxCreateNumericArray(2,size,mxUINT64_CLASS,mxREAL);
      for(k=0; k<size[1]; k++){ ((uint64_t*)mxGetData(value))[k]=MAT(A,i,j,LDA)->_mpfr_d[k]; }
      mxSetField(ret,j*m+i,"digits",value);
    }
  }
  return ret;
}
