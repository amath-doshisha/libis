#include "mex.h"
#include <isys.h>

void mx2rmat(int m, int n, rmulti **A, int LDA, const mxArray *src)
{
  mwSize size[2]={1,1};
  mxArray *value=NULL;
  int i,j,k;

  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      // prec
      value=mxGetField(src,j*m+i,"prec");
      if(value!=NULL && mxIsInt64(value)){ rround(MAT(A,i,j,LDA),(*(int64_t*)mxGetData(value))); }
      else{ mexErrMsgIdAndTxt("MATLAB:rmat2mx","The arg should be Struct with the feild 'prec'."); }
      // sign
      value=mxGetField(src,j*m+i,"sign");
      if(value!=NULL && mxIsInt32(value)){ MAT(A,i,j,LDA)->_mpfr_sign=(*(int32_t*)mxGetData(value)); }
      else{ mexErrMsgIdAndTxt("MATLAB:rmat2mx","The arg should be Struct with the feild 'sign'."); }
      // exp
      value=mxGetField(src,j*m+i,"exp");
      if(value!=NULL && mxIsInt64(value)){ MAT(A,i,j,LDA)->_mpfr_exp=(*(int64_t*)mxGetData(value)); }
      else{ mexErrMsgIdAndTxt("MATLAB:rmat2mx","The arg should be Struct with the feild 'exp'."); }
      // digits
      value=mxGetField(src,j*m+i,"digits");
      if(value!=NULL && mxIsUint64(value)){
	for(k=0; k<rget_size(MAT(A,i,j,LDA)); k++){
	  MAT(A,i,j,LDA)->_mpfr_d[k]=((uint64_t*)mxGetData(value))[k];
	}
      }
      else{ mexErrMsgIdAndTxt("MATLAB:rmat2mx","The arg should be Struct with the feild 'digits'."); }
    }
  }
  return;
}

