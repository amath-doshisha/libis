#include "mex.h"
#include <isys.h>

void mx2cmat(int m, int n, cmulti **A, int LDA, const mxArray *src)
{
  mwSize size[2]={1,1};
  mxArray *value=NULL;
  int i,j,k;

  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      // real part
      // prec
      value=mxGetField(src,j*m+i,"r_prec");
      if(value!=NULL && mxIsInt64(value)){ rround(C_R(MAT(A,i,j,LDA)),(*(int64_t*)mxGetData(value))); }
      else{ mexErrMsgIdAndTxt("MATLAB:mx2cmat","The arg should be Struct with the feild 'prec'."); }
      // sign
      value=mxGetField(src,j*m+i,"r_sign");
      if(value!=NULL && mxIsInt32(value)){ C_R(MAT(A,i,j,LDA))->_mpfr_sign=(*(int32_t*)mxGetData(value)); }
      else{ mexErrMsgIdAndTxt("MATLAB:mx2cmat","The arg should be Struct with the feild 'sign'."); }
      // exp
      value=mxGetField(src,j*m+i,"r_exp");
      if(value!=NULL && mxIsInt64(value)){ C_R(MAT(A,i,j,LDA))->_mpfr_exp=(*(int64_t*)mxGetData(value)); }
      else{ mexErrMsgIdAndTxt("MATLAB:mx2cmat","The arg should be Struct with the feild 'exp'."); }
      // digits
      value=mxGetField(src,j*m+i,"r_digits");
      if(value!=NULL && mxIsUint64(value)){
	for(k=0; k<rget_size(C_R(MAT(A,i,j,LDA))); k++){
	  C_R(MAT(A,i,j,LDA))->_mpfr_d[k]=((uint64_t*)mxGetData(value))[k];
	}
      }
      else{ mexErrMsgIdAndTxt("MATLAB:mx2cmat","The arg should be Struct with the feild 'digits'."); }

      // imaginary part
      // prec
      value=mxGetField(src,j*m+i,"i_prec");
      if(value!=NULL && mxIsInt64(value)){ rround(C_I(MAT(A,i,j,LDA)),(*(int64_t*)mxGetData(value))); }
      else{ mexErrMsgIdAndTxt("MATLAB:mx2cmat","The arg should be Struct with the feild 'prec'."); }
      // sign
      value=mxGetField(src,j*m+i,"i_sign");
      if(value!=NULL && mxIsInt32(value)){ C_I(MAT(A,i,j,LDA))->_mpfr_sign=(*(int32_t*)mxGetData(value)); }
      else{ mexErrMsgIdAndTxt("MATLAB:mx2cmat","The arg should be Struct with the feild 'sign'."); }
      // exp
      value=mxGetField(src,j*m+i,"i_exp");
      if(value!=NULL && mxIsInt64(value)){ C_I(MAT(A,i,j,LDA))->_mpfr_exp=(*(int64_t*)mxGetData(value)); }
      else{ mexErrMsgIdAndTxt("MATLAB:mx2cmat","The arg should be Struct with the feild 'exp'."); }
      // digits
      value=mxGetField(src,j*m+i,"i_digits");
      if(value!=NULL && mxIsUint64(value)){
	for(k=0; k<rget_size(C_I(MAT(A,i,j,LDA))); k++){
	  C_I(MAT(A,i,j,LDA))->_mpfr_d[k]=((uint64_t*)mxGetData(value))[k];
	}
      }
      else{ mexErrMsgIdAndTxt("MATLAB:mx2cmat","The arg should be Struct with the feild 'digits'."); }
    }
  }
  return;
}

