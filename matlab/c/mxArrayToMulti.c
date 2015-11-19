/**
 * MATLABの構造体から多倍長型へ変換
 */
void mxArrayToMulti(rmulti *y, const mxArray *x, int i, const char *field_names[])
{
  mxArray *value=NULL;
  int t;
  // prec
  value=mxGetField(x,i,field_names[0]);
  if(value!=NULL && mxIsInt64(value)){ rround(y,(*(int64_t*)mxGetData(value))); }
  else{ mexErrMsgIdAndTxt("MATLAB:rmat2mx","The arg should be Struct with the feild 'prec'."); }
  // sign
  value=mxGetField(x,i,field_names[1]);
  if(value!=NULL && mxIsInt32(value)){ y->_mpfr_sign=(*(int32_t*)mxGetData(value)); }
  else{ mexErrMsgIdAndTxt("MATLAB:rmat2mx","The arg should be Struct with the feild 'sign'."); }
  // exp
  value=mxGetField(x,i,field_names[2]);
  if(value!=NULL && mxIsInt64(value)){ y->_mpfr_exp=(*(int64_t*)mxGetData(value)); }
  else{ mexErrMsgIdAndTxt("MATLAB:rmat2mx","The arg should be Struct with the feild 'exp'."); }
  // digits
  value=mxGetField(x,i,field_names[3]);
  if(value!=NULL && mxIsUint64(value)){
    for(t=0; t<rget_size(y); t++){ y->_mpfr_d[t]=((uint64_t*)mxGetData(value))[t]; }
  }
  else{ mexErrMsgIdAndTxt("MATLAB:rmat2mx","The arg should be Struct with the feild 'digits'."); }  
}

//EOF
