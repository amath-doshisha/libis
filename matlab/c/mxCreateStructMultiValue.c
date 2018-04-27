/**
 * @brief 多倍長型配列からmxArrayの構造体配列への変換の1要素のみの分担
 */
void mxCreateStructMultiValue(mxArray *y, int i, rmulti *x, const char *field_names[])
{
  int t;
  mwSize scalar[2]={1,1},size[2]={1,1};
  mxArray *value=NULL;
  // prec
  value=mxCreateNumericArray(2,scalar,mxINT64_CLASS,mxREAL);
  (*(int64_t*)mxGetData(value))=x->_mpfr_prec;
  mxSetField(y,i,field_names[0],value);
  // sign
  value=mxCreateNumericArray(2,scalar,mxINT32_CLASS,mxREAL);
  (*(int32_t*)mxGetData(value))=x->_mpfr_sign;
  mxSetField(y,i,field_names[1],value);
  // exp
  value=mxCreateNumericArray(2,scalar,mxINT64_CLASS,mxREAL);
  (*(int64_t*)mxGetData(value))=x->_mpfr_exp;
  mxSetField(y,i,field_names[2],value);
  // digits
  size[1]=rget_size(x);
  value=mxCreateNumericArray(2,size,mxUINT64_CLASS,mxREAL);
  for(t=0; t<size[1]; t++){ ((uint64_t*)mxGetData(value))[t]=x->_mpfr_d[t]; }
  mxSetField(y,i,field_names[3],value);
}

//EOF
