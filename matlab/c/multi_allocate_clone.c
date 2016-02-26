/**
 * @breif y=multi_allocate_clone(x)
 */
multi *multi_allocate_clone(multi *x)
{
  multi *y=NULL;
  // allocate by default precision
  y=multi_allocate(_T(x),_M(x),_N(x),_L(x));
  // copy
       if(_T(y)=='r'){ rmat3_copy(_M(y),_N(y),_L(y),_R(y),_LD1(y),_LD2(y),_R(x),_LD1(x),_LD2(x)); }
  else if(_T(y)=='c'){ cmat3_copy(_M(y),_N(y),_L(y),_C(y),_LD1(y),_LD2(y),_C(x),_LD1(x),_LD2(x)); }
  else{ MATLAB_ERROR("multi_copy: Unkown type"); }
  return y;
}

//EOF
