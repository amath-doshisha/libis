/**
 * @breif A=matgen_dhToda(...)
void riep_EXTdhToda_TN(int m, int N, int M, rmulti **A, int LDA, rmulti **Q, int LDQ, rmulti **E, int LDE, rmulti **lambda, rmulti **c, int debug)
 */
void multi_matgen_dhToda(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  int M,N,m;
  multi *A=NULL,*lambda=NULL,*c;
  if(nlhs!=1)                   { mexErrMsgIdAndTxt("MATLAB:multi_mex:maxlhs","Too many output arguments."); }
  if(!(IS_NUMR(nrhs,prhs,N0)))  { MATLAB_ERROR("multi_matgen_dhToda: The 1st-arg should be scalar."); }
  if(!(IS_NUMR(nrhs,prhs,N0+1))){ MATLAB_ERROR("multi_matgen_dhToda: The 2nd-arg should be scalar."); }
  if(!(IS_STRT(nrhs,prhs,N0+2))){ MATLAB_ERROR("multi_matgen_dhToda: The 3rd-arg should be Struct."); }
  if(!(IS_STRT(nrhs,prhs,N0+3))){ MATLAB_ERROR("multi_matgen_dhToda: The 4th-arg should be Struct."); }
  // allocate by clone
  N=GET_DOUBLE(prhs[N0])[0];
  M=GET_DOUBLE(prhs[N0+1])[0];
  lambda=multi_allocate_mxArray(prhs[N0+2]);
  c=multi_allocate_mxArray(prhs[N0+3]);
  // check size
  m=_M(lambda)*_N(lambda)*_L(lambda);
  if(m<0){ MATLAB_ERROR("multi_matgen_dhToda: Illegal size."); }  
  // allocate
  A=multi_allocate('r',m,m,1);
  // operations
  riep_EXTdhToda_TN(m,N,M,_R(A),_LD1(A),NULL,0,NULL,0,_R(lambda),_R(c),0);
  plhs[0]=mxCreateStructMulti(A);
  // free
  A=multi_free(A);
  lambda=multi_free(lambda);
  c=multi_free(c);
  return;
}

//EOF
