
/**
 * @breif 
 */
int subs_index_allocate(int *dims, int **index, const mxArray *s, int M, int N, int L, char *col_end)
{
#define SUBS_2ND_ERROR {MATLAB_ERROR("multi_subsref: The 2nd-arg shoud be Struct of subsref.")}
  char *buf=NULL;
  int ndim=1,i,j;
  mxArray *subs=NULL,*type=NULL;
  // check s
  if(mxGetNumberOfDimensions(s)!=2){ SUBS_2ND_ERROR; }
  if(mxGetDimensions(s)[0]!=1){ SUBS_2ND_ERROR; }
  if(mxGetDimensions(s)[1]!=1){ SUBS_2ND_ERROR; }
  if(!STR_EQ(mxGetFieldNameByNumber(s,0),"type")){ SUBS_2ND_ERROR; }
  if(!STR_EQ(mxGetFieldNameByNumber(s,1),"subs")){ SUBS_2ND_ERROR; }
  // check s.type=="()"
  type=mxGetField(s,0,"type");
  if(type!=NULL && mxIsChar(type)){ buf=mxArrayToString(type); }else{ SUBS_2ND_ERROR; }
  if(!STR_EQ(buf,"()")){ SUBS_2ND_ERROR; }
  if(buf!=NULL){ mxFree(buf); }
  // check s.subs
  subs=mxGetField(s,0,"subs");
  if(subs==NULL || !mxIsCell(subs)){ SUBS_2ND_ERROR; }
  if(mxGetNumberOfDimensions(subs)!=2){ SUBS_2ND_ERROR; }
  if(mxGetDimensions(subs)[0]!=1){ SUBS_2ND_ERROR; }
  // get number of dimensions
  ndim=mxGetDimensions(subs)[1];
  if(ndim>3){ MATLAB_ERROR("multi_subsref: The number of dimensions is bigger than 3.");}
  // get indecies
  for(i=0; i<ndim && i<3; i++){
    // check s.subs(i)
    if(mxGetCell(subs,i)==NULL){ SUBS_2ND_ERROR; }
    // if double
    if(mxIsDouble(mxGetCell(subs,i))){    
      // chek row vector
      if(mxGetNumberOfDimensions(mxGetCell(subs,i))!=2){ SUBS_2ND_ERROR; }
      if(mxGetDimensions(mxGetCell(subs,i))[0]!=1){ SUBS_2ND_ERROR; }
      // get size of indecis of i-th
      dims[i]=(int)mxGetDimensions(mxGetCell(subs,i))[1];
      if(mxGetPr(mxGetCell(subs,i))==NULL){ SUBS_2ND_ERROR; }
      index[i]=ivec_allocate(dims[i]);
      for(j=0; j<dims[i]; j++){
	index[i][j]=(int)(mxGetPr(mxGetCell(subs,i))[j])-1;
      }
      *col_end='e';
    }else if(mxIsChar(mxGetCell(subs,i))){
      buf=mxArrayToString(mxGetCell(subs,i));
      if(!STR_EQ(buf,":")){ SUBS_2ND_ERROR; }
      if(buf!=NULL){ mxFree(buf); }
      if(ndim==1){ dims[i]=M*N*L; }
      else{ if(i==0){ dims[i]=M; }else if(i==1){ dims[i]=N; }else if(i==2){ dims[i]=L; }else { SUBS_2ND_ERROR; } }
      index[i]=ivec_allocate(dims[i]);
      ivec_grid(dims[i],index[i]);

      *col_end='c';
    }else{ SUBS_2ND_ERROR; }
  }
  for(i=ndim; i<3; i++){
    index[i]=ivec_allocate(1);
    index[i][0]=0;
  }
  if(*col_end=='e'){
    // check indecies
    for(i=0; i<ndim; i++){
      for(j=0; j<dims[i]; j++){
	if(ndim==1 && (index[i][j]<0 || index[i][j]>=M*N*L)){ MATLAB_ERROR("multi_subsref: Index is out of range.");}
	if(ndim>1 && i==0 && (index[i][j]<0 || index[i][j]>=M)){ MATLAB_ERROR("multi_subsref: Index is out of range.");}
	if(ndim>1 && i==1 && (index[i][j]<0 || index[i][j]>=N)){ MATLAB_ERROR("multi_subsref: Index is out of range.");}
	if(ndim>1 && i==2 && (index[i][j]<0 || index[i][j]>=L)){ MATLAB_ERROR("multi_subsref: Index is out of range.");}
      }
    }
  }
  return ndim;
}

/**
 * @breif 
 */
void subs_index_free(int** index){
  int i;
  for(i=0; i<3; i++){ index[i]=ivec_free(index[i]); }
  return;
}


/**
 * @breif y=x(s)
 */
void multi_subsref(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  int m=1,n=1,l=1;
  multi *x=NULL,*y=NULL;
  int ndim=3,dims[]={1,1,1},*index[3];
  char col_end='e';//'c'or'e'//colon,end
  if(nlhs>1){ mexErrMsgIdAndTxt("MATLAB:multi_mex:maxlhs","Too many output arguments."); }
  if(!(IS_STRT(nrhs,prhs,N0))){ MATLAB_ERROR("multi_subsref: The 1st-arg should be Struct."); }
  if(!(IS_STRT(nrhs,prhs,N0+1))){ MATLAB_ERROR("multi_subsref: The 2nd-arg should be Struct."); }
  // load the 1st-arg
  x=multi_allocate_mxArray(prhs[N0]);
  // load the 2nd-arg
  ndim=subs_index_allocate(dims,index,prhs[N0+1],_M(x),_N(x),_L(x),&col_end);
  // copy by indecies
  if(col_end=='e'&&_M(x)==1&&_N(x)==1&&_L(x)!=1){    
    y=multi_allocate(_T(x),_LD1(x),_LD2(x),_L(x));
         if(_T(x)=='r'){ rmat3_clone(_LD1(x),_LD2(x),_L(x),_R(y),_LD1(y),_LD2(y),_R(x),_LD1(x),_LD2(x)); }//y=x
    else if(_T(x)=='c'){ cmat3_clone(_LD1(x),_LD2(x),_L(x),_C(y),_LD1(y),_LD2(y),_C(x),_LD1(x),_LD2(x)); }//y=x
    else{ MATLAB_ERROR("multi_subsref: Unkown type"); }
  }else if((col_end=='c')||(col_end=='e'&&_M(x)!=1&&_N(x)==1&&_L(x)==1)){
    y=multi_allocate(_T(x),dims[0],1,1);
         if(_T(x)=='r'){ rvec_clone_index(_M(y),_R(y),_R(x),index[0]); }
    else if(_T(x)=='c'){ cvec_clone_index(_M(y),_C(y),_C(x),index[0]); }
    else{ MATLAB_ERROR("multi_subsref: Unkown type"); }
  }else if(ndim==1){
    y=multi_allocate(_T(x),dims[0],1,1);
         if(_T(x)=='r'){ rvec_clone_index(_M(y),_R(y),_R(x),index[0]); _N(y)=_M(y);_LD2(y)=_LD1(y);_M(y)=1;_LD1(y)=1; }
    else if(_T(x)=='c'){ cvec_clone_index(_M(y),_C(y),_C(x),index[0]); _N(y)=_M(y);_LD2(y)=_LD1(y);_M(y)=1;_LD1(y)=1; }
    else{ MATLAB_ERROR("multi_subsref: Unkown type"); }
  }else{
    y=multi_allocate(_T(x),dims[0],dims[1],dims[2]);
         if(_T(x)=='r'){ rmat3_clone_index(_M(y),_N(y),_L(y),_R(y),_LD1(y),_LD2(y),_R(x),_LD1(x),_LD2(x),index[0],index[1],index[2]); }
    else if(_T(x)=='c'){ cmat3_clone_index(_M(y),_N(y),_L(y),_C(y),_LD1(y),_LD2(y),_C(x),_LD1(x),_LD2(x),index[0],index[1],index[2]); }
    else{ MATLAB_ERROR("multi_subsref: Unkown type"); }
  }
  // done
  plhs[0]=mxCreateStructMulti(y);  
  subs_index_free(index);
  x=multi_free(x);
  y=multi_free(y);
  return;
}

//EOF
