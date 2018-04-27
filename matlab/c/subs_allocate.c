
/**
 * @brief 
 */
subs_index_t *subs_index_allocate(const mxArray *s)
{
#define SUBS_2ND_ERROR {MATLAB_ERROR("multi_subsref: The 2nd-arg shoud be Struct of subsref.")}
  char *buf=NULL;
  int i,j;
  mxArray *subs=NULL,*type=NULL;
  subs_index_t *ret=NULL;
  // check s
  if(!mxIsStruct(s)){ SUBS_2ND_ERROR; }  
  if(mxGetNumberOfDimensions(s)!=2){ SUBS_2ND_ERROR; }
  if(mxGetDimensions(s)[0]!=1){ SUBS_2ND_ERROR; }
  if(mxGetDimensions(s)[1]!=1){ SUBS_2ND_ERROR; }
  if(!STR_EQ(mxGetFieldNameByNumber(s,0),"type")){ SUBS_2ND_ERROR; }
  if(!STR_EQ(mxGetFieldNameByNumber(s,1),"subs")){ SUBS_2ND_ERROR; }
  // allocate
  ret=malloc(sizeof(subs_index_t)*1);
  ret->ndim=0;
  ret->dims=ivec_allocate(3);
  ivec_zeros(3,ret->dims);
  ret->index=malloc(sizeof(int*)*3);
  for(i=0; i<3; i++){ ret->index[i]=NULL; }
  //  s = 
  //    type: '()'
  //    subs: {':'} or {[1,2,3],':'} etc..
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
  ret->ndim=mxGetDimensions(subs)[1];
  if(ret->ndim>3){ MATLAB_ERROR("multi_subsref: The number of dimensions is bigger than 3.");}
  // get indecies
  for(i=0; i<ret->ndim && i<3; i++){
    // check s.subs(i)
    if(mxGetCell(subs,i)==NULL){ SUBS_2ND_ERROR; }
    // if double
    if(mxIsDouble(mxGetCell(subs,i))){    
      // chek row vector
      if(mxGetNumberOfDimensions(mxGetCell(subs,i))!=2){ SUBS_2ND_ERROR; }
      if(mxGetDimensions(mxGetCell(subs,i))[0]!=1){ SUBS_2ND_ERROR; }
      // get size of indecis
      ret->dims[i]=(int)mxGetDimensions(mxGetCell(subs,i))[1];
      // get indecies
      if(mxGetPr(mxGetCell(subs,i))==NULL){ SUBS_2ND_ERROR; }
      ret->index[i]=ivec_allocate(ret->dims[i]);      
      for(j=0; j<ret->dims[i]; j++){ ret->index[i][j]=(int)(mxGetPr(mxGetCell(subs,i))[j])-1; }
    }else if(mxIsChar(mxGetCell(subs,i))){
      // check ':'
      buf=mxArrayToString(mxGetCell(subs,i));
      if(!STR_EQ(buf,":")){ SUBS_2ND_ERROR; }
      if(buf!=NULL){ mxFree(buf); }
      ret->dims[i]=0;
      ret->index[i]=NULL;
    }else{ SUBS_2ND_ERROR; }
  }
  // done
  return ret;
}

/**
 * @breif 
 */
void subs_index_free(subs_index_t *s){
  int i;
  if(s==NULL){ return; }
  s->dims=ivec_free(s->dims);
  for(i=0; i<3; i++){ s->index[i]=ivec_free(s->index[i]); }
  if(s->index!=NULL){ free(s->index); }
  if(s!=NULL){ free(s); }
  return;
}
