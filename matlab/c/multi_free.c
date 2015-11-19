/**
 * @breif 多倍長型のメモリの解放
 */
multi *multi_free(multi *A)
{
  if(A!=NULL && A->type=='r'){ A->p0=(void*)rmat3_free(A->LD1,A->LD2,A->L,(rmulti**)(A->p0)); free(A); A=NULL; }
  if(A!=NULL && A->type=='c'){ A->p0=(void*)cmat3_free(A->LD1,A->LD2,A->L,(cmulti**)(A->p0)); free(A); A=NULL; }
  if(A!=NULL && A->type=='R'){ A->p0=(void*)rmat3_free(A->LD1,A->LD2,A->L,(rmulti**)(A->p0)); A->p1=(void*)rmat3_free(A->LD1,A->LD2,A->L,(rmulti**)(A->p1)); free(A); A=NULL; }
  if(A!=NULL && A->type=='C'){ A->p0=(void*)cmat3_free(A->LD1,A->LD2,A->L,(cmulti**)(A->p0)); A->p1=(void*)cmat3_free(A->LD1,A->LD2,A->L,(cmulti**)(A->p1)); free(A); A=NULL; }
  return A;
}

//EOF
