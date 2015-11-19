/**
 * @breif 多倍長型のメモリの割当
 */
multi *multi_allocate(char type, int m, int n, int l)
{
  multi *A=NULL;
  A=(multi*)malloc(sizeof(multi_struct)*1);
  A->type=type;
  A->M=m;  A->LD1=m;
  A->N=n;  A->LD2=n;
  A->L=l;
       if(type=='r'){ A->p0=(void*)rmat3_allocate(_LD1(A),_LD2(A),_L(A)); A->p1=NULL;  }
  else if(type=='c'){ A->p0=(void*)cmat3_allocate(_LD1(A),_LD2(A),_L(A)); A->p1=NULL;  }
  else if(type=='R'){ A->p0=(void*)rmat3_allocate(_LD1(A),_LD2(A),_L(A)); A->p1=(void*)rmat3_allocate(_LD1(A),_LD2(A),_L(A));  }
  else if(type=='C'){ A->p0=(void*)cmat3_allocate(_LD1(A),_LD2(A),_L(A)); A->p1=(void*)cmat3_allocate(_LD1(A),_LD2(A),_L(A));  }
  else              { MATLAB_ERROR("multi_allocate: Unknow type."); }
  return A;
}

//EOF
