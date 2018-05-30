#include<isys.h>
#include <time.h>


#define PROG             "mateig_krawczyk"
#define BIN_SUFFIX       "dat"
#define TXT_SUFFIX       "txt"
#define DEFAULT_PREC     128

#define RMA(A,M,N) { A=rmat_allocate(M,N); }
#define RMF(A,M,N) { A=rmat_free(M,N,A); }
#define RVA(X,N)   { X=rvec_allocate(N); }
#define RVF(X,N)   { X=rvec_free(N,X); }
#define RA(X)      { X=rallocate(); }
#define RF(X)      { X=rfree(X); }
#define CMA(A,M,N) { A=cmat_allocate(M,N); }
#define CMF(A,M,N) { A=cmat_free(M,N,A); }
#define CMR(M,N,A,LDA,B) { cmat_round((M),(N),(A),(LDA),(B)); }
#define CVA(X,N)   { X=cvec_allocate(N); }
#define CVF(X,N)   { X=cvec_free(N,X); }
#define CVR(N,X,B) { cvec_round((N),(X),(B)); }
#define CA(X)      { X=callocate(); }
#define CF(X)      { X=cfree(X); }
#define RR(X,B)    { rround((X),(B)); }



void usage()
{
  printf("Usage: %s input1.dat input2.dat input3.dat [options]\n",PROG);
  printf("  input1.dat         Binary data file of the eigenvalues.\n");
  printf("  input2.dat         Binary data file of the eigenvectors.\n");
  printf("  input3.dat         Binary data file of the target matrix.\n");
  printf("\n");
  printf("  -prec              Set precition as compute.\n");
  printf("  -real              Set real mode.\n");
  printf("  -complex           Set complex.\n");
  printf("  -o-result fname    The result is saved in the file 'fname'.\n");
  printf("  -o-eigval fname    The error bounds of the eigenvalues are saved in the file 'fname'.\n");
  printf("  -o-eigvec fname    The error bounds of the eigenvectors are saved in the file 'fname'.\n");
  printf("  -O                 The matrix is saved in the file whose name is automatically determined.\n");
  printf("  -d dir             Set the directory of output file as 'dir'.\n");
  printf("  -f num             Set output format as 'f' with num digits.\n");
  printf("  -e num             Set output format as 'e' with num digits.\n");
  printf("  -nv                No verbose mode.\n");
  printf("  -v                 Verbose mode, level=1.\n");
  printf("  -vv                Verbose mode, level=2.\n");
  printf("  -vvv               Verbose mode, level=3.\n");
  printf("  -help              This message.\n");
  printf("\n");
  printf("Examples:\n");
  printf("     # %s A_eigval.dat A_eigvec.dat A.dat\n",PROG);
  printf("\n");
  exit(0);
}

int main(int argc, char *argv[])
{
  FILE *fid=NULL;
  strings *path=NULL,*path2=NULL;
  char *in_fname_eigval=NULL,*in_fname_eigvec=NULL,*in_fname_matrix=NULL;
  char *out_fname_result=NULL,*out_fname_eigval=NULL,*out_fname_eigvec=NULL,*out_dir=NULL,buf[10000];
  char format[]="f",mode='r';
  int debug=0,LDA=0,LDX=0,LDE=0,n_A=0,n_L=0,n_X=0,ret=-1,i,n=0,k=0,ie_max=0;
  int prec=DEFAULT_PREC,autoname=0,digits=2;
  rmulti **rA=NULL,**rLambda=NULL,**rX=NULL,**rE_eigval=NULL,**rE_eigvec=NULL,*E_max=NULL;
  cmulti **cA=NULL,**cLambda=NULL,**cX=NULL,**cE_eigval=NULL,**cE_eigvec=NULL;
  int *iE_eigval=NULL,*iE_eigvec=NULL,*iLambda=NULL,*iX=NULL,iE_eigval_max=0;
  int iLambda_max=0,iX_max=0;
  rmulti  *Lambda_max=NULL,*X_max=NULL;
  clock_t t1,t2;


  t1 = clock();
  // init
  in_fname_eigval=char_new("",NULL);
  in_fname_eigvec=char_new("",NULL);
  in_fname_matrix=char_new("",NULL);
  out_fname_result=char_new("",NULL);
  out_fname_eigval=char_new("",NULL);
  out_fname_eigvec=char_new("",NULL);
  out_dir=char_new("",NULL);

  // options
  i=1;
  while(i<argc){
    if     (STR_EQ(argv[i],"-help"))                 { usage(); }
    else if(STR_EQ(argv[i],"-prec") && i+1<argc)     { prec=atoi(argv[++i]); }
    else if(STR_EQ(argv[i],"-o-result") && i+1<argc) { out_fname_result=char_renew(out_fname_result,argv[++i],NULL); }
    else if(STR_EQ(argv[i],"-o-eigval") && i+1<argc) { out_fname_eigval=char_renew(out_fname_eigval,argv[++i],NULL); }
    else if(STR_EQ(argv[i],"-o-eigvec") && i+1<argc) { out_fname_eigvec=char_renew(out_fname_eigvec,argv[++i],NULL); }
    else if(STR_EQ(argv[i],"-O"))                    { autoname=1; }
    else if(STR_EQ(argv[i],"-d") && i+1<argc)        { out_dir=char_renew(out_dir,argv[++i],NULL); if(out_dir[strlen(out_dir)-1]!='/'){ out_dir=char_cat(out_dir,"/"); } }
    else if(STR_EQ(argv[i],"-f") && i+1<argc)        { strcpy(format,"f"); digits=atoi(argv[++i]); }
    else if(STR_EQ(argv[i],"-e") && i+1<argc)        { strcpy(format,"e"); digits=atoi(argv[++i]); }
    else if(STR_EQ(argv[i],"-real"))                 { mode='r'; }
    else if(STR_EQ(argv[i],"-complex"))              { mode='c'; }
    else if(STR_EQ(argv[i],"-nv"))                   { debug=0; }
    else if(STR_EQ(argv[i],"-v"))                    { debug=1; }
    else if(STR_EQ(argv[i],"-vv"))                   { debug=2; }
    else if(STR_EQ(argv[i],"-vvv"))                  { debug=3; }
    else if(STR_EQ(argv[i],"-vvvv"))                 { debug=4; }
    else if(STR_EQ(argv[i],"-v5"))                   { debug=-1; }
    else if(strlen(in_fname_eigval)<=0)              { in_fname_eigval=char_renew(in_fname_eigval,argv[i],NULL); }
    else if(strlen(in_fname_eigvec)<=0)              { in_fname_eigvec=char_renew(in_fname_eigvec,argv[i],NULL); }
    else if(strlen(in_fname_matrix)<=0)              { in_fname_matrix=char_renew(in_fname_matrix,argv[i],NULL); }
    else if(STR_EQ_N(argv[i],"-",1))                 { usage(); }
    else                                             { usage(); }
    i++;
  }

  // check
  if(strlen(in_fname_eigval)<=0 || strlen(in_fname_eigvec)<=0 || strlen(in_fname_matrix)<=0){ usage(); exit(0); }

  // set default precision
  set_default_prec(prec);

  // autoname
  if(autoname){
    // split path
    path=strings_split_path(in_fname_matrix,NULL);
    // dir
    if(strlen(out_dir)<=0 && strlen(path->str[0])>0){ out_dir=char_renew(out_dir,path->str[0],NULL); out_dir=char_cat(out_dir,"/"); }
    // basename
    out_fname_result=char_renew(out_fname_result,path->str[1],NULL);
    out_fname_eigval=char_renew(out_fname_eigval,path->str[1],NULL);
    out_fname_eigvec=char_renew(out_fname_eigvec,path->str[1],NULL);
    if(!STR_EQ(path->str[2],BIN_SUFFIX) && strlen(path->str[2])>0){
      out_fname_result=char_cat(out_fname_result,"."); out_fname_result=char_cat(out_fname_result,path->str[2]);
      out_fname_eigval=char_cat(out_fname_eigval,"."); out_fname_eigval=char_cat(out_fname_eigval,path->str[2]);
      out_fname_eigvec=char_cat(out_fname_eigvec,"."); out_fname_eigvec=char_cat(out_fname_eigvec,path->str[2]);
    }
    // split path-eigval
    path=strings_del(path);
    path=strings_split_path(in_fname_eigval,NULL);
    path2=strings_split(path->str[1],"_","{","}",NULL);
    if((path2->n)>0){
      // basename
      out_fname_result=char_cat(out_fname_result,",");
      out_fname_eigval=char_cat(out_fname_eigval,",");
      out_fname_eigvec=char_cat(out_fname_eigvec,",");
      out_fname_result=char_cat(out_fname_result,path2->str[(path2->n)-1]);
      out_fname_eigval=char_cat(out_fname_eigval,path2->str[(path2->n)-1]);
      out_fname_eigvec=char_cat(out_fname_eigvec,path2->str[(path2->n)-1]);
    }
    path=strings_del(path);
    path2=strings_del(path2);
    // split path-eigvec
    path=strings_del(path);
    path=strings_split_path(in_fname_eigvec,NULL);
    path2=strings_split(path->str[1],"_","{","}",NULL);
    if((path2->n)>0){
      // basename
      out_fname_result=char_cat(out_fname_result,",");
      out_fname_eigval=char_cat(out_fname_eigval,",");
      out_fname_eigvec=char_cat(out_fname_eigvec,",");
      out_fname_result=char_cat(out_fname_result,path2->str[(path2->n)-1]);
      out_fname_eigval=char_cat(out_fname_eigval,path2->str[(path2->n)-1]);
      out_fname_eigvec=char_cat(out_fname_eigvec,path2->str[(path2->n)-1]);
    }
    path=strings_del(path);
    path2=strings_del(path2);
    // append
    sprintf(buf,"_text{%dbits|%ckrawczyk|result}.%s",prec,mode,TXT_SUFFIX); out_fname_result=char_cat(out_fname_result,buf);
    sprintf(buf,"_rvec{%dbits|%ckrawczyk|eigval}.%s",prec,mode,BIN_SUFFIX); out_fname_eigval=char_cat(out_fname_eigval,buf);
    sprintf(buf,"_rmat{%dbits|%ckrawczyk|eigvec}.%s",prec,mode,BIN_SUFFIX); out_fname_eigvec=char_cat(out_fname_eigvec,buf);
  }

  // load input-matrix
  rA=rmat_bin_load(&LDA,&n_A,"%s",in_fname_matrix);
  if(rA==NULL){ cA=cmat_bin_load(&LDA,&n_A,"%s",in_fname_matrix); }
  if(rA==NULL && cA==NULL){ printf("\nError! The input file '%s' can NOT be opened, or is NOT supported type.\n\n",in_fname_matrix); usage(); exit(0); }
  if(LDA!=n_A){ printf("\nError! The input matrix should be square. The size of '%s' is %d-by-%d.\n\n",in_fname_matrix,LDA,n_A); usage(); exit(0); }
  if(debug!=-1){
    if(strlen(out_fname_result)<=0 || debug>0){ printf("Input_Matrix_File: %s\n",in_fname_matrix); }
  }
  if(strlen(out_fname_result)<=0 || debug>1){
    if(debug!=-1){
      printf("Input_Matrix_Size: %d-by-%d\n",LDA,n_A);
    }
    if(rA!=NULL){
      if(debug!=-1){
      printf("Input_Matrix_Type: rmat\n");
      printf("Input_Matrix:\n");
      rmat_print(LDA,n_A,rA,LDA,NULL,format[0],digits);
      }
    }
    if(cA!=NULL){
      if(debug!=-1){
      printf("Input_Matrix_Type: cmat\n");
      printf("Input_Matrix:\n");
      cmat_print(LDA,n_A,cA,LDA,NULL,format[0],digits);
      }
    }
  }

  // load input-eigenvalues
  if(cLambda==NULL && mode=='r') { rLambda=rvec_bin_load(&n_L,"%s",in_fname_eigval); }
  if(rLambda==NULL && mode=='c') { cLambda=cvec_bin_load(&n_L,"%s",in_fname_eigval); }
  if(rLambda==NULL && cLambda==NULL){ printf("\nError! The input file '%s' can NOT be opened, or is NOT supported type.\n\n",in_fname_eigvec); usage(); exit(0); }
  if(debug!=-1){
  if(strlen(out_fname_result)<=0 || debug>0){ printf("Input_Eigenvalues_File: %s\n",in_fname_eigval); }
  }
  if(strlen(out_fname_result)<=0 || debug>1){
    if(debug!=-1){
      printf("Input_Eigenvalues_Size: %d\n",n_L);
    }
    if(rLambda!=NULL){
      if(debug!=-1){
      printf("Input_Eigenvalues_Type: rvec\n");
      printf("Input_Eigenvalues:\n");
      rvec_print(n_L,rLambda,NULL,format[0],digits);
      }
    }
    if(cLambda!=NULL){
      if(debug!=-1){
      printf("Input_Eigenvalues_Type: cvec\n");
      printf("Input_Eigenvalues:\n");
      cvec_print_prec(n_L,cLambda,NULL,format,digits);
      }
    }
  }

  // load input-eigenvectors
  rX=rmat_bin_load(&LDX,&n_X,"%s",in_fname_eigvec);
  if(rX==NULL){ cX=cmat_bin_load(&LDX,&n_X,"%s",in_fname_eigvec); }
  if(rX==NULL && cX==NULL){ printf("\nError! The input file '%s' can NOT be opened, or is NOT supported type.\n\n",in_fname_eigval); usage(); exit(0); }
  if(debug!=-1){
    if(strlen(out_fname_result)<=0 || debug>0){
      printf("Input_Eigenvectors_File: %s\n",in_fname_eigvec);
    }
  }
  if(strlen(out_fname_result)<=0 || debug>1){
    if(debug!=-1){
    printf("Input_Eigenvectors_Size: %d-by-%d\n",LDX,n_X);
    }
    if(rX!=NULL){
      if(debug!=-1){
	printf("Input_Eigenvectors_Type: rmat\n");
	printf("Input_Eigenvectors:\n");
	rmat_print(LDX,n_X,rX,LDX,NULL,format[0],digits);
      }
    }
    if(cX!=NULL){
      if(debug!=-1){
	printf("Input_Eigenvectors_Type: cmat\n");
      printf("Input_Eigenvectors:\n");
      cmat_print(LDX,n_X,cX,LDX,NULL,format[0],digits);
      }
    }
  }

  // check sizes
  if(LDA!=n_A || n_A!=LDX || n_L!=n_X){ printf("\nError! The sizes of the input files are illeagal.\n\n"); exit(0); }
  n=n_A;
  k=n_X;
  if(debug!=-1){
    if(strlen(out_fname_result)<=0 || debug>1){
      printf("Size: %d\n",n);
      printf("Eigenpairs: %d\n",k);
    }
  }

  // allocate
  E_max=rallocate();
  rset_zero(E_max);


  // compute
  if(mode=='r'){
    // real mode
    LDE=n;
    rE_eigvec=rmat_allocate(LDE,k);
    iE_eigvec=imat_allocate(LDE,k);
    iX=imat_allocate(LDX,n_X);
    rE_eigval=rvec_allocate(k);
    iE_eigval=ivec_allocate(k);
    iLambda=ivec_allocate(k);
    Lambda_max=rallocate();
    X_max=rallocate();
    if(cA!=NULL && rA==NULL){
      rA=rmat_allocate(LDA,n_A);
      cmat_real_clone(LDA,n_A,rA,LDA,cA,LDA);
      cA=cmat_free(LDA,n_A,cA);      
    }
    if(cX!=NULL && rX==NULL){
      rX=rmat_allocate(LDX,n_X);
      cmat_real_clone(LDX,n_X,rX,LDX,cX,LDX);
      cX=cmat_free(LDX,n_X,cX);
    }
    if(cLambda!=NULL && rLambda==NULL){
      rLambda=rvec_allocate(n_L);
      rvec_real_cvec_clone(n_L,rLambda,cLambda);
      cLambda=cvec_free(n_L,cLambda);
    }
    // compute
    if(rA!=NULL && rX!=NULL && rLambda!=NULL){
      ret=ireig_krawczyk(n,k,rE_eigvec,LDE,rE_eigval,rA,LDA,rX,LDX,rLambda,debug-2);
      rmax_abs_rvec(Lambda_max,k,rLambda);
      iLambda_max=rget_exp(Lambda_max);
      ivec_set_all(k,iLambda,iLambda_max);
      rvec_get_exp(k,iE_eigval,rE_eigval);
      ivec_sub(k,iE_eigval,iLambda);
      ivec_add_scalar(k,iE_eigval,1);
      rmat_max_abs(X_max,LDX,n_X,rX,LDX);
      iX_max=rget_exp(X_max);
      imat_set(LDX,n_X,iX,LDX,iX_max);
      rmat_get_exp(LDE,k,iE_eigvec,LDE,rE_eigvec,LDE);
      imat_sub(n,k,iE_eigvec,LDE,iX,LDX);
      imat_add_scalar(n,k,iE_eigvec,LDE,1);
      rmax_rvec(E_max,n,rE_eigval);       
      /*     rvec_get_exp(k,iLambda,rLambda);
      rmat_get_exp(LDX,n_X,iX,LDX,rX,LDX);
      */
      iE_eigval_max=ivec_max(n,iE_eigval);       
    }    
  }else{
    // complex mode
    LDE=n;
    cE_eigvec=cmat_allocate(LDE,k);
    iE_eigvec=imat_allocate(LDE,k);
    iX=imat_allocate(LDX,n_X);
    cE_eigval=cvec_allocate(k);
    iE_eigval=ivec_allocate(k);
    iLambda=ivec_allocate(k);
    Lambda_max=rallocate();
    X_max=rallocate();
    if(cA==NULL && rA!=NULL){
      cA=cmat_allocate(LDA,n_A);
      cmat_clone_r(LDA,n_A,cA,LDA,rA,LDA);
      rA=rmat_free(LDA,n_A,rA);      
    }
    if(cX==NULL && rX!=NULL){
      cX=cmat_allocate(LDX,n_X);
      cmat_clone_r(LDX,n_X,cX,LDX,rX,LDX);
      rX=rmat_free(LDX,n_X,rX);      
    }
    if(cLambda==NULL && rLambda!=NULL){
      cLambda=cvec_allocate(n_L);
      cvec_clone_rvec(n_L,cLambda,rLambda);
      rLambda=rvec_free(n_L,rLambda);      
    }
    // compute
    if(cA!=NULL && cX!=NULL && cLambda!=NULL){
      ret=iceig_krawczyk(n,k,cE_eigvec,LDE,cE_eigval,cA,LDA,cX,LDX,cLambda,debug-2);
      rmax_abs_cvec(Lambda_max,k,cLambda);
      iLambda_max=rget_exp(Lambda_max);
      ivec_set_all(k,iLambda,iLambda_max);
      cvec_get_exp(k,iE_eigval,cE_eigval);
      ivec_sub(k,iE_eigval,iLambda);
      ivec_add_scalar(k,iE_eigval,1);
      cmat_max_abs(X_max,LDX,n_X,cX,LDX);
      iX_max=rget_exp(X_max);
      imat_set(LDX,n_X,iX,LDX,iX_max);
      cmat_get_exp(LDE,k,iE_eigvec,LDE,cE_eigvec,LDE);
      imat_sub(n,k,iE_eigvec,LDE,iX,LDX);
      imat_add_scalar(n,k,iE_eigvec,LDE,1);
      
      /*cvec_get_exp(k,iLambda,cLambda);
      cvec_get_exp(k, iE_eigval,cE_eigval);
      ivec_sub(k,iE_eigval,iLambda);
      ivec_add_scalar(k,iE_eigval,1);
      cmat_get_exp(LDX,n_X,iX,LDX,cX,LDX);
      cmat_get_exp(LDE,k,iE_eigvec,LDE,cE_eigvec,LDE);
      imat_sub(n,k,iE_eigvec,LDE,iX,LDX);
      imat_add_scalar(n,k,iE_eigvec,LDE,1);
      */
      }
  }
  // output 
  if(mode=='r'){
    //real mode
    if(debug!=-1){
      if(strlen(out_fname_result)>0 && debug>1){ 
      printf("Absolute_Error_Bounds_Eigenvalues:\n");
      rvec_print(k,rE_eigval,NULL,'e',3);
      printf("Relative_Error_Bounds_Log2_Eigenvalues:\n");
      imat_print(k,1,iE_eigval,k,NULL);
      printf("Relative_Error_Bounds_Log2_Eigenvalues_max:\t");
      printf("%d\n",iE_eigval_max);
      printf("Absolute_Error_Bounds_Eigenvectors:\n");
      rmat_print(LDE,k,rE_eigvec,LDE,NULL,'e',3);
      printf("Relative_Error_Bounds_Log2_Eigenvectors:\n");
      imat_print(LDE,k,iE_eigvec,LDE,NULL);
    }
      if(strlen(out_fname_result)<=0 || debug>0){
	printf("MPFR_Precision: %d\n",prec);
	printf("Result: ");
	if(ret==0){ print_green(); printf("succeeded"); }
      else      { print_red();   printf("failed"); }
	print_reset(); printf("\n");
      }
      if(strlen(out_fname_result)>0 && debug>0){ printf("Output_File_Result:       %s%s\n",out_dir,out_fname_result); }
      if(strlen(out_fname_eigval)>0 && debug>0){ printf("Output_File_Eigenvalues:  %s%s\n",out_dir,out_fname_eigval); }
      if(strlen(out_fname_eigvec)>0 && debug>0){ printf("Output_File_Eigenvectors: %s%s\n",out_dir,out_fname_eigvec); }
    }else{
      mpfr_printf("%d\t%.2Re\n",n,E_max);
    }
  }else{
    //complex mode
    if(strlen(out_fname_result)>0 && debug>1){ 
      printf("Absolute_Error_Bounds_Eigenvalues:\n");
      cvec_print(k,cE_eigval,NULL,'e',3);
      printf("Relative_Error_Bounds_Log2_Eigenvalues:\n");
      imat_print(k,1,iE_eigval,k,NULL);
      printf("Absolute_Error_Bounds_Eigenvectors:\n");
      cmat_print(LDE,k,cE_eigvec,LDE,NULL,'e',3);
      printf("Relative_Error_Bounds_Log2_Eigenvectors:\n");
      imat_print(LDE,k,iE_eigvec,LDE,NULL);
    }
    if(strlen(out_fname_result)<=0 || debug>0){
      printf("MPFR_Precision: %d\n",prec);
      printf("Result: ");
      if(ret==0){ print_green(); printf("succeeded"); }
      else      { print_red();   printf("failed"); }
      print_reset(); printf("\n");
    }
    if(strlen(out_fname_result)>0 && debug>0){ printf("Output_File_Result:       %s%s\n",out_dir,out_fname_result); }
    if(strlen(out_fname_eigval)>0 && debug>0){ printf("Output_File_Eigenvalues:  %s%s\n",out_dir,out_fname_eigval); }
    if(strlen(out_fname_eigvec)>0 && debug>0){ printf("Output_File_Eigenvectors: %s%s\n",out_dir,out_fname_eigvec); }
    if(debug==-1){
      ie_max=ivec_max(k,iE_eigval);
      printf("%d\t%d\n",prec,ie_max);
    }
 }

  // save
  if(strlen(out_fname_result)>0){
    sprintf(buf,"%s%s",out_dir,out_fname_result);
    fid=fopen(buf,"w");
    if(ret==0){ fprintf(fid,"succeeded\n"); }
    else      { fprintf(fid,"failed\n"); }
    fclose(fid);
  }
  if(mode=='r'){
    if(strlen(out_fname_eigval)>0){
      rvec_bin_save(k,rE_eigval,"%s%s",out_dir,out_fname_eigval);
    }
    if(strlen(out_fname_eigvec)>0){
      rmat_bin_save(LDE,k,rE_eigvec,LDE,"%s%s",out_dir,out_fname_eigvec);
    }
  }else{
    if(strlen(out_fname_eigval)>0){
      cvec_bin_save(k,cE_eigval,"%s%s",out_dir,out_fname_eigval);
    }
    if(strlen(out_fname_eigvec)>0){
      cmat_bin_save(LDE,k,cE_eigvec,LDE,"%s%s",out_dir,out_fname_eigvec);
    }
  }
  // done
  rA=rmat_free(LDA,n_A,rA);
  cA=cmat_free(LDA,n_A,cA);
  cLambda=cvec_free(n_L,cLambda);
  rLambda=rvec_free(n_L,rLambda);
  rX=rmat_free(n,n,rX);
  cX=cmat_free(n,n,cX);
  rE_eigvec=rmat_free(LDE,k,rE_eigvec);
  cE_eigvec=cmat_free(LDE,k,cE_eigvec);
  iE_eigvec=imat_free(iE_eigvec);
  rE_eigval=rvec_free(k,rE_eigval);
  cE_eigval=cvec_free(k,cE_eigval);
  iE_eigval=ivec_free(iE_eigval);
  iLambda=ivec_free(iLambda);
  iX=imat_free(iX);
  Lambda_max=rmfree(Lambda_max);
  X_max=rmfree(X_max);
  E_max=rmfree(E_max);
  in_fname_eigval=char_del(in_fname_eigval);
  in_fname_eigvec=char_del(in_fname_eigvec);
  in_fname_matrix=char_del(in_fname_matrix);
  out_fname_result=char_del(out_fname_result);
  out_fname_eigval=char_del(out_fname_eigval);
  out_fname_eigvec=char_del(out_fname_eigvec);
  out_dir=char_del(out_dir);

  t2 = clock();
  printf("実行時間:%f\n",(double)(t2 - t1) / CLOCKS_PER_SEC);

  return 0;
}
