#include<isys.h>
#include <time.h>


#define PROG                 "mateig"
#define TXT_SUFFIX           "txt"
#define BIN_SUFFIX           "dat"
#define METHOD_NAME_LENGTH   128
#define METHOD_RHPEIG        "rhpeig"
#define METHOD_CHPEIG        "chpeig"
#define METHOD_CHPEIG_VERIFY "chpeig-verify"
#define METHOD_RHQR          "rhqr"
#define METHOD_CHQR          "chqr"
#define METHOD_RHQR_MT       "rhqr-mt"
#define METHOD_CHQR_MT       "chqr-mt"
//#define DEFAULT_PREC       1024
#define DEFAULT_PREC         53
#define KRAWCZYK_PREC_RATIO  2

void usage()
{
  printf("Usage: %s input.dat [options]\n",PROG);
  printf("  -i-qr               Set the initial vectors by computing QR method.\n");
  printf("  -i-rand             Set the initial vectors randomly.\n");
  printf("  -i-vec fname        Set the initial vectors from the data file of eigenvectors.\n");
  printf("  -i-val fname        Set the initial vectors from the data file of eigenvalues.\n");
  printf("  -rhpeig             Set the method as hyperplane constrained method on real field.\n");
  printf("  -chpeig             Set the method as hyperplane constrained method on complex field.\n");
  printf("  -chpeig-verify prec Set the method as hyperplane constrained method on complex field with verification.\n");
  printf("  -rhqr               Set the method as Hessenberg-type QR method on real field.\n");
  printf("  -chqr               Set the method as Hessenberg-type QR method on complex field.\n");
  printf("  -prec num           Set the precision used in output mode.\n");
  printf("  -k-yes              Set the error bounds mode as ON.\n");
  printf("  -k-no               Set the error bounds mode as OFF.\n");
  printf("  -no                 No output.\n");
  printf("  -o-log fname        The log is saved in the file 'fname'.\n");
  printf("  -o-val fname        The matrix of eigenvectors is saved in the file 'fname'.\n");
  printf("  -o-vec fname        The vector of eigenvalues is saved in the file 'fname'.\n");
  printf("  -o-val-e fname      The error bounds for matrix of eigenvectors is saved in the file 'fname'.\n");
  printf("  -o-vec-e fname      The error bounds for vector of eigenvalues is saved in the file 'fname'.\n");
  printf("  -O                  The matrix is saved in the file whose name is automatically determined.\n");
  printf("  -d dir              Set the directory of output file as 'dir'.\n");
  printf("  -f num              Set output format as 'f' with num digits.\n");
  printf("  -e num              Set output format as 'e' with num digits.\n");
  printf("  -nv                 No verbose mode.\n");
  printf("  -v                  Verbose mode, level=1.\n");
  printf("  -vv                 Verbose mode, level=2.\n");
  printf("  -vvv                Verbose mode, level=3.\n");
  printf("  -vvvv               Verbose mode, level=4.\n");
  printf("  -help               This message.\n");
  printf("\n");
  printf("Examples:\n");
  printf("     # %s input.dat\n",PROG);
  printf("\n");
  exit(0);
}

int main(int argc, char *argv[])
{
  FILE *fid=NULL;
  strings *path=NULL;
  char *in_fname=NULL,*init_fname=NULL,mode_init='Q',mode='?';
  char *out_fname_log=NULL,*out_fname_val=NULL,*out_fname_vec=NULL,*out_fname_mat=NULL,*out_fname_val_e=NULL,*out_fname_vec_e=NULL,*out_dir=NULL;
  char format[]="e",method[METHOD_NAME_LENGTH+1];
  int n=0,k=0,debug=2,i,mX=0,nX=0,mA=0,nA=0,prec=DEFAULT_PREC,autoname=0,digits=6,kmode=1,kprec=-1,kret=1,prec_verify=128;
  cmulti **cA=NULL,**cX=NULL,**cL=NULL,**cXe=NULL,**cLe=NULL;
  rmulti **rA=NULL,**rX=NULL,**rL=NULL,**rXe=NULL,**rLe=NULL,*rXe_max=NULL,*rXe_max_log2=NULL,*rLe_max=NULL,*rLe_max_log2=NULL,*Lmax=NULL;
  cmulti **cB=NULL;
  rmulti **rB=NULL;
  clock_t clock_start,clock_end;
  double time_hpeig=0,time_qr=0,time_ii=0,time_total=0,time_k=0;

  // init
  in_fname=char_new("",NULL);
  init_fname=char_new("",NULL);
  out_fname_log=char_new("",NULL);
  out_fname_val=char_new("",NULL);
  out_fname_val_e=char_new("",NULL);
  out_fname_vec=char_new("",NULL);
  out_fname_vec_e=char_new("",NULL);
  out_fname_mat=char_new("",NULL);
  out_dir=char_new("",NULL);
  strcpy(method,METHOD_CHPEIG_VERIFY);

  // options
  i=1;
  while(i<argc){
    if     (STR_EQ(argv[i],"-help"))                { usage(); }
    else if(STR_EQ(argv[i],"-i-qr"))                { mode_init='Q'; }
    else if(STR_EQ(argv[i],"-i-rand"))              { mode_init='R';  }
    else if(STR_EQ(argv[i],"-i-vec") && i+1<argc)   { mode_init='X'; init_fname=char_renew(init_fname,argv[++i],NULL); }
    else if(STR_EQ(argv[i],"-i-val") && i+1<argc)   { mode_init='L'; init_fname=char_renew(init_fname,argv[++i],NULL); }
    else if(STR_EQ(argv[i],"-rhpeig"))              { strcpy(method,METHOD_RHPEIG); }
    else if(STR_EQ(argv[i],"-chpeig"))              { strcpy(method,METHOD_CHPEIG); }
    else if(STR_EQ(argv[i],"-chpeig-verify") && i+1<argc){ strcpy(method,METHOD_CHPEIG_VERIFY); prec_verify=atoi(argv[++i]); }
    else if(STR_EQ(argv[i],"-rhqr"))                { strcpy(method,METHOD_RHQR); }
    else if(STR_EQ(argv[i],"-chqr"))                { strcpy(method,METHOD_CHQR); }
    else if(STR_EQ(argv[i],"-rhqr-mt"))             { strcpy(method,METHOD_RHQR_MT); }
    else if(STR_EQ(argv[i],"-chqr-mt"))             { strcpy(method,METHOD_CHQR_MT); }
    else if(STR_EQ(argv[i],"-prec") && i+1<argc)    { prec=atoi(argv[++i]); }
    else if(STR_EQ(argv[i],"-k-yes"))               { kmode=1; }
    else if(STR_EQ(argv[i],"-k-no"))                { kmode=0; }
    else if(STR_EQ(argv[i],"-k-prec") && i+1<argc)  { kprec=atoi(argv[++i]); }
    else if(STR_EQ(argv[i],"-o-log") && i+1<argc)   { out_fname_log=char_renew(out_fname_log,argv[++i],NULL); }
    else if(STR_EQ(argv[i],"-o-val") && i+1<argc)   { out_fname_val=char_renew(out_fname_val,argv[++i],NULL); }
    else if(STR_EQ(argv[i],"-o-vec") && i+1<argc)   { out_fname_vec=char_renew(out_fname_vec,argv[++i],NULL); }
    else if(STR_EQ(argv[i],"-o-val-e") && i+1<argc) { out_fname_val_e=char_renew(out_fname_val_e,argv[++i],NULL); }
    else if(STR_EQ(argv[i],"-o-vec-e") && i+1<argc) { out_fname_vec_e=char_renew(out_fname_vec_e,argv[++i],NULL); }
    else if(STR_EQ(argv[i],"-O"))                   { autoname=1; }
    else if(STR_EQ(argv[i],"-no"))                  { autoname=0; }
    else if(STR_EQ(argv[i],"-d") && i+1<argc)       { out_dir=char_renew(out_dir,argv[++i],NULL); if(out_dir[strlen(out_dir)-1]!='/'){ out_dir=char_cat(out_dir,"/"); } }
    else if(STR_EQ(argv[i],"-f") && i+1<argc)       { strcpy(format,"f"); digits=atoi(argv[++i]); }
    else if(STR_EQ(argv[i],"-e") && i+1<argc)       { strcpy(format,"e"); digits=atoi(argv[++i]); }
    else if(STR_EQ(argv[i],"-nv"))                  { debug=0; }
    else if(STR_EQ(argv[i],"-v"))                   { debug=1; }
    else if(STR_EQ(argv[i],"-vv"))                  { debug=2; }
    else if(STR_EQ(argv[i],"-vvv"))                 { debug=3; }
    else if(STR_EQ(argv[i],"-vvvv"))                { debug=4; }
    else if(STR_EQ_N(argv[i],"-",1))                { printf("Error!\nIllegal option: %s\n\n",argv[i]);  usage(); }
    else if(strlen(in_fname)<=0)                    { in_fname=char_renew(in_fname,argv[i],NULL); }
    else                                            { printf("Error!\nIllegal option: %s\n\n",argv[i]); usage(); }
    i++;
  }
  if(strlen(in_fname)<=0){ printf("Error!\nNo input file.\n\n"); usage(); }

  // mode check
  if     (char_eq(method,METHOD_CHPEIG))       { mode='c'; }
  else if(char_eq(method,METHOD_CHPEIG_VERIFY)){ mode='c'; }
  else if(char_eq(method,METHOD_RHPEIG))       { mode='r'; }
  else if(char_eq(method,METHOD_RHQR))         { mode='r'; }
  else if(char_eq(method,METHOD_CHQR))         { mode='c'; }
  else if(char_eq(method,METHOD_RHQR_MT))      { mode='r'; }
  else if(char_eq(method,METHOD_CHQR_MT))      { mode='c'; }
  else{ ERROR_AT; exit(0); }

 // output
  if(kprec<0){ kprec=prec*KRAWCZYK_PREC_RATIO; }

  // load A
  cA=cmat_bin_load(&mA,&nA,"%s",in_fname);
  if(cA==NULL){ printf("Error! The input file '%s' can NOT be opened, or is NOT supported type.\n\n",in_fname); usage(); }
  if(mA!=nA){ printf("Error! The input matrix should be square. The size of the input matrix is %d-by-%d.\n\n",mA,nA); usage(); }
  n=nA;
  if(mode=='r'){
    rA=rmat_allocate_prec(mA,nA,prec);
    cmat_real_clone(mA,nA,rA,mA,cA,mA);
    cA=cmat_free(mA,nA,cA);    
  }

  // autoname
  if(autoname){
    path=strings_split_path(in_fname,BIN_SUFFIX);
    // get the directory name from path if it exists.
    if(strlen(out_dir)<=0 && strlen(path->str[0])>0){ out_dir=char_renew_sprintf(out_dir,NULL,"%s/",path->str[0]); }
    // set filename
    out_fname_log=char_renew_sprintf(out_fname_log,NULL,"%s_text{%dbits|%s|log}.%s",path->str[1],prec,method,TXT_SUFFIX);
    path=strings_del(path);
  }
  if(strlen(out_fname_log)>0){
    out_fname_log=char_renew_sprintf(out_fname_log,NULL,"%s%s",out_dir,out_fname_log);
    fid=fopen(out_fname_log,"w");
  }
  if(debug>=1){
    printf("Debug_Level: %d\n",debug);                                  if(fid!=NULL){ fprintf(fid,"Debug_Level: %d\n",debug); }
    printf("Precision: %d\n",prec);                                     if(fid!=NULL){ fprintf(fid,"Precision: %d\n",prec); }
    if(prec_verify>0){ printf("Precision_Verify: %d\n",prec_verify);    if(fid!=NULL){ fprintf(fid,"Precision_Verify: %d\n",prec_verify); } }
    printf("Method: %s\n",method);                                      if(fid!=NULL){ fprintf(fid,"Method: %s\n",method); }
    printf("Input_File: %s\n",in_fname);                                if(fid!=NULL){ fprintf(fid,"Input_File: %s\n",in_fname); }
    printf("Mode: %s\n",(mode=='r'?"real":(mode=='c'?"complex":"?")));  if(fid!=NULL){ fprintf(fid,"Mode: %s\n",(mode=='r'?"real":(mode=='c'?"complex":"?"))); }
    printf("Input_Matrix_Size: %d-by-%d\n",mA,nA);                      if(fid!=NULL){ fprintf(fid,"Input_Matrix_Size: %d -by- %d\n",mA,nA); }
  }
  if(debug>=3){
    printf("Input_Matrix_Entries:\n");
    if     (cA!=NULL){ cmat_print(mA,nA,cA,mA,NULL,format[0],digits); }
    else if(rA!=NULL){ rmat_print(mA,nA,rA,mA,NULL,format[0],digits); }
    else{ ERROR_AT; exit(0); }
  }
  if(debug>=1){
    printf("Output_File_Log:  %s\n",(strlen(out_fname_log)>0?out_fname_log:"NULL")); 
    if(fid!=NULL){ fprintf(fid,"Output_File_Log:  %s\n",(strlen(out_fname_log)>0?out_fname_log:"NULL")); }
  }

  // initail vectors
  if(char_eq(method,METHOD_RHQR) || char_eq(method,METHOD_CHQR) || char_eq(method,METHOD_RHQR_MT) || char_eq(method,METHOD_RHQR_MT)){ mode_init='Q'; }
  if(mode_init=='R'){
    if(debug>=1){
      printf("Initial_Vectors: random\n");
      if(fid!=NULL){ fprintf(fid,"Initial_Vectors: random\n"); }
    }
    mX=n; nX=n;
    cX=cmat_allocate_prec(mX,nX,prec); cmat_set_nan(n,n,cX,mX);
    cL=cvec_allocate_prec(n,prec);     cvec_set_nan(n,cL);
  }else if(mode_init=='X'){
    if(debug>=1){
      printf("Initial_Vectors: loaded from %s\n",(strlen(init_fname)>0?init_fname:"NULL"));
      if(fid!=NULL){ fprintf(fid,"Initial_Vectors: loaded from %s\n",(strlen(init_fname)>0?init_fname:"NULL")); }
    }
    if(strlen(init_fname)<=0){ mX=n; nX=n; cX=cmat_allocate_prec(mX,nX,prec); cmat_set_nan(n,n,cX,mX); }
    else{
      cX=cmat_bin_load(&mX,&nX,"%s",init_fname);
      cmat_round(mX,nX,cX,mX,prec);
      if(cX==NULL || mX!=n || nX!=n){
	printf("Error! The input file '%s' for initial vectors can NOT be opened, or is NOT supported type.\n\n",in_fname);
	usage();
      }
    }
    cL=cvec_allocate_prec(n,prec); cvec_set_nan(n,cL);
  }else if(mode_init=='L'){
    if(debug>=1){
      printf("Initial_Vectors: computed by II with %s\n",(strlen(init_fname)>0?init_fname:"NULL"));
      if(fid!=NULL){ fprintf(fid,"Initial_Vectors: computed by II with %s\n",(strlen(init_fname)>0?init_fname:"NULL")); }
    }
    mX=n; nX=n; cX=cmat_allocate_prec(mX,nX,prec); cmat_set_nan(n,n,cX,mX);
    cL=cvec_allocate_prec(n,prec); cvec_set_nan(n,cL);
  }else if(mode_init=='Q'){
    if(debug>=1){
      printf("Initial_Vectors: computed by QR\n");
      if(fid!=NULL){ fprintf(fid,"Initial_Vectors: computed by QR\n"); }
    }
    if(rA!=NULL){
      // computing eigenvalues by QR-method
      cL=cvec_allocate_prec(n,prec);
      clock_start=clock();
      if(char_eq(method,METHOD_RHQR))        { reig_hqr(n,cL,rA,mA,debug-1); }
      else if(char_eq(method,METHOD_RHQR_MT)){ rB=rmat_allocate_prec(mA,nA,prec); reig_hqr_mt(mA,nA,rB,mA,cL,rA,mA,debug-1); }
      clock_end=clock();
      time_qr=(double)(clock_end-clock_start)/CLOCKS_PER_SEC;
      time_total+=time_qr;
      // convert complex to real
      rL=rvec_allocate_prec(n,prec);
      rvec_real_cvec_clone(n,rL,cL);
      cL=cvec_free(n,cL);
      // computing eigenvectors by inverse iteration method
      mX=n; nX=n; rX=rmat_allocate_prec(mX,nX,prec);
      clock_start=clock();
      reig_ii(n,rX,mX,rA,mA,rL,debug-1);
      clock_end=clock();
      time_ii=(double)(clock_end-clock_start)/CLOCKS_PER_SEC;
      time_total+=time_ii;
      printf("Computing_Time_QR: %.15f sec\n",time_qr);
      printf("Computing_Time_II: %.15f sec\n",time_ii);
    }else if(cA!=NULL){
      // computing eigenvalues by QR-method
      cL=cvec_allocate_prec(n,prec);
      clock_start=clock();
      if(char_eq(method,METHOD_CHQR))        { ceig_hqr(n,cL,cA,mA,debug-1); }
      else if(char_eq(method,METHOD_CHQR_MT)){ cB=cmat_allocate_prec(mA,nA,prec); ceig_hqr_mt(mA,nA,cB,mA,cL,cA,mA,debug-1); }
      clock_end=clock();
      time_qr=(double)(clock_end-clock_start)/CLOCKS_PER_SEC;
      time_total+=time_qr;
      // computing eigenvectors by inverse iteration method
      mX=n; nX=n; cX=cmat_allocate_prec(mX,nX,prec);
      clock_start=clock();
      ceig_ii(n,cX,mX,cA,mA,cL,debug-1);
      clock_end=clock();
      time_ii=(double)(clock_end-clock_start)/CLOCKS_PER_SEC;
      time_total+=time_ii;
      printf("Computing_Time_QR: %.15f sec\n",time_qr);
      printf("Computing_Time_II: %.15f sec\n",time_ii);
    }else{ ERROR_AT; exit(0); }
  }else{ ERROR_AT; exit(0); }
  if(cX!=NULL && mode=='r'){
    rX=rmat_allocate_prec(mX,nX,prec);
    cmat_real_clone(mX,nX,rX,mX,cX,mX);
    cX=cmat_free(mX,nX,cX);
  }
  if(cL!=NULL && mode=='r'){
    rL=rvec_allocate_prec(n,prec);
    rvec_real_cvec_clone(n,rL,cL);
    cL=cvec_free(n,cL);
  }
  if(debug>=3){
    printf("Initail_Vectors_Entries:");
    if     (mode=='c' && cX!=NULL){ if(cmat_has_nan(mX,nX,cX,mX)){ printf(" NaN\n"); }else{ printf("\n"); cmat_print(mX,nX,cX,mX,NULL,format[0],digits); }}
    else if(mode=='r' && rX!=NULL){ if(rmat_has_nan(mX,nX,rX,mX)){ printf(" NaN\n"); }else{ printf("\n"); rmat_print(mX,nX,rX,mX,NULL,format[0],digits); }}
    else                          { printf("NULL\n"); }
  }

  // compute eigendecomposition
  clock_start=clock();
  if     (char_eq(method,METHOD_RHPEIG) && rA!=NULL && rX!=NULL && rL!=NULL){ k=rhpeig(n,rX,n,rL,rA,mA,debug-1); reig_sort(n,k,rL,rX,n); }
  else if(char_eq(method,METHOD_CHPEIG) && cA!=NULL && cX!=NULL && cL!=NULL){ k=chpeig(n,cX,n,cL,cA,mA,debug-1); ceig_sort(n,k,cL,cX,n); }
  else if(char_eq(method,METHOD_CHPEIG_VERIFY) && cA!=NULL && cX!=NULL && cL!=NULL){
    cXe=cmat_allocate_prec(mX,nX,prec_verify*2);
    cLe=cvec_allocate_prec(n,prec_verify*2);
    k=chpeig_verify(n,cX,n,cL,cXe,mX,cLe,cA,mA,prec_verify,&prec,&kprec,debug-1);
    ceig_sort(n,k,cL,cX,n);
  }
  else if(char_eq(method,METHOD_RHQR) && rA!=NULL && rX!=NULL && rL!=NULL){ k=n; }
  else if(char_eq(method,METHOD_CHQR) && cA!=NULL && cX!=NULL && cL!=NULL){ k=n; }
  else if(char_eq(method,METHOD_RHQR_MT) && rA!=NULL && rX!=NULL && rL!=NULL){ k=n; }
  else if(char_eq(method,METHOD_CHQR_MT) && cA!=NULL && cX!=NULL && cL!=NULL){ k=n; }
  else{ ERROR_AT; exit(0); }
  clock_end=clock();
  time_hpeig=(double)(clock_end-clock_start)/CLOCKS_PER_SEC;
  time_total+=time_hpeig;

  // autoname
  if(autoname){
    path=strings_split_path(in_fname,BIN_SUFFIX);
    // get the directory name from path if it exists.
    if(strlen(out_dir)<=0 && strlen(path->str[0])>0){ out_dir=char_renew_sprintf(out_dir,NULL,"%s/",path->str[0]); }
    // set filename
    out_fname_val=char_renew_sprintf(out_fname_val,NULL,"%s_%cvec{%dbits|%s|eigval}.%s",path->str[1],mode,prec,method,BIN_SUFFIX);
    out_fname_vec=char_renew_sprintf(out_fname_vec,NULL,"%s_%cmat{%dbits|%s|eigvec}.%s",path->str[1],mode,prec,method,BIN_SUFFIX);
    if(char_eq(method,METHOD_RHQR_MT) || char_eq(method,METHOD_CHQR_MT)){ out_fname_mat=char_renew_sprintf(out_fname_mat,NULL,"%s_%cmat{%dbits|%s|QRmat}.%s",path->str[1],mode,prec,method,BIN_SUFFIX); }
    path=strings_del(path);
  }

  if(debug>=1){
    printf("Output_File_Eigenvalues:  %s%s\n",out_dir,(strlen(out_fname_val)>0?out_fname_val:"NULL"));   
    if(fid!=NULL){ fprintf(fid,"Output_File_Eigenvalues:  %s%s\n",out_dir,(strlen(out_fname_val)>0?out_fname_val:"NULL")); }
    printf("Output_File_Eigenvectors: %s%s\n",out_dir,(strlen(out_fname_vec)>0?out_fname_vec:"NULL"));
    if(fid!=NULL){ fprintf(fid,"Output_File_Eigenvectors: %s%s\n",out_dir,(strlen(out_fname_vec)>0?out_fname_vec:"NULL")); }
    if(char_eq(method,METHOD_RHQR_MT) || char_eq(method,METHOD_CHQR_MT)){
      printf("Output_File_Matrix: %s%s\n",out_dir,(strlen(out_fname_mat)>0?out_fname_mat:"NULL"));
      if(fid!=NULL){ fprintf(fid,"Output_File_Matrix: %s%s\n",out_dir,(strlen(out_fname_mat)>0?out_fname_mat:"NULL")); }    
    }
  }

  // output results
  if(debug>=1){
    printf("Num_Obtained_Vectos: %d\n",k);
    if(fid!=NULL){ fprintf(fid,"Num_Obtained_Vectos: %d\n",k); }
    printf("Num_Not_Obtained_Vectos: %d\n",n-k);
    if(fid!=NULL){ fprintf(fid,"Num_Not_Obtained_Vectos: %d\n",n-k); }
    printf("Computing_Time_HPEIG: %.15f sec\n",time_hpeig);
    if(fid!=NULL){ fprintf(fid,"Computing_Time_HPEIG: %.15f sec\n",time_hpeig); }
    printf("Computing_Time_Total: %.15f sec\n",time_total);
    if(fid!=NULL){ fprintf(fid,"Computing_Time_Total: %.15f sec\n",time_total); }
  }
  if(debug>=3){
    printf("Eigenvectors:\n");
    if     (rX!=NULL){ rmat_print(n,k,rX,n,NULL,format[0],digits); }
    else if(cX!=NULL){ cmat_print(n,k,cX,n,NULL,format[0],digits); }
    else { ERROR_AT; exit(0); }
  }
  Lmax=rallocate_prec(prec);
  if(rL!=NULL){ rmax_abs_rvec(Lmax,k,rL); }
  if(cL!=NULL){ rmax_abs_cvec(Lmax,k,cL); }
  if(debug>=1){
    mpfr_printf("Eigenvalues_Max: %.7Re\n",Lmax);
    if(fid!=NULL){ mpfr_fprintf(fid,"Eigenvalues_Max: %.7Re\n",Lmax); }
  }
  if(debug>=2){
    printf("Eigenvalues:\n");
    if     (rL!=NULL){ rvec_print(k,rL,NULL,format[0],digits); }
    else if(cL!=NULL){ cvec_print(k,cL,NULL,format[0],digits); }
    else { ERROR_AT; exit(0); }
  }
  if(debug>=3 && char_eq(method,METHOD_RHQR_MT)){
    printf("Matrix after QR:\n");
    if     (rB!=NULL){ rmat_print(n,k,rB,n,NULL,format[0],digits); }
    else { ERROR_AT; exit(0); }
  }
  if(debug>=3 && char_eq(method,METHOD_CHQR_MT)){
    printf("Matrix after QR:\n");
    if(cB!=NULL){ cmat_print(n,k,cB,n,NULL,format[0],digits); }
    else { ERROR_AT; exit(0); }
  }
  if(strlen(out_fname_val)>0){
    if     (rL!=NULL){ rvec_bin_save(n,rL,"%s%s",out_dir,out_fname_val); }
    else if(cL!=NULL){ cvec_bin_save(n,cL,"%s%s",out_dir,out_fname_val); }
    else { ERROR_AT; exit(0); }
    if(debug>=1 && fid!=NULL){ fprintf(stderr,"saved: %s%s\n",out_dir,out_fname_val); }
  }
  if(strlen(out_fname_vec)>0){
    if     (rX!=NULL){ rmat_bin_save(n,n,rX,n,"%s%s",out_dir,out_fname_vec); }
    else if(cX!=NULL){ cmat_bin_save(n,n,cX,n,"%s%s",out_dir,out_fname_vec); }
    else { ERROR_AT; exit(0); }
    if(debug>=1 && fid!=NULL){ fprintf(stderr,"saved: %s%s\n",out_dir,out_fname_vec); }
  }
  if(strlen(out_fname_mat)>0){
    if     (rB!=NULL){ rmat_bin_save(n,n,rB,n,"%s%s",out_dir,out_fname_mat); }
    else if(cB!=NULL){ cmat_bin_save(n,n,cB,n,"%s%s",out_dir,out_fname_mat); }
    else { ERROR_AT; exit(0); }
    if(debug>=1 && fid!=NULL){ fprintf(stderr,"saved: %s%s\n",out_dir,out_fname_mat); }
  }

  // error bounds
  if(debug>=1){
    printf("Error_Bounds_Mode: %s\n",kmode?"Yes":"No");
    if(fid!=NULL){ fprintf(fid,"Error_Bounds_Mode: %s\n",kmode?"Yes":"No"); }
  }
  if(kmode){
    if(debug>=1){
      printf("Error_Bounds_Prec: %d\n",kprec);
      if(fid!=NULL){ fprintf(fid,"Error_Bounds_Prec: %d\n",kprec); }
    }
    rXe_max=rallocate_prec(kprec);
    rLe_max=rallocate_prec(kprec);
    rXe_max_log2=rallocate_prec(kprec);
    rLe_max_log2=rallocate_prec(kprec);
    if(rA!=NULL && rX!=NULL && rL!=NULL){
      rXe=rmat_allocate_prec(mX,nX,kprec);
      rLe=rvec_allocate_prec(n,kprec);
      // computing
      clock_start=clock();
      kret=ireig_krawczyk(n,k,rXe,mX,rLe,rA,mA,rX,mX,rL,debug-1);
      clock_end=clock();
      time_k=(double)(clock_end-clock_start)/CLOCKS_PER_SEC;
    }else if(cA!=NULL && cX!=NULL && cL!=NULL){
      cXe=cmat_allocate_prec(mX,nX,kprec);
      cLe=cvec_allocate_prec(n,kprec);
      // computing
      clock_start=clock();
      kret=iceig_krawczyk(n,k,cXe,mX,cLe,cA,mA,cX,mX,cL,debug-1);
      clock_end=clock();
      time_k=(double)(clock_end-clock_start)/CLOCKS_PER_SEC;
    }else{ ERROR_AT; exit(0); }

    // output
    if(debug>=1){
      printf("Error_Bounds_Result: %s\n",kret?"No":"Yes");
      if(fid!=NULL){ fprintf(fid,"Error_Bounds_Result: %s\n",kret?"No":"Yes"); }
      printf("Computing_Time_Error_Bounds: %.15f sec\n",time_k);
      if(fid!=NULL){ fprintf(fid,"Computing_Time_Error_Bounds: %.15f sec\n",time_k); }
    }

    // autoname
    if(kmode && autoname){
      path=strings_split_path(in_fname,BIN_SUFFIX);
      // get the directory name from path if it exists.
      if(strlen(out_dir)<=0 && strlen(path->str[0])>0){ out_dir=char_renew_sprintf(out_dir,NULL,"%s/",path->str[0]); }
      out_fname_val_e=char_renew_sprintf(out_fname_val_e,NULL,"%s_%cvec{%dbits|%s|eigval}_%cvec{krawczyk|%dbits}.%s",path->str[1],mode,prec,method,mode,kprec,BIN_SUFFIX);
      out_fname_vec_e=char_renew_sprintf(out_fname_vec_e,NULL,"%s_%cmat{%dbits|%s|eigvec}_%cmat{krawczyk|%dbits}.%s",path->str[1],mode,prec,method,mode,kprec,BIN_SUFFIX);
      path=strings_del(path);
    }

    if(kmode && debug>=1){
      printf("Output_File_Error_Bounds_Eigenvalues:  %s%s\n",out_dir,(strlen(out_fname_val_e)>0?out_fname_val_e:"NULL")); 
      if(fid!=NULL){ fprintf(fid,"Output_File_Error_Bounds_Eigenvalues:  %s%s\n",out_dir,(strlen(out_fname_val_e)>0?out_fname_val_e:"NULL")); }
      printf("Output_File_Error_Bounds_Eigenvectors: %s%s\n",out_dir,(strlen(out_fname_vec_e)>0?out_fname_vec_e:"NULL")); 
      if(fid!=NULL){ fprintf(fid,"Output_File_Error_Bounds_Eigenvectors: %s%s\n",out_dir,(strlen(out_fname_vec_e)>0?out_fname_vec_e:"NULL")); }
    }

    if(kret){
      // failed
      out_fname_vec_e=char_renew(out_fname_vec_e,"",NULL);
      out_fname_val_e=char_renew(out_fname_val_e,"",NULL);
      if(rX!=NULL){ rmat_set_nan(mX,nX,rX,mX); } else if(cX!=NULL){ cmat_set_nan(mX,nX,cX,mX); } else{ ERROR_AT; exit(0); }
      if(rL!=NULL){ rvec_set_nan(n,rLe); }       else if(cL!=NULL){ cvec_set_nan(n,cLe); }       else{ ERROR_AT; exit(0); }
    }else{
      // seccess
      if(rXe!=NULL){ rmat_max_abs(rXe_max,n,k,rXe,n); }
      if(cXe!=NULL){ cmat_max_absc(rXe_max,n,k,cXe,n); }
      if(rLe!=NULL){ rmax_abs_rvec(rLe_max,k,rLe); }
      if(cLe!=NULL){ rmax_abs_cvec(rLe_max,k,cLe); }
      rdiv_rr(rLe_max,rLe_max,Lmax); rlog2_r(rLe_max_log2,rLe_max); rlog2_r(rXe_max_log2,rXe_max);
      if(debug>=1){
	mpfr_printf("Error_Bounds_Eigenvectors_Max: %.1Re 2^%d\n",rXe_max,(int)ceil(rget_d(rXe_max_log2)));
	if(fid!=NULL){ mpfr_fprintf(fid,"Error_Bounds_Eigenvectors_Max: %.1Re 2^%d\n",rXe_max,(int)ceil(rget_d(rXe_max_log2))); }
	mpfr_printf("Error_Bounds_Eigenvalues_Max_Per_Max: %.1Re 2^%d\n",rLe_max,(int)ceil(rget_d(rLe_max_log2)));
	if(fid!=NULL){ mpfr_fprintf(fid,"Error_Bounds_Eigenvalues_Max_Per_Max: %.1Re 2^%d\n",rLe_max,(int)ceil(rget_d(rLe_max_log2))); }
      }
      if(debug>=3){
	printf("Error_Bounds_Eigenvectors:\n");
	if     (rXe!=NULL){ rmat_print(n,k,rXe,n,NULL,'e',0); }
	else if(cXe!=NULL){ cmat_print(n,k,cXe,n,NULL,'e',0); }
	else{ ERROR_AT; exit(0); }
	printf("Error_Bounds_Eigenvalues:\n");
	if     (rLe!=NULL){ rvec_print(k,rLe,NULL,'e',0); }
	else if(cLe!=NULL){ cvec_print(k,cLe,NULL,'e',0); }
	else{ ERROR_AT; exit(0); }
      }
      if(strlen(out_fname_val_e)>0){
	if     (rLe!=NULL){ rvec_bin_save(n,rLe,"%s%s",out_dir,out_fname_val_e); }
	else if(cLe!=NULL){ cvec_bin_save(n,cLe,"%s%s",out_dir,out_fname_val_e); }
	else{ ERROR_AT; exit(0); }
	if(debug>=1 && fid!=NULL){ fprintf(stderr,"saved: %s%s\n",out_dir,out_fname_val_e); }
      }
      if(strlen(out_fname_vec_e)>0){
	if     (rXe!=NULL){ rmat_bin_save(n,n,rXe,n,"%s%s",out_dir,out_fname_vec_e); }
	else if(cXe!=NULL){ cmat_bin_save(n,n,cXe,n,"%s%s",out_dir,out_fname_vec_e); }
	else{ ERROR_AT; exit(0); }
	if(debug>=1 && fid!=NULL){ fprintf(stderr,"saved: %s%s\n",out_dir,out_fname_vec_e); }
      } 
    }
  }

  // done
  if(fid!=NULL){ fclose(fid); fid=NULL; }
  cA=cmat_free(mA,nA,cA);
  cB=cmat_free(mA,nA,cB);
  rA=rmat_free(mA,nA,rA);
  rB=rmat_free(mA,nA,rB);
  cX=cmat_free(nX,nX,cX);
  rX=rmat_free(nX,nX,rX);
  cXe=cmat_free(nX,nX,cXe);
  rXe=rmat_free(nX,nX,rXe);
  cL=cvec_free(n,cL);
  rL=rvec_free(n,rL);
  cLe=cvec_free(n,cLe);
  rLe=rvec_free(n,rLe);
  rXe_max=rmfree(rXe_max);
  rLe_max=rmfree(rLe_max);
  rXe_max_log2=rmfree(rXe_max_log2);
  rLe_max_log2=rmfree(rLe_max_log2);
  Lmax=rmfree(Lmax);
  in_fname=char_del(in_fname);
  out_fname_val=char_del(out_fname_val);
  out_fname_vec=char_del(out_fname_vec);
  out_dir=char_del(out_dir);
  return 0;
}
