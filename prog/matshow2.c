#include<isys.h>

#define PROG "matshow2"
#define MAX_LENGTH 10000
#define BIN_SUFFIX "dat"
#define TXT_SUFFIX "txt"
#define METHOD_NAME_LENGTH  128
#define METHOD_RSUB       "rsub_coeff"
#define METHOD_CSUB       "csub_coeff"
#define DEFAULT_PREC     4096

#define RA(X,P)    ((X)=rallocate_prec(P))
#define RAP(X,Y)   ((X)=rallocate_prec(rget_prec(Y)))
#define RF(X)      ((X)=rmfree(X))
#define CA(X,P)    ((X)=callocate_prec(P))
#define CF(X)      ((X)=cmfree(X))
#define RVA(X,N,P)  { X=rvec_allocate_prec(N,P); }
#define RVAr(X,Y,N) { X=rvec_allocate_prec(N,rget_prec(Y)); }
#define RVArv(X,Y,N){ X=rvec_allocate_prec(N,rvec_get_prec_max(N,Y)); }
#define RVF(X,N)    { X=rvec_free(N,X); }
#define CVA(X,N,P)  { X=cvec_allocate_prec(N,P); }
#define CVF(X,N)    { X=cvec_free(N,X); }

enum { NOTYPE=-1, TEXT=0, DVEC, ZVEC, RVEC, CVEC, DMAT, ZMAT, RMAT, CMAT };
char *type_name[]={ "text", "dvec", "zvec", "rvec", "cvec", "dmat", "zmat", "rmat", "cmat" };

int rmat_max_abs_coeff(int *v_m, int *v_n, rmulti *value, int m, int n, rmulti **A, int LDA)
{
  int i,j;
  rmulti *a=NULL;
  a=rallocate_prec(rget_prec(value));
  rabs_r(value,MAT(A,0,0,LDA));  // value=abs(x[0])
  *v_m=0; *v_n=0;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      rabs_r(a,MAT(A,i,j,LDA));  // a=abs(x[i])
      if(gt_rr(a,value)){            // a>value
	rset_r(value,a);         // value=a
	*v_m=i; *v_n=j;
      }
    }
  }
  a=rmfree(a);
  return 0;
}

int cmat_max_abs_coeff(int *v_m, int *v_n, rmulti *value, int m, int n, cmulti **A, int LDA)
{
  int i,j;
  rmulti *a=NULL;
  a=rallocate_prec(rget_prec(value));
  rabs2_c(value,MAT(A,0,0,LDA));  // value=abs(x[0])^2
  *v_m=0; *v_n=0;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      rabs2_c(a,MAT(A,i,j,LDA));  // a=abs(x[i])^2
      if(gt_rr(a,value)){            // a>value
	rset_r(value,a);         // value=a
	*v_m=i; *v_n=j;
      }
    }
  }
  rsqrt_r(value,value); // value=sqrt(value)
  a=rmfree(a);            // free
  return 0;
}

void usage()
{
  printf("Usage: %s input.dat [options]\n",PROG);
  printf("  -input1            Input first file.\n");
  printf("  -input2            Input second file.\n");
  printf("  -rsub-coeff        Do subtraction on real matrix coefficient.\n");
  printf("  -csub-coeff        Do subtraction on complex matrix coefficient.\n");
  printf("  -o-mt fname        The matrix is saved in the file 'fname'.\n");
  printf("  -o-rlog fname      The matrix is saved in the file 'fname'.\n");
  printf("  -O                 The matrix is saved in the file whose name is automatically determined.\n");
  printf("  -d dir             Set the directory of output file as 'dir'.\n");
  printf("  -f num             Set print format as 'f' with num digits.\n");
  printf("  -e num             Set print format as 'e' with num digits.\n");
  printf("  -prec num          Set the precision used in output mode.\n");
  printf("  -nv                Set no verbose mode.\n");
  printf("  -v                 Set verbose mode, level 1.\n");
  printf("  -vv                Set verbose mode, level 2.\n");
  printf("  -vvv               Set verbose mode, leval 3.\n");
  printf("  -help              This message.\n");
  printf("\n");
  printf("Examples:\n");
  printf("     # %s -input1 _____.dat -input2 _____.dat \n",PROG);
  printf("\n");
  exit(0);
}


//まだmatsub_coeffのみ
int main(int argc, char *argv[])
{
  FILE *fid=NULL;
  strings *path=NULL;
  char *in_fname1=NULL,*in_fname2=NULL;
  char *out_fname_log=NULL,*out_fname_mt=NULL,*out_dir=NULL;
  char format[]="f",method[METHOD_NAME_LENGTH+1],buf[MAX_LENGTH],mode='?';
  int m1=0,m2=0,max_m=0,Emax_m=0,n1=0,n2=0,max_n=0,Emax_n=0,debug=1,i,prec=DEFAULT_PREC,flag=0,pflag=0;
  int loaded=0,autoname=0,type_in1=NOTYPE,type_in2=NOTYPE,digits=3,*index=NULL;
  double *dx1=NULL,*dx2=NULL,*dA1=NULL,*dA2=NULL;
  dcomplex *zx1=NULL,*zx2=NULL,*zA1=NULL,*zA2=NULL;
  dcomplex *x1_z=NULL;
  rmulti **rx1=NULL,**rx2=NULL,**rA1=NULL,**rA2=NULL,**rE=NULL;
  cmulti **cx1=NULL,**cx2=NULL,**cA1=NULL,**cA2=NULL,**cE=NULL;
  rmulti *rEmax=NULL,*rmax=NULL;

  // init
  in_fname1=char_new("",NULL);
  in_fname2=char_new("",NULL);
  out_fname_log=char_new("",NULL);
  out_fname_mt=char_new("",NULL);
  out_dir=char_new("",NULL);
  strcpy(method,METHOD_RSUB);

  // options
  i=1;
  while(i<argc){
    if     (STR_EQ(argv[i],"-help"))                                 { usage(); }
    else if(STR_EQ(argv[i],"-input1") && i+1<argc)                   { in_fname1=char_renew(in_fname1,argv[++i],NULL); }
    else if(STR_EQ(argv[i],"-input2") && i+1<argc)                   { in_fname2=char_renew(in_fname2,argv[++i],NULL); }
    else if(STR_EQ(argv[i],"-rsub-coeff"))                           { flag=1; strcpy(method,METHOD_RSUB); }
    else if(STR_EQ(argv[i],"-csub-coeff"))                           { flag=1; strcpy(method,METHOD_CSUB); }
    else if(STR_EQ(argv[i],"-o-mt") && i+1<argc)                     { ++i; out_fname_mt=char_renew(out_fname_mt,argv[i],NULL); }
    else if(STR_EQ(argv[i],"-o-log") && i+1<argc)                    { ++i; out_fname_log=char_renew(out_fname_log,argv[i],NULL); }
    else if(STR_EQ(argv[i],"-O"))                                    { autoname=1; }
    else if(STR_EQ(argv[i],"-d") && i+1<argc)                        { ++i; out_dir=char_renew(out_dir,argv[i],NULL); if(out_dir[strlen(out_dir)-1]!='/'){ out_dir=char_cat(out_dir,"/"); } }
    else if(STR_EQ(argv[i],"-f") && i+1<argc)                        { strcpy(format,"f"); ++i; digits=atoi(argv[i]); }
    else if(STR_EQ(argv[i],"-e") && i+1<argc)                        { strcpy(format,"e"); ++i; digits=atoi(argv[i]); }
    else if(STR_EQ(argv[i],"-prec") && i+1<argc)                     { strcpy(format,"e"); ++i; prec=atoi(argv[i]); pflag=1; }
    else if(STR_EQ(argv[i],"-nv"))                                   { debug=0; }
    else if(STR_EQ(argv[i],"-v"))                                    { debug=1; }
    else if(STR_EQ(argv[i],"-vv"))                                   { debug=2; }
    else if(STR_EQ(argv[i],"-vvv"))                                  { debug=3; }
    else if(STR_EQ_N(argv[i],"-",1))                                 { usage(); }
    else                                                             { usage(); }
    i++;
  }

  // mode check
  if     (char_eq(method,METHOD_RSUB)){ mode='r'; }
  else if(char_eq(method,METHOD_CSUB)){ mode='c'; }
  else{ ERROR_AT; exit(0); }

  // load input file1
  if(strlen(in_fname1)<=0){ usage(); }
  //if(!loaded){ dx1=dvec_bin_load(&n1,"%s",in_fname1); }     if(dx1!=NULL){ type_in1=DVEC; loaded=1; }
  //if(!loaded){ zx1=zvec_bin_load(&n1,"%s",in_fname1); }     if(zx1!=NULL){ type_in1=ZVEC; loaded=1; }
  if(!loaded && mode=='r'){ rx1=rvec_bin_load(&n1,"%s",in_fname1); }     if(rx1!=NULL){ type_in1=RVEC; loaded=1; }
  if(!loaded && mode=='c'){ cx1=cvec_bin_load(&n1,"%s",in_fname1); }     if(cx1!=NULL){ type_in1=CVEC; loaded=1; }
  //if(!loaded){ dA1=dmat_bin_load(&m1,&n1,"%s",in_fname1); } if(dA1!=NULL){ type_in1=DMAT; loaded=1; }
  //if(!loaded){ zA1=zmat_bin_load(&m1,&n1,"%s",in_fname1); } if(zA1!=NULL){ type_in1=ZMAT; loaded=1; }
  if(!loaded && mode=='r'){ rA1=rmat_bin_load(&m1,&n1,"%s",in_fname1); } if(rA1!=NULL){ type_in1=RMAT; loaded=1; }
  if(!loaded && mode=='c'){ cA1=cmat_bin_load(&m1,&n1,"%s",in_fname1); } if(cA1!=NULL){ type_in1=CMAT; loaded=1; }
  if(!loaded){ printf("Error! The input file '%s%s' can NOT be opened, or is NOT supported type.\n\n",out_dir,in_fname1); usage(); }
  if(debug>=2){
    printf("Input1_File: %s\n",in_fname1);
    printf("Input1_Type: %s\n",type_name[type_in1]);
    if(DVEC<=type_in1 && type_in1<=CVEC){ printf("Input1_Size: %d\n",n1); }
    if(DMAT<=type_in1 && type_in1<=CMAT){ printf("Input1_Size: %d-by-%d\n",m1,n1); }
    if(type_in1==DVEC || type_in1==ZVEC || type_in1==DMAT || type_in1==ZMAT){ printf("Input1_Precision: double\n"); }
    else if(type_in1==RVEC){ printf("Input1_Precision_Max: %d\n",rvec_get_prec_max(n1,rx1)); }
    else if(type_in1==CVEC){ printf("Input1_Precision_Max: %d\n",cvec_get_prec_max(n1,cx1)); }
    else if(type_in1==RMAT){ printf("Input1_Precision_Max: %d\n",rmat_get_prec_max(m1,n1,rA1,m1)); }
    else if(type_in1==CMAT){ printf("Input1_Precision_Max: %d\n",cmat_get_prec_max(m1,n1,cA1,m1)); }
  }
  if(debug>=3){
    if(DVEC<=type_in1 && type_in1<=CVEC){ printf("Input1_Vector:\n"); }
    if(DMAT<=type_in1 && type_in1<=CMAT){ printf("Input1_Matrix:\n"); }
    if(dx1!=NULL){ dvec_print(n1,dx1,NULL,format[0],digits); }
    if(zx1!=NULL){ zvec_print(n1,zx1,NULL,format[0],digits); }
    if(rx1!=NULL){ rvec_print(n1,rx1,NULL,format[0],digits); }
    if(cx1!=NULL){ cvec_print(n1,cx1,NULL,format[0],digits); }
    if(dA1!=NULL){ dmat_print(m1,n1,dA1,m1,NULL,format[0],digits); }
    if(zA1!=NULL){ zmat_print(m1,n1,zA1,m1,NULL,format[0],digits); }
    if(rA1!=NULL){ rmat_print(m1,n1,rA1,m1,NULL,format[0],digits); }
    if(cA1!=NULL){ cmat_print(m1,n1,cA1,m1,NULL,format[0],digits); }
  }

  // load input file2
  loaded=0;
  if(strlen(in_fname2)<=0){ usage(); }
  //if(!loaded){ dx2=dvec_bin_load(&n2,"%s",in_fname2); } if(dx2!=NULL){ type_in2=DVEC; loaded=1; }
  //if(!loaded){ zx2=zvec_bin_load(&n2,"%s",in_fname2); } if(zx2!=NULL){ type_in2=ZVEC; loaded=1; }
  if(!loaded && mode=='r'){ rx2=rvec_bin_load(&n2,"%s",in_fname2); } if(rx2!=NULL){ type_in2=RVEC; loaded=1; }
  if(!loaded && mode=='c'){ cx2=cvec_bin_load(&n2,"%s",in_fname2); } if(cx2!=NULL){ type_in2=CVEC; loaded=1; }
  //if(!loaded){ dA2=dmat_bin_load(&m2,&n2,"%s",in_fname2); } if(dA2!=NULL){ type_in2=DMAT; loaded=1; }
  //if(!loaded){ zA2=zmat_bin_load(&m2,&n2,"%s",in_fname2); } if(zA2!=NULL){ type_in2=ZMAT; loaded=1; }
  if(!loaded && mode=='r'){ rA2=rmat_bin_load(&m2,&n2,"%s",in_fname2); } if(rA2!=NULL){ type_in2=RMAT; loaded=1; }
  if(!loaded && mode=='c'){ cA2=cmat_bin_load(&m2,&n2,"%s",in_fname2); } if(cA2!=NULL){ type_in2=CMAT; loaded=1; }
  if(!loaded){ printf("Error! The input file '%s%s' can NOT be opened, or is NOT supported type.\n\n",out_dir,in_fname2); usage(); }
  if(debug>=2){
    printf("Input2_File: %s\n",in_fname2);
    printf("Input2_Type: %s\n",type_name[type_in2]);
    if(DVEC<=type_in2 && type_in2<=CVEC){ printf("Input2_Size: %d\n",n2); }
    if(DMAT<=type_in2 && type_in2<=CMAT){ printf("Input2_Size: %d-by-%d\n",m2,n2); }
    if(type_in2==DVEC || type_in2==ZVEC || type_in2==DMAT || type_in2==ZMAT){ printf("Input2_Precision: double\n"); }
    else if(type_in2==RVEC){ printf("Input2_Precision_Max: %d\n",rvec_get_prec_max(n2,rx2)); }
    else if(type_in2==CVEC){ printf("Input2_Precision_Max: %d\n",cvec_get_prec_max(n2,cx2)); }
    else if(type_in2==RMAT){ printf("Input2_Precision_Max: %d\n",rmat_get_prec_max(m2,n2,rA2,m2)); }
    else if(type_in2==CMAT){ printf("Input2_Precision_Max: %d\n",cmat_get_prec_max(m2,n2,cA2,m2)); }
  }
  if(debug>=3){
    if(DVEC<=type_in2 && type_in2<=CVEC){ printf("Input2_Vector:\n"); }
    if(DMAT<=type_in2 && type_in2<=CMAT){ printf("Input2_Matrix:\n"); }
    if(dx2!=NULL){ dvec_print(n2,dx2,NULL,format[0],digits); }
    if(zx2!=NULL){ zvec_print(n2,zx2,NULL,format[0],digits); }
    if(rx2!=NULL){ rvec_print(n2,rx2,NULL,format[0],digits); }
    if(cx2!=NULL){ cvec_print(n2,cx2,NULL,format[0],digits); }
    if(dA2!=NULL){ dmat_print(m2,n2,dA2,m2,NULL,format[0],digits); }
    if(zA2!=NULL){ zmat_print(m2,n2,zA2,m2,NULL,format[0],digits); }
    if(rA2!=NULL){ rmat_print(m2,n2,rA2,m2,NULL,format[0],digits); }
    if(cA2!=NULL){ cmat_print(m2,n2,cA2,m2,NULL,format[0],digits); }
  }

  // check
  if(n1!=n2){ printf("Error! The input file's size is NOT same!\n\n"); exit(0); }
  if(m1!=n1){ printf("Error! The input file is NOT square!\n\n"); exit(0); }
  if(m2!=n2){ printf("Error! The input file is NOT square!\n\n"); exit(0); }

  // get prec
  if(pflag==0 && DVEC<=type_in1 && type_in1<=CVEC && DVEC<=type_in2 && type_in2<=CVEC){
    if(mode=='r'){ prec=MAX2(rvec_get_prec_max(n1,rx1),rvec_get_prec_max(n2,rx2)); }
    if(mode=='c'){ prec=MAX2(cvec_get_prec_max(n1,cx1),cvec_get_prec_max(n2,cx2)); }
  }
  if(pflag==0 && DMAT<=type_in1 && type_in1<=CMAT && DMAT<=type_in2 && type_in2<=CMAT){
    if(mode=='r'){ prec=MAX2(rmat_get_prec_max(m1,n1,rA1,m1),rmat_get_prec_max(m2,n2,rA2,m2)); }
    if(mode=='c'){ prec=MAX2(cmat_get_prec_max(m1,n1,cA1,m1),cmat_get_prec_max(m2,n2,cA2,m2)); }
  }

  // autoname
  if(autoname){
    path=strings_split_path(in_fname2,BIN_SUFFIX);
    if(strlen(out_dir)<=0 && strlen(path->str[0])>0){ out_dir=char_renew(out_dir,path->str[0],NULL); out_dir=char_cat(out_dir,"/"); }
    if(flag==1){
      sprintf(buf,"%s_text{%dbits|%s|log}.%s",path->str[1],prec,method,TXT_SUFFIX); out_fname_log=char_cat(out_fname_log,buf);
      sprintf(buf,"%s_%cmat{%dbits|%s}.%s",path->str[1],mode,prec,method,BIN_SUFFIX); out_fname_mt=char_cat(out_fname_mt,buf);
    }
    path=strings_del(path);
  }
  if(strlen(out_fname_log)>0){
    strcpy(buf,out_dir); strcat(buf,"/"); strcat(buf,out_fname_log);
    fid=fopen(buf,"w");
  }

  // log
  if(debug>=1){
    printf("Debug_Level: %d\n",debug);     if(fid!=NULL){ fprintf(fid,"Debug_Level: %d\n",debug);     }
    printf("Precision: %d\n",prec);        if(fid!=NULL){ fprintf(fid,"Precision: %d\n",prec);        }
    printf("Method: %s\n",method);         if(fid!=NULL){ fprintf(fid,"Method: %s\n",method);         }
    printf("Input_Vector_Size: %d\n",n1);  if(fid!=NULL){ fprintf(fid,"Input_Vector_Size: %d\n",n1);  }
    printf("Input1_File: %s\n",in_fname1); if(fid!=NULL){ fprintf(fid,"Input_File1: %s\n",in_fname1); }
    printf("Input2_File: %s\n",in_fname2); if(fid!=NULL){ fprintf(fid,"Input_File2: %s\n",in_fname2); }
    printf("Mode: %s\n",(mode=='r'?"real":(mode=='c'?"complex":"?")));  
    if(fid!=NULL){ fprintf(fid,"Mode: %s\n",(mode=='r'?"real":(mode=='c'?"complex":"?"))); }
    printf("Output_File_Log:  %s%s\n",out_dir,(strlen(out_fname_log)>0?out_fname_log:"NULL"));
    if(fid!=NULL){ fprintf(fid,"Output_File_Log:  %s%s\n",out_dir,(strlen(out_fname_log)>0?out_fname_log:"NULL")); }
    printf("Output_File_mt:  %s%s\n",out_dir,(strlen(out_fname_mt)>0?out_fname_mt:"NULL"));
    if(fid!=NULL){ fprintf(fid,"Output_File_m:  %s%s\n",out_dir,(strlen(out_fname_mt)>0?out_fname_mt:"NULL")); }
  }

  // vec sorting
  if(type_in1==RVEC && type_in2==RVEC){
    index=ivec_allocate(n1);
    x1_z=zvec_allocate(n1);
    rvec_get_zvec(n1,x1_z,rx1);
    zvec_sort(n1,x1_z,index); zvec_reverse(n1,x1_z); ivec_reverse(n1,index);
    rvec_swap_index(n1,rx1,index);
    rvec_swap_index(n2,rx2,index);
    ivec_set_grid(n1,index);
    reig_sort_value_guide(n1,n1,rx2,NULL,0,rx1);
  }
  if(type_in1==CVEC && type_in2==CVEC){
    index=ivec_allocate(n1);
    x1_z=zvec_allocate(n1);
    cvec_get_zvec(n1,x1_z,cx1);
    zvec_sort(n1,x1_z,index); zvec_reverse(n1,x1_z); ivec_reverse(n1,index);
    cvec_swap_index(n1,cx1,index);
    cvec_swap_index(n2,cx2,index);
    ivec_set_grid(n1,index);
    ceig_sort_value_guide(n1,n1,cx2,NULL,0,cx1);
  }

  //compute mode [sub_coeff]
  if(flag==1){
    RA(rEmax,prec);
    RA(rmax,prec);
    if(type_in1==RMAT && type_in2==RMAT){
      rE=rmat_allocate_prec(m1,n1,prec);
      rmat_sub(m1,n1,rE,m1,rA1,m1,rA2,m2);
      rmat_max_abs_coeff(&Emax_m,&Emax_n,rEmax,m1,n1,rE,m1);
      if(rmat_get_prec_max(m1,n1,rA1,m1)>rmat_get_prec_max(m2,n2,rA2,m2)){ 
	rmat_max_abs_coeff(&max_m,&max_n,rmax,m1,n1,rA1,m1);
      }else{
	rmat_max_abs_coeff(&max_m,&max_n,rmax,m2,n2,rA2,m2);
      }
      rdiv_rr(rEmax,rEmax,rmax);
    }
    if(type_in1==CMAT && type_in2==CMAT){
      cE=cmat_allocate_prec(m1,n1,prec);
      cmat_sub(m1,n1,cE,m1,cA1,m1,cA2,m2);
      cmat_max_abs_coeff(&Emax_m,&Emax_n,rEmax,m1,n1,cE,m1);
      if(cmat_get_prec_max(m1,n1,cA1,m1)>cmat_get_prec_max(m2,n2,cA2,m2)){ 
	cmat_max_abs_coeff(&max_m,&max_n,rmax,m1,n1,cA1,m1); 
      }else{
	cmat_max_abs_coeff(&max_m,&max_n,rmax,m2,n2,cA2,m2);
      }
      rdiv_rr(rEmax,rEmax,rmax);
    }
  }

  //output results
  if(debug>=2){
    if(rE!=NULL){ printf("Output_coeff_Error=\n"); rmat_print(m1,n1,rE,m2,NULL,format[0],digits); }
    if(cE!=NULL){ printf("Output_coeff_Error=\n"); cmat_print(m1,n1,cE,m2,NULL,format[0],digits); }
  }
  if(debug>=1){
    if(rEmax!=NULL){ mpfr_printf("Output_coeff_Error_max: [%d,%d] = %.7Re\n",Emax_m,Emax_n,rEmax); }
  }

  //save
  if(strlen(out_fname_mt)>0){
    if(flag==1 && rE!=NULL){ rmat_bin_save(m1,n1,rE,m1,"%s%s",out_dir,out_fname_mt); }
    if(flag==1 && cE!=NULL){ cmat_bin_save(m1,n1,cE,m1,"%s%s",out_dir,out_fname_mt); }
    if(debug>=1){ fprintf(stderr,"saved: %s%s\n",out_dir,out_fname_mt); }
  }
  if(strlen(out_fname_log)>0){
    if(flag==1 && rEmax!=NULL && fid!=NULL){ mpfr_fprintf(fid,"Output_coeff_Error_max: [%d,%d] = %.7Re\n",Emax_m,Emax_n,rEmax); }
    if(debug>=1){ fprintf(stderr,"saved: %s%s\n",out_dir,out_fname_log); }
  }

  // done
  if(fid!=NULL){ fclose(fid); fid=NULL; }
  in_fname1=char_del(in_fname1);
  in_fname2=char_del(in_fname2);
  out_fname_log=char_del(out_fname_log);
  out_fname_mt=char_del(out_fname_mt);
  out_dir=char_del(out_dir);
  index=ivec_free(index);
  dx1=dvec_free(dx1);
  dx2=dvec_free(dx2);
  zx1=zvec_free(zx1);
  zx2=zvec_free(zx2);
  x1_z=zvec_free(x1_z);
  RVF(rx1,n1);
  RVF(rx2,n2);
  CVF(cx1,n1);
  CVF(cx2,n2);
  dA1=dmat_free(dA1);
  dA2=dmat_free(dA2);
  zA1=zmat_free(zA1);
  zA2=zmat_free(zA2);
  rA1=rmat_free(m1,n1,rA1);
  rA2=rmat_free(m2,n2,rA2);
  rE=rmat_free(m2,n2,rE);
  cA1=cmat_free(m1,n1,cA1);
  cA2=cmat_free(m2,n2,cA2);
  cE=cmat_free(m2,n2,cE);
  RF(rEmax);
  RF(rmax);
  return 0;
}
