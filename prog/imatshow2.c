#include<isys.h>

#define PROG "imatshow2"
#define MAX_LENGTH 10000
#define BIN_SUFFIX "dat"
#define TXT_SUFFIX "txt"
#define METHOD_NAME_LENGTH  128
#define METHOD_RSUB       "rsub"
#define METHOD_CSUB       "csub"
#define METHOD_RREL       "rrelerror"
#define METHOD_CREL       "crelerror"
#define METHOD_RABS       "rabserror"
#define METHOD_CABS       "cabserror"
#define DEFAULT_PREC     1024

#define RA(X,P)    ((X)=rallocate_prec(P))
#define RAP(X,Y)   ((X)=rallocate_prec(rget_prec(Y)))
#define RF(X)      ((X)=rfree(X))
#define CA(X,P)    ((X)=callocate_prec(P))
#define CF(X)      ((X)=cfree(X))
#define RVA(X,N,P)  { X=rvec_allocate_prec(N,P); }
#define RVAr(X,Y,N) { X=rvec_allocate_prec(N,rget_prec(Y)); }
#define RVArv(X,Y,N){ X=rvec_allocate_prec(N,rvec_get_prec_max(N,Y)); }
#define RVF(X,N)    { X=rvec_free(N,X); }
#define CVA(X,N,P)  { X=cvec_allocate_prec(N,P); }
#define CVF(X,N)    { X=cvec_free(N,X); }


enum { NOTYPE=-1, TEXT=0, DVEC, ZVEC, RVEC, CVEC, DMAT, ZMAT, RMAT, CMAT };
char *type_name[]={ "text", "dvec", "zvec", "rvec", "cvec", "dmat", "zmat", "rmat", "cmat" };

int irvec_relerror(int n, rmulti **z0, rmulti **z1, rmulti **x0, rmulti **x1, rmulti **y0, rmulti **y1)
{
  int e=0;
  rmulti **ry0=NULL,**ry1=NULL;
  RVArv(ry0,z0,n); RVArv(ry1,z1,n);
  irvec_abs_sub(n,z0,z1,x0,x1,y0,y1);
  irvec_abs(n,ry0,ry1,y0,y1);
  irvec_div(n,z0,z1,z0,z1,ry0,ry1);
  RVF(ry0,n); RVF(ry1,n);
  return e;
}

int irvec_abserror(rmulti *z0, rmulti *z1, int n, rmulti **x0, rmulti **x1, rmulti **y0, rmulti **y1)
{
  int e=0;
  rmulti *ry0=NULL,*ry1=NULL;
  rmulti **e0=NULL,**e1=NULL;
  RAP(ry0,z0); RAP(ry1,z1);
  RVAr(e0,z0,n); RVAr(e1,z1,n);
  irvec_abs_sub(n,e0,e1,x0,x1,y0,y1);
  irvec_umax(z0,z1,n,e0,e1);
  irvec_umax_abs(ry0,ry1,n,y0,y1);
  irdiv(z0,z1,z0,z1,ry0,ry1);
  RF(ry0); RF(ry1);
  RVF(e0,n); RVF(e1,n);
  return e;
}

int icvec_relerror(int n, rmulti **z0, rmulti **z1, cmulti **x0, cmulti **x1, cmulti **y0, cmulti **y1)
{
  int e=0;
  rmulti **ry0=NULL,**ry1=NULL;
  RVArv(ry0,z0,n); RVArv(ry1,z1,n);
  icvec_abs_sub(n,z0,z1,x0,x1,y0,y1);
  icvec_abs(n,ry0,ry1,y0,y1);
  irvec_div(n,z0,z1,z0,z1,ry0,ry1);
  RVF(ry0,n); RVF(ry1,n);
  return e;
}
 
int icvec_abserror(rmulti *z0, rmulti *z1, int n, cmulti **x0, cmulti **x1, cmulti **y0, cmulti **y1)
{
  int e=0;
  rmulti *ry0=NULL,*ry1=NULL;
  rmulti **e0=NULL,**e1=NULL;
  RAP(ry0,z0); RAP(ry1,z1);
  RVAr(e0,z0,n); RVAr(e1,z1,n);
  icvec_abs_sub(n,e0,e1,x0,x1,y0,y1);
  irvec_umax(z0,z1,n,e0,e1);
  icvec_umax_abs(ry0,ry1,n,y0,y1);
  irdiv(z0,z1,z0,z1,ry0,ry1);
  RF(ry0); RF(ry1);
  RVF(e0,n); RVF(e1,n);
  return e;
}

void usage()
{
  printf("Usage: %s input.dat [options]\n",PROG);
  printf("  -input1            Input first file.\n");
  printf("  -input2            Input second file.\n");
  printf("  -interval1         Input first file's interval.\n");
  printf("  -interval2         Input second file's interval.\n");
  printf("  -rsub              Do subtraction on real field.\n");
  printf("  -csub              Do subtraction on complex field.\n");
  printf("  -rrel              Calculate relative error on real field.\n");
  printf("  -crel              Calculate relative error on complex field.\n");
  printf("  -rabs              Calculate absolute error on real field.\n");
  printf("  -cabs              Calculate absolute error on complex field.\n");
  printf("  -o-mid fname       The matrix is saved in the file 'fname'.\n");
  printf("  -o-rad fname       The matrix is saved in the file 'fname'.\n");
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
  printf("     # %s -input1 _____.dat -input2 _____.dat -interval1 _________.dat \n",PROG);
  printf("     # %s -input1 _____.dat -input2 _____.dat -interval1 _________.dat -interval2 _________.dat \n",PROG);
  printf("\n");
  exit(0);
}

int main(int argc, char *argv[])
{
  FILE *fid=NULL;
  //FILE *fid=NULL,*fid1=NULL,*fid2=NULL;
  strings *path=NULL;
  char *in_fname1=NULL,*in_fname2=NULL,*in_fname3=NULL,*in_fname4=NULL;
  char *out_fname_log=NULL,*out_fname_m=NULL,*out_fname_r=NULL,*out_dir=NULL;
  char format[]="e",method[METHOD_NAME_LENGTH+1],buf[MAX_LENGTH],mode='?';
  int m1=0,m2=0,m3=0,m4=0,n1=0,n2=0,n3=0,n4=0,debug=1,i,prec=DEFAULT_PREC,flag=0;
  int loaded=0,autoname=0,type_in1=NOTYPE,type_in2=NOTYPE,type_in3=NOTYPE,type_in4=NOTYPE,digits=7,*index=NULL;
  double *dx1=NULL,*dx2=NULL,*dx3=NULL,*dx4=NULL,*dA1=NULL,*dA2=NULL,*dA3=NULL,*dA4=NULL;
  dcomplex *zx1=NULL,*zx2=NULL,*zx3=NULL,*zx4=NULL,*zA1=NULL,*zA2=NULL,*zA3=NULL,*zA4=NULL;
  dcomplex *x1_z=NULL;
  rmulti **rx1=NULL,**rx1_1=NULL,**rx1_2=NULL,**rx2=NULL,**rx3=NULL,**rx3_1=NULL,**rx3_2=NULL,**rx4=NULL,**rA1=NULL,**rA2=NULL,**rA3=NULL,**rA4=NULL;
  rmulti **mid_vec=NULL,**rad_vec=NULL,*mid=NULL,*rad=NULL;
  cmulti **cx1=NULL,**cx1_1=NULL,**cx1_2=NULL,**cx2=NULL,**cx3=NULL,**cx3_1=NULL,**cx3_2=NULL,**cx4=NULL,**cA1=NULL,**cA2=NULL,**cA3=NULL,**cA4=NULL;

  // init
  in_fname1=char_new("",NULL);
  in_fname2=char_new("",NULL);
  in_fname3=char_new("",NULL);
  in_fname4=char_new("",NULL);
  out_fname_log=char_new("",NULL);
  out_fname_m=char_new("",NULL);
  out_fname_r=char_new("",NULL);
  out_dir=char_new("",NULL);
  strcpy(method,METHOD_RSUB);

  // options
  i=1;
  while(i<argc){
    if     (STR_EQ(argv[i],"-help"))                                 { usage(); }
    else if(STR_EQ(argv[i],"-input1") && i+1<argc)                   { in_fname1=char_renew(in_fname1,argv[++i],NULL); }
    else if(STR_EQ(argv[i],"-input2") && i+1<argc)                   { in_fname3=char_renew(in_fname3,argv[++i],NULL); }
    else if(STR_EQ(argv[i],"-interval1") && i+1<argc)                { in_fname2=char_renew(in_fname2,argv[++i],NULL); }
    else if(STR_EQ(argv[i],"-interval2") && i+1<argc)                { in_fname4=char_renew(in_fname4,argv[++i],NULL); }
    else if(STR_EQ(argv[i],"-rsub"))                                 { flag=1; strcpy(method,METHOD_RSUB); }
    else if(STR_EQ(argv[i],"-csub"))                                 { flag=1; strcpy(method,METHOD_CSUB); }
    else if(STR_EQ(argv[i],"-rrel"))                                 { flag=2; strcpy(method,METHOD_RREL);}
    else if(STR_EQ(argv[i],"-crel"))                                 { flag=2; strcpy(method,METHOD_CREL);}
    else if(STR_EQ(argv[i],"-rabs"))                                 { flag=3; strcpy(method,METHOD_RABS);}
    else if(STR_EQ(argv[i],"-cabs"))                                 { flag=3; strcpy(method,METHOD_CABS);}
    else if(STR_EQ(argv[i],"-o-mid") && i+1<argc)                    { ++i; out_fname_m=char_renew(out_fname_m,argv[i],NULL); }
    else if(STR_EQ(argv[i],"-o-rad") && i+1<argc)                    { ++i; out_fname_r=char_renew(out_fname_r,argv[i],NULL); }
    else if(STR_EQ(argv[i],"-O"))                                    { autoname=1; }
    else if(STR_EQ(argv[i],"-d") && i+1<argc)                        { ++i; out_dir=char_renew(out_dir,argv[i],NULL); if(out_dir[strlen(out_dir)-1]!='/'){ out_dir=char_cat(out_dir,"/"); } }
    else if(STR_EQ(argv[i],"-f") && i+1<argc)                        { strcpy(format,"f"); ++i; digits=atoi(argv[i]); }
    else if(STR_EQ(argv[i],"-e") && i+1<argc)                        { strcpy(format,"e"); ++i; digits=atoi(argv[i]); }
    else if(STR_EQ(argv[i],"-prec") && i+1<argc)                     { strcpy(format,"e"); ++i; prec=atoi(argv[i]); }
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
  else if(char_eq(method,METHOD_RREL)){ mode='r'; }
  else if(char_eq(method,METHOD_CREL)){ mode='c'; }
  else if(char_eq(method,METHOD_RABS)){ mode='r'; }
  else if(char_eq(method,METHOD_CABS)){ mode='c'; }
  else{ ERROR_AT; exit(0); }

  // load input file1
  if(strlen(in_fname1)<=0){ usage(); }
  if(!loaded){ dx1=dvec_bin_load(&n1,"%s",in_fname1); }     if(dx1!=NULL){ type_in1=DVEC; loaded=1; }
  if(!loaded){ zx1=zvec_bin_load(&n1,"%s",in_fname1); }     if(zx1!=NULL){ type_in1=ZVEC; loaded=1; }
  if(!loaded){ rx1=rvec_bin_load(&n1,"%s",in_fname1); }     if(rx1!=NULL){ type_in1=RVEC; loaded=1; }
  if(!loaded){ cx1=cvec_bin_load(&n1,"%s",in_fname1); }     if(cx1!=NULL){ type_in1=CVEC; loaded=1; }
  if(!loaded){ dA1=dmat_bin_load(&m1,&n1,"%s",in_fname1); } if(dA1!=NULL){ type_in1=DMAT; loaded=1; }
  if(!loaded){ zA1=zmat_bin_load(&m1,&n1,"%s",in_fname1); } if(zA1!=NULL){ type_in1=ZMAT; loaded=1; }
  if(!loaded){ rA1=rmat_bin_load(&m1,&n1,"%s",in_fname1); } if(rA1!=NULL){ type_in1=RMAT; loaded=1; }
  if(!loaded){ cA1=cmat_bin_load(&m1,&n1,"%s",in_fname1); } if(cA1!=NULL){ type_in1=CMAT; loaded=1; }
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
    if(dx1!=NULL){ dvec_print(n1,dx1,NULL,format,digits); }
    if(zx1!=NULL){ zvec_print(n1,zx1,NULL,format,digits); }
    if(rx1!=NULL){ rvec_print(n1,rx1,NULL,format,digits); }
    if(cx1!=NULL){ cvec_print(n1,cx1,NULL,format,digits); }
    if(dA1!=NULL){ dmat_print(m1,n1,dA1,m2,NULL,format,digits); }
    if(zA1!=NULL){ zmat_print(m1,n1,zA1,m2,NULL,format,digits); }
    if(rA1!=NULL){ rmat_print(m1,n1,rA1,m2,NULL,format,digits); }
    if(cA1!=NULL){ cmat_print(m1,n1,cA1,m2,NULL,format,digits); }
  }

  // load input file2(interval file)
  loaded=0;
  if(strlen(in_fname2)<=0){ usage(); }
  if(!loaded){ dx2=dvec_bin_load(&n2,"%s",in_fname2); }     if(dx2!=NULL){ type_in2=DVEC; loaded=1; }
  if(!loaded){ zx2=zvec_bin_load(&n2,"%s",in_fname2); }     if(zx2!=NULL){ type_in2=ZVEC; loaded=1; }
  if(!loaded){ rx2=rvec_bin_load(&n2,"%s",in_fname2); }     if(rx2!=NULL){ type_in2=RVEC; loaded=1; }
  if(!loaded){ cx2=cvec_bin_load(&n2,"%s",in_fname2); }     if(cx2!=NULL){ type_in2=CVEC; loaded=1; }
  if(!loaded){ dA2=dmat_bin_load(&m2,&n2,"%s",in_fname2); } if(dA2!=NULL){ type_in2=DMAT; loaded=1; }
  if(!loaded){ zA2=zmat_bin_load(&m2,&n2,"%s",in_fname2); } if(zA2!=NULL){ type_in2=ZMAT; loaded=1; }
  if(!loaded){ rA2=rmat_bin_load(&m2,&n2,"%s",in_fname2); } if(rA2!=NULL){ type_in2=RMAT; loaded=1; }
  if(!loaded){ cA2=cmat_bin_load(&m2,&n2,"%s",in_fname2); } if(cA2!=NULL){ type_in2=CMAT; loaded=1; }
  if(!loaded){ printf("Error! The input file '%s%s' can NOT be opened, or is NOT supported type.\nPlease use matshow2.\n",out_dir,in_fname2); usage(); }
  if(debug>=2){
    printf("Input1_interval_File: %s\n",in_fname2);
    printf("Input1_interval_Type: %s\n",type_name[type_in2]);
    if(DVEC<=type_in2 && type_in2<=CVEC){ printf("Input1_interval_Size: %d\n",n2); }
    if(DMAT<=type_in2 && type_in2<=CMAT){ printf("Input1_interval_Size: %d-by-%d\n",m2,n2); }
    if(type_in2==DVEC || type_in2==ZVEC || type_in2==DMAT || type_in2==ZMAT){ printf("Input1_interval_Precision: double\n"); }
    else if(type_in2==RVEC){ printf("Input1_interval_Precision_Max: %d\n",rvec_get_prec_max(n2,rx2)); }
    else if(type_in2==CVEC){ printf("Input1_interval_Precision_Max: %d\n",cvec_get_prec_max(n2,cx2)); }
    else if(type_in2==RMAT){ printf("Input1_interval_Precision_Max: %d\n",rmat_get_prec_max(m2,n2,rA2,m2)); }
    else if(type_in2==CMAT){ printf("Input1_interval_Precision_Max: %d\n",cmat_get_prec_max(m2,n2,cA2,m2)); }
  }
  if(debug>=3){
    if(DVEC<=type_in2 && type_in2<=CVEC){ printf("Input1_interval_Vector:\n"); }
    if(DMAT<=type_in2 && type_in2<=CMAT){ printf("Input1_interval_Matrix:\n"); }
    if(dx2!=NULL){ dvec_print(n2,dx2,NULL,format,digits); }
    if(zx2!=NULL){ zvec_print(n2,zx2,NULL,format,digits); }
    if(rx2!=NULL){ rvec_print(n2,rx2,NULL,format,digits); }
    if(cx2!=NULL){ cvec_print(n2,cx2,NULL,format,digits); }
    if(dA2!=NULL){ dmat_print(m2,n2,dA2,m2,NULL,format,digits); }
    if(zA2!=NULL){ zmat_print(m2,n2,zA2,m2,NULL,format,digits); }
    if(rA2!=NULL){ rmat_print(m2,n2,rA2,m2,NULL,format,digits); }
    if(cA2!=NULL){ cmat_print(m2,n2,cA2,m2,NULL,format,digits); }
  }

  // load input file3
  loaded=0;
  if(strlen(in_fname3)<=0){ usage(); }
  //if(!loaded){ dx3=dvec_bin_load(&n3,"%s",in_fname3); } if(dx3!=NULL){ type_in3=DVEC; loaded=1; }
  //if(!loaded){ zx3=zvec_bin_load(&n3,"%s",in_fname3); } if(zx3!=NULL){ type_in3=ZVEC; loaded=1; }
  if(!loaded && mode=='r'){ rx3=rvec_bin_load(&n3,"%s",in_fname3); } if(rx3!=NULL){ type_in3=RVEC; loaded=1; }
  if(!loaded && mode=='c'){ cx3=cvec_bin_load(&n3,"%s",in_fname3); } if(cx3!=NULL){ type_in3=CVEC; loaded=1; }
  //if(!loaded){ dA3=dmat_bin_load(&m3,&n3,"%s",in_fname3); } if(dA3!=NULL){ type_in3=DMAT; loaded=1; }
  //if(!loaded){ zA3=zmat_bin_load(&m3,&n3,"%s",in_fname3); } if(zA3!=NULL){ type_in3=ZMAT; loaded=1; }
  //if(!loaded){ rA3=rmat_bin_load(&m3,&n3,"%s",in_fname3); } if(rA3!=NULL){ type_in3=RMAT; loaded=1; }
  //if(!loaded){ cA3=cmat_bin_load(&m3,&n3,"%s",in_fname3); } if(cA3!=NULL){ type_in3=CMAT; loaded=1; }
  if(!loaded){ printf("Error! The input file '%s%s' can NOT be opened, or is NOT supported type.\n\n",out_dir,in_fname3); usage(); }
  if(debug>=2){
    printf("Input2_File: %s\n",in_fname3);
    printf("Input2_Type: %s\n",type_name[type_in3]);
    if(DVEC<=type_in3 && type_in3<=CVEC){ printf("Input2_Size: %d\n",n3); }
    if(DMAT<=type_in3 && type_in3<=CMAT){ printf("Input2_Size: %d-by-%d\n",m3,n3); }
    if(type_in3==DVEC || type_in3==ZVEC || type_in3==DMAT || type_in3==ZMAT){ printf("Input2_Precision: double\n"); }
    else if(type_in3==RVEC){ printf("Input2_Precision_Max: %d\n",rvec_get_prec_max(n3,rx3)); }
    else if(type_in3==CVEC){ printf("Input2_Precision_Max: %d\n",cvec_get_prec_max(n3,cx3)); }
    else if(type_in3==RMAT){ printf("Input2_Precision_Max: %d\n",rmat_get_prec_max(m3,n3,rA3,m3)); }
    else if(type_in3==CMAT){ printf("Input2_Precision_Max: %d\n",cmat_get_prec_max(m3,n3,cA3,m3)); }
  }
  if(debug>=3){
    if(DVEC<=type_in3 && type_in3<=CVEC){ printf("Input2_Vector:\n"); }
    if(DMAT<=type_in3 && type_in3<=CMAT){ printf("Input2_Matrix:\n"); }
    if(dx3!=NULL){ dvec_print(n3,dx3,NULL,format,digits); }
    if(zx3!=NULL){ zvec_print(n3,zx3,NULL,format,digits); }
    if(rx3!=NULL){ rvec_print(n3,rx3,NULL,format,digits); }
    if(cx3!=NULL){ cvec_print(n3,cx3,NULL,format,digits); }
    if(dA3!=NULL){ dmat_print(m3,n3,dA3,m3,NULL,format,digits); }
    if(zA3!=NULL){ zmat_print(m3,n3,zA3,m3,NULL,format,digits); }
    if(rA3!=NULL){ rmat_print(m3,n3,rA3,m3,NULL,format,digits); }
    if(cA3!=NULL){ cmat_print(m3,n3,cA3,m3,NULL,format,digits); }
  }

  // load input file4(interval file)
  loaded=0;
  if(!loaded){ dx4=dvec_bin_load(&n4,"%s",in_fname4); }     if(dx4!=NULL){ type_in4=DVEC; loaded=1; }
  if(!loaded){ zx4=zvec_bin_load(&n4,"%s",in_fname4); }     if(zx4!=NULL){ type_in4=ZVEC; loaded=1; }
  if(!loaded){ rx4=rvec_bin_load(&n4,"%s",in_fname4); }     if(rx4!=NULL){ type_in4=RVEC; loaded=1; }
  if(!loaded){ cx4=cvec_bin_load(&n4,"%s",in_fname4); }     if(cx4!=NULL){ type_in4=CVEC; loaded=1; }
  if(!loaded){ dA4=dmat_bin_load(&m4,&n4,"%s",in_fname4); } if(dA4!=NULL){ type_in4=DMAT; loaded=1; }
  if(!loaded){ zA4=zmat_bin_load(&m4,&n4,"%s",in_fname4); } if(zA4!=NULL){ type_in4=ZMAT; loaded=1; }
  if(!loaded){ rA4=rmat_bin_load(&m4,&n4,"%s",in_fname4); } if(rA4!=NULL){ type_in4=RMAT; loaded=1; }
  if(!loaded){ cA4=cmat_bin_load(&m4,&n4,"%s",in_fname4); } if(cA4!=NULL){ type_in4=CMAT; loaded=1; }
  if(debug>=2){
    printf("Input2_interval_File: %s\n",in_fname4);
    printf("Input2_interval_Type: %s\n",type_name[type_in4]);
    if(DVEC<=type_in4 && type_in4<=CVEC){ printf("Input2_interval_Size: %d\n",n4); }
    if(DMAT<=type_in4 && type_in4<=CMAT){ printf("Input2_interval_Size: %d-by-%d\n",m4,n4); }
    if(type_in4==DVEC || type_in4==ZVEC || type_in4==DMAT || type_in4==ZMAT){ printf("Input2_interval_Precision: double\n"); }
    else if(type_in4==RVEC){ printf("Input2_interval_Precision_Max: %d\n",rvec_get_prec_max(n4,rx4)); }
    else if(type_in4==CVEC){ printf("Input2_interval_Precision_Max: %d\n",cvec_get_prec_max(n4,cx4)); }
    else if(type_in4==RMAT){ printf("Input2_interval_Precision_Max: %d\n",rmat_get_prec_max(m4,n4,rA4,m4)); }
    else if(type_in4==CMAT){ printf("Input2_interval_Precision_Max: %d\n",cmat_get_prec_max(m4,n4,cA4,m4)); }
  }
  if(debug>=3){
    if(DVEC<=type_in4 && type_in4<=CVEC){ printf("Input2_interval_Vector:\n"); }
    if(DMAT<=type_in4 && type_in4<=CMAT){ printf("Input2_interval_Matrix:\n"); }
    if(dx4!=NULL){ dvec_print(n4,dx4,NULL,format,digits); }
    if(zx4!=NULL){ zvec_print(n4,zx4,NULL,format,digits); }
    if(rx4!=NULL){ rvec_print(n4,rx4,NULL,format,digits); }
    if(cx4!=NULL){ cvec_print(n4,cx4,NULL,format,digits); }
    if(dA4!=NULL){ dmat_print(m4,n4,dA4,m4,NULL,format,digits); }
    if(zA4!=NULL){ zmat_print(m4,n4,zA4,m4,NULL,format,digits); }
    if(rA4!=NULL){ rmat_print(m4,n4,rA4,m4,NULL,format,digits); }
    if(cA4!=NULL){ cmat_print(m4,n4,cA4,m4,NULL,format,digits); }
  }

  // check
  if(n1!=n2){ printf("Error! The input file's size is NOT same!\n\n"); exit(0); }

  // get prec
  if(mode=='r'){
    prec=MAX2(rvec_get_prec_max(n1,rx1),rvec_get_prec_max(n3,rx3));
    if(rvec_get_prec_max(n2,rx2)>prec){ prec=rvec_get_prec_max(n2,rx2); }
  }
  if(mode=='c'){
    prec=MAX2(cvec_get_prec_max(n1,cx1),cvec_get_prec_max(n3,cx3));
    if(cvec_get_prec_max(n2,cx2)>prec){ prec=cvec_get_prec_max(n2,cx2); }
  }

  // autoname
  if(autoname){
    path=strings_split_path(in_fname3,BIN_SUFFIX);
    if(strlen(out_dir)<=0 && strlen(path->str[0])>0){ out_dir=char_renew(out_dir,path->str[0],NULL); out_dir=char_cat(out_dir,"/"); }
    if(flag==1){
      sprintf(buf,"%s_text{%dbits|%s|log}.%s",path->str[1],prec,method,TXT_SUFFIX); out_fname_log=char_cat(out_fname_log,buf);
      sprintf(buf,"%s_%cvec{%dbits|%s|mid}.%s",path->str[1],mode,prec,method,BIN_SUFFIX); out_fname_m=char_cat(out_fname_m,buf);
      sprintf(buf,"%s_%cvec{%dbits|%s|rad}.%s",path->str[1],mode,prec,method,BIN_SUFFIX); out_fname_r=char_cat(out_fname_r,buf);
    }
    else if(flag==2){
      sprintf(buf,"%s_text{%dbits|%s|log}.%s",path->str[1],prec,method,TXT_SUFFIX); out_fname_log=char_cat(out_fname_log,buf);
      sprintf(buf,"%s_rvec{%dbits|%s|mid}.%s",path->str[1],prec,method,BIN_SUFFIX); out_fname_m=char_cat(out_fname_m,buf); 
      sprintf(buf,"%s_rvec{%dbits|%s|rad}.%s",path->str[1],prec,method,BIN_SUFFIX); out_fname_r=char_cat(out_fname_r,buf); 
    }
    else if(flag==3){
      sprintf(buf,"%s_text{%dbits|%s|log}.%s",path->str[1],prec,method,TXT_SUFFIX); out_fname_log=char_cat(out_fname_log,buf);
      sprintf(buf,"%s_rmulti{%dbits|%s|mid}.%s",path->str[1],prec,method,BIN_SUFFIX); out_fname_m=char_cat(out_fname_m,buf);
      sprintf(buf,"%s_rmulti{%dbits|%s|rad}.%s",path->str[1],prec,method,BIN_SUFFIX); out_fname_r=char_cat(out_fname_r,buf);
    }
    path=strings_del(path);
  }
  if(strlen(out_fname_log)>0){
    strcpy(buf,out_dir); strcat(buf,"/"); strcat(buf,out_fname_log);
    fid=fopen(buf,"w");
  }

  // log
  if(debug>=1){
    printf("Debug_Level: %d\n",debug);          if(fid!=NULL){ fprintf(fid,"Debug_Level: %d\n",debug); }
    printf("Precision: %d\n",prec);             if(fid!=NULL){ fprintf(fid,"Precision: %d\n",prec); }
    printf("Method: %s\n",method);              if(fid!=NULL){ fprintf(fid,"Method: %s\n",method); }
    printf("Input_Vector_Size: %d\n",n1);       if(fid!=NULL){ fprintf(fid,"Input_Vector_Size: %d\n",n1); }
    printf("Input1_File: %s\n",in_fname1);      if(fid!=NULL){ fprintf(fid,"Input_File1: %s\n",in_fname1); }
    printf("Input1_Interval: %s\n",in_fname2);  if(fid!=NULL){ fprintf(fid,"Input_interval1: %s\n",in_fname2); }
    printf("Input2_File: %s\n",in_fname3);      if(fid!=NULL){ fprintf(fid,"Input_File2: %s\n",in_fname3); }
    printf("Input2_Interval: %s\n",in_fname4);  if(fid!=NULL){ fprintf(fid,"Input_interval2: %s\n",in_fname4); }
    printf("Mode: %s\n",(mode=='r'?"real":(mode=='c'?"complex":"?")));  if(fid!=NULL){ fprintf(fid,"Mode: %s\n",(mode=='r'?"real":(mode=='c'?"complex":"?"))); }
  }
  if(debug>=1){
    printf("Output_File_Log:  %s%s\n",out_dir,(strlen(out_fname_log)>0?out_fname_log:"NULL"));
    if(fid!=NULL){ fprintf(fid,"Output_File_Log:  %s%s\n",out_dir,(strlen(out_fname_log)>0?out_fname_log:"NULL")); }
    printf("Output_File_m:  %s%s\n",out_dir,(strlen(out_fname_m)>0?out_fname_m:"NULL"));
    if(fid!=NULL){ fprintf(fid,"Output_File_m:  %s%s\n",out_dir,(strlen(out_fname_m)>0?out_fname_m:"NULL")); }
    printf("Output_File_r: %s%s\n",out_dir,(strlen(out_fname_r)>0?out_fname_r:"NULL"));
    if(fid!=NULL){ fprintf(fid,"Output_File_r: %s%s\n",out_dir,(strlen(out_fname_r)>0?out_fname_r:"NULL")); }
  }

  // sorting
  if(type_in1==RVEC && type_in3==RVEC){
    index=ivec_allocate(n1);
    x1_z=zvec_allocate(n1);
    rvec_get_z(n1,x1_z,rx1);
    zvec_sort(n1,x1_z,index); zvec_reverse(n1,x1_z); ivec_reverse(n1,index);
    rvec_swap_index(n1,rx1,index);
    rvec_swap_index(n2,rx2,index);
    ivec_grid(n1,index);
    reig_sort_value_guide(n1,n1,rx3,NULL,0,rx1);
  }
  if(type_in1==CVEC && type_in3==CVEC){
    index=ivec_allocate(n1);
    x1_z=zvec_allocate(n1);
    cvec_get_z(n1,x1_z,cx1);
    zvec_sort(n1,x1_z,index); zvec_reverse(n1,x1_z); ivec_reverse(n1,index);
    cvec_swap_index(n1,cx1,index);
    cvec_swap_index(n2,cx2,index);
    ivec_grid(n1,index);
    ceig_sort_value_guide(n1,n1,cx3,NULL,0,cx1);
  }

  //add pm
  if(type_in1==RVEC && type_in3==RVEC){
    RVA(rx1_1,n1,prec); RVA(rx1_2,n1,prec);
    RVA(rx3_1,n1,prec); RVA(rx3_2,n1,prec);
    if(type_in2==RVEC){ irvec_add_pm(n1,rx1_1,rx1_2,rx1,rx1,rx2); }
    else              { irvec_copy(n1,rx1_1,rx1_2,rx1,rx1); }
    if(type_in4==RVEC){ irvec_add_pm(n1,rx3_1,rx3_2,rx3,rx3,rx4); }
    else              { irvec_copy(n1,rx3_1,rx3_2,rx3,rx3); }
  }
  if(type_in1==CVEC && type_in3==CVEC){
    CVA(cx1_1,n1,prec); CVA(cx1_2,n1,prec);
    CVA(cx3_1,n1,prec); CVA(cx3_2,n1,prec);
    if(type_in2==CVEC){ icvec_add_pm(n1,cx1_1,cx1_2,cx1,cx1,cx2); }
    else              { icvec_copy(n1,cx1_1,cx1_2,cx1,cx1); }
    if(type_in4==CVEC){ icvec_add_pm(n1,cx3_1,cx3_2,cx3,cx3,cx4); }
    else              { icvec_copy(n1,cx3_1,cx3_2,cx3,cx3); }
  }

  //compute mode [sub] 
  if(flag==1 && type_in1==RVEC && type_in3==RVEC){ irvec_sub(n1,rx1_1,rx1_2,rx3_1,rx3_2,rx1_1,rx1_2); }
  if(flag==1 && type_in1==CVEC && type_in3==CVEC){ icvec_sub(n1,cx1_1,cx1_2,cx3_1,cx3_2,cx1_1,cx1_2); }

  //compute mode [rel]
  if(flag==2){
    RVA(mid_vec,n1,prec);  RVA(rad_vec,n1,prec);
    RA(mid,prec); RA(rad,prec);
    if(type_in1==RVEC && type_in3==RVEC){
      irvec_relerror(n1,mid_vec,rad_vec,rx3_1,rx3_2,rx1_1,rx1_2);
      irvec_umax(mid,rad,n1,mid_vec,rad_vec);
      irvec_mr(n1,mid_vec,rad_vec,mid_vec,rad_vec);
      irmr(mid,rad,mid,rad);
    }
    if(type_in1==CVEC && type_in3==CVEC){
      icvec_relerror(n1,mid_vec,rad_vec,cx3_1,cx3_2,cx1_1,cx1_2);
      irvec_umax(mid,rad,n1,mid_vec,rad_vec);
      irvec_mr(n1,mid_vec,rad_vec,mid_vec,rad_vec);
      irmr(mid,rad,mid,rad);
    }
  }

  //compute mode [abs]
  if(flag==3){
    RA(mid,prec); RA(rad,prec);
    if(type_in1==RVEC && type_in3==RVEC){
      irvec_abserror(mid,rad,n1,rx3_1,rx3_2,rx1_1,rx1_2);
      irmr(mid,rad,mid,rad);
    }
    if(type_in1==CVEC && type_in3==CVEC){
      icvec_abserror(mid,rad,n1,cx3_1,cx3_2,cx1_1,cx1_2);
      irmr(mid,rad,mid,rad);
    }
  }
  
  //output results
  if(debug>=2){
    if(flag==1 && rx1_1!=NULL && rx1_2!=NULL){ printf("Output_Error=\n"); irvec_print(n1,rx1_1,rx1_2,NULL,"e",10); }
    if(flag==1 && cx1_1!=NULL && cx1_2!=NULL){ printf("Output_Error=\n"); icvec_print(n1,cx1_1,cx1_2,NULL,"e",10); }
    if(flag==2){ printf("Output_Error_mid_rad=\n"); irvec_print(n1,mid_vec,rad_vec,NULL,"e",10); }
  }
  if(debug>=1){
    if(flag==1){}
    if(flag==2 || flag==3){
      mpfr_printf("Output_Error_MAX_mid_rad: %.7Re  %.7Re\n",mid,rad);
      if(fid!=NULL){
	mpfr_fprintf(fid,"Output_Error_size_and_MAX_mid: %d\t%.7Re\n",n1,mid); 
	mpfr_fprintf(fid,"Output_Error_size_and_MAX_rad: %d\t%.7Re\n",n1,rad); 
      }
    }
  }

  //save
  if(strlen(out_fname_m)>0){
    if(flag==1 && mid_vec!=NULL){ }
    if(flag==2 && mid_vec!=NULL){ rvec_bin_save(n1,mid_vec,"%s%s",out_dir,out_fname_m); }
    if(flag==3 && mid!=NULL){ 
      //strcpy(buf,out_dir); strcat(buf,"/"); strcat(buf,out_fname_m);
      //fid1=fopen(buf,"w");
      //if(fid1!=NULL){ rbin_save(mid,fid1); }
    }
    if(debug>=1){ fprintf(stderr,"saved: %s%s\n",out_dir,out_fname_m); }
  }
  if(strlen(out_fname_r)>0){
    if(flag==1 && rad_vec!=NULL){ }
    if(flag==2 && rad_vec!=NULL){ rvec_bin_save(n1,rad_vec,"%s%s",out_dir,out_fname_r); }
    if(flag==3 && rad!=NULL){ 
      //strcpy(buf,out_dir); strcat(buf,"/"); strcat(buf,out_fname_r);
      //fid2=fopen(buf,"w");
      //if(fid2!=NULL){ rbin_save(rad,fid2); }
    }
    if(debug>=1){ fprintf(stderr,"saved: %s%s\n",out_dir,out_fname_r); }
  }

  // done
  if(fid!=NULL){ fclose(fid); fid=NULL; }
  //if(fid1!=NULL){ fclose(fid1); fid1=NULL; }
  //if(fid2!=NULL){ fclose(fid2); fid2=NULL; }
  in_fname1=char_del(in_fname1);
  in_fname2=char_del(in_fname2);
  in_fname3=char_del(in_fname3);
  in_fname4=char_del(in_fname4);
  out_fname_log=char_del(out_fname_log);
  out_fname_m=char_del(out_fname_m);
  out_fname_r=char_del(out_fname_r);
  out_dir=char_del(out_dir);
  index=ivec_free(index);
  dx1=dvec_free(dx1);
  dx2=dvec_free(dx2);
  dx3=dvec_free(dx3);
  dx4=dvec_free(dx4);
  zx1=zvec_free(zx1);
  zx2=zvec_free(zx2);
  zx3=zvec_free(zx3);
  zx4=zvec_free(zx4);
  x1_z=zvec_free(x1_z);
  RVF(rx1,n1);
  RVF(rx1_1,n1);
  RVF(rx1_2,n1);
  RVF(rx2,n2);
  RVF(rx3,n3);
  RVF(rx3_1,n1);
  RVF(rx3_2,n1);
  RVF(rx4,n4);
  CVF(cx1,n1);
  CVF(cx1_1,n1);
  CVF(cx1_2,n1);
  CVF(cx2,n2);
  CVF(cx3,n3);
  CVF(cx3_1,n3);
  CVF(cx3_2,n3);
  CVF(cx4,n4);
  dA1=dmat_free(dA1);
  dA2=dmat_free(dA2);
  dA3=dmat_free(dA3);
  dA4=dmat_free(dA4);
  zA1=zmat_free(zA1);
  zA2=zmat_free(zA2);
  zA3=zmat_free(zA3);
  zA4=zmat_free(zA4);
  rA1=rmat_free(m1,n1,rA1);
  rA2=rmat_free(m2,n2,rA2);
  rA3=rmat_free(m3,n3,rA3);
  rA4=rmat_free(m4,n4,rA4);
  cA1=cmat_free(m1,n1,cA1);
  cA2=cmat_free(m2,n2,cA2);
  cA3=cmat_free(m3,n3,cA3);
  cA4=cmat_free(m4,n4,cA4);
  RVF(mid_vec,n1);
  RVF(rad_vec,n1);
  RF(mid);
  RF(rad);
  return 0;
}
