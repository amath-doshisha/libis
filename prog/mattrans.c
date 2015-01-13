#include<isys.h>

#define PROG "mattrans"
#define BASENAME_APPEND "_T"
#define BIN_SUFFIX "dat"

enum { NOTYPE=-1, DVEC=0, ZVEC, RVEC, CVEC, DMAT, ZMAT, RMAT, CMAT };
char *type_name[]={ "dvec", "zvec", "rvec", "cvec", "dmat", "zmat", "rmat", "cmat" };

void usage()
{
  printf("Usage: %s input.dat [options]\n",PROG);
  printf("  -o fname           The matrix is saved in the file 'fname'.\n");
  printf("  -O                 The matrix is saved in the file whose name is automatically determined.\n");
  printf("  -d dir             Set the directory of output file as 'dir'.\n");
  printf("  -f num             Set output format as 'f' with num digits.\n");
  printf("  -e num             Set output format as 'e' with num digits.\n");
  printf("  -nv                No verbose mode.\n");
  printf("  -v                 Verbose mode.\n");
  printf("  -help              This message.\n");
  printf("\n");
  printf("Examples:\n");
  printf("     # %s input.dat\n",PROG);
  printf("\n");
  exit(0);
}

int main(int argc, char *argv[])
{
  strings *path=NULL;
  char *in_fname=NULL,*out_fname=NULL,*out_dir=NULL;
  char format[]="f";
  int m=0,n=0,debug=0,i;
  int loaded=0,autoname=0,type=NOTYPE,digits=2;
  double *dA=NULL,*dB=NULL;
  dcomplex *zA=NULL,*zB=NULL;
  rmulti **rA=NULL,**rB=NULL;
  cmulti **cA=NULL,**cB=NULL;

  // init
  in_fname=char_new("",NULL);
  out_fname=char_new("",NULL);
  out_dir=char_new("",NULL);

  // options
  i=1;
  while(i<argc){
    if     (STR_EQ(argv[i],"-help"))              { usage(); }
    else if(STR_EQ(argv[i],"-o") && i+1<argc)     { ++i; out_fname=char_renew(out_fname,argv[i],NULL); }
    else if(STR_EQ(argv[i],"-O"))                 { autoname=1; }
    else if(STR_EQ(argv[i],"-d") && i+1<argc)     { ++i; out_dir=char_renew(out_dir,argv[i],NULL); if(out_dir[strlen(out_dir)-1]!='/'){ out_dir=char_cat(out_dir,"/"); } }
    else if(STR_EQ(argv[i],"-f") && i+1<argc)     { strcpy(format,"f"); ++i; digits=atoi(argv[i]); }
    else if(STR_EQ(argv[i],"-e") && i+1<argc)     { strcpy(format,"e"); ++i; digits=atoi(argv[i]); }
    else if(STR_EQ(argv[i],"-nv"))                { debug=0; }
    else if(STR_EQ(argv[i],"-v"))                 { debug=1; }
    else if(STR_EQ(argv[i],"-vv"))                { debug=2; }
    else if(STR_EQ_N(argv[i],"-",1))              { usage(); }
    else if(strlen(in_fname)<=0)                  { in_fname=char_renew(in_fname,argv[i],NULL); }
    else                                          { usage(); }
    i++;
  }

  // check
  if(strlen(in_fname)<=0){ usage(); exit(0); }

  // autoname
  if(autoname){
    path=strings_split_path(in_fname,NULL);
    if(strlen(out_dir)<=0 && strlen(path->str[0])>0){
      out_dir=char_renew(out_dir,path->str[0],NULL);
      out_dir=char_cat(out_dir,"/");
    }
    out_fname=char_renew(out_fname,path->str[1],NULL);
    if(!STR_EQ(path->str[2],BIN_SUFFIX) && strlen(path->str[2])>0){
      out_fname=char_cat(out_fname,".");
      out_fname=char_cat(out_fname,path->str[2]);
    }
    out_fname=char_cat(out_fname,BASENAME_APPEND);
    out_fname=char_cat(out_fname,".");
    out_fname=char_cat(out_fname,BIN_SUFFIX);
    path=strings_del(path);
  }

  // load
  if(!loaded){ dA=dmat_bin_load(&m,&n,"%s",in_fname); } if(dA!=NULL){ type=DMAT; loaded=1; }
  if(!loaded){ zA=zmat_bin_load(&m,&n,"%s",in_fname); } if(zA!=NULL){ type=ZMAT; loaded=1; }
  if(!loaded){ rA=rmat_bin_load(&m,&n,"%s",in_fname); } if(rA!=NULL){ type=RMAT; loaded=1; }
  if(!loaded){ cA=cmat_bin_load(&m,&n,"%s",in_fname); } if(cA!=NULL){ type=CMAT; loaded=1; }
  if(!loaded){ printf("Error! The input file '%s' can NOT be opened, or is NOT supported type.\n\n",in_fname); usage(); exit(0); }

  // transpose
  if(dA!=NULL){ dB=dmat_allocate(n,m); dmat_copy_t(m,n,dB,n,dA,m); }
  if(zA!=NULL){ zB=zmat_allocate(n,m); zmat_copy_t(m,n,zB,n,zA,m); }
  if(rA!=NULL){ rB=rmat_allocate(n,m); rmat_clone_t(m,n,rB,n,rA,m); }
  if(cA!=NULL){ cB=cmat_allocate(n,m); cmat_clone_t(m,n,cB,n,cA,m); }

  // show
  if(strlen(out_fname)<=0 || debug>0){
    printf("Input_File: %s\n",in_fname);
    printf("Input_Type: %s\n",type_name[type]);
    if(DMAT<=type && type<=CMAT){ printf("Input_Size: %d-by-%d\n",m,n); }
  }
  if(strlen(out_fname)<=0 || debug>1){
    printf("Input_Matrix:\n");
    if(dA!=NULL){ dmat_print(m,n,dA,m,NULL,format,digits); }
    if(zA!=NULL){ zmat_print(m,n,zA,m,NULL,format,digits); }
    if(rA!=NULL){ rmat_print(m,n,rA,m,NULL,format,digits); }
    if(cA!=NULL){ cmat_print(m,n,cA,m,NULL,format,digits); }
  }
  if(strlen(out_fname)<=0 || debug>0){
    printf("Output_Type: %s\n",type_name[type]);
    if(DMAT<=type && type<=CMAT){ printf("Output_Size: %d-by-%d\n",n,m); }
  }
  if(strlen(out_fname)<=0 || debug>1){
    printf("Output_Matrix:\n");
    if(dB!=NULL){ dmat_print(n,m,dB,n,NULL,format,digits); }
    if(zB!=NULL){ zmat_print(n,m,zB,n,NULL,format,digits); }
    if(rB!=NULL){ rmat_print(n,m,rB,n,NULL,format,digits); }
    if(cB!=NULL){ cmat_print(n,m,cB,n,NULL,format,digits); }
  }

  // save
  if(strlen(out_fname)>0){
    if(debug>0){ printf("Output_File: %s%s\n",out_dir,out_fname); }
    if(dA!=NULL){ dmat_bin_save(n,m,dB,n,"%s%s",out_dir,out_fname); }
    if(zA!=NULL){ zmat_bin_save(n,m,zB,n,"%s%s",out_dir,out_fname); }
    if(rA!=NULL){ rmat_bin_save(n,m,rB,n,"%s%s",out_dir,out_fname); }
    if(cA!=NULL){ cmat_bin_save(n,m,cB,n,"%s%s",out_dir,out_fname); }
  }

  // done
  in_fname=char_del(in_fname);
  out_fname=char_del(out_fname);
  out_dir=char_del(out_dir);
  dA=dmat_free(dA);
  zA=zmat_free(zA);
  rA=rmat_free(m,n,rA);
  cA=cmat_free(m,n,cA);
  dB=dmat_free(dB);
  zB=zmat_free(zB);
  rB=rmat_free(n,m,rB);
  cB=cmat_free(n,m,cB);
  return 0;
}
