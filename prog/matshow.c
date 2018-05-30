#include<isys.h>

#define PROG "matshow"
#define MAX_LENGTH 10000
#define BIN_SUFFIX "dat"
#define TXT_SUFFIX "txt"
#define CRN(C,S) { C=char_renew(C,S,NULL); }
#define CCT(C,S) { C=char_cat(C,S); }

enum { NOTYPE=-1, TEXT=0, DVEC, ZVEC, RVEC, CVEC, DMAT, ZMAT, RMAT, CMAT };
char *type_name[]={ "text", "dvec", "zvec", "rvec", "cvec", "dmat", "zmat", "rmat", "cmat" };

void usage()
{
  printf("Usage: %s input.dat [options]\n",PROG);
  printf("  -f num             Set print format as 'f' with num digits.\n");
  printf("  -e num             Set print format as 'e' with num digits.\n");
  printf("  -o-t               Set save mode as text.\n");
  printf("  -cplane            Set save mode as text and complex plane.\n");
  printf("  -o-d               Set save mode as double.\n");
  printf("  -o-z               Set save mode as double complex.\n");
  printf("  -o-r               Set save mode as multiple real.\n");
  printf("  -o-c               Set save mode as multiple complex.\n");
  printf("  -o fname           The matrix is saved in the file 'fname'.\n");
  printf("  -O                 The matrix is saved in the file whose name is automatically determined.\n");
  printf("  -d dir             Set the directory of output file as 'dir'.\n");
  printf("  -prec num          Set the precision used in output mode.\n");
  printf("  -nv                Set no verbose mode.\n");
  printf("  -v                 Set verbose mode, level 1.\n");
  printf("  -vv                Set verbose mode, level 2.\n");
  printf("  -vvv               Set verbose mode, leval 3.\n");
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
  char *fname_in=NULL,*fname_out=NULL,*dir_out=NULL;
  char format[]="f",mode='t',buf[MAX_LENGTH];
  int m=0,n=0,debug=2,i,prec=53;
  int loaded=0,autoname=0,type_in=NOTYPE,type_out=TEXT,digits=7;
  int cplane=0;
  double *dx=NULL,*dA=NULL;
  dcomplex *zx=NULL,*zA=NULL;
  rmulti **rx=NULL,**rA=NULL;
  cmulti **cx=NULL,**cA=NULL;

  // init
  fname_in=char_new("",NULL);
  fname_out=char_new("",NULL);
  dir_out=char_new("",NULL);

  // options
  i=1;
  while(i<argc){
    if     (STR_EQ(argv[i],"-help"))              { usage(); }
    else if(STR_EQ(argv[i],"-nv"))                { debug=0; }
    else if(STR_EQ(argv[i],"-v"))                 { debug=1; }
    else if(STR_EQ(argv[i],"-vv"))                { debug=2; }
    else if(STR_EQ(argv[i],"-vvv"))               { debug=3; }
    else if(STR_EQ(argv[i],"-v0"))                { debug=-1; }
    else if(STR_EQ(argv[i],"-o") && i+1<argc)     { ++i; fname_out=char_renew(fname_out,argv[i],NULL); }
    else if(STR_EQ(argv[i],"-o-t") || STR_EQ(argv[i],"-o-txt"))     { mode='t'; }
    else if(STR_EQ(argv[i],"-o-d") || STR_EQ(argv[i],"-o-double"))  { mode='d'; }
    else if(STR_EQ(argv[i],"-o-z") || STR_EQ(argv[i],"-o-dcomplex")){ mode='z'; }
    else if(STR_EQ(argv[i],"-o-r") || STR_EQ(argv[i],"-o-real"))    { mode='r'; }
    else if(STR_EQ(argv[i],"-o-c") || STR_EQ(argv[i],"-o-complex")) { mode='c'; }
    else if(STR_EQ(argv[i],"-O"))                 { autoname=1; }
    else if(STR_EQ(argv[i],"-d") && i+1<argc)     { ++i; dir_out=char_renew(dir_out,argv[i],NULL); if(dir_out[strlen(dir_out)-1]!='/'){ dir_out=char_cat(dir_out,"/"); } }
    else if(STR_EQ(argv[i],"-f") && i+1<argc)     { strcpy(format,"f"); ++i; digits=atoi(argv[i]); }
    else if(STR_EQ(argv[i],"-e") && i+1<argc)     { strcpy(format,"e"); ++i; digits=atoi(argv[i]); }
    else if(STR_EQ(argv[i],"-prec") && i+1<argc)  { strcpy(format,"e"); ++i; prec=atoi(argv[i]); }
    else if(STR_EQ(argv[i],"-cplane"))            { cplane=1; }
    else if(STR_EQ_N(argv[i],"-",1))              { usage(); }
    else if(strlen(fname_in)<=0)                  { fname_in=char_renew(fname_in,argv[i],NULL); }
    else                                          { usage(); }
    i++;
  }

  // load input file
  if(strlen(fname_in)<=0){ usage(); }
  if(!loaded){ dx=dvec_bin_load(&n,"%s",fname_in); } if(dx!=NULL){ type_in=DVEC; loaded=1; }
  if(!loaded){ zx=zvec_bin_load(&n,"%s",fname_in); } if(zx!=NULL){ type_in=ZVEC; loaded=1; }
  if(!loaded){ rx=rvec_bin_load(&n,"%s",fname_in); } if(rx!=NULL){ type_in=RVEC; loaded=1; }
  if(!loaded){ cx=cvec_bin_load(&n,"%s",fname_in); } if(cx!=NULL){ type_in=CVEC; loaded=1; }
  if(!loaded){ dA=dmat_bin_load(&m,&n,"%s",fname_in); } if(dA!=NULL){ type_in=DMAT; loaded=1; }
  if(!loaded){ zA=zmat_bin_load(&m,&n,"%s",fname_in); } if(zA!=NULL){ type_in=ZMAT; loaded=1; }
  if(!loaded){ rA=rmat_bin_load(&m,&n,"%s",fname_in); } if(rA!=NULL){ type_in=RMAT; loaded=1; }
  if(!loaded){ cA=cmat_bin_load(&m,&n,"%s",fname_in); } if(cA!=NULL){ type_in=CMAT; loaded=1; }
  if(!loaded){ printf("Error! The input file '%s%s' can NOT be opened, or is NOT supported type.\n\n",dir_out,fname_in); usage(); }
  if(debug>=1){
    printf("Input_File: %s\n",fname_in);
    printf("Input_Type: %s\n",type_name[type_in]);
    if(DVEC<=type_in && type_in<=CVEC){ printf("Input_Size: %d\n",n); }
    if(DMAT<=type_in && type_in<=CMAT){ printf("Input_Size: %d-by-%d\n",m,n); }
    if(type_in==DVEC || type_in==ZVEC || type_in==DMAT || type_in==ZMAT){ printf("Input_Precision: double\n"); }
    else if(type_in==RVEC){ printf("Input_Precision_Max: %d\n",rvec_get_prec_max(n,rx)); }
    else if(type_in==CVEC){ printf("Input_Precision_Max: %d\n",cvec_get_prec_max(n,cx)); }
    else if(type_in==RMAT){ printf("Input_Precision_Max: %d\n",rmat_get_prec_max(m,n,rA,m)); }
    else if(type_in==CMAT){ printf("Input_Precision_Max: %d\n",cmat_get_prec_max(m,n,cA,m)); }
  }
  if(debug>=2){
    if(DVEC<=type_in && type_in<=CVEC){ printf("Input_Vector:\n"); }
    if(DMAT<=type_in && type_in<=CMAT){ printf("Input_Matrix:\n"); }
    if(dx!=NULL){ dvec_print(n,dx,NULL,format[0],digits); }
    if(zx!=NULL){ zvec_print(n,zx,NULL,format[0],digits); }
    if(rx!=NULL){ rvec_print(n,rx,NULL,format[0],digits); }
    if(cx!=NULL){ cvec_print(n,cx,NULL,format[0],digits); }
    if(dA!=NULL){ dmat_print(m,n,dA,m,NULL,format[0],digits); }
    if(zA!=NULL){ zmat_print(m,n,zA,m,NULL,format[0],digits); }
    if(rA!=NULL){ rmat_print(m,n,rA,m,NULL,format[0],digits); }
    if(cA!=NULL){ cmat_print(m,n,cA,m,NULL,format[0],digits); }
  }
  if(debug==-1){
    if(dx!=NULL){ dvec_print(n,dx,NULL,format[0],digits); }
    if(zx!=NULL){ zvec_print(n,zx,NULL,format[0],digits); }
    if(rx!=NULL){ rvec_print(n,rx,NULL,format[0],digits); }
    if(cx!=NULL){ cvec_print(n,cx,NULL,format[0],digits); }
    if(dA!=NULL){ dmat_print(m,n,dA,m,NULL,format[0],digits); }
    if(zA!=NULL){ zmat_print(m,n,zA,m,NULL,format[0],digits); }
    if(rA!=NULL){ rmat_print(m,n,rA,m,NULL,format[0],digits); }
    if(cA!=NULL){ cmat_print(m,n,cA,m,NULL,format[0],digits); }
  }

  // output type
  if     (mode=='t')                                  { type_out=TEXT; }
  else if(mode=='d' && DVEC<=type_in && type_in<=CVEC){ type_out=DVEC; }
  else if(mode=='d' && DMAT<=type_in && type_in<=CMAT){ type_out=DMAT; }
  else if(mode=='z' && DVEC<=type_in && type_in<=CVEC){ type_out=ZVEC; }
  else if(mode=='z' && DMAT<=type_in && type_in<=CMAT){ type_out=ZMAT; }
  else if(mode=='r' && DVEC<=type_in && type_in<=CVEC){ type_out=RVEC; }
  else if(mode=='r' && DMAT<=type_in && type_in<=CMAT){ type_out=RMAT; }
  else if(mode=='c' && DVEC<=type_in && type_in<=CVEC){ type_out=CVEC; }
  else if(mode=='c' && DMAT<=type_in && type_in<=CMAT){ type_out=CMAT; }
  else                                                { usage(); }  
  if(debug>=1){
    printf("Output_Type: %s\n",type_name[type_out]);
    printf("Output_Precision: %d\n",prec);
  }

  // generate output filename
  if(autoname){
    path=strings_split_path(fname_in,BIN_SUFFIX);
    if(strlen(dir_out)<=0 && strlen(path->str[0])>0){ CRN(dir_out,path->str[0]); CCT(dir_out,"/"); }
    CRN(fname_out,"");
    if(type_out==TEXT){
      sprintf(buf,"%s_%s.%s",path->str[1],type_name[type_out],TXT_SUFFIX); CCT(fname_out,buf);
    }else if(DVEC==type_out || ZVEC==type_out || DMAT==type_out || ZMAT==type_out){
      sprintf(buf,"%s_%s.%s",path->str[1],type_name[type_out],BIN_SUFFIX); CCT(fname_out,buf);
    }else if(RVEC==type_out || CVEC==type_out || RMAT==type_out || CMAT==type_out){
      sprintf(buf,"%s_%s{%dbits}.%s",path->str[1],type_name[type_out],prec,BIN_SUFFIX); CCT(fname_out,buf);
    }
    path=strings_del(path);
  }

  // save
  if(debug>=1){
    if(strlen(fname_out)>0){ printf("Output_File: %s%s\n",dir_out,fname_out);  }
    else                   { printf("Output_File: NULL\n"); }
  }
  if(strlen(fname_out)>0){
    if(type_out==TEXT){
      if(dx!=NULL){ dvec_save(n,dx,0,"%s%s",dir_out,fname_out); }
      if(zx!=NULL){ zvec_save(n,zx,0,"%s%s",dir_out,fname_out); }
      if(rx!=NULL){ rvec_save(n,rx,0,(int)(prec*0.6),"%s%s",dir_out,fname_out); }
      if(cx!=NULL && !cplane){ cvec_save(n,cx,0,(int)(prec*0.6),"%s%s",dir_out,fname_out); }
      if(cx!=NULL &&  cplane){ cvec_save_cplane(n,cx,(int)(prec*0.6),"%s%s",dir_out,fname_out); }
      if(dA!=NULL){ dmat_save(m,n,dA,m,"%s%s",dir_out,fname_out); }
      if(zA!=NULL){ zmat_save(m,n,zA,m,"%s%s",dir_out,fname_out); }
      if(rA!=NULL){ rmat_save(m,n,rA,m,(int)(prec*0.6),"%s%s",dir_out,fname_out); }
      if(cA!=NULL){ cmat_save(m,n,cA,m,(int)(prec*0.6),"%s%s",dir_out,fname_out); }
    }
    else if(type_out==DMAT && type_in==DMAT){                                                              dmat_bin_save(m,n,dA,m,"%s%s",dir_out,fname_out); }
    else if(type_out==DMAT && type_in==ZMAT){ dA=dmat_allocate(m,n);           zmat_real(m,n,dA,m,zA,m);   dmat_bin_save(m,n,dA,m,"%s%s",dir_out,fname_out); dA=dmat_free(dA); }
    else if(type_out==DMAT && type_in==RMAT){ dA=dmat_allocate(m,n);           rmat_get_d(m,n,dA,m,rA,m);  dmat_bin_save(m,n,dA,m,"%s%s",dir_out,fname_out); dA=dmat_free(dA); }
    else if(type_out==DMAT && type_in==CMAT){ dA=dmat_allocate(m,n);           cmat_get_d(m,n,dA,m,cA,m);  dmat_bin_save(m,n,dA,m,"%s%s",dir_out,fname_out); dA=dmat_free(dA); }
    else if(type_out==ZMAT && type_in==DMAT){ zA=zmat_allocate(m,n);           zmat_copy_d(m,n,zA,m,dA,m); zmat_bin_save(m,n,zA,m,"%s%s",dir_out,fname_out); zA=zmat_free(zA); }
    else if(type_out==ZMAT && type_in==ZMAT){                                                              zmat_bin_save(m,n,zA,m,"%s%s",dir_out,fname_out); }
    else if(type_out==ZMAT && type_in==RMAT){ zA=zmat_allocate(m,n);           rmat_get_z(m,n,zA,m,rA,m);  zmat_bin_save(m,n,zA,m,"%s%s",dir_out,fname_out); zA=zmat_free(zA); }
    else if(type_out==ZMAT && type_in==CMAT){ zA=zmat_allocate(m,n);           cmat_get_z(m,n,zA,m,cA,m);  zmat_bin_save(m,n,zA,m,"%s%s",dir_out,fname_out); zA=zmat_free(zA); }
    else if(type_out==RMAT && type_in==DMAT){ rA=rmat_allocate_prec(m,n,prec); rmat_set_d(m,n,rA,m,dA,m);  rmat_bin_save(m,n,rA,m,"%s%s",dir_out,fname_out); rA=rmat_free(m,n,rA); }
    else if(type_out==RMAT && type_in==ZMAT){ rA=rmat_allocate_prec(m,n,prec); rmat_set_z(m,n,rA,m,zA,m);  rmat_bin_save(m,n,rA,m,"%s%s",dir_out,fname_out); rA=rmat_free(m,n,rA); }
    else if(type_out==RMAT && type_in==RMAT){                                  rmat_round(m,n,rA,m,prec);  rmat_bin_save(m,n,rA,m,"%s%s",dir_out,fname_out); }
    else if(type_out==RMAT && type_in==CMAT){ rA=rmat_allocate_prec(m,n,prec); cmat_real(m,n,rA,m,cA,m);   rmat_bin_save(m,n,rA,m,"%s%s",dir_out,fname_out); rA=rmat_free(m,n,rA); }
    else if(type_out==CMAT && type_in==DMAT){ cA=cmat_allocate_prec(m,n,prec); cmat_set_d(m,n,cA,m,dA,m);  cmat_bin_save(m,n,cA,m,"%s%s",dir_out,fname_out); cA=cmat_free(m,n,cA); }
    else if(type_out==CMAT && type_in==ZMAT){ cA=cmat_allocate_prec(m,n,prec); cmat_set_z(m,n,cA,m,zA,m);  cmat_bin_save(m,n,cA,m,"%s%s",dir_out,fname_out); cA=cmat_free(m,n,cA); }
    else if(type_out==CMAT && type_in==RMAT){ cA=cmat_allocate_prec(m,n,prec); cmat_copy_rmat(m,n,cA,m,rA,m); cmat_bin_save(m,n,cA,m,"%s%s",dir_out,fname_out); cA=cmat_free(m,n,cA); }
    else if(type_out==CMAT && type_in==CMAT){                                  cmat_round(m,n,cA,m,prec);  cmat_bin_save(m,n,cA,m,"%s%s",dir_out,fname_out); }
    else if(type_out==DVEC && type_in==DVEC){                                                              dvec_bin_save(n,dx,"%s%s",dir_out,fname_out); }
    else if(type_out==DVEC && type_in==ZVEC){ dx=dvec_allocate(n);             dvec_real_zvec(n,dx,zx);         dvec_bin_save(n,dx,"%s%s",dir_out,fname_out); dx=dvec_free(dx); }
    else if(type_out==DVEC && type_in==RVEC){ dx=dvec_allocate(n);             rvec_get_dvec(n,dx,rx);        dvec_bin_save(n,dx,"%s%s",dir_out,fname_out); dx=dvec_free(dx); }
    else if(type_out==DVEC && type_in==CVEC){ dx=dvec_allocate(n);             cvec_get_dvec(n,dx,cx);        dvec_bin_save(n,dx,"%s%s",dir_out,fname_out); dx=dvec_free(dx); }
    else if(type_out==ZVEC && type_in==DVEC){ zx=zvec_allocate(n);             zvec_set_dvec(n,zx,dx);       zvec_bin_save(n,zx,"%s%s",dir_out,fname_out); zx=zvec_free(zx); }
    else if(type_out==ZVEC && type_in==ZVEC){                                                              zvec_bin_save(n,zx,"%s%s",dir_out,fname_out); }
    else if(type_out==ZVEC && type_in==RVEC){ zx=zvec_allocate(n);             rvec_get_zvec(n,zx,rx);        zvec_bin_save(n,zx,"%s%s",dir_out,fname_out); zx=zvec_free(zx); }
    else if(type_out==ZVEC && type_in==CVEC){ zx=zvec_allocate(n);             cvec_get_zvec(n,zx,cx);        zvec_bin_save(n,zx,"%s%s",dir_out,fname_out); zx=zvec_free(zx); }
    else if(type_out==RVEC && type_in==DVEC){ rx=rvec_allocate_prec(n,prec);   rvec_set_dvec(n,rx,dx);        rvec_bin_save(n,rx,"%s%s",dir_out,fname_out); rx=rvec_free(n,rx); }
    else if(type_out==RVEC && type_in==ZVEC){ rx=rvec_allocate_prec(n,prec);   rvec_set_zvec(n,rx,zx);        rvec_bin_save(n,rx,"%s%s",dir_out,fname_out); rx=rvec_free(n,rx); }
    else if(type_out==RVEC && type_in==RVEC){                                  rvec_round_rvec(n,rx,prec);      rvec_bin_save(n,rx,"%s%s",dir_out,fname_out); }
    else if(type_out==RVEC && type_in==CVEC){ rx=rvec_allocate_prec(n,prec);   rvec_real_cvec(n,rx,cx);         rvec_bin_save(n,rx,"%s%s",dir_out,fname_out); rx=rvec_free(n,rx); }
    else if(type_out==CVEC && type_in==DVEC){ cx=cvec_allocate_prec(n,prec);   cvec_set_dvec(n,cx,dx);        cvec_bin_save(n,cx,"%s%s",dir_out,fname_out); cx=cvec_free(n,cx); }
    else if(type_out==CVEC && type_in==ZVEC){ cx=cvec_allocate_prec(n,prec);   cvec_set_zvec(n,cx,zx);        cvec_bin_save(n,cx,"%s%s",dir_out,fname_out); cx=cvec_free(n,cx); }
    else if(type_out==CVEC && type_in==RVEC){ cx=cvec_allocate_prec(n,prec);   cvec_set_rvec(n,cx,rx);       cvec_bin_save(n,cx,"%s%s",dir_out,fname_out); cx=cvec_free(n,cx); }
    else if(type_out==CVEC && type_in==CVEC){                                  cvec_round_cvec(n,cx,prec);      cvec_bin_save(n,cx,"%s%s",dir_out,fname_out); }
    else{ ERROR_AT; exit(0); } 
  }

  // done
  fname_in=char_del(fname_in);
  fname_out=char_del(fname_out);
  dir_out=char_del(dir_out);
  dx=dvec_free(dx);
  zx=zvec_free(zx);
  rx=rvec_free(n,rx);
  cx=cvec_free(n,cx);
  dA=dmat_free(dA);
  zA=zmat_free(zA);
  rA=rmat_free(m,n,rA);
  cA=cmat_free(m,n,cA);
  return 0;
}
