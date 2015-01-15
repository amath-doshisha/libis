#include<isys.h>
#include <time.h>

#define MATNAME "dhToda"
#define MATTYPE_DEFAULT "TN"
#define TXT_SUFFIX  "txt"
#define BIN_SUFFIX "dat"
#define PROG "matgen_dhToda"
#define N 1000  
#define FILE_NAME_LENGTH_MAX 100000


void usage()
{
  printf("Usage: %s [options]\n",PROG);
  printf("  -o fname           The matrix is saved in the file 'fname'.\n");
  printf("  -O                 The matrix is saved in the file whose name is automatically determined.\n");
  printf("  -d dir             The matrix is saved in the directory 'dir'.\n");
  printf("  -n num             Set the size of matrix as num.\n");
  printf("  -m num             Set the size of matrix as num.\n");
  printf("  -size num          Set the size of matrix as num.\n");
  printf("  -f digits          Set the print format 'f' with digits.\n");
  printf("  -e digits          Set the print format 'e' with digits.\n");
  printf("  -nv                Set no verbose mode.\n");
  printf("  -v                 Set verbose mode, level 1.\n");
  printf("  -vv                Set verbose mode, level 2.\n");
  printf("  -vvv               Set verbose mode, level 3.\n");
  printf("  -vvvv              Set verbose mode, level 4.\n");
  printf("  -lambda-min num    Set the lambda_min as num.\n");
  printf("  -lambda-min-e num  Set the lambda_min as 2^num.\n");
  printf("  -lambda-a          Set the lambda as arithmetic interval.\n");
  printf("  -lambda-g          Set the lambda as geometric interval.\n");
  printf("  -lambda-i          Set the lambda as grid integers.\n");
  printf("  -lambda-I          Set the lambda as grid integers powered to M.\n");
  printf("  -lambda-file fname Set the lambda as input file.\n");
  printf("  -help              This message.\n");
  printf("\n");
  printf("Examples:\n");
  printf("     # %s -size 10 -M 2 -O\n",PROG);
  printf("\n");
  exit(0);
}


/**
 @brief rmulti型のベクトルの読み込み.
*/
void dvec_load(int n, double *x, char* fmt, ...)
{
  int i;
  char fname[FILE_NAME_LENGTH_MAX];
  va_list argp;
  FILE *fid;
  // file name
  va_start(argp,fmt);
  vsprintf(fname,fmt,argp);
  // open file
  if((fid=fopen(fname,"r"))==0){ ERROR_AT; printf("Cant't open file: %s\n",fname); exit(0); }
  // real
  for(i=0;i<n;i++){ if(fscanf(fid,"%lf",&x[i])==EOF){ ERROR_EXIT("Error in dvec_load(), fscanf(), fname=%s\n",fname); } }
  // close
  fclose(fid);
}


void set_set_lambda_arith(int m, double *lambda, double parameter)
{
  int i,n;
  n=m+1;
  for(i=1;i<n;i++){ lambda[i-1]=1-(1-parameter)*(double)(i-1)/(m-1); }
}

void set_set_lambda_geome(int m,  double *lambda, double parameter)
{
  int i,n;
  n=m+1;
  for(i=1;i<n;i++){ lambda[i-1]=pow(parameter,(i-1)/(double)(m-1)); }
}

void set_set_lambda_arith2(int m, double *lambda, double parameter)
{
  int i,n,k;
  n=m+1;
  for(i=1;i<n;i++){
    k=i%2;
    lambda[i-1]=1-(1-parameter)*(double)(i-1)/(m-1);
    if(k==0){ lambda[i-1]=-lambda[i-1]; }
  }
}

void set_set_lambda_geome2(int m,  double *lambda, double parameter)
{
  int i,n,k;
  n=m+1;
  for(i=1;i<n;i++){
    k=i%2;
    lambda[i-1]=pow(parameter,(i-1)/(double)(m-1));
    if(k==0){ lambda[i-1]=-lambda[i-1]; }
  }
}

int main(int argc, char *argv[])
{
  int i,LDA,LDQ,debug=1,M=1,prec=1024,autoname=0,m=0,digits=15,matgen=1;
  rmulti **A=NULL,**E=NULL,**Q=NULL, **rLambda=NULL,**c=NULL;
  char *fname=NULL,*fname_lambda=NULL,*fname_Q=NULL,*fname_E=NULL,*in_fname=NULL,*dir=NULL,format[]="f",*mattype=NULL,*l_mode=NULL,*c_mode=NULL;
  double *dLambda=NULL,lambda_min=0.5,time=-1;
  clock_t t1,t2;

  // init
  fname=char_new("",NULL);
  fname_lambda=char_new("",NULL);
  fname_E=char_new("",NULL);
  fname_Q=char_new("",NULL);
  in_fname=char_new("",NULL);
  dir=char_new("./",NULL);
  l_mode=char_new("i",NULL);
  c_mode=char_new("ones",NULL);
  mattype=char_new(MATTYPE_DEFAULT,NULL);

  // options
  i=1;
  while(i<argc){
    if     (STR_EQ(argv[i],"-help"))                    { usage(); }
    else if(STR_EQ(argv[i],"-nv"))                      { debug=0; }
    else if(STR_EQ(argv[i],"-v"))                       { debug=1; }
    else if(STR_EQ(argv[i],"-vv"))                      { debug=2; }
    else if(STR_EQ(argv[i],"-vvv"))                     { debug=3; }
    else if(STR_EQ(argv[i],"-vvvv"))                    { debug=4; }
    else if(STR_EQ(argv[i],"-prec") && i+1<argc)        { ++i; prec=atoi(argv[i]); }
    else if(STR_EQ(argv[i],"-n") && i+1<argc)           { ++i; m=atoi(argv[i]); }
    else if(STR_EQ(argv[i],"-m") && i+1<argc)           { ++i; m=atoi(argv[i]); }
    else if(STR_EQ(argv[i],"-M") && i+1<argc)           { ++i; M=atoi(argv[i]); }
    else if(STR_EQ(argv[i],"-size") && i+1<argc)        { ++i; m=atoi(argv[i]); }
    else if(STR_EQ(argv[i],"-e") && i+1<argc)           { ++i; format[0]='e'; digits=atoi(argv[i]); }
    else if(STR_EQ(argv[i],"-f") && i+1<argc)           { ++i; format[0]='f'; digits=atoi(argv[i]); }
    else if(STR_EQ(argv[i],"-o") && i+1<argc)           { ++i; fname=char_renew(fname,argv[i],NULL); }
    else if(STR_EQ(argv[i],"-O"))                       { autoname=1; }
    else if(STR_EQ(argv[i],"-d") && i+1<argc)           { ++i; dir=char_renew_sprintf(dir,NULL,"%s/",argv[i]); }
    else if(STR_EQ(argv[i],"-lambda-min") && i+1<argc)  { ++i; lambda_min=atof(argv[i]); }
    else if(STR_EQ(argv[i],"-lambda-min-e") && i+1<argc){ ++i; lambda_min=pow(2,atoi(argv[i])); }
    else if(STR_EQ(argv[i],"-lambda-a"))                { l_mode=char_renew(l_mode,"a",NULL); }
    else if(STR_EQ(argv[i],"-lambda-g"))                { l_mode=char_renew(l_mode,"g",NULL); }
    else if(STR_EQ(argv[i],"-lambda-i"))                { l_mode=char_renew(l_mode,"i",NULL); }
    else if(STR_EQ(argv[i],"-lambda-I"))                { l_mode=char_renew(l_mode,"I",NULL); }
    else if(STR_EQ(argv[i],"-lambda-file") && i+2<argc) { ++i; in_fname=char_renew(in_fname,argv[i],NULL); ++i; l_mode=char_renew(l_mode,argv[i]," \t\n"); }
    else if(STR_EQ(argv[i],"-c-ones"))                  { c_mode=char_renew(c_mode,"ones",NULL); }
    else if(STR_EQ(argv[i],"-c-grid"))                  { c_mode=char_renew(c_mode,"grid",NULL); }
    else if(STR_EQ(argv[i],"-A-yes"))                   { matgen=1; }
    else if(STR_EQ(argv[i],"-A-no"))                    { matgen=0; }
    else                                                { usage(); }
    i++;
  }

  // check
  if(prec<26 || m<2 || M<=0){ usage(); }

  // precision
  if(debug>=1){
    printf("MPFR_Precision: %d [bits]\n",prec);
    printf("Eigenvalues_Mode: ");
    if     (char_eq(l_mode,"a")){ printf("lambda=a; arithmetic interval with [1,%g)\n",lambda_min); }
    else if(char_eq(l_mode,"g")){ printf("lambda=a; geometric interval with [1,%g)\n",lambda_min); }
    else if(char_eq(l_mode,"i")){ printf("lambda=i; [1,2,..,%d]\n",m); }
    else if(char_eq(l_mode,"I")){ printf("lambda=I; [1^%d,2^%d,..,%d^%d]\n",M,M,M,m); }
    else if(strlen(in_fname)>0) { printf("lambda=%s; input=%s\n",l_mode,in_fname); }
    else                        { ERROR_AT; exit(0); }
    printf("Constants_Mode: ");
    if     (char_eq(c_mode,"ones")){ printf("c=oens; [1,1,..,1]\n"); } 
    else if(char_eq(c_mode,"grid")){ printf("c=grid; [1,2,..,%d]\n",m); }
    else                           { ERROR_AT; exit(0); }
  }

  // autoname
  if(autoname){
    if(char_eq(l_mode,"i") || strlen(in_fname)>0){
      fname=char_renew_sprintf(fname,NULL,"%srmat{%dbits|%s|%s|size=%04d|M=%04d|lambda=%s|c=%s}.%s",dir,prec,MATNAME,mattype,m,M,l_mode,c_mode,BIN_SUFFIX);
      fname_E=char_renew_sprintf(fname_E,NULL,"%srvec{%dbits|%s|%s|size=%04d|M=%04d|lambda=%s|c=%s}_E.%s",dir,prec,MATNAME,mattype,m,M,l_mode,c_mode,BIN_SUFFIX);
      fname_Q=char_renew_sprintf(fname_Q,NULL,"%srmat{%dbits|%s|%s|size=%04d|M=%04d|lambda=%s|c=%s}_Q.%s",dir,prec,MATNAME,mattype,m,M,l_mode,c_mode,BIN_SUFFIX);
      fname_lambda=char_renew_sprintf(fname_lambda,NULL,"%sdvec{%dbits|%s|%s|size=%04d|M=%04d|lambda=%s|c=%s}_lambda.%s",dir,prec,MATNAME,mattype,m,M,l_mode,c_mode,BIN_SUFFIX);
    }else{
      fname=char_renew_sprintf(fname,NULL,"%srmat{%dbits|%s|%s|size=%04d|M=%04d|lambda=%s|min=%g|c=%s}.%s",dir,prec,MATNAME,mattype,m,M,l_mode,lambda_min,c_mode,BIN_SUFFIX);
      fname_E=char_renew_sprintf(fname_E,NULL,"%srvec{%dbits|%s|%s|size=%04d|M=%04d|lambda=%s|min=%g|c=%s}_E.%s",dir,prec,MATNAME,mattype,m,M,l_mode,lambda_min,c_mode,BIN_SUFFIX);
      fname_Q=char_renew_sprintf(fname_Q,NULL,"%srmat{%dbits|%s|%s|size=%04d|M=%04d|lambda=%s|min=%g|c=%s}_Q.%s",dir,prec,MATNAME,mattype,m,M,l_mode,lambda_min,c_mode,BIN_SUFFIX);
      fname_lambda=char_renew_sprintf(fname_lambda,NULL,"%sdvec{%dbits|%s|%s|size=%04d|M=%04d|lambda=%s|min=%g|c=%s}_lambda.%s",dir,prec,MATNAME,mattype,m,M,l_mode,lambda_min,c_mode,BIN_SUFFIX);
    }
  }

  // allocate
  LDA=m;
  LDQ=m;
  A=rmat_allocate_prec(LDA,m,prec);
  rLambda=rvec_allocate_prec(m,prec);
  c=rvec_allocate_prec(m,prec);
  dLambda=dvec_allocate(m);
  E=rvec_allocate_prec(m,prec);
  Q=rmat_allocate_prec(LDQ,M,prec);

  //matrix mode
  if(matgen==1){ rmat_set_ones(m,m,A,LDA); }

  // init lambda
  if     (char_eq(l_mode,"a")){ set_set_lambda_arith(m,dLambda,lambda_min); }
  else if(char_eq(l_mode,"g")){ set_set_lambda_geome(m,dLambda,lambda_min); }
  else if(char_eq(l_mode,"i")){ dvec_set_grid(m,dLambda); dvec_add_d(m,dLambda,dLambda,1); }
  else if(char_eq(l_mode,"I")){ dvec_set_grid(m,dLambda); dvec_add_d(m,dLambda,dLambda,1); dvec_pow_si(m,dLambda,dLambda,M); }
  else if(strlen(in_fname)>0) { dvec_load(m,dLambda,in_fname); }
  else                        { ERROR_AT; exit(0); }
  rvec_set_d(m,rLambda,dLambda);

  // init c
  if(char_eq(c_mode,"grid")){ rvec_set_grid(m,c); rvec_add_d(m,c,c,1);  }
  else                      { rvec_set_ones(m,c); }

  // show
  if(debug>=1){ dvec_print(m,dLambda,"Specified_Eigenvalues:",format,digits); }

  // generate matrix
  t1=clock();
  riep_dhToda_TN(m,M,A,LDA,Q,LDQ,E,rLambda,c,debug-1);
  t2=clock();
  time=(double)(t2-t1)/CLOCKS_PER_SEC;

  // show
  if(debug>=1 || strlen(fname)<=0){ rmat_print(m,m,A,m,"Gerated_Matrix:",format,digits); }
  if(debug>=1){ printf("Computing_Time: %f [sec]\n",time); }

  // save
  if(strlen(fname_lambda)>0){ dvec_bin_save(m,dLambda,"%s",fname_lambda); printf("Output_File_Lambda: %s\n",fname_lambda); }
  if(strlen(fname)>0)       { rmat_bin_save(m,m,A,m,"%s",fname);          printf("Output_File_Matrix: %s\n",fname); }
  if(strlen(fname_E)>0)     { rvec_bin_save(m-1,E,"%s",fname_E);          printf("Output_File_E:      %s\n",fname_E); }
  if(strlen(fname_Q)>0)     { rmat_bin_save(m,M,Q,LDQ,"%s",fname_Q);      printf("Output_File_Q:      %s\n",fname_Q); }

  // done
  fname=char_del(fname);
  fname_lambda=char_del(fname_lambda);
  fname_Q=char_del(fname_Q);
  fname_E=char_del(fname_E);
  in_fname=char_del(in_fname);
  l_mode=char_del(l_mode);
  c_mode=char_del(c_mode);
  dir=char_del(dir);
  c=rvec_free(m,c);
  A=rmat_free(LDA,m,A);
  Q=rmat_free(LDQ,M,Q);
  E=rvec_free(m,E);
  rLambda=rvec_free(m,rLambda);
  dLambda=dvec_free(dLambda); 
  return 0;
}

//EOF
