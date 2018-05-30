#include<isys.h>

//////////////////////////////////////////////

#define MATTYPE "dmat"
#define MATNAME_DEFAULT "toeplitz"
#define DATA_SUFFIX ".dat"
#define PROG "matgen"
#define MAX_LENGTH 100000

double ENTRY_DEFAULT[]={1.0, 2.0, 0.0, 1.5};
double ENTRY_DEFAULT_SIZE=4;
#define ENTRY_DEFAULT_OFFSET 1

void usage()
{
  printf("Usage: %s [options]\n",PROG);
  printf("  -o fname           The matrix is saved in the file 'fname'.\n");
  printf("  -O                 The matrix is saved in the file whose name is automatically determined.\n");
  printf("  -d dir             The matrix is saved in the directory 'dir'.\n");
  printf("  -n num             Set the size of matrix as num.\n");
  printf("  -toeplitz          Generate Toeplitz matrix.\n");
  printf("  -frank             Generate Frank matrix.\n");
  printf("  -ipjfact           Generate ipjfact matrix.\n");
  printf("  -param num         Set the parameter of matrix as num.\n");
  printf("  -entry n0,n1,..,nN Set the entries of Toeplitz matrix as n0,n1,..,nN.\n");
  printf("  -offset num        Set the offset for the entries of Toeplitz as num.\n");
  printf("  -nv                No verbose mode.\n");
  printf("  -v                 Verbose mode.\n");
  printf("  -vv                More verbose mode.\n");
  printf("  -help              This message.\n");
  printf("\n");
  printf("Examples:\n");
  printf("     # %s -n 10 -entry 1,2,0,1.5 -offset 1 -o toeplitz4.dat\n",PROG);
  printf("     # %s -n 10 -entry 1,2,0,1.5 -offset 1 -O\n",PROG);
  printf("\n");
  exit(0);
}

char *set_filename(char *matname, int m, int n, int offset, int k, double *entry, int param)
{
  int i;
  char *fname=NULL;
  if(char_eq(matname,"toeplitz")){
    fname=char_renew_sprintf(NULL,NULL,"%s{%s|offset=%d|entry=",MATTYPE,matname,offset);
    for(i=0; i<k; i++){
      if(i!=0){ fname=char_renew_sprintf(fname,NULL,"%s,%g",fname,entry[i]); }
      else    { fname=char_renew_sprintf(fname,NULL,"%s%g",fname,entry[i]); }
    }
    fname=char_renew_sprintf(fname,NULL,"%s|size=%04dx%04d}%s",fname,m,n,DATA_SUFFIX);
  }else if(char_eq(matname,"frank")){
    fname=char_renew_sprintf(NULL,NULL,"%s{%s|param=%1d|size=%04dx%04d}%s",MATTYPE,matname,param,m,n,DATA_SUFFIX);
  }else if(char_eq(matname,"ipjfact")){
    fname=char_renew_sprintf(NULL,NULL,"%s{%s|param=%1d|size=%04dx%04d}%s",MATTYPE,matname,param,m,n,DATA_SUFFIX);
  }else{ printf("Error in set_filename()\n"); exit(0); }
  return fname;
}

void begin(char *matname, int m, int n, int k, double *entry, int offset, int param)
{
  int i;
  if(char_eq(matname,"toeplitz")){
    printf("Generating the %d-by-%d Toeplitz matrix with offset=%d, entry={ ",m,n,offset);
    if(entry!=NULL){
      for(i=0; i<k; i++){
	if(i!=0){ printf(", "); }
	printf("%g",entry[i]);
      }
    }
    printf(" }...\n");
  }
  else if(char_eq(matname,"frank"))  { printf("Generating the %d-by-%d Frank matrix with pram=%d...\n",m,n,param); }
  else if(char_eq(matname,"ipjfact")){ printf("Generating the %d-by-%d ipjfact matrix with param=%d...\n",m,n,param); }
  else{ printf("Error in begin()\n"); exit(0); }
}

int main(int argc, char *argv[])
{
  int i,m=0,n=0,k=0,debug=0,offset=ENTRY_DEFAULT_OFFSET,param=0;
  int autoname=0;
  double *A=NULL;
  double *entry=NULL;
  char dir[MAX_LENGTH+1];
  char *fname=NULL,*matname=NULL;

  // init
  fname=char_new("",NULL);
  matname=char_new(MATNAME_DEFAULT,NULL);
  strcpy(dir,"");
  entry=dvec_allocate(ENTRY_DEFAULT_SIZE);
  dvec_set_dvec(ENTRY_DEFAULT_SIZE,entry,ENTRY_DEFAULT);

  // options
  i=1;
  while(i<argc){
    if     (STR_EQ(argv[i],"-help"))              { usage(); }
    else if(STR_EQ(argv[i],"-nv"))                { debug=0; }
    else if(STR_EQ(argv[i],"-v"))                 { debug=1; }
    else if(STR_EQ(argv[i],"-vv"))                { debug=2; }
    else if(STR_EQ(argv[i],"-n") && i+1<argc)     { ++i; n=atoi(argv[i]); m=n; }
    else if(STR_EQ(argv[i],"-size") && i+2<argc)  { ++i; m=atoi(argv[i]); ++i; n=atoi(argv[i]); }
    else if(STR_EQ(argv[i],"-o") && i+1<argc)     { ++i; fname=char_renew(fname,argv[i],NULL); }
    else if(STR_EQ(argv[i],"-O"))                 { autoname=1; }
    else if(STR_EQ(argv[i],"-d") && i+1<argc)     { ++i; strcpy(dir,argv[i]); if(dir[strlen(dir)-1]!='/'){ strcat(dir,"/"); } }
    else if(STR_EQ(argv[i],"-toeplitz"))          { matname=char_renew(matname,"toeplitz",NULL); }
    else if(STR_EQ(argv[i],"-frank"))             { matname=char_renew(matname,"frank",NULL); }
    else if(STR_EQ(argv[i],"-ipjfact"))           { matname=char_renew(matname,"ipjfact",NULL); }
    //    else if(STR_EQ(argv[i],"-entry") && i+1<argc) { ++i; entry=dvec_free(entry); entry=dvec_allocate_s(argv[i],&k); }
    else if(STR_EQ(argv[i],"-offset") && i+1<argc){ ++i; offset=atoi(argv[i]); }
    else if(STR_EQ(argv[i],"-param") && i+1<argc) { ++i; param=atoi(argv[i]); }
    else                                          { usage(); }
    i++;
  }

  // check
  if(k<=0 || m<=0 || n<=0 || entry==NULL){ begin(matname,m,n,k,entry,offset,param); print_red(); printf("Error!\n\n"); print_reset(); usage(); }
  if(debug>0){ begin(matname,m,n,k,entry,offset,param); }

  // autoname
  if(autoname){ fname=set_filename(matname,m,n,offset,k,entry,param); }

  // generator
  A=dmat_allocate(m,n);
  if     (char_eq(matname,"toeplitz")){ dmat_toeplitz(m,n,A,m,k,entry,offset); }
  else if(char_eq(matname,"frank"))   { dmat_frank(m,n,A,m,param); }
  else if(char_eq(matname,"ipjfact")) { dmat_ipjfact(m,n,A,m,param); }
  else{ printf("Error!\n"); exit(0); }

  // show
  if(debug>1 || strlen(fname)<=0){ dmat_print(m,n,A,m,"A=",'f',2); }

  // save
  if(strlen(fname)>0){
    dmat_bin_save(m,n,A,m,"%s%s",dir,fname);
    if(debug>0){ printf("Saved in %s%s\n",dir,fname); }
  }

  // done
  A=dmat_free(A);
  entry=dvec_free(entry);
  fname=char_del(fname);
  matname=char_del(matname);
  return 0;
}

//EOF
