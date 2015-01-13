#include<isys.h>
#include <time.h>

#define MATNAME "dhToda"
#define TXT_SUFFIX  "txt"
#define BIN_SUFFIX "dat"
#define PROG "matgen_dhToda"
#define N 1000  

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
  printf("  -p num             Set the parameter=num.\n");
  printf("  -p-e num           Set the parameter=2^(num).\n");
  printf("  -mode-a            Set the mode=arithmetic.\n");
  printf("  -mode-g            Set the mode=geometric.\n");
  printf("  -mode-i            Set the mode=integers.\n");
  printf("  -help              This message.\n");
  printf("\n");
  printf("Examples:\n");
  printf("     # %s -size 10 -M 2 -O\n",PROG);
  printf("\n");
  exit(0);
}

void rmat_gen_dhToda_set_lambda_arith(int m, double *lambda, double parameter)
{
  int i,n;
  n=m+1;
  for(i=1;i<n;i++){
    //    lambda[i-1]=1-(1-parameter)*(double)(i-1)/(m);
    lambda[i-1]=1-(1-parameter)*(double)(i-1)/(m-1);
  }
}

void rmat_gen_dhToda_set_lambda_geome(int m,  double *lambda, double parameter)
{
  int i,n;
  n=m+1;
  for(i=1;i<n;i++){
    //    lambda[i-1]=pow(parameter,(i-1)/(double)(m));
    lambda[i-1]=pow(parameter,(i-1)/(double)(m-1));
  }
}

void rmat_gen_dhToda_set_lambda_arith2(int m, double *lambda, double parameter)
{
  int i,n,k;
  n=m+1;
  for(i=1;i<n;i++){
    k=i%2;
    //    lambda[i-1]=1-(1-parameter)*(double)(i-1)/(m);
    lambda[i-1]=1-(1-parameter)*(double)(i-1)/(m-1);
    if(k==0){
      lambda[i-1]=-lambda[i-1];
    }
  }
}

void rmat_gen_dhToda_set_lambda_geome2(int m,  double *lambda, double parameter)
{
  int i,n,k;
  n=m+1;
  for(i=1;i<n;i++){
    k=i%2;
    //    lambda[i-1]=pow(parameter,(i-1)/(double)(m));
    lambda[i-1]=pow(parameter,(i-1)/(double)(m-1));
    if(k==0){
      lambda[i-1]=-lambda[i-1];
    }
  }
}

void rmat_gen_dhToda(int m, rmulti **A, int LDA, int M, rmulti **lambda, rmulti **c, rmulti *E_max, int debug,char mode)
{
  int prec=0,k=0,n=0,l=0,f_size=0,n_size=0,*Q_size=NULL,*E_size=NULL,LDR=0,i=0;
  rmulti *a=NULL,**f=NULL,***Q=NULL,***E=NULL,**sigma=NULL,**R=NULL;
  // init
  n_size=(M+1)*(m-1)+2*M;
  // precision
  prec=rmat_get_prec_max(m,m,A,LDA);
  // allocate
  sigma=rvec_allocate_prec(m,prec);
  a=rallocate_prec(prec);
  LDR=m; R=rmat_allocate_prec(LDR,m,prec);
  Q_size=ivec_allocate(m);
  Q=malloc(m*sizeof(rmulti**));
  for(k=0; k<m; k++){
    Q_size[k]=n_size-k*(M+1);
    Q[k]=rvec_allocate_prec(Q_size[k],prec);
  }
  E_size=ivec_allocate(m);
  E=malloc(m*sizeof(rmulti**));
  for(k=0; k<m; k++){
    E_size[k]=n_size-(k+1)*(M+1)+1;
    E[k]=rvec_allocate_prec(E_size[k],prec);
  }
  f_size=Q_size[0]+1; f=rvec_allocate_prec(f_size,prec);
  // set sigma
  rinv_d(a,M);
  if(mode=='A' || mode=='G'){
    rvec_pow_r(m,sigma,lambda,a);
    if(M%2==1){
      for(i=0;i<m;i++){
	if(i%2==1){
	  rneg(sigma[i],sigma[i]);
	}
      }
    }
  }else{
    rvec_pow_r(m,sigma,lambda,a);
  }
  // set f
  for(n=0; n<f_size; n++){
    rset_d(f[n],0);
    for(k=0; k<m; k++){
      rpow_si(a,sigma[k],n); // a=(sigma[i])^n
      radd_mul(f[n],c[k],a); // f[i]=f[i]+c[i]*(sigma[i])^n
    }
  }
  // set Q[0]
  for(n=0; n<Q_size[0]; n++){
    if(n+1<f_size){
      rdiv(Q[0][n],f[n+1],f[n]); // Q[n][1]=f[n+1]/f[n]
    }else{ rset_nan(Q[0][n]); }
  }
  // set E[0]
  k=0;
  for(n=0; n<E_size[k]; n++){
    if(n+M<Q_size[k] && n<Q_size[k]){
      rsub(E[k][n],Q[k][n+M],Q[k][n]);       //    E[k][n]=Q[k][n+M]-Q[k][n];
    }else{ rset_nan(E[k][n]); }
  }  
  // set QE-table
  for(k=1; k<m; k++){
    for(n=0; n<Q_size[k]; n++){
      rdiv(a,E[k-1][n+1],E[k-1][n]);
      rmul(Q[k][n],a,Q[k-1][n+M]);     // Q[k][n]=(E[k-1][n+1]*Q[k-1][n+M])/E[k-1][n];
    }    
    for(n=0; n<E_size[k]; n++){
      //radd(a,Q[k][n+M],E[k-1][n+1]);
      //rsub(E[k][n],a,Q[k][n]);
      rsub(a,Q[k][n+M],Q[k][n]);
      radd(E[k][n],a,E[k-1][n+1]);     // E[k][n]=Q[k][n+M]-Q[k][n]+E[k-1][n+1];
    }
  }
  rvec_max_abs(E_max,E_size[m-1],E[m-1]);
  // debug
  if(debug>0){
    printf("Q=\n");
    for(n=0; n<n_size; n++){ for(k=0; k<m; k++){ if(n<Q_size[k]){ mpfr_printf("%.3Re ",Q[k][n]); } } printf("\n"); }
    printf("E=\n");
    for(n=0; n<n_size; n++){ for(k=0; k<m; k++){ if(n<E_size[k]){ mpfr_printf("%.3Re ",E[k][n]); } } printf("\n"); }
  }
  // set A=L
  for(n=0; n<m; n++){
    for(k=0; k<m; k++){
      if     (n==k)    { rset_one(MAT(A,n,k,LDA)); }      // A[n][k]=1;
      else if((n-k)==1){ rcopy(MAT(A,n,k,LDA),E[k][0]); } // A[n][k]=E[k][0];
      else             { rset_zero(MAT(A,n,k,LDA)); }     // A[n][k]=0;
    }
  }
  // set A=A*R
  for(l=M-1;l>=0;l--){
    for(n=0;n<m;n++){
      for(k=0;k<m;k++){
	if     ((k-n)==1){ rset_one(MAT(R,n,k,LDR)); }      // R[n][k]=1;
	else if(n==k)    { rcopy(MAT(R,n,k,LDR),Q[k][l]); } // R[n][k]=Q[k][l];
	else             { rset_zero(MAT(R,n,k,LDR)); }     // R[n][k]=0;
      }
    }
    rmat_prod(m,m,m,A,LDA,A,LDA,R,LDR); // A=A*R
  }
  // done
  a=rfree(a);
  f=rvec_free(f_size,f);
  sigma=rvec_free(m,sigma);
  R=rmat_free(LDR,m,R);
  for(k=0; k<m; k++){ Q[k]=rvec_free(Q_size[k],Q[k]); } free(Q); Q=NULL;
  for(k=0; k<m; k++){ E[k]=rvec_free(E_size[k],E[k]); } free(E); E=NULL;
  Q_size=ivec_free(Q_size);
  E_size=ivec_free(E_size);
  return;
}

int main(int argc, char *argv[])
{
  int i,LDA,debug=2,M=1,prec=1024,autoname=0,m=0,digits=15;
  rmulti **A=NULL,*E_max=NULL, **rLambda=NULL,**c=NULL;
  char buf[10000],*fname=NULL,*fname_lambda=NULL,*dir=NULL,format[]="f";
  char mode='i',*c_mode=NULL;
  double *dLambda=NULL,p=0.5;
  clock_t t1,t2;


  t1 = clock();
  // init
  fname=char_new("",NULL);
  fname_lambda=char_new("",NULL);
  dir=char_new("./",NULL);
  c_mode=char_new("ones",NULL);

  // options
  i=1;
  while(i<argc){
    if     (STR_EQ(argv[i],"-help"))              { usage(); }
    else if(STR_EQ(argv[i],"-nv"))                { debug=0; }
    else if(STR_EQ(argv[i],"-v"))                 { debug=1; }
    else if(STR_EQ(argv[i],"-vv"))                { debug=2; }
    else if(STR_EQ(argv[i],"-vvv"))               { debug=3; }
    else if(STR_EQ(argv[i],"-vvvv"))              { debug=4; }
    else if(STR_EQ(argv[i],"-prec") && i+1<argc)  { ++i; prec=atoi(argv[i]); }
    else if(STR_EQ(argv[i],"-n") && i+1<argc)     { ++i; m=atoi(argv[i]); }
    else if(STR_EQ(argv[i],"-m") && i+1<argc)     { ++i; m=atoi(argv[i]); }
    else if(STR_EQ(argv[i],"-m") && i+1<argc)     { ++i; m=atoi(argv[i]); }
    else if(STR_EQ(argv[i],"-size") && i+1<argc)  { ++i; m=atoi(argv[i]); }
    else if(STR_EQ(argv[i],"-p") && i+1<argc)     { ++i; p=atof(argv[i]); }
    else if(STR_EQ(argv[i],"-p-e") && i+1<argc)   { ++i; p=pow(2,atoi(argv[i])); }
    else if(STR_EQ(argv[i],"-e") && i+1<argc)     { ++i; format[0]='e'; digits=atoi(argv[i]); }
    else if(STR_EQ(argv[i],"-f") && i+1<argc)     { ++i; format[0]='f'; digits=atoi(argv[i]); }
    else if(STR_EQ(argv[i],"-o") && i+1<argc)     { ++i; strcpy(fname,argv[i]); }
    else if(STR_EQ(argv[i],"-O"))                 { autoname=1; }
    else if(STR_EQ(argv[i],"-d") && i+1<argc)     { dir=char_renew(dir,argv[++i],NULL); if(dir[strlen(dir)-1]!='/'){ dir=char_cat(dir,"/"); } }
    else if(STR_EQ(argv[i],"-M") && i+1<argc)     { ++i; M=atoi(argv[i]); }
    else if(STR_EQ(argv[i],"-mode-a"))            { mode='a'; }
    else if(STR_EQ(argv[i],"-mode-g"))            { mode='g'; }
    else if(STR_EQ(argv[i],"-mode-A"))           { mode='A'; }
    else if(STR_EQ(argv[i],"-mode-G"))           { mode='G'; }
    else if(STR_EQ(argv[i],"-mode-i"))            { mode='i'; }
    else if(STR_EQ(argv[i],"-mode-I"))            { mode='I'; }
    else if(STR_EQ(argv[i],"-c-ones"))            { c_mode=char_renew(c_mode,"ones",NULL); }
    else if(STR_EQ(argv[i],"-c-grid"))            { c_mode=char_renew(c_mode,"grid",NULL); }
    else                                          { usage(); }
    i++;
  }

  // check
  if(prec<26 || m<2 || M<=0){ usage(); }

  // precision
  if(debug>=1){
    printf("MPFR_Precision: %d\n",prec);
    printf("Eigenvalues_Mode: ");
    if     (mode=='a'){ printf("arithmetic interval with [1,%g)\n",p); }
    else if(mode=='g'){ printf("geometric interval with [1,%g)\n",p); }
    else if(mode=='A'){ printf("arithmetic interval with [1,%g)\n",p); }
    else if(mode=='G'){ printf("geometric interval with [1,%g)\n",p); }
    else if(mode=='i'){ printf("integers with [1,%d]\n",m); }
    else              { printf("NULL\n"); }
  }
  //  set_default_prec(prec);

  // autoname
  if(autoname){
    if(mode=='i'){
      sprintf(buf,"rmat{%dbits|%s|size=%04d|M=%04d|c=%s|mode=%c}.%s",prec,MATNAME,m,M,c_mode,mode,BIN_SUFFIX);        fname=char_renew(fname,buf,NULL);
      sprintf(buf,"dvec{%dbits|%s|size=%04d|M=%04d|c=%s|mode=%c}_lambda.%s",prec,MATNAME,m,M,c_mode,mode,BIN_SUFFIX); fname_lambda=char_renew(fname_lambda,buf,NULL);
    }else{
      sprintf(buf,"rmat{%dbits|%s|size=%04d|M=%04d|c=%s|mode=%c|p=%g}.%s",prec,MATNAME,m,M,c_mode,mode,p,BIN_SUFFIX);        fname=char_renew(fname,buf,NULL);
      sprintf(buf,"dvec{%dbits|%s|size=%04d|M=%04d|c=%s|mode=%c|p=%g}_lambda.%s",prec,MATNAME,m,M,c_mode,mode,p,BIN_SUFFIX); fname_lambda=char_renew(fname_lambda,buf,NULL);
    }
  }

  // allocate
  LDA=m;
  A=rmat_allocate_prec(LDA,m,prec);
  rLambda=rvec_allocate_prec(m,prec);
  c=rvec_allocate_prec(m,prec);
  dLambda=dvec_allocate(m);
  E_max=rallocate_prec(prec);

  // init lambda
  if     (mode=='a'){ rmat_gen_dhToda_set_lambda_arith(m,dLambda,p); }
  else if(mode=='g'){ rmat_gen_dhToda_set_lambda_geome(m,dLambda,p); }
  else if(mode=='A'){ rmat_gen_dhToda_set_lambda_arith2(m,dLambda,p); }
  else if(mode=='G'){ rmat_gen_dhToda_set_lambda_geome2(m,dLambda,p); }
  else if(mode=='i'){ dvec_set_grid(m,dLambda); dvec_add_d(m,dLambda,dLambda,1); }
  else if(mode=='I'){ dvec_set_grid(m,dLambda); dvec_add_d(m,dLambda,dLambda,1); dvec_pow_si(m,dLambda,dLambda,M); }
  else              { dvec_set_grid(m,dLambda); dvec_add_d(m,dLambda,dLambda,1); }
  rvec_set_d(m,rLambda,dLambda);
  rvec_abs(m,rLambda,rLambda);

  // init c
  if(char_eq(c_mode,"grid")){
    rvec_set_grid(m,c); rvec_add_d(m,c,c,1);
  }else{
    rvec_set_ones(m,c);
  }

  // generate matrix
  rmat_gen_dhToda(m,A,LDA,M,rLambda,c,E_max,debug-2,mode);
  if(debug>=1){ mpfr_printf("Max_E_At_Bound: %.2Re\n",E_max); }

  // show
  if(debug>=2){
    dvec_print(m,dLambda,"lambda=",format,digits);
    rmat_print(m,m,A,m,"A=",format,digits);
  }

  // save
  if(debug>=1){
    if(strlen(fname)>0){ printf("Output_File_Matrix: %s%s\n",dir,fname); }
    else               { printf("Output_File_Matrix: NULL\n"); }
    if(strlen(fname_lambda)>0){ printf("Output_File_Lambda: %s%s\n",dir,fname_lambda); }
    else                      { printf("Output_File_Lambda: NULL\n"); }
  }
  if(strlen(fname)>0){
    rmat_bin_save(m,m,A,m,"%s%s",dir,fname);
    dvec_bin_save(m,dLambda,"%s%s",dir,fname_lambda);
  }

  // done
  fname=char_del(fname);
  fname_lambda=char_del(fname_lambda);
  dir=char_del(dir);
  c=rvec_free(m,c);
  E_max=rfree(E_max);
  A=rmat_free(LDA,m,A);
  rLambda=rvec_free(m,rLambda);
  dLambda=dvec_free(dLambda);



  t2 = clock();
  printf("実行時間:\t%f\n",(double)(t2 - t1) / CLOCKS_PER_SEC);
  
  return 0;
}

//EOF
