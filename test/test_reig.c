#include<isys.h>

int main(int argc, char *argv[])
{
  int i,n=10,LDA,LDX,debug=0,k,prec=256;
  rmulti **A=NULL,**X=NULL,**lambda=NULL,**E=NULL,**Ebin=NULL,*a=NULL;
  double tp[]={1.0, 2.0, 1.0, 0.0};

  // options
  i=1;
  while(i<argc){
    if(STR_EQ(argv[i],"-d0")) debug=0;
    else if(STR_EQ(argv[i],"-d1")) debug=1;
    else if(STR_EQ(argv[i],"-d2")) debug=2;
    else if(STR_EQ(argv[i],"-d3")) debug=3;
    else if(STR_EQ(argv[i],"-sym")) { tp[0]=1; tp[1]=2; tp[2]=1; tp[3]=0; }
    else if(STR_EQ(argv[i],"-t3"))  { tp[0]=1; tp[1]=2; tp[2]=3; tp[3]=0; }
    else if(STR_EQ(argv[i],"-t4"))  { tp[0]=1; tp[1]=2; tp[2]=0; tp[3]=1.5; }
    else if(STR_EQ(argv[i],"-prec") && i+1<argc) { prec=atoi(argv[++i]); }
    else if(STR_EQ(argv[i],"-n") && i+1<argc) { n=atoi(argv[++i]); }
    else{ printf("Error!\n"); exit(0); }
    i++;
  }
  k=n;

  // set default precision
  set_default_prec(prec);
  printf("prec=%d\n",prec);
  
  // allocation
  LDA=n; A=rmat_allocate(LDA,n);
  LDX=n; X=rmat_allocate(LDX,n);
  lambda=rvec_allocate(n);
  E=rvec_allocate(n);
  Ebin=rvec_allocate(n);
  a=rallocate();
  
  // initialize
  rmat_toeplitz(n,n,A,LDA,4,tp,1);
  rmat_set_eye(n,n,X,LDX);
  rvec_set_zeros(n,lambda);

  // print
  rmat_print(n,n,A,LDA,"A=",'f',0);

  // compute
  k=rhpeig(n,X,LDX,lambda,A,LDA,debug);
  reig_sort(n,k,lambda,X,LDX);

  // output  
  reig_max_abs_residuals(k,E,A,LDA,X,LDX,lambda);
  rvec_log2_abs(k,Ebin,E);
  printf("Number of computed eigenpairs = %d\n",k);
  if(k<n){
    printf("Failed to compute.\n");
  }else{
    printf("Succeeded to compute.\n");
  }

  rmat_print(n,k,X,LDX,"X=",'f',1);
  rvec_print(k,lambda,"lambda=",'e',2);
  rvec_print(k,E,"E=",'e',0);
  rvec_print(k,Ebin,"log2(E)=",'f',1);
  rvec_max_abs(a,k,E); mpfr_printf("max(E)=%.0Re\n",a);
  rvec_max(a,k,Ebin);  mpfr_printf("max(log2(E))=%.2Rf\n",a);

  // free
  A=rmat_free(LDA,n,A);
  X=rmat_free(LDX,n,X);
  lambda=rvec_free(n,lambda);
  E=rvec_free(n,E);
  Ebin=rvec_free(n,Ebin);
  a=rfree(a);

  // done
  return 0;
}
