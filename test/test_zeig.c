#include<isys.h>

int main(int argc, char *argv[])
{
  int i,n=10,LDA,LDX,debug=0,k;
  double *E=NULL,*Ebin=NULL,tp[]={1.0, 2.0, 1.0, 0.0};
  dcomplex *A=NULL,*X=NULL,*lambda=NULL;

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
    else if(STR_EQ(argv[i],"-n") && i+1<argc) { n=atoi(argv[++i]); }
    else{ printf("Error!\n"); exit(0); }
    i++;
  }
  k=n;

  // allocation
  LDA=n; A=zmat_allocate(LDA,n);
  LDX=n; X=zmat_allocate(LDX,n);
  lambda=zvec_allocate(n);
  E=dvec_allocate(n);
  Ebin=dvec_allocate(n);

  // initialize
  zmat_toeplitz(n,n,A,LDA,4,tp,1);
  zmat_set_eye(n,n,X,LDX);
  zvec_set_zeros(n,lambda);

  // print
  zmat_print(n,n,A,LDA,"A=",'f',0);

  // compute
  k=zhpeig(n,A,LDA,X,LDX,lambda,debug);
  zeig_sort(n,k,lambda,X,LDX);

  // output  
  zeig_residual_norm_max(k,E,A,LDA,X,LDX,lambda);
  dvec_log2_abs(k,Ebin,E);
  printf("Number of computed eigenpairs = %d\n",k);
  if(k<n){
    printf("Failed to compute.\n");
  }else{
    printf("Succeeded to compute.\n");
  }

  zmat_print(n,k,X,LDX,"X=",'f',1);
  zvec_print(k,lambda,"lambda=",'f',2);
  dvec_print(k,E,"E=",'e',0);
  dvec_print(k,Ebin,"log2(E)=",'f',1);
  printf("max(E)=%.0e\n",dvec_norm_max(k,E));
  printf("max(log2(E))=%.2g\n",dvec_max(k,Ebin));

  // free
  A=zmat_free(A);
  X=zmat_free(X);
  lambda=zvec_free(lambda);
  E=dvec_free(E);
  Ebin=dvec_free(Ebin);

  // done
  return 0;
}
