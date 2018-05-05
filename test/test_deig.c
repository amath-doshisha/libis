#include<isys.h>

int main(int argc, char *argv[])
{
  int i,n=10,LDA,LDX,debug=0,k;
  double *A=NULL,*X=NULL,*lambda=NULL,*E=NULL,*Ebin=NULL,tp[]={1.0, 2.0, 1.0, 0.0};

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
  LDA=n; A=dmat_allocate(LDA,n);
  LDX=n; X=dmat_allocate(LDX,n);
  lambda=dvec_allocate(n);
  E=dvec_allocate(n);
  Ebin=dvec_allocate(n);
  
  // initialize
  dmat_toeplitz(n,n,A,LDA,4,tp,1);
  dmat_set_eye(n,n,X,LDX);
  dvec_set_zeros(n,lambda);

  // print
  dmat_print(n,n,A,LDA,"A=",'f',0);

  // compute
  k=dhpeig(n,A,LDA,X,LDX,lambda,debug);
  deig_sort(n,k,lambda,X,LDX);

  // output  
  deig_residual_norm_max(k,E,A,LDA,X,LDX,lambda);
  dvec_log2_abs(k,Ebin,E);
  printf("Number of computed eigenpairs = %d\n",k);
  if(k<n){
    printf("Failed to compute.\n");
  }else{
    printf("Succeeded to compute.\n");
  }

  dmat_print(n,k,X,LDX,"X=",'f',1);
  dvec_print(k,lambda,"lambda=",'e',2);
  dvec_print(k,E,"E=",'e',0);
  dvec_print(k,Ebin,"log2(E)=",'f',1);
  printf("max(E)=%.0e\n",dvec_norm_max(k,E));
  printf("max(log2(E))=%.2g\n",dvec_max(k,Ebin));

  // free
  A=dmat_free(A);
  A=dmat_free(X);
  lambda=dvec_free(lambda);
  E=dvec_free(E);
  Ebin=dvec_free(Ebin);

  // done
  return 0;
}
