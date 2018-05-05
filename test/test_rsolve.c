// 日本語
#include<isys.h>

enum { MAT_RAND=0, MAT_T3, MAT_T4, MAT_CAUCHY, MAT_SINGULAR };

int main(int argc, char *argv[])
{
  int i,mat=MAT_RAND,n=10,LDA,prec=53,prec2=512,info;
  rmulti **A0=NULL,**A=NULL,**b=NULL,**x=NULL,**y=NULL,**r=NULL,*E=NULL;
  char *str[]={"1","2","1","2",NULL};
  double tp3[]={1.0,-1.0,2.0},tp4[]={1.0,-1.0,1.0,1.5};

  // options
  for(i=0; i<argc; i++){
    if(STR_EQ(argv[i],"-n") && i+1<argc) { n=atoi(argv[++i]); }
    else if(STR_EQ(argv[i],"-p") && i+1<argc) { prec2=atoi(argv[++i]); }
    else if(STR_EQ(argv[i],"-r"))  { mat=MAT_RAND; }
    else if(STR_EQ(argv[i],"-t3")) { mat=MAT_T3; }
    else if(STR_EQ(argv[i],"-t4")) { mat=MAT_T4; }
    else if(STR_EQ(argv[i],"-c"))  { mat=MAT_CAUCHY; }
    else if(STR_EQ(argv[i],"-s"))  { mat=MAT_SINGULAR; n=2; }
  }

  
  // set default precision
  set_default_prec(prec);
  printf("default precision is %d bits.\n",get_default_prec());

  // allocation
  LDA=n; A=rmat_allocate(LDA,n);
  LDA=n; A0=rmat_allocate(LDA,n);
  b=rvec_allocate(n);
  x=rvec_allocate(n);
  y=rvec_allocate(n);
  r=rvec_allocate(n);
  E=rallocate();


  printf("================\n");

  // initialize
  if(mat==MAT_RAND)        { init_genrand(0); rmat_set_rand(n,n,A0,LDA,2,-1); }
  else if(mat==MAT_T3)     { rmat_toeplitz(n,n,A0,LDA,3,tp3,1); }
  else if(mat==MAT_T4)     { rmat_toeplitz(n,n,A0,LDA,4,tp4,1); }
  else if(mat==MAT_CAUCHY) { rmat_cauchy(n,n,A0,LDA); }
  else if(mat==MAT_SINGULAR) { rmat_set_s(n,n,A0,LDA,str,2); }
  rmat_print(n,n,A0,LDA,"A=",'f',2);
  rvec_set_rand(n,b,2,-1);
  rvec_print_prec(n,b,"b=","f",20);

  printf("================\n");
  printf("[デフォルト精度で計算]\n");
  rmat_copy(n,n,A,LDA,A0,LDA);  rvec_copy(n,x,b); rsolve(n,1,x,n,A,LDA,&info); rvec_copy(n,y,x);
  if(info==0){
    rvec_print_prec(n,x,"x=","e",20);
    printf("[残差を計算]\n");
    rsolve_residual(n,1,r,n,A0,LDA,x,n,b,n);
    rvec_print_prec(n,r,"r=b-A0*x=","e",2);
    rvec_max_abs(E,n,r);
    printf("||r||=%.2e\n",rget_d(E));
  }else{
    printf("rank(A)=%d\n",info);
  }


  printf("================\n");
  printf("[精度%dbitsで計算]\n",prec2);
  set_default_prec(prec2);
  rvec_round(n,x,prec2);
  rmat_copy(n,n,A,LDA,A0,LDA);
  rvec_copy(n,x,b);
  rsolve(n,1,x,n,A,LDA,&info);
  if(info==0){
    rvec_print_prec(n,x,"x=","e",20);
    printf("[残差を計算]\n");
    rvec_round(n,r,prec2);
    rsolve_residual(n,1,r,n,A0,LDA,x,n,b,n);
    rvec_print_prec(n,r,"r=b-A0*x=","e",2);
    rvec_max_abs(E,n,r);
    printf("||r||=%.2e\n",rget_d(E));
    printf("================\n");
    printf("[精度%dbitsで計算した結果を丸める]\n",prec2);
    rvec_round(n,r,prec);
    rvec_round(n,x,prec);
    rvec_print_prec(n,x,"x=","e",20);
    printf("[残差を計算]\n");
    rsolve_residual(n,1,r,n,A0,LDA,x,n,b,n);
    rvec_print_prec(n,r,"r=b-A0*x=","e",2);
    rvec_max_abs(E,n,r);
    printf("||r||=%.2e\n",rget_d(E));
    printf("[デフォルト精度の誤差を評価]\n");
    rvec_round(n,r,prec);
    rvec_sub(n,r,x,y);
    rvec_div_2exp(n,r,r,rvec_get_prec_max(n,x));
    rvec_print_prec(n,r,"r=(x-y)/log2(x)=","e",2);
    rvec_max_abs(E,n,r);
    printf("||r||=%.2e\n",rget_d(E));
    rvec_print_exp(n,r,"r=log2(x-y)-log2(x)=");

  }else{
    printf("rank(A)=%d\n",info);
  }
  rvec_save_log2_abs(n,x,-3,5,"foo.txt");

  // free
  E=rfree(E);
  A=rmat_free(LDA,n,A);
  A0=rmat_free(LDA,n,A0);
  b=rvec_free(n,b);
  x=rvec_free(n,x);
  y=rvec_free(n,x);
  r=rvec_free(n,r);

  // done
  return 0;
}
