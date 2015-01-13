// 日本語
#include<isys.h>

enum { MAT_RAND=0, MAT_T3, MAT_T4, MAT_CAUCHY };

int main(int argc, char *argv[])
{
  int i,mat=MAT_RAND,n=10,LDA,prec=53,prec2=512,info,e;
  cmulti **A0=NULL,**A=NULL,**b=NULL,**x=NULL,**y=NULL,**r=NULL;
  rmulti *E=NULL;
  double tp3[]={1.0,-1.0,2.0},tp4[]={1.0,-1.0,1.0,1.5};

  // options
  for(i=0; i<argc; i++){
    if(STR_EQ(argv[i],"-n") && i+1<argc) { n=atoi(argv[++i]); }
    else if(STR_EQ(argv[i],"-p") && i+1<argc) { prec2=atoi(argv[++i]); }
    else if(STR_EQ(argv[i],"-r"))  { mat=MAT_RAND; }
    else if(STR_EQ(argv[i],"-t3")) { mat=MAT_T3; }
    else if(STR_EQ(argv[i],"-t4")) { mat=MAT_T4; }
    else if(STR_EQ(argv[i],"-c"))  { mat=MAT_CAUCHY; }
  }

  // set default precision
  set_default_prec(prec);
  printf("default precision is %d bits.\n",get_default_prec());
  set_auto_prec_disabled();

  // allocation
  LDA=n; A=cmat_allocate(LDA,n);
  LDA=n; A0=cmat_allocate(LDA,n);
  b=cvec_allocate(n);
  x=cvec_allocate(n);
  y=cvec_allocate(n);
  r=cvec_allocate(n);
  E=rallocate();

  printf("================\n");

  // initialize
  if(mat==MAT_RAND)        { init_genrand(0); cmat_set_rand(n,n,A0,LDA,2,-1); }
  else if(mat==MAT_T3)     { cmat_toeplitz(n,n,A0,LDA,3,tp3,1); }
  else if(mat==MAT_T4)     { cmat_toeplitz(n,n,A0,LDA,4,tp4,1); }
  else if(mat==MAT_CAUCHY) { cmat_cauchy(n,n,A0,LDA); }
  cmat_print(n,n,A0,LDA,"A=","f",2);
  cvec_set_rand(n,b,2,-1);
  cvec_print_prec(n,b,"b=","f",20);

  printf("================\n");
  printf("[デフォルト精度で計算]\n");
  cmat_copy(n,n,A,LDA,A0,LDA);  cvec_copy(n,x,b);  e=csolve(n,1,x,n,A,LDA,&info); cvec_copy(n,y,x);
  if(info==0){
    printf("(R=%d) ",e);  cvec_print_prec(n,x,"x=","e",20);
    printf("[無誤差で残差を計算]\n");
    set_auto_prec_enabled();
    e=csolve_residual(n,1,r,n,A0,LDA,x,n,b,n);
    set_auto_prec_disabled();
    printf("(R=%d) ",e);
    cvec_print_prec(n,r,"r=b-A0*x=","e",2);
    cvec_max_abs(E,n,r);
    printf("||r||=%.2e\n",rget_d(E));
  }else{
    printf("rank(A)=%d\n",info);
  }

  printf("================\n");
  printf("[精度%dbitsで計算]\n",prec2);
  set_default_prec(prec2); cmat_round(n,n,A,LDA,prec2); cvec_round(n,x,prec2);
  cmat_copy(n,n,A,LDA,A0,LDA);  cvec_copy(n,x,b);  e=csolve(n,1,x,n,A,LDA,&info);
  set_default_prec(prec);
  if(info==0){
    printf("(R=%d) ",e);  cvec_print_prec(n,x,"x=","e",20);
    printf("[無誤差で残差を計算]\n");
    set_auto_prec_enabled();
    e=csolve_residual(n,1,r,n,A0,LDA,x,n,b,n);
    set_auto_prec_disabled();
    printf("(R=%d) ",e);
    cvec_print_prec(n,r,"r=b-A0*x=","e",2);
    cvec_max_abs(E,n,r);
    printf("||r||=%.2e\n",rget_d(E));

    printf("================\n");
    printf("[精度%dbitsで計算した結果を丸める]\n",prec2);
    cvec_round(n,r,prec);
    e=cvec_round(n,x,prec);
    printf("(R=%d) ",e);  cvec_print_prec(n,x,"x=","e",20);
    printf("[無誤差で残差を計算]\n");
    set_auto_prec_enabled();
    e=csolve_residual(n,1,r,n,A0,LDA,x,n,b,n);
    set_auto_prec_disabled();
    printf("(R=%d) ",e);
    cvec_print_prec(n,r,"r=b-A0*x=","e",2);
    cvec_max_abs(E,n,r);
    printf("||r||=%.2e\n",rget_d(E));

    printf("[デフォルト精度の誤差を評価]\n");
    set_auto_prec_enabled();
    cvec_round(n,r,prec);
    e=cvec_sub(n,r,x,y);
    e+=cvec_div_2exp(n,r,r,cvec_get_exp_max(n,x),cvec_get_exp_max(n,x));
    printf("(R=%d) ",e);
    cvec_print_prec(n,r,"r=(x-y)/log2(x)=","e",2);
    cvec_max_abs(E,n,r);
    printf("||r||=%.2e\n",rget_d(E));
    cvec_print_exp(n,r,"r=(x-y)/log2(x)=");
  }else{
    printf("rank(A)=%d\n",info);
  }


  // free
  E=rfree(E);
  A=cmat_free(LDA,n,A);
  A0=cmat_free(LDA,n,A0);
  b=cvec_free(n,b);
  x=cvec_free(n,x);
  y=cvec_free(n,y);
  r=cvec_free(n,r);

  // done
  return 0;
}
