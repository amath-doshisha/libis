#include<stdio.h>
#include<stdlib.h>
#include<isys.h>

#define RMA(A,M,N,P){ A=rmat_allocate_prec(M,N,P); }
#define RMF(A,M,N) { A=rmat_free(M,N,A); }
#define RVA(X,N,P) { X=rvec_allocate_prec(N,P); }
#define RVF(X,N)   { X=rvec_free(N,X); }
#define RA(X,P)    { X=rallocate_prec(P); }
#define RF(X)      { X=rfree(X); }

void usage()
{
  printf("Usage:\n");
  exit(0);
}

int main(int argc, const char *argv[])
{
  int debug=0,i,m,step_max=-1,prec=53,kprec=128;
  double x0[]={1,0.5};
  double x1[]={1,1};
  func_t *fF=NULL;
  rmulti **x=NULL,**e=NULL;

  // option
  i=1;
  while(i<argc){
    if(STR_EQ(argv[i],"--help")){ usage(); }
    else if(STR_EQ(argv[i],"-v")){ debug=1; }
    else if(STR_EQ(argv[i],"-vv")){ debug=2; }
    else if(i+1<argc && STR_EQ(argv[i],"-n")){ step_max=atoi(argv[++i]); }
    else if(i+1<argc && STR_EQ(argv[i],"-prec")){ prec=atoi(argv[++i]); }
    else if(i+1<argc && STR_EQ(argv[i],"-kprec")){ kprec=atoi(argv[++i]); }
    else{ usage(); }
    i++;
  }

  // init
  printf("step_max=%d\n",step_max);
  func_eval(func_script("begin(x,y,z,w,v,u)"));
  fF=func_eval(func_script("{x-sin(3/2*(x+y)), y-cos(3/2*(x+y))}"));
  printf("fF="); func_print(fF); printf("\n");
  m=func_asize(fF);
  printf("m=%d\n",m);
  RVA(x,m,prec);
  RVA(e,m,kprec);
  rvec_set_nan(m,x);
  rvec_set_d(m,x,x1);
  rvec_set_d(m,x,x0);  
  rvec_print(m,x,"x0=","f",6);

  // solve
  rsolve_newton(m,x,fF,step_max,debug);
  rvec_print(m,x,"x=","f",20);
  
  // error
  printf("prec=%d\n",prec);
  if(rsolve_krawczyk(m,e,x,fF,debug)){ print_red(); printf("failed.\n"); }
  else{ print_green(); printf("succeeded.\n"); }
  print_reset();
  rvec_print(m,e,"e=","e",1);


  // done
  fF=func_del(fF);
  func_clear();
  RVF(x,m); RVF(e,m);
  printf("func_new_del_check_sum=%d\n",func_new_del_check_sum());  
  return 0;
}

//EOF
