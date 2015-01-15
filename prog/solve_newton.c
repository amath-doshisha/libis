#include<stdio.h>
#include<stdlib.h>
#include<isys.h>

#define PROG_NAME "solve_newton"
#define DEFAULT_EQS "{x-sin(3/2*(x+y)), y-cos(3/2*(x+y))}"

#define RMA(A,M,N,P){ A=rmat_allocate_prec(M,N,P); }
#define RMF(A,M,N) { A=rmat_free(M,N,A); }
#define RVA(X,N,P) { X=rvec_allocate_prec(N,P); }
#define RVF(X,N)   { X=rvec_free(N,X); }
#define RA(X,P)    { X=rallocate_prec(P); }
#define RF(X)      { X=rfree(X); }

void usage()
{
  printf("\n");
  printf(" Usage: %s [options]\n",PROG_NAME);
  printf("\n");
  printf("   -eqs str    Set the equations as str.\n");
  printf("   -f num      Set the print format as 'f' with num ditits.\n");
  printf("   -e num      Set the print format as 'e' with num ditits.\n");
  printf("   -prec num   Set the precision in computing as num.\n");
  printf("   -kprec num  Set the precision in error bound computing as num.\n");
  printf("   -n num      Set the max of steps as num.\n");
  printf("   -o-x fname  Save the solution to fname.\n");
  printf("   -o-e fname  Save the error bound to fname.\n");
  printf("   -nv         Silent mode.\n");
  printf("   -v          Verbose mode, leve=1.\n");
  printf("   -vv         Verbose mode, leve=2.\n");
  printf("\n");
  exit(0);
}

int main(int argc, const char *argv[])
{
  int debug=1,i,m,step_max=-1,prec=53,kprec=-1,digits=15,info;
  double x0[]={1,0.5};
  double x1[]={1,1};
  func_t *fF=NULL;
  rmulti **x=NULL,**e=NULL;
  char *eqs=NULL,*out_x=NULL,*out_e=NULL,format[]="f";

  // init
  out_x=char_new(""," \t\n");
  out_e=char_new(""," \t\n");
  eqs=char_new(DEFAULT_EQS," \t\n");

  // option
  i=1;
  while(i<argc){
    if     (STR_EQ(argv[i],"--help"))            { usage(); }
    else if(STR_EQ(argv[i],"-nv"))               { debug=0; }
    else if(STR_EQ(argv[i],"-v"))                { debug=1; }
    else if(STR_EQ(argv[i],"-vv"))               { debug=2; }
    else if(i+1<argc && STR_EQ(argv[i],"-n"))    { step_max=atoi(argv[++i]); }
    else if(i+1<argc && STR_EQ(argv[i],"-f"))    { digits=atoi(argv[++i]); format[0]='f'; }
    else if(i+1<argc && STR_EQ(argv[i],"-e"))    { digits=atoi(argv[++i]); format[0]='e'; }
    else if(i+1<argc && STR_EQ(argv[i],"-prec")) { prec=atoi(argv[++i]); }
    else if(i+1<argc && STR_EQ(argv[i],"-kprec")){ kprec=atoi(argv[++i]); }
    else if(i+1<argc && STR_EQ(argv[i],"-eqs"))  { eqs=char_renew(eqs,argv[++i]," \t\n"); }
    else if(i+1<argc && STR_EQ(argv[i],"-o-x"))  { out_x=char_renew(out_x,argv[++i]," \t\n"); }
    else if(i+1<argc && STR_EQ(argv[i],"-o-e"))  { out_e=char_renew(out_e,argv[++i]," \t\n"); }
    else{ usage(); }
    i++;
  }

  // init
  if(kprec<=0){ if(prec<64){ kprec=128; }else{ kprec=2*prec; } }
  func_eval(func_script("begin(x,y,z,w,v,u)"));
  fF=func_eval(func_script(eqs));
  m=func_asize(fF);
  RVA(x,m,prec);
  RVA(e,m,kprec);
  rvec_set_nan(m,x);
  rvec_set_d(m,x,x1);
  rvec_set_d(m,x,x0);

  // message
  if(debug>=1){
    printf("Equations: "); func_print(fF); printf("\n");
    printf("Eequsiots_Dimension: %d\n",m);
    printf("Newton_Precision: %d\n",prec);
    if(step_max<0){
      printf("Newton_Steps_Max: infinite_loop\n");
    }else{
      printf("Newton_Steps_Max: %d\n",step_max);
    }
    rvec_print(m,x,"Initial_Vector:",format,digits);
    printf("Error_Estimate_Precision: %d\n",kprec);
  }

  // solve
  rsolve_newton(m,x,fF,step_max,debug-1);

  // message
  if(debug>=1 || strlen(out_x)<=0){
    rvec_print(m,x,"Solution:",format,digits);
  }
  
  // error
  info=rsolve_krawczyk(m,e,x,fF,debug-1);

  // message
  if(debug>=1){
    printf("Error_Estimate_Result: ");
    if(info){ print_red(); printf("failed"); }  else{ print_green(); printf("succeeded"); }
    print_reset();
    printf("\n");
    rvec_print(m,e,"Error_Bound:","e",3);
  }

  // save
  if(debug>=1){
    printf("Output_Solution: %s\n",(strlen(out_x)<=0?"null":out_x));
    printf("Output_Error_Bound: %s\n",(strlen(out_e)<=0?"null":out_e));
  }
  if(strlen(out_x)>0){
    rvec_bin_save(m,x,out_x);
    printf("Solution is saved as '%s'\n",out_x);
  }
  if(strlen(out_e)>0){
    rvec_bin_save(m,e,out_e);
    printf("Error bounds is saved as '%s'\n",out_e);
  }


  // done
  eqs=char_del(eqs);
  fF=func_del(fF);
  func_clear();
  RVF(x,m); RVF(e,m);
  return 0;
}

//EOF
