#include<isys.h>

int main(int argc, const char *argv[])
{
  int i,b;
  rmulti *x=NULL;
  
  for(i=1; i<argc; i++){
    if(strcmp(argv[i],"--h")){ printf("Usage: test_rmulti\n"); exit(0); }
  }
  
  set_default_prec(128);
  x=rallocate();

  // nan
  rset_nan(x);
  b=rget_length(x);
  mpfr_printf("x=%0.40Rf\tlen=%d\t",x,b);
  rput(x);

  // zero
  rset_zero(x);
  b=rget_length(x);
  mpfr_printf("x=%0.40Rf\tlen=%d\t",x,b);
  rput(x);

  // 1
  rset_one(x);
  b=rget_length(x);
  mpfr_printf("x=%0.40Rf\tlen=%d\t",x,b);
  rput(x);

  // -1
  rset_one(x); rneg_r(x,x);
  b=rget_length(x);
  mpfr_printf("x=%0.40Rf\tlen=%d\t",x,b);
  rput(x);

  // 1/2
  rset_d(x,0.5);
  b=rget_length(x);
  mpfr_printf("x=%0.40Rf\tlen=%d\t",x,b);
  rput(x);
 
  // loop
  rset_d(x,2);
  b=rget_length(x);
  mpfr_printf("x=%0.40Rf\tlen=%d\t",x,b);
  rput(x);
  for(i=0; i<get_default_prec(); i++){
    rsub_rd(x,x,1); // x=x-1
    rdiv_rd(x,x,2); // x=x/2
    radd_rd(x,x,1);  // x=x+1
    b=rget_length(x);
    mpfr_printf("x=%0.40Rf\tlen=%d\t",x,b);
    rput(x);
  }

  x=rmfree(x);
  
  /*
  cmulti *x=NULL;
  set_default_prec(128);
  x=callocate();
  if(argc>=2){
    printf("str='%s'\n",argv[1]);
    cset_s(x,argv[1]);
  }else{
    cset_s(x,"1234");
  }
  mpfr_printf("x=%.10Rg%+.10Rgi\n",C_R(x),C_I(x));
  x=cfree(x);
  */
  
  /*
  rmulti *x=NULL;
  set_default_prec(128);
  x=rallocate();
  if(argc>=2){
    printf("str='%s'\n",argv[1]);
    rset_s(x,argv[1]);
  }else{
    rset_s(x,"1234");
  }
  mpfr_printf("x=%.40Rg\n",x);  
  x=rfree(x);
  */

  /*
  int n=6;
  func_t *f=NULL;
  rmulti *x=NULL,**rx=NULL,**ry=NULL,*rz=NULL;
  cmulti **cx=NULL,**cy=NULL,*cz=NULL;

  set_default_prec(53);

  x=rallocate();
  rset_d(x,1);
  rput_info(x,"x");

  rx=rvec_allocate(n);
  ry=rvec_allocate(2);
  rz=rallocate();
  cx=cvec_allocate(n);
  cy=cvec_allocate(2);
  cz=callocate();

  func_eval(func_script("begin(x0,x1,x2,x3,x4,x5)"));
  //f=func_script("x0+x1+x2+x3+x4+x5");
  //  f=func_script("x0*x1*x2*x3*x4*x5");
  //f=func_script("x0*x1^2*x2^3*x3^4*x4^5*x5^6");
  // f=func_script("x0^2*x1^3*x2^(-2)*x3^(-3)*x4*x5");
  //f=func_script("(((2/3)*x0^2*x1^3*x2^(-2)*x3^(-3)*x4*x5))");
  f=func_eval(func_script("atan(((2/3)*x0^2*x1^3*x2^(-2)*x3^(-3)*x4*x5))"));
  printf("f="); func_print(f); printf("\n");

  rvec_set_rand(n,rx,2,-1);
  rvec_print(n,rx,"rx=","e",20);
  irvec_func(ry[0],ry[1],f,n,rx,rx);
  rvec_func(rz,f,n,rx);
  rvec_print(2,ry,"ry=","e",20);
  mpfr_printf("rz=%.20Re\n",rz);

  cvec_set_rand(n,cx,2,-1);
  cvec_print(n,cx,"cx=","e",20);
  icvec_func(cy[0],cy[1],f,n,cx,cx);
  cvec_func(cz,f,n,cx);
  cvec_print(2,cy,"cy=","e",20);
  mpfr_printf("cz=%.20Re    %.20Re\n",C_R(cz),C_I(cz));


  rx=rvec_free(n,rx);
  ry=rvec_free(2,ry);
  rz=rfree(rz);
  cx=cvec_free(n,cx);
  cy=cvec_free(2,cy);
  cz=cfree(cz);
  f=func_del(f);
  func_clear();
  printf("check sum: %d\n",func_new_del_check_sum()); 
  */
  return 0;
}
