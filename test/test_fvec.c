#include<stdio.h>
#include<stdlib.h>
#include"isys.h"

int main()
{
  int n,LDJ,digits=5;
  char *begin="begin(x,y,z,w,v,u)";
  char *str="{-x+x^2+2*y^2+2*z^2+2*w^2,-y+2*x*y+2*y*z+2*z*w,-z+y^2+2*x*z+2*y*w,-1+x+2*y+2*z+2*w}";
  //  const char *sx[]={ "2", "-1", "3", "-4", NULL };
  func_t *fF=NULL,*fJ=NULL;
  rmulti **rF=NULL,**rF0=NULL,**rF1=NULL,**rJ=NULL,**rJ0=NULL,**rJ1=NULL,**rx=NULL;
  cmulti **cF=NULL,**cF0=NULL,**cF1=NULL,**cFc=NULL,**cFr=NULL,**cJ=NULL,**cJ0=NULL,**cJ1=NULL,**cJc=NULL,**cJr=NULL,**cx=NULL;

  // vars
  func_eval(func_script(begin));

  // F
  fF=func_eval(func_script(str));
  n=func_asize(fF);
  printf("n=%d\n",n);
  rF=rvec_allocate(n);
  rF0=rvec_allocate(n);
  rF1=rvec_allocate(n);
  cF=cvec_allocate(n);
  cF0=cvec_allocate(n);
  cF1=cvec_allocate(n);
  cFc=cvec_allocate(n);
  cFr=cvec_allocate(n);

  // x
  rx=rvec_allocate(n);
  cx=cvec_allocate(n);
  //  rvec_set_s(n,rx,sx);
  rvec_set_rand(n,rx,2,-1);
  //  cvec_set_s(n,cx,sx);
  cvec_set_rand(n,cx,2,-1);

  // J
  fJ=func_eval(func_grad(func_retain(fF),func_var1_list(n)));
  LDJ=n;
  rJ=rmat_allocate(LDJ,n);
  rJ0=rmat_allocate(LDJ,n);
  rJ1=rmat_allocate(LDJ,n);
  cJ=cmat_allocate(LDJ,n);
  cJ0=cmat_allocate(LDJ,n);
  cJ1=cmat_allocate(LDJ,n);
  cJc=cmat_allocate(LDJ,n);
  cJr=cmat_allocate(LDJ,n);

  // map
  rvec_func_list(n,rF,fF,n,rx);
  irvec_func_list(n,rF0,rF1,fF,n,rx,rx);
  rmat_func_list2(n,n,rJ,LDJ,fJ,n,rx);
  irmat_func_list2(n,n,rJ0,rJ1,LDJ,fJ,n,rx,rx);
  cvec_func_list(n,cF,fF,n,cx);
  icvec_func_list(n,cF0,cF1,fF,n,cx,cx);
  icvec_center_radius(n,cFc,cFr,cF0,cF1);
  cmat_func_list2(n,n,cJ,LDJ,fJ,n,cx);
  icmat_func_list2(n,n,cJ0,cJ1,LDJ,fJ,n,cx,cx);
  icmat_center_radius(n,n,cJc,LDJ,cJr,LDJ,cJ0,LDJ,cJ1,LDJ);


  // output
  printf("fF="); func_print(fF); printf("\n");
  printf("fJ="); func_print(fJ); printf("\n");
  printf("------- real \n");
  rvec_print(n,rx,"x=",'f',digits);
  rvec_print(n,rF,"F=",'f',digits);
  irvec_print(n,rF0,rF1,"[F]=",'f',digits);
  rmat_print(n,n,rJ,LDJ,"J=",'f',digits);
  irmat_print(n,n,rJ0,LDJ,rJ1,LDJ,"[J]=",'f',digits);
  printf("------- complex \n");
  cvec_print(n,cx,"x=",'f',digits);
  cvec_print(n,cF,"F=",'f',digits);
  icvec_print(n,cF0,cF1,"[F]=",'f',digits);
  cvec_print(n,cFc,"Fc=",'f',3);
  cvec_print(n,cFr,"Fr=",'e',3);
  cmat_print(n,n,cJ,LDJ,"J=",'f',digits);
  icmat_print(n,n,cJ0,LDJ,cJ1,LDJ,"[J]=",'f',digits);
  cmat_print(n,n,cJc,LDJ,"Jc=",'f',3);
  cmat_print(n,n,cJr,LDJ,"Jr=",'e',3);


  // free
  fF=func_del(fF);
  fJ=func_del(fJ);
  rx=rvec_free(n,rx);
  rF=rvec_free(n,rF);
  rF0=rvec_free(n,rF0);
  rF1=rvec_free(n,rF1);
  rJ=rmat_free(LDJ,n,rJ);
  rJ0=rmat_free(LDJ,n,rJ0);
  rJ1=rmat_free(LDJ,n,rJ1);
  cx=cvec_free(n,cx);
  cF=cvec_free(n,cF);
  cF0=cvec_free(n,cF0);
  cF1=cvec_free(n,cF1);
  cJ=cmat_free(LDJ,n,cJ);
  cJ0=cmat_free(LDJ,n,cJ0);
  cJ1=cmat_free(LDJ,n,cJ1);
  func_clear();
  printf("check sum: %d\n",func_new_del_check_sum());  

  // done
  return 0;
}
