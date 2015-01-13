#include<stdio.h>
#include<stdlib.h>
#include"isys.h"

#define FR(f) func_retain(f)

  const char *str[]={
    "nan","nan^2","inf","inf^2","0","1","3","-3","(-3)","123","(123)","-123","(-123)","171/13","(171/13)","-171/13",
    "(-171/13)","1/2+1/3", "(1/2)*(1/3)", "(1/2)/(1/3)","123.", "12.34", "0.123", "-0.1","-.1","0^0","0^1","0^2",
    "0^(-2)","1^0","1^2","1^(-2)","2^0","2^2","2^(-2)","(2/3)^2","(2/3)^1","(2/3)^0","(2/3)^(-1)","(2/3)^(-2)","1-1",
    "0.", "1.", "2.-1.", "2.*0.5", "2.^0.","5^(1/2)+2^(1/2)+3^(1/2)+4^(1/2+3/2)","2^(0.5)","2.^(1/2)",
    "2.^(0.5)","#0","#1","#2", "#0^2","#1^(-2)",
    "x","y","x^0","x^12","y^(12)","y^(-12)","(x+y)^(1/2)","(x*y)^2","(x*y)^1","(x*y)^0","(2*x*y)^2",
    "(2*x*y)^1","(2*x*y)^0","(x+x*y)^3","(x*2*y)^2+x*y+(2/3*x*y)^2","(x+x*y)","(x+x*y)^2","(x+x*y)^3","x^x",
    "2*y*(2)*x^2*(-2)*y^3*(-2)*x^5","(-2/3)*x^(-2)*(-3/2)*y^(-3)*x^2*y^3","2/3*x*3*y*(1/2)","2*x*3*sqrt(2)",
    "sqrt(0.)", "sqrt(1.)", "sqrt(2.)", "sqrt(-1.)", "sqrt(-2.)","sqrt(-1)+5","(sqrt(-1)+5)*(3*sqrt(-1)-4)",
    "sqrt(-1.)+5", "(sqrt(-1.)+5)*(3*sqrt(-1)-4)","sqrt(2)+sqrt(-2)+3*sqrt(2)-2*sqrt(-2)",
    "sqrt(2.)+sqrt(-2.)+3*sqrt(2)-2*sqrt(-2)",
    "1+sqrt(2)+1.1+3.3+3/2",
    "1+sqrt(2.)+1.1+3.3+3/2","(3/2)*1.1*1*3.3*sqrt(2)","(3/2)*1.1*1*3.3*sqrt(2.)","sqrt(2)^4*3","sqrt(2)^4*3.",
    "sqrt(2.)^4*3","sqrt(2)^5*3","sqrt(2)^5*3.","sqrt(2.)^5*3","sqrt(2)^(-4)*3","sqrt(2)^(-4)*3.",
    "sqrt(2.)^(-4)*3","sqrt(2)^(-5)*3","sqrt(2)^(-5)*3.","sqrt(2.)^(-5)*3","sqrt(2)*sqrt(2)","sqrt(2)^5*sqrt(2)^3",
    "sqrt(2)^sqrt(2)+sqrt(3)^sqrt(3)","sqrt(2)^sqrt(2)+sqrt(3)^sqrt(3.)",
    "sqrt(2)*sqrt(3)+sqrt(2)*sqrt(4)+sqrt(4)*sqrt(3)+sqrt(5)*sqrt(6)",
    "sqrt(2)*sqrt(3)+sqrt(2)*sqrt(4.)+sqrt(4)*sqrt(3)+sqrt(5)*sqrt(6)",
    "sqrt(2.)^sqrt(2)","sqrt(2)^sqrt(2.)","sqrt(2.)^sqrt(2.)",
    "(1+sqrt(5))*(1+sqrt(3))*(1+sqrt(2))*(1+sqrt(4))*(1+sqrt(1))",
    "(1+sqrt(5))*(1+sqrt(3))*(1+sqrt(2))*(1+sqrt(4.))*(1+sqrt(1))",
    "x+2","x+2.","0.1*x^2", "y+0.5", "0.1*x^2", "0.5*y","sqrt(2)*x+2*x+sqrt(3)*y+3*y",
    "sqrt(2)*x+2.*x+sqrt(3)*y+3.*y","sqrt(-2)*x+2*x+sqrt(3)*y+3*y+sqrt(2)*sqrt(3)*z","sqrt(-2)*x+2.*x+sqrt(3)*y+3.*y",
    "sqrt(-2.)*x+2*x+sqrt(3.)*y+3*y","2*sqrt(2)*x+x+3*sqrt(3)*x+4*sqrt(2)*x+y+5*sqrt(2)*y",
    "x+y", "(x+y)^2", "(x+y)^3", "(x+y)^4", "(x+y)^5", "(x+y)^6","x+1+0+1/2+(-3/2)-x","x-1",
    "expand((x+y)^3)",
    "x-y-x+y+2*x-2*y+z-2*x-3*y-4*z","2*x^2*y+(-3*x^2*y+x)",
    "2*x^2*y+2*x+(-x+x*y^2)-x^2*y","x*y^2*z^3+x^2*y^3*z+x^3*y*z^2",
    "x*y^3*z+y^5+x^2*y^2","x^3*y*z^2+x^2*y^3*z+x*y^2*z^3","x^2*z^2+x*y^2*z+x^3+z^2","3*(x+x*y)","(x+x*y)*(x+y)",
    "x+(y+x^2)*x","(x*(z*(w+v)+(y+x)*(x+y))+1)*(x+y)","(x*2*y*(z+w))^2","((x*2*y*(z+w))^2)^3","x+2*y^2+3*z^3",
    "sqrt(x+y*z)","exp(x+y*z)","log(x+y*z)","sin(x+y*z)","(sin(x+y*z))^2","cos(x+y*z)","tan(x+y*z)","sinh(x+y*z)",
    "cosh(x+y*z)","tanh(x+y*z)","asin(x+y*z)","acos(x+y*z)","atan(x+y*z)","asinh(x+y*z)","acosh(x+y*z)","atanh(x+y*z)",
    "diff(x^x,x)",
    "(x^x)^(-2)",
    "1+x+sin(x*y)+sin(x^2)+cos(x)+cos(y)+(sin(x))^2+(cos(x))^2+2*cos(x)-cos(y)+sqrt(2)*sqrt(3)*(sin(x))^2",
    "(x+y)^2+(x+y)^3+(x+y+z)^2+x^2+(x+y+z)^3+(y+z)^3+(x+y)^2+sqrt(2)*(x+y)^2+sqrt(2)*sqrt(3)*(x+y)^2",
    "(x+2*y)^2+(x+2*y)^3+2*(x+2*y)^2+3*(x+2*y)^3",
    "sqrt(2)+3*sqrt(2)",
    "sqrt(2)*x+sqrt(3)*x+sqrt(2)*x+sqrt(2)*sqrt(3)*x+sqrt(2)*sqrt(3)*y+sqrt(2)*sqrt(3)*z*sin(x)+x^2",
    "(x+y)^2*(x+y)^3","(x+y)^2*(x+y)^3*(x+y)^(x)*(x+y)^5","sin(x)*cos(x)*sin(x)*sin(y)*cos(x)",
    "5+x+y+2*x+3*y+4*x+5*y+x^2+y^2-7*x-9*y+(sqrt(2)*x)^2",
    "5+x+y+x+y", "eval(5+x+y+x+y)",
    "x*sqrt(2)*y*sin(x)*cos(y)",
    "2*x*(x+2*y)","2*x*(x+2*y)+3*y*(y+x)","5./10.","(x^2*y^3)/(x*y^5)","(5*x^2*y^3)/(3*x*y^2)^3",
    "(-1.5e-8)*x", "+1.2e2", "(5e-8)*x", "+2e2", "+2e21", "+2e-128", "2e12", "-2.e-12", "-1e-12",
    "{x,y,z}",
    "{x*y,y*z,z*x}","{x+2*y+3*z,-x-2*y-3*z,x+(1/2)*y+(1/3)*z}","{x+2*y,-x-2*y}",
    "{1/2, 1/3, 1/4, 1/5}","{1/2, 1/sqrt(y), 1/sqrt(x), 1/sqrt(2)}",
    "[1 sqrt(2) 1/3 5]", "[1/sqrt(2) 1/sqrt(3) 1/sqrt(4); 1/sqrt(3) 1 3]",
    "[1 sqrt(-2) 1/3 5]","[1/sqrt(-2) 1/sqrt(3) 1/sqrt(4); 1/sqrt(3) 1 3]",
    "{x,y,z}","{{x,y,z},{2*x,2*y,2*z},{-x,-y,-z}}",
    "[1 2 3]+[4 5 6]","-[1 2 3]","[1 2 3]-[4 5 6]","-[1 2 3]-[4 5 6]",
    "x/y","x/y/z","x*y/z*u/v",
    "grad(x*y^2*z^3*w^4,{x,y,z,w})","grad(x*y^2*z^3*w^4,4)",
    "grad({x+y+z,x*y+y*z+z*x,x*y*z},{x,y,z})",
    "-x","x-y","-x+y","-x-y","1/2-1/3","0.5-0.3","1.5-I","(1/2+I)-(1/3+2*I)","I","i","PI",
    "f(x,y)", "f(u,v)", "f(1,2)",
    NULL};


int main()
{
  int i,n;
  rmulti *ry=NULL;
  cmulti *cy=NULL;
  char key[1024];
  func_t *f=NULL;

  print_dark_gray(); func_print(func_scope(0)); print_reset(); printf("\n"); 

  for(n=0; str[n]!=NULL; n++) ;
  printf("n=%d\n",n);

  ry=rallocate();
  cy=callocate();

  print_red();  printf("> "); f=func_scriptf_echo("begin(x,y,z,w,v,u)"); print_reset(); printf("\n");
  //print_red();  printf("> "); f=func_scriptf_echo("begin"); print_reset(); printf("\n");
  print_yellow(); printf("> "); func_print(f); print_reset(); printf("\n");
  f=func_eval(f); func_print(f); printf("\n\n"); f=func_del(f);

  print_red();  printf("> "); f=func_scriptf_echo("f=NULL"); print_reset(); printf("\n");
  //  print_red();  printf("> "); f=func_scriptf_echo("f(a,b)=a^2+3*b^2"); print_reset(); printf("\n");
  //print_red();  printf("> "); f=func_scriptf_echo("f(u,v)=u^2+3*v^2"); print_reset(); printf("\n");
  //print_red();  printf("> "); f=func_scriptf_echo("begin"); print_reset(); printf("\n");
  print_yellow(); printf("> "); func_print(f); print_reset(); printf("\n");
  f=func_eval(f); func_print(f); printf("\n\n"); f=func_del(f);

  //  print_red(); printf("> "); f=func_scriptf_echo("lx={-2, 1,-1, -2, 3}"); print_reset(); printf("\n"); printf("lx="); func_scriptf("print(lx)"); printf("\n"); if(f!=NULL){ func_print(f); printf("\n"); f=func_del(f);  }
  //  print_red(); printf("> "); f=func_scriptf_echo("rx=[-2. 1. -1. -2. 3.]"); print_reset(); printf("\n"); printf("rx="); func_scriptf("print(rx)"); printf("\n"); if(f!=NULL){ func_print(f); printf("\n"); f=func_del(f);  }
  //  print_red(); printf("> "); f=func_scriptf_echo("cx=[-2. I -1. -2. 3.]"); print_reset(); printf("\n"); printf("cx="); func_scriptf("print(cx)"); printf("\n"); if(f!=NULL){ func_print(f); printf("\n"); f=func_del(f);  }


  for(i=0; i<n; i++){
    sprintf(key,"f%d",i);
    print_red(); printf("> "); f=func_scriptf_echo("%s",str[i]); print_reset(); printf("\n");
    print_yellow(); printf("> "); func_print(f); print_reset(); printf("\n");
    f=func_eval(f); func_print(f); printf("\n\n"); f=func_del(f);

    print_red(); printf("> "); f=func_scriptf_echo("f%d=%s",i,str[i]); print_reset(); printf("\n");
    print_yellow(); printf("> "); func_print(f); print_reset(); printf("\n");
    f=func_eval(f); printf("> "); func_print(f); printf("\n"); f=func_del(f);
    printf("f%d=",i); func_eval(func_scriptf("print(f%d)",i)); printf("\n\n");

    //    print_red(); printf("> "); f=func_scriptf_echo("g%d=expand(f%d)",i,i); print_reset(); printf("\n");
    //    print_yellow(); printf("> "); func_print(f); print_reset(); printf("\n");
    //    f=func_eval(f); printf("> "); func_print(f); printf("\n"); f=func_del(f);
    //    printf("g%d=",i); func_eval(func_scriptf("print(g%d)",i)); printf("\n\n");

    //    print_red(); printf("> "); f=func_scriptf_echo("hx%d=diff(f%d,x)",i,i); print_reset(); printf("\n");
    //    print_yellow(); printf("> "); func_print(f); print_reset(); printf("\n");
    //    f=func_eval(f); printf("> "); func_print(f); printf("\n"); f=func_del(f);
    //    printf("hx%d=",i); func_eval(func_scriptf("print(hx%d)",i)); printf("\n\n");


    //    print_purple(); printf("> "); f=func_scriptf_echo("hy%d=diff(f%d,y)",i,i); print_reset(); printf("\n"); printf("hy%d=",i); func_scriptf("print(hy%d)",i); printf("\n"); if(f!=NULL){ func_print(f); printf("\n"); f=func_del(f);  }
    //    print_purple(); printf("> "); f=func_scriptf_echo("hz%d=diff(f%d,z)",i,i); print_reset(); printf("\n"); printf("hz%d=",i); func_scriptf("print(hz%d)",i); printf("\n"); if(f!=NULL){ func_print(f); printf("\n"); f=func_del(f);  }
    //    print_cyan(); printf("> "); f=func_scriptf_echo("kl%d=f%d(lx)",i,i); print_reset(); printf("\n"); printf("kl%d=",i); func_scriptf("print(kl%d)",i); printf("\n"); if(f!=NULL){ func_print(f); printf("\n"); f=func_del(f);  }
    //    print_green(); printf("> "); f=func_scriptf_echo("kr%d=evalf(f%d(rx))",i,i); print_reset(); printf("\n"); printf("kr%d=",i); func_scriptf("print(kr%d)",i); printf("\n"); if(f!=NULL){ func_print(f); printf("\n"); f=func_del(f);  }
    //    rvec_mapped_rmulti(ry,func_find(key),5,func_rvec_p(func_find("rx"))); print_yellow(); mpfr_printf("kr%d=(%.5Rg)",i,ry); print_reset(); printf("\n");
    //    print_green(); printf("> "); f=func_scriptf_echo("kc%d=evalf(f%d(cx))",i,i); print_reset(); printf("\n"); printf("kc%d=",i); func_scriptf("print(kc%d)",i); printf("\n"); if(f!=NULL){ func_print(f); printf("\n"); f=func_del(f);  }
    //    cvec_mapped_cmulti(cy,func_find(key),5,func_cvec_p(func_find("cx"))); print_yellow(); mpfr_printf("kc%d=(%.5Rg%+.5Rg*I)",i,C_R(cy),C_I(cy)); print_reset(); printf("\n");
  }

  print_dark_gray(); func_print(func_scope(0)); print_reset(); printf("\n"); 
  print_red(); printf("> "); f=func_eval(func_scriptf_echo("end")); print_reset(); printf("\n"); if(f!=NULL){ func_print(f); printf("\n"); f=func_del(f); }
  print_dark_gray(); func_print(func_scope(0)); print_reset(); printf("\n"); 

  // free
  ry=rfree(ry);
  cy=cfree(cy);

  print_red(); printf("> clear\n"); print_reset(); func_clear();
  printf("check sum: %d\n",func_new_del_check_sum());  
  return 0;
}
