#include<stdio.h>
#include<stdlib.h>
#include"isys.h"


#define FR(F) func_retain(F)


void output(char *name)
{
  printf("%s=",name); func_eval(func_scriptf("print(%s)",name)); printf("\n");
}


int main()
{
  //set_default_prec(1024);
  int i,n;
  const char *str[]={
    "{x^2+y^2-5,x*y-2,}",
    //    "{2*x^2+y^2+z^2-1,3*x^2+4*z^2-y,9*x-5*z}",                                    //maple比較1
    //    "{2*x*y*z-5*x*z^2+3*y,x^2+z^2-4,3*x^2-4*z^2-2}",                              //maple比較2
    //    "{0.2*x^3-2.4*x*y,0.93*x^2*y-2*y^2+x}",                                       //maple比較3
    //    "{(2/10)*x^3-(24/10)*x*y,(93/100)*x^2*y-2*y^2+x}",                            //maple比較3
    //    "{2*x*y^2*z-x^2*z^2+x^3,2*x^2+z^2,3*x^2-4*z^2-2}",                            //maple比較4
    //    "{0.1*x*y^2*z-x^2*z^2+x^3,0.03*x^2+z^2,0.16*x^2-0.4*z^2-2}",                  //maple比較5
    //    "{(1/10)*x*y^2*z-x^2*z^2+x^3,(3/100)*x^2+z^2,(16/100)*x^2-(4/10)*z^2-2}",     //maple比較5a
    NULL};

  // vars
  func_eval(func_script("begin(x,y,z,w,v,u)"));

  //allocation   

  for(n=0; str[n]!=NULL; n++) ;
  printf("n=%d\n",n);

  for(i=0; i<n; i++){
    printf("=======================\n");
    func_eval(func_scriptf("f=%s",str[i])); output("f");
    func_set_mono_order(FUNC_MONO_LEX);     print_light_red();   printf("--------lex\n");     print_reset();
    func_eval(func_script("f=expand(f)")); printf("original "); output("f");
    func_eval(func_script("g=gbasis(f)")); printf("original "); output("g");

    func_set_mono_order(FUNC_MONO_GRLEX);   print_light_green(); printf("--------grlex\n");   print_reset();
    func_eval(func_script("f=expand(f)")); printf("original "); output("f");
    func_eval(func_script("g=gbasis(f)")); printf("original "); output("g");

    func_set_mono_order(FUNC_MONO_GREVLEX); print_purple();      printf("--------grevlex\n"); print_reset();
    func_eval(func_script("f=expand(f)")); printf("original "); output("f");
    func_eval(func_script("g=gbasis(f)")); printf("original "); output("g");
  }
  
  
//////////////////////////////////////////////////////////////////////////////////////     

/*
 printf("==================================== solution \n\n");
  for(i=0; i<n; i++){
    N=func_asize(f[i]);
    for(count=0; count<3; count++){
      if(count==0)      { func_set_mono_order(FUNC_MONO_LEX);     print_light_red(); printf("--------LEX\n"); }
      else if(count==1) { func_set_mono_order(FUNC_MONO_GRLEX);   print_light_green(); printf("--------GRLEX\n"); }
      else if(count==2) { func_set_mono_order(FUNC_MONO_GREVLEX); printf("\x1b[35m"); printf("--------GREVLEX\n"); }
      f[i]=func_expand(f[i]);
      print_reset();
      set_default_prec(53);
      gsolution=func_poly_list_reduced_groebner(FR(f[i]),0);
      printf("---------grobner solve\n");
      gsolution=func_poly_list_solve_and_arrange(gsolution);
      printf("solution = \n"); func_print(gsolution); printf("\n");
      printf("\n");
      if(gsolution!=NULL){
	if(func_asize(gsolution)!=0){
	  size1=func_asize(gsolution);
	  size2=func_asize(gsolution->a[0]->a[1]);
	  x1=cvec_allocate(size2);
	  for(i1=0; i1<size1; i1++){
	    for(i2=0; i2<size2; i2++){
	      func_change_cmulti(x1[i2],gsolution->a[i1]->a[1]->a[i2]);
	    }
	    printf("----- grobner solution\n");
	    cvec_print(size2,x1,"x1=","f",40);
	  }
	  x1=cvec_free(size2,x1);
	}
      }
      gsolution=func_del(gsolution);
    }
  }
  printf("======================================== end \n");    
  print_reset(); printf("\n");    
*/

//////////////////////////////////////////////////////////////////////////////////////     
/*
  printf("======================================== test all program \n\n");
  for(i=0; i<n; i++){
    printf("original f[%d]=",i); func_print(f[i]); printf("\n");
    N=func_list_size(f[i]);
    for(count=0; count<3; count++){
      if(count==0)      { func_set_mono_order(FUNC_MONO_LEX);     print_light_red(); printf("--------LEX\n"); }
      else if(count==1) { func_set_mono_order(FUNC_MONO_GRLEX);   print_light_green(); printf("--------GRLEX\n"); }
      else if(count==2) { func_set_mono_order(FUNC_MONO_GREVLEX); printf("\x1b[35m"); printf("--------GREVLEX\n"); }
      f[i]=func_expand(f[i]);
      print_reset();
      printf("f[%d]=",i); func_print(f[i]); printf("\n");
      set_default_prec(53);
      gsolution=func_poly_list_reduced_groebner(FR(f[i]),debug);
      printf("---------grobner \n");
      printf("g[%d]=",i); func_print(gsolution); printf("\n"); print_reset();
      printf("---------grobner solution\n");
      gsolution=func_poly_list_solve_and_arrange(gsolution);
      printf("solution = \n"); func_print(gsolution); printf("\n");
      printf("\n");
      if(gsolution!=NULL){
	if(func_asize(gsolution)!=0){
	  fnewton=func_list(0);
	  error=func_list(0);
	  size1=func_asize(gsolution);
	  size2=func_asize(FUNC_ARG1(FUNC_ARG(gsolution,0)));
	  x1=cvec_allocate(size2);
	  for(i1=0; i1<size1; i1++){
	    for(i2=0; i2<size2; i2++){
	      func_change_cmulti(x1[i2],FUNC_ARG(FUNC_ARG1(FUNC_ARG(gsolution,i1)),i2));
	    }
	    set_default_prec(1024);
	    x2=cvec_allocate(size2);
	    cvec_copy(size2,x2,x1);
	    cnvergence(N,f[i],x2,kmax,eta,1024,0);
	    f2=func_change_cmulti_list(size2,x2);
	    func_args_append(fnewton,FR(f2));
	    f2=func_del(f2);
	    cvec_sub(N,x2,x2,x1);
	    f2=func_change_cmulti_list(size2,x2);
	    func_args_append(error,FR(f2));
	    f2=func_del(f2);
	  }
	  x1=cvec_free(size2,x1);
	  x2=cvec_free(size2,x2);
	}
	printf("----- cnewton solution\n");
	printf("newton=\n"); func_print(fnewton); printf("\n\n");
	printf("----- error\n");
	printf("error=\n"); func_print(error); printf("\n\n");
      }
      gsolution=func_del(gsolution);
      fnewton=func_del(fnewton);
      error=func_del(error);
    }
  }
  printf("======================================== end \n");    
  print_reset(); printf("\n");    
*/

  // free
  func_clear();
  printf("check sum: %d\n",func_new_del_check_sum());
  return 0;
}

//EOF
