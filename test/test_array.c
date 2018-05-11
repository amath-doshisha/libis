#include<math.h>
#include<isys.h>

int main(int argc, const char *argv[])
{
  int i,k,l;
  char x_type='r',y_type='r';
  int x_ndim=5,x_dim[]={2,2,2,1,2};
  int y_ndim=5,y_dim[]={1,2,1,2,1};
  array *x=NULL,*y=NULL,*z=NULL;

  for(i=0; i<argc; i++){
    if(STR_EQ(argv[i],"-h")){ exit(0); }
    if(STR_EQ(argv[i],"-x-type") && i+1<argc){ x_type=argv[++i][0]; }
    if(STR_EQ(argv[i],"-y-type") && i+1<argc){ y_type=argv[++i][0]; }
    if(STR_EQ(argv[i],"-x-ndim") && i+1<argc){ x_ndim=atoi(argv[++i]); }
    if(STR_EQ(argv[i],"-y-ndim") && i+1<argc){ y_ndim=atoi(argv[++i]); }
  }

  x=array_allocate(x_type,x_ndim,x_dim);
  y=array_allocate(y_type,y_ndim,y_dim);
  
  //  init_genrand(0);
  //  array_set_rand(x,-1000,1000);
  //  array_set_rand(y,-1000,1000);
  
  array_set_grid(x);
  array_set_grid(y);
  z=array_add(x,y);
  //  z=array_add(y,x);
//  z=array_sub(x,y);


  printf("---\n");  array_print(x,"x",'g',3);
  printf("---\n");  array_print(y,"y",'g',3);
  printf("---\n");  array_print(z,"z",'g',3);
  printf("---\n");


  printf("size(x)="); ivec_put(ARRAY_NDIM(x),ARRAY_DIM_P(x),"x"); printf("\n");
  printf("size(y)="); ivec_put(ARRAY_NDIM(y),ARRAY_DIM_P(y),"x"); printf("\n");
  printf("size(z)="); ivec_put(ARRAY_NDIM(z),ARRAY_DIM_P(z),"x"); printf("\n");

  if(array_same_dim_check(x,y)){ printf("dim(x)=dim(y)\n"); }else{ printf("dim(x)!=dim(y)\n"); }
  k=array_get_subdim(x,y);
  l=array_get_subdim(y,x);
       if(k>0){ printf("dim(x)<dim(y), k=%d\n",k); }
  else if(l>0){ printf("dim(x)<dim(y), l=%d\n",l); }
  else        { printf("dim(x)<>dim(y)\n"); }
  printf("compatible(x,y)=%d\n",array_compatible_dim_check(x,y));
  

  x=array_free(x);
  y=array_free(y);
  z=array_free(z);


  /*
  int s_digits=3,digits=6,i,x_type='d',y_type='m',z_type='s',ndim=2,dim[]={5,3,2,2,3},adim[]={3,3};
  char f='f',mode[]="          ";
  array *x=NULL,*y=NULL,*z=NULL,*a=NULL;

  // options
  for(i=0; i<argc; i++){
    if(STR_EQ(argv[i],"-ndim") && i+1<argc){ ndim=atoi(argv[++i]); }
    else if(STR_EQ(argv[i],"-zeros")){ strcpy(mode,"zeros"); }
    else if(STR_EQ(argv[i],"-ones")){ strcpy(mode,"ones"); }
    else if(STR_EQ(argv[i],"-nan")){ strcpy(mode,"nan"); }
    else if(STR_EQ(argv[i],"-inf")){ strcpy(mode,"inf"); }
    else if(STR_EQ(argv[i],"-grid")){ strcpy(mode,"grid"); }
    else if(STR_EQ(argv[i],"-rand")){ strcpy(mode,"rand"); }
    else if(STR_EQ(argv[i],"-x-s")){ x_type='s'; }
    else if(STR_EQ(argv[i],"-x-i")){ x_type='i'; }
    else if(STR_EQ(argv[i],"-x-d")){ x_type='d'; }
    else if(STR_EQ(argv[i],"-x-z")){ x_type='z'; }
    else if(STR_EQ(argv[i],"-x-r")){ x_type='r'; }
    else if(STR_EQ(argv[i],"-x-c")){ x_type='c'; }
    else if(STR_EQ(argv[i],"-x-R")){ x_type='R'; }
    else if(STR_EQ(argv[i],"-x-C")){ x_type='C'; }
    else if(STR_EQ(argv[i],"-y-copy"))  { y_type='a'; }
    else if(STR_EQ(argv[i],"-y-clone")) { y_type='b'; }
    else if(STR_EQ(argv[i],"-y-char"))  { y_type='c'; }
    else if(STR_EQ(argv[i],"-y-int"))   { y_type='i'; }
    else if(STR_EQ(argv[i],"-y-double")){ y_type='d'; }
    else if(STR_EQ(argv[i],"-y-multi")) { y_type='m'; }
    else if(STR_EQ(argv[i],"-y-imulti")){ y_type='I'; }
    else if(STR_EQ(argv[i],"-y-real"))  { y_type='r'; }
    else if(STR_EQ(argv[i],"-y-complex")){ y_type='C'; }   
    else if(STR_EQ(argv[i],"-z-copy"))  { z_type='a'; }
    else if(STR_EQ(argv[i],"-z-clone")) { z_type='b'; }
    else if(STR_EQ(argv[i],"-z-char"))  { z_type='c'; }
    else if(STR_EQ(argv[i],"-z-int"))   { z_type='i'; }
    else if(STR_EQ(argv[i],"-z-double")){ z_type='d'; }
    else if(STR_EQ(argv[i],"-z-multi")) { z_type='m'; }
    else if(STR_EQ(argv[i],"-z-imulti")){ z_type='I'; }
    else if(STR_EQ(argv[i],"-z-real"))  { z_type='r'; }
    else if(STR_EQ(argv[i],"-z-complex")){ z_type='C'; }
    else if(STR_EQ(argv[i],"-char-digits") && i+1<argc){ s_digits=atoi(argv[++i]); }
    else if(STR_EQ(argv[i],"-digits") && i+1<argc){ digits=atoi(argv[++i]); }
    else if(STR_EQ(argv[i],"-f")){ f='f'; }
    else if(STR_EQ(argv[i],"-e")){ f='e'; }
    else if(STR_EQ(argv[i],"-g")){ f='g'; }
  }

  // allocate of x
  x=array_allocate(x_type,ndim,dim);
  init_genrand(0);
  array_set_rand(x,-1000,1000);

  // init of x
  if(STR_EQ(mode,"zeros")){ array_set_zeros(x); }
  else if(STR_EQ(mode,"ones")){ array_set_ones(x); }
  else if(STR_EQ(mode,"nan")){ array_set_nan(x); }
  else if(STR_EQ(mode,"inf")){ array_set_inf(x); }
  else if(STR_EQ(mode,"inf")){ array_set_all_d(x,-M_PI); }
  else if(STR_EQ(mode,"rand")){ init_genrand(0); array_set_rand(x,-1000,1000); }
  else if(STR_EQ(mode,"grid")){ array_set_grid(x); }
  
  // output of x
  printf("type(x)=%c\n",ARRAY_TYPE(x));
  printf("size(x)=["); ivec_put(ARRAY_NDIM(x),ARRAY_DIM_P(x),","); printf("]\n");
  array_print(x,"x",f,digits);
  //  array_put(x);

  printf("==============\n");
  
  // y=func(x)
       if(y_type=='c'){ printf("y=char(x)\n");    y=array_get_char   (x,f,s_digits); }
  else if(y_type=='a'){ printf("y=round(x)\n");   y=array_copy       (x); }
  else if(y_type=='b'){ printf("y=x\n");          y=array_clone      (x); }
  else if(y_type=='i'){ printf("y=int(x)\n");     y=array_get_int    (x); }
  else if(y_type=='d'){ printf("y=double(x)\n");  y=array_get_double (x); }
  else if(y_type=='m'){ printf("y=multi(x)\n");   y=array_get_multi  (x); }
  else if(y_type=='I'){ printf("y=imulti(x)\n");  y=array_get_imulti (x); }
  else if(y_type=='r'){ printf("y=real(x)\n");    y=array_get_real   (x); }
  else if(y_type=='C'){ printf("y=complex(x)\n"); y=array_get_complex(x); }
  else{ y=array_clone(x); }

  // output of y
  if(y!=NULL){
    printf("type(y)=%c\n",ARRAY_TYPE(y));
    printf("size(y)=["); ivec_put(ARRAY_NDIM(y),ARRAY_DIM_P(y),","); printf("]\n");
    array_print(y,"y",f,digits);
    //  array_put(y);
  }else{
    printf("y=NULL\n");
  }
  
  printf("==============\n");
  
  // z=func(y)
       if(z_type=='c'){ printf("z=char(y)\n");    z=array_get_char   (y,f,s_digits); }
  else if(z_type=='a'){ printf("z=round(y)\n");   z=array_copy       (y); }
  else if(z_type=='b'){ printf("z=y\n");          z=array_clone      (y); }
  else if(z_type=='i'){ printf("z=int(y)\n");     z=array_get_int    (y); }
  else if(z_type=='d'){ printf("z=double(y)\n");  z=array_get_double (y); }
  else if(z_type=='m'){ printf("z=multi(y)\n");   z=array_get_multi  (y); }
  else if(z_type=='I'){ printf("z=imulti(y)\n");  z=array_get_imulti (y); }
  else if(z_type=='r'){ printf("z=real(y)\n");    z=array_get_real   (y); }
  else if(z_type=='C'){ printf("z=complex(y)\n"); z=array_get_complex(y); }
  else                { printf("z=y\n");          z=array_clone      (y); }

  // output of y
  if(z!=NULL){    
    printf("type(z)=%c\n",ARRAY_TYPE(z));
    printf("size(z)=["); ivec_put(ARRAY_NDIM(z),ARRAY_DIM_P(z),","); printf("]\n");
    array_print(z,"z",f,digits);
    //  array_put(z);
  }else{
    printf("z=NULL\n");
  }

  printf("==============\n");
  // allocate of a
  a=array_allocate('a',1,adim);
  ARRAY_AVEC(a,0)=x;
  ARRAY_AVEC(a,1)=y;
  ARRAY_AVEC(a,2)=z;
  array_put(a);
  array_print(a,"a",f,digits);

  //  x=array_free(x);
  //  y=array_free(y);
  //  z=array_free(z);
  a=array_free(a);
  */
}
