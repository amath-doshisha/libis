#include<isys.h>

#define SPACE " \t\n"
#define SEP   ","
#define MASK0 "{[("
#define MASK1 "}])"

int main(int argc, char *argv[])
{
  char *s;
  strings *list=NULL;
    
  if(argc>=2){
    s=malloc(sizeof(char)*(strlen(argv[1])+1));
    strcpy(s,argv[1]);
  }else{
    s=malloc(sizeof(char)*(strlen("1234567890")+1));
    strcpy(s,"1234567890");
  }
  printf("s=%s\n",s);
  list=strings_split_number(s);
  strings_print(list);
  printf("\n");
    
  /*
  int k=0;
  const char *str[]={"  a ","b","c","d",NULL};
  const char *Str[]={"A","   B ","C",NULL};
  const char *ex1="sin(x+y*(a+b)+c)";
  const char *ex2="f(x,y,z)";
  const char *ex3="(x,y,z)";
  char *s=NULL,*p=NULL;
  strings *list=NULL;

  if(argc>=2){
    k=atoi(argv[1]);
  }
  printf("k=%d\n",k);

  list=strings_new_str(str,SPACE);  strings_print(list); printf("\n");
  strings_item_replace(list,k,strings_new_str(Str,SPACE));  strings_print(list); printf("\n");
  list=strings_del(list);

  printf("%s=",ex1);
  //  strings_item_replace(list,-1,strings_split_mask(ex,MASK0,MASK1,SPACE));
  list=strings_split_mask(ex1,MASK0,MASK1,SPACE);
  strings_print(list); printf("\n");
  list=strings_del(list);

  printf("%s=",ex2);
  list=strings_split_mask(ex2,MASK0,MASK1,SPACE);
  strings_item_replace(list,1,strings_split(list->str[1],SEP,MASK0,MASK1,SPACE));
  strings_print(list); printf("\n");
  list=strings_del(list);

  printf("%s=",ex3);
  list=strings_split_mask(ex3,MASK0,MASK1,SPACE);
  strings_item_replace(list,0,strings_split(list->str[0],SEP,MASK0,MASK1,SPACE));
  strings_print(list); printf("\n");
  list=strings_del(list);

  s="abc1234567890";  p="abc";  printf("match('%s',/%s/)=%d\n",s,p,str_match(s,p,NULL,NULL));
  s="1234abc567890";  p="abc";  printf("match('%s',/%s/)=%d\n",s,p,str_match(s,p,NULL,NULL));
  s="1234567890abc";  p="abc";  printf("match('%s',/%s/)=%d\n",s,p,str_match(s,p,NULL,NULL));
  s="1234567890ab";  p="abc";  printf("match('%s',/%s/)=%d\n",s,p,str_match(s,p,NULL,NULL));
  s="1234ab567890";  p="abc";  printf("match('%s',/%s/)=%d\n",s,p,str_match(s,p,NULL,NULL));

  s="abc1234567890";  p="^abc";  printf("match('%s',/%s/)=%d\n",s,p,str_match(s,p,NULL,NULL));
  s="ab1234567890";  p="^abc";  printf("match('%s',/%s/)=%d\n",s,p,str_match(s,p,NULL,NULL));
  s="1234abc567890";  p="^abc";  printf("match('%s',/%s/)=%d\n",s,p,str_match(s,p,NULL,NULL));
  s="1234567890abc";  p="^abc";  printf("match('%s',/%s/)=%d\n",s,p,str_match(s,p,NULL,NULL));
  s="1234567890ab";  p="^abc";  printf("match('%s',/%s/)=%d\n",s,p,str_match(s,p,NULL,NULL));

  s="abc1234567890";  p="abc$";  printf("match('%s',/%s/)=%d\n",s,p,str_match(s,p,NULL,NULL));
  s="ab1234567890";  p="abc$";  printf("match('%s',/%s/)=%d\n",s,p,str_match(s,p,NULL,NULL));
  s="1234abc567890";  p="abc$";  printf("match('%s',/%s/)=%d\n",s,p,str_match(s,p,NULL,NULL));
  s="1234567890abc";  p="abc$";  printf("match('%s',/%s/)=%d\n",s,p,str_match(s,p,NULL,NULL));
  s="1234567890ab";  p="abc$";  printf("match('%s',/%s/)=%d\n",s,p,str_match(s,p,NULL,NULL));

  s="ab";  p="^abc$";  printf("match('%s',/%s/)=%d\n",s,p,str_match(s,p,NULL,NULL));
  s="abc";  p="^abc$";  printf("match('%s',/%s/)=%d\n",s,p,str_match(s,p,NULL,NULL));
  s="abc0";  p="^abc$";  printf("match('%s',/%s/)=%d\n",s,p,str_match(s,p,NULL,NULL));
  s="0abc";  p="^abc$";  printf("match('%s',/%s/)=%d\n",s,p,str_match(s,p,NULL,NULL));

  s="1234567890";  p="%a";  printf("match('%s',/%s/)=%d\n",s,p,str_match(s,p,NULL,NULL));
  s="1234a567890";  p="%a";  printf("match('%s',/%s/)=%d\n",s,p,str_match(s,p,NULL,NULL));
  s="1234a567890";  p="^%a";  printf("match('%s',/%s/)=%d\n",s,p,str_match(s,p,NULL,NULL));
  s="-1234567890";  p="^-";  printf("match('%s',/%s/)=%d\n",s,p,str_match(s,p,NULL,NULL));
  s="-1234567890";  p="^-?";  printf("match('%s',/%s/)=%d\n",s,p,str_match(s,p,NULL,NULL));
  s="1234567890";  p="^-?";  printf("match('%s',/%s/)=%d\n",s,p,str_match(s,p,NULL,NULL));

  s="1234567890.";  p=".%d?$";  printf("match('%s',/%s/)=%d\n",s,p,str_match(s,p,NULL,NULL));
  s="1234567890.";  p=".%d+$";  printf("match('%s',/%s/)=%d\n",s,p,str_match(s,p,NULL,NULL));
  s="1234567890.1";  p=".%d+$";  printf("match('%s',/%s/)=%d\n",s,p,str_match(s,p,NULL,NULL));
  s="1234567890.12";  p=".%d+$";  printf("match('%s',/%s/)=%d\n",s,p,str_match(s,p,NULL,NULL));
  s="1234567890.";  p=".%d*$";  printf("match('%s',/%s/)=%d\n",s,p,str_match(s,p,NULL,NULL));
  s="1234567890.1";  p=".%d*$";  printf("match('%s',/%s/)=%d\n",s,p,str_match(s,p,NULL,NULL));
  s="1234567890.12";  p=".%d*$";  printf("match('%s',/%s/)=%d\n",s,p,str_match(s,p,NULL,NULL));
  s="1234567890abcde";  p="^%d+%a+$";  printf("match('%s',/%s/)=%d\n",s,p,str_match(s,p,MASK0,MASK1));


  s="(1234567890)";  p="^%(%m+%)$";  printf("match('%s',/%s/)=%d\n",s,p,str_match(s,p,MASK0,MASK1));
  s="(1234(56)7890)";  p="^%(%m+%)$";  printf("match('%s',/%s/)=%d\n",s,p,str_match(s,p,MASK0,MASK1));
  s="abc(1234567890)def";  p="^%(%m+%)$";  printf("match('%s',/%s/)=%d\n",s,p,str_match(s,p,MASK0,MASK1));
  s="abc(1234{5678)90}def";  p="^%(%m+%)$";  printf("match('%s',/%s/)=%d\n",s,p,str_match(s,p,MASK0,MASK1));
  s="abc(1234{5678}90)def";  p="^%(%m+%)$";  printf("match('%s',/%s/)=%d\n",s,p,str_match(s,p,MASK0,MASK1));
  s="abc(12[34{5678}9]0)def";  p="^%(%m+%)$";  printf("match('%s',/%s/)=%d\n",s,p,str_match(s,p,MASK0,MASK1));
  s="(12[34{5678}9]0)";  p="^%(%m+%)$";  printf("match('%s',/%s/)=%d\n",s,p,str_match(s,p,MASK0,MASK1));
  */

  return 0;
}
