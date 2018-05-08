#include<stdio.h>
#include<stdlib.h>
//#include<strings.h>
#include<string.h>
#include<limits.h>
#include"is_macros.h"
#include"is_strings.h"
#include"is_print.h"
#include"is_ivec.h"

static int __strings_obj_id=0;

/////////////////////////////////////////////

char *char_new(char *str, char *skip)
{
  int len;
  char *p=NULL;
  if(str==NULL){
    p=NULL;
  }else{
    while(char_match_any(*str,skip,NULL)){ str++; }
    len=(int)strlen(str);
    while(len>0 && char_match_any(*(str+len-1),skip,NULL)){ len--; }
    p=(char*)malloc(sizeof(char)*(len+1));
    strncpy(p,str,len);
    p[len]='\0';
  }
  return p;
}

char *char_del(char *p)
{
  if(p!=NULL){ free(p); }
  return NULL;
}

////////////////////////////////////////////////

char *char_clone(char *str)
{
  char *p=NULL;
  p=(char*)malloc(sizeof(char)*(strlen(str)+1));
  strcpy(p,str);
  return p;
}

////////////////////////////////////////////////

char *char_renew(char *old_str, char *str, char *skip)
{
  old_str=char_del(old_str);
  return char_new(str,skip);
}

////////////////////////////////////////////////

#define CHAR_RENEW_SPRINTF_BUF 10000

char *char_renew_sprintf(char *old, char *skip, char* fmt, ...)
{
  char buf[CHAR_RENEW_SPRINTF_BUF+1];
  va_list argp;
  va_start(argp,fmt);
  vsprintf(buf,fmt,argp);
  old=char_del(old);
  return char_renew(old,buf,skip);
}

////////////////////////////////////////////////

char *char_cat(char *str_old, char *str_append)
{
  char *str_new=NULL;
  str_new=(char*)malloc(sizeof(char)*(strlen(str_old)+strlen(str_append)+1));
  strcpy(str_new,str_old);
  strcat(str_new,str_append);
  str_old=char_del(str_old);
  return str_new;
}

////////////////////////////////////////////////


int char_cmp(char *a, char *b)
{
  if(a==NULL && b==NULL){ return 1; } // a=b
  if(a==NULL){ return -1; }           // a<b
  if(b==NULL){ return +1; }           // a>b
  return strcmp(a,b);
}

int char_eq(char *a, char *b)
{
  return char_cmp(a,b)==0;
}


////////////////////////////////////////////////////////

int char_match_any(char c, char *str, int *pos)
{
  int value=0;
  if(str==NULL) return 0;
  if(pos!=NULL) (*pos)=0;
  while(!value && (*str)!='\0'){
    if(c==(*str)){ value=1; }
    else{ if(pos!=NULL) (*pos)++; }
    str++;
  }
  return value;
}

////////////////////////////////////////////////////////

int str_has_any_char(char *s, char *c)
{
  int ret=0;
  if(s==NULL || c==NULL){ return 0; }
  while(ret==0 && (*s)!='\0'){
    if(char_match_any((*s),c,NULL)){ ret=1; }
    s++;
  }
  return ret;
}

////////////////////////////////////////////////////////

char *str_create_mask(char *s, char *mask_begin, char *mask_end)
{
  int n,count,i,k;
  char *mask=NULL,*m=NULL;
  n=(int)strlen(s);
  mask=(char*)malloc(sizeof(char)*(n+1));
  m=(char*)malloc(sizeof(char)*(2*(n+1)+1));
  for(count=0, i=0; i<n; i++){
    if(char_match_any(s[i],mask_begin,&k)){
      m[2*count]=mask_end[k]; m[2*count+1]='\0';
      if(count==0) mask[i]='('; else mask[i]='*';
      count++;
    }else if(count>0 && char_match_any(s[i],m+2*(count-1),&k)){
      count--;
      if(count==0) mask[i]=')'; else mask[i]='*';
    }else if(count>0){
      mask[i]='*';
    }else{
      mask[i]='.';
    }
  }
  mask[n]='\0';
  m=char_del(m);
  return mask;
}

#define STR_DIGITS    "0123456789"
#define STR_ALPHA     "abcdefghijklmnopqrstuvwxyz"
#define STR_ALPHA_ALL "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ"
#define STR_IDENT     "0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ_"
#define STR_SIGN      "+-"
#define STR_POINT     "."
#define STR_EXP       "eE"

int str_match(char *s, char *p0, char *mask_begin, char *mask_end)
{
  char *p=NULL;
  char *m0=NULL,*m=NULL;
  int done=0,matched=0,dp=1,ds=1,count=0,doing=0;
  // マスクの抜き出し
  m0=str_create_mask(s,mask_begin,mask_end);
  // 不一致になるか，パターンの最後まで到達するまでループを継続
  count=0; p=p0; m=m0;
  while(!done){
    // 文字どうしのチェック
    if     ((*p)=='^')                  { dp=1; ds=0; matched=1; } // 先頭指定^がある場合はマッチングを開始する
    else if((*p)=='%' && (*(p+1))=='^') { dp=2; ds=1; if((*m)!='*' && (*s)=='^') matched=1; else matched=0; }     // %^
    else if((*p)=='%' && (*(p+1))=='?') { dp=2; ds=1; if((*m)!='*' && (*s)=='?') matched=1; else matched=0; }     // %?
    else if((*p)=='%' && (*(p+1))=='+') { dp=2; ds=1; if((*m)!='*' && (*s)=='+') matched=1; else matched=0; }     // %+
    else if((*p)=='%' && (*(p+1))=='*') { dp=2; ds=1; if((*m)!='*' && (*s)=='*') matched=1; else matched=0; }     // %*
    else if((*p)=='%' && (*(p+1))=='s') { dp=2; ds=1; if((*m)!='*') matched=char_match_any((*s)," \t\n",NULL); else matched=0; }       // %s 空白
    else if((*p)=='%' && (*(p+1))=='p') { dp=2; ds=1; if((*m)!='*') matched=char_match_any((*s),"+-",NULL); else matched=0; }          // %p +-
    else if((*p)=='%' && (*(p+1))=='d') { dp=2; ds=1; if((*m)!='*') matched=char_match_any((*s),STR_DIGITS,NULL); else matched=0; }    // %d 数字
    else if((*p)=='%' && (*(p+1))=='a') { dp=2; ds=1; if((*m)!='*') matched=char_match_any((*s),STR_ALPHA,NULL); else matched=0; }     // %a アルファベット小文字
    else if((*p)=='%' && (*(p+1))=='A') { dp=2; ds=1; if((*m)!='*') matched=char_match_any((*s),STR_ALPHA_ALL,NULL); else matched=0; } // %A アルファベット
    else if((*p)=='%' && (*(p+1))=='I') { dp=2; ds=1; if((*m)!='*') matched=char_match_any((*s),STR_IDENT,NULL); else matched=0; }     // %I 識別子
    else if((*p)=='%' && (*(p+1))=='(') { dp=2; ds=1; if((*m)=='(') matched=1; else matched=0; }                                       // %( マスク開始
    else if((*p)=='%' && (*(p+1))==')') { dp=2; ds=1; if((*m)==')') matched=1; else matched=0; }                                       // %) マスク終了
    else if((*p)=='%' && (*(p+1))=='m') { dp=2; ds=1; if((*m)=='*') matched=1; else matched=0; }                                       // %m マスク途中
    else                                { dp=1; ds=1; if((*m)!='*' && (*s)==(*p)) { matched=1; } else { matched=0; } }                 // 一文字
    // 0回または1回の繰り返し
    if((*(p+dp))=='?'){
      if(!matched){ ds=0; }
      matched=1; dp++;
    }
    // 1回以上の繰り返し
    if((*(p+dp))=='+'){
      if(matched){ count++; dp=0; }
      else{
	if(count){ matched=1; } else{ matched=0; }
	count=0; ds=0; dp++;
      }
    }
    // 0回以上の繰り返し
    if((*(p+dp))=='*'){
      if(matched){ count++; dp=0; }
      else       { count=0; dp++; ds=0; matched=1; }
    }
    // 次のパターンへ移動
    if     (doing  &&  matched)              { doing=1; matched=1; p+=dp; }
    else if(doing  && !matched && (*p0)=='^'){ doing=1; matched=0; while((*p)!='\0') p++; } // 先頭指定^がある場合は，マッチングを終了
    else if(doing  && !matched && (*p0)!='^'){ doing=0; matched=0; p=p0; }                  // 先頭指定^がない場合は，マッチングを最初からやり直す
    else if(!doing &&  matched)              { doing=1; matched=1; p+=dp; }                 // マッチングを開始する
    else if(!doing && !matched)              { doing=0; matched=0; }                        // マッチングを未開始のまま継続する
    // 次の文字へ移動
    if((*s)=='\0'){ ds=0; }
    s+=ds;
    m+=ds;
    // 終了？
    if((*p)=='\0'){
      if(doing && matched){ done=1; matched=1; }
      else                { done=1; matched=0; }
    }else if((*p)=='$' && (*(p+1))=='\0'){
      if((*s)=='\0')      { done=1; matched=1; }
      else                { done=1; matched=0; }
    }else if((*s)=='\0'){
      if(!doing)          { done=1; matched=0; }
    }
  }
  // done
  m0=char_del(m0);
  return matched;
}

////////////////////////////////////////////////////////

strings* strings_new(int n)
{
  int i;
  strings *list;
  list=(strings*)malloc(sizeof(strings));
  list->id=__strings_obj_id++;
  if(n>0){
    list->n=n;
    list->str=(char**)malloc(sizeof(char*)*(list->n));
    for(i=0; i<list->n; i++){ list->str[i]=NULL; }
  }else{
    list->n=0;
    list->str=NULL;
  }
#ifdef DEBUG
  print_dark_gray();
  printf("strings_new id=%010d n=%d",list->id,list->n);
  print_reset();
  printf("\n");
#endif
  return list;
}

strings* strings_new_str(char *str[], char *skip)
{
  int i;
  strings *list;
  for(i=0; str[i]!=NULL; i++) ;
  list=strings_new(i);
  for(i=0; i<strings_size(list); i++){ strings_item_set(list,i,str[i],skip); }
  return list;
}

strings* strings_new_clone(strings *src)
{
  int i;
  strings *list;
  if(src==NULL){ return NULL;; }
  list=strings_new(strings_size(src));
  for(i=0; i<strings_size(list); i++){ strings_item_set(list,i,strings_at(src,i),NULL); }
  return list;
}

////////////////////////////////////////////////////////

int strings_size(strings *str)
{
  if(str==NULL){ return 0; }
  else         { return str->n; }
}

char *strings_at(strings *str, int i)
{
  if(str==NULL){ return NULL; }
  else         { return str->str[i]; }
}


/////////////////////////////////////////////////////////////

void strings_del_str(strings *list)
{
  int i;
  if(list==NULL) return;
  if(list->str!=NULL){
    for(i=0; i<strings_size(list); i++){ list->str[i]=char_del(list->str[i]); }
    free(list->str);
    list->str=NULL;
    list->n=0;
  }
}

strings* strings_del(strings *list)
{
  if(list==NULL) return NULL;
#ifdef DEBUG
  print_dark_gray();
  printf("strings_del id=%010d n=%d",list->id,list->n);
  print_reset();
  printf("\n");
#endif
  strings_del_str(list);
  free(list); list=NULL;
  return NULL;
}

///////////////////////////////////////////////////////

void strings_reset(strings *list, int n, char **str)
{
  if(list==NULL) return;
  strings_del_str(list);
  list->str=str;
  list->n=n; 
}

void strings_resize(strings *list, int n)
{
  int i;
  char **str=NULL;
  if(list==NULL) return;
  if(list->n==n) return;
  str=(char**)malloc(sizeof(char*)*n);
  if(list->n<n){
    for(i=0; i<n; i++){
      if(i<(list->n)){ str[i]=list->str[i]; list->str[i]=NULL; }
      else           { str[i]=NULL; }
    }
  }else{
    for(i=0; i<list->n; i++){
      if(i<n){ str[i]=list->str[i]; list->str[i]=NULL; }
      else   { list->str[i]=char_del(list->str[i]); }
    }
  }
  strings_reset(list,n,str);
}

/////////////////////////////////////////////////////////

int strings_index(strings *list, char *str)
{
  int i;
  if(list==NULL || (list->n)<=0 || (list->str)==NULL){ return -1; }
  for(i=0; i<list->n; i++){
    if(STR_EQ(list->str[i],str)){
      return i;
    }
  }
  return -1;
}

/////////////////////////////////////////////////////////////

void strings_print(strings *list)
{
  int i;
  if(list==NULL) {
    printf("NULL\n");
    return;
  }
  printf("{");
  for(i=0; i<list->n; i++){
    if(list->str[i]==NULL){ printf("NULL"); }
    else                  { printf("'%s'",list->str[i]); }
    if(i<list->n-1){ printf(", "); }
  }
  printf("}");  
}

void strings_put(strings *list, char *name)
{
  int i;
  if(list==NULL || list->str==NULL || list->n<=0) {
    if(name!=NULL) printf("%s=",name);
    printf("NULL\n");
    return;
  }
  for(i=0; i<list->n; i++){
    if(list->str[i]==NULL){
      //      if(name==NULL) printf("NULL\n");
      //      else printf("%s[%d]=NULL\n",name,i);
    }else{
      if(name==NULL) printf("%s\n",list->str[i]);
      else printf("%s[%d]=%s\n",name,i,list->str[i]);
    }
  }
}

/////////////////////////////////////////////////////////////

void strings_item_del(strings *list, int n)
{
  if(list==NULL) return;
  if(list->str==NULL) return;
  if(n<0 || n>=(list->n)) return;
  if(list->str[n]==NULL) return;
  list->str[n]=char_del(list->str[n]);
}

void strings_item_set(strings *list, int n, char *str, char *skip)
{
  if(list==NULL) return;
  if(n<0) return;
  strings_item_del(list,n);
  if((n+1)>=(list->n)){
    strings_resize(list,n+1);
  }
  list->str[n]=char_del(list->str[n]);
  list->str[n]=char_new(str,skip);
}

void strings_item_replace(strings *list, int k, strings *list2)
{
  int n,i,j,l;
  char **str=NULL;
  if(list==NULL) return;
  if(list2==NULL) return;
  n=(list->n)+(list2->n);
  if(k>=0 && k<list->n){ n--; }
  str=(char**)malloc(sizeof(char*)*n);
  if(k<0){ k=-1; }
  if(k>=list->n){ k=list->n; }
  if(k<0) j=-1; else j=0;
  for(i=0,l=0; i<n; i++){
    if(j==k){
      str[i]=list2->str[l];
      list2->str[l]=NULL;
      l++;
      if(l>=list2->n && j<list->n){
	if(j>=0){ strings_item_del(list,j); }
	j++;
      }
    }else{
      str[i]=list->str[j];
      list->str[j]=NULL;
      j++;
    }
  }
  strings_reset(list,n,str);
  list2=strings_del(list2);
}

/////////////////////////////////////////////////////////////

strings* strings_split_path(char *str, char *suffix)
{
  strings *list=NULL;
  char *s=NULL,*p=NULL;
  // allocate
  list=strings_new(3);
  s=char_new(str,NULL);
  // jump to the end
  p=s; while((*p)!='\0'){ p++; } p--;
  // search
  while((*p)!='.' && (*p)!='/' && p!=s){ p--; }
  // found suffix
  if((*p)=='.'){
    (*p)='\0'; list->str[2]=char_new(p+1,NULL);   // suffix
    // search
    while((*p)!='/' && p!=s){ p--; }
    if((*p)=='/'){
      (*p)='\0'; list->str[1]=char_new(p+1,NULL); // basename
      list->str[0]=char_new(s,NULL);              // dir
    }else{
      list->str[1]=char_new(s,NULL);   // basename
      list->str[0]=char_new("",NULL);  // dir
    }
  }else if((*p)=='/'){
    list->str[2]=char_new("",NULL);              // suffix
    (*p)='\0'; list->str[1]=char_new(p+1,NULL);  // basename
    list->str[0]=char_new(s,NULL);               // dir
  }else{
    list->str[2]=char_new("",NULL); // suffix
    list->str[1]=char_new(s,NULL);  // basename
    list->str[0]=char_new("",NULL); // dir
  }
  // check suffix
  if(suffix!=NULL && !char_eq(list->str[2],suffix)){
    list->str[1]=char_cat(list->str[1],".");
    list->str[1]=char_cat(list->str[1],list->str[2]);
    list->str[2]=char_renew(list->str[2],"",NULL);
  }
  // done
  s=char_del(s);
  return list;
}

/////////////////////////////////////////////////////////////

strings* strings_split(char *str, char *sep, char *mask_begin, char *mask_end, char *skip)
{
  int i,len,n,flag;
  char *s,*t,*m=NULL;
  strings *list;
  // copy
  len=(int)strlen(str);
  s=char_new(str,skip);
  m=str_create_mask(str,mask_begin,mask_end);
  // count tokens
  n=0; flag=0;
  for(i=0; i<len; i++){
    if(m[i]=='.' && char_match_any(s[i],sep,NULL)){ s[i]='\0'; flag=0; }
    else{ if(flag==0){ flag=1; n++; } }
  }
  // allocate list
  list=strings_new(n);
  // copy strings
  t=s;
  for(i=0; i<list->n; i++){
    while((*t)=='\0') t++;
    strings_item_set(list,i,t,skip);
    while((*t)!='\0') t++;
  }
  // free
  s=char_del(s);
  m=char_del(m);
  // done
  return list;
}

strings* strings_split_mask(char *str, char *mask_begin, char *mask_end, char *skip)
{
  int i,len,n,flag;
  char *s,*t,*m=NULL;
  strings *list;
  // copy
  len=(int)strlen(str);
  s=char_new(str,skip);
  m=str_create_mask(str,mask_begin,mask_end);
  // count tokens
  n=0; flag=0;
  for(i=0; i<len; i++){
    if(m[i]=='(' || m[i]==')'){ s[i]='\0'; flag=0; }
    else{ if(flag==0){ flag=1; n++; } }
  }
  // allocate list
  list=strings_new(n);
  // copy strings
  t=s;
  for(i=0; i<list->n; i++){
    while((*t)=='\0') t++;
    strings_item_set(list,i,t,skip);
    while((*t)!='\0') t++;
  }
  // free
  s=char_del(s);
  m=char_del(m);
  // done
  return list;
}

///////////////////////////////////

strings* strings_split_number(char *s)
{
  strings *list=NULL;
  int done,mode,*p,i,j,n,c,w,q;
  char *t;
  // check mode
  // メモリの長さ
  n=(int)strlen(s)+1;
  // モードの保存用
  p=ivec_allocate(n);
  // 浮動小数点数=[+-][数字][.][数字][e][+-][数字]
  // mode=0: なし
  // mode=1: 符号
  // mode=2: 整数部の数字
  // mode=3: 小数点'.'
  // mode=4: 小数部の数字
  // mode=5: 指数部の'e'
  // mode=6: 指数部の+-
  // mode=7: 指数部の数字
  for(done=0, mode=0, i=0; !done; i++){
         if(mode==0 && char_match_any(s[i],STR_SIGN,  NULL)){ mode=1; }
    else if(mode==0 && char_match_any(s[i],STR_DIGITS,NULL)){ mode=2; }
    else if(mode==0 && char_match_any(s[i],STR_POINT, NULL)){ mode=3; }
    else if(mode==1 && char_match_any(s[i],STR_DIGITS,NULL)){ mode=2; }
    else if(mode==1 && char_match_any(s[i],STR_POINT, NULL)){ mode=3; }
    else if(mode==2 && char_match_any(s[i],STR_DIGITS,NULL)){ mode=2; }
    else if(mode==2 && char_match_any(s[i],STR_POINT, NULL)){ mode=3; }
    else if(mode==2 && char_match_any(s[i],STR_EXP,   NULL)){ mode=5; }
    else if(mode==3 && char_match_any(s[i],STR_DIGITS,NULL)){ mode=4; }
    else if(mode==3 && char_match_any(s[i],STR_EXP,   NULL)){ mode=5; }
    else if(mode==4 && char_match_any(s[i],STR_DIGITS,NULL)){ mode=4; }
    else if(mode==4 && char_match_any(s[i],STR_EXP,   NULL)){ mode=5; }
    else if(mode==5 && char_match_any(s[i],STR_SIGN,  NULL)){ mode=6; }
    else if(mode==5 && char_match_any(s[i],STR_DIGITS,NULL)){ mode=7; }
    else if(mode==6 && char_match_any(s[i],STR_DIGITS,NULL)){ mode=7; }
    else if(mode==7 && char_match_any(s[i],STR_DIGITS,NULL)){ mode=7; }
    else                                                    { mode=0; }
    if(mode==0 && char_match_any(s[i],STR_SIGN,  NULL)){ mode=1; }
    else if(mode==0 && char_match_any(s[i],STR_DIGITS,NULL)){ mode=2; }
    else if(mode==0 && char_match_any(s[i],STR_POINT, NULL)){ mode=3; }
    if(s[i]=='\0'){ done=1; }
    p[i]=mode;
  }
  // count letters and words
  for(q=0, c=0, w=0, i=0; i<n; i++){
    // 文字数のカウント mode=1〜7
    if(p[i]>0){ c++; }
         if(p[i]>q)          { q=p[i];      }
    else if(p[i]>0 && p[i]<q){ q=p[i]; w++; }
  }
  w++;
  // split strings
  t=malloc(sizeof(char)*(c+w));
  for(i=0, j=0, q=0; i<n; i++){
         if(p[i]>q)          { q=p[i];              }
    else if(p[i]>0 && p[i]<q){ q=p[i]; t[j++]='\0'; }
         if(p[i]>0)    { t[j++]=s[i]; }
    else if(s[i]=='\0'){ t[j++]='\0'; }
  }
  // allocate list
  list=strings_new(w);
  i=0; j=0;
  if(t[i]!='\0'){
    strings_item_set(list,j,&t[i],NULL);
    j++;
  }
  for(i=0, q=0; i<c+w; i++){
    if(t[i]=='\0'){ q=1; }
    if(t[i]!='\0' && q){
      strings_item_set(list,j,&t[i],NULL);
      j++;
      q=0;
    }
  }
  // free
  p=ivec_free(p);
  free(t);
  // done
  return list;
}

///////////////////////////////////


int strings_cmp(strings *f, strings *g)
{
  int i,value;
  for(i=0; i<MIN2(f->n,g->n); i++){
    value=char_cmp(f->str[i],g->str[i]);
    if(value){ return value; }
  }
  if     (f->n < g->n){ return -1; }
  else if(f->n > g->n){ return +1; }
  else                { return  0; }
}

///////////////////////////////////

strings *strings_split_str_to_irmulti(char *str)
{
  int i;
  strings *a=NULL,*b=NULL,*c=NULL;
  c=strings_new(4);
  a=strings_split(str,", ",NULL,NULL,"()[] \n\t");  
  if(strings_size(a)>=1 && strings_at(a,0)!=NULL){    
    if(STR_EQ(strings_at(a,0),"i")){ strings_item_set(a,0,"1i","()[] \n\t"); }
    b=strings_split_number(strings_at(a,0));
    if(strings_size(b)>=1 && strings_at(b,0)!=NULL){ strings_item_set(c,0,strings_at(b,0),"()[] \n\t"); }
    if(strings_size(b)>=2 && strings_at(b,1)!=NULL){ strings_item_set(c,1,strings_at(b,1),"()[] \n\t"); }
    b=strings_del(b);
    if(str_has_any_char(strings_at(a,0),"iIjJ") && strings_at(c,0)!=NULL && strings_at(c,1)==NULL){
      strings_item_set(c,1,strings_at(c,0),"()[] \n\t");
      strings_item_del(c,0);
    }    
  }
  if(strings_size(a)>=2 && strings_at(a,1)!=NULL){
    if(STR_EQ(strings_at(a,1),"i")){ strings_item_set(a,1,"1i","()[] \n\t"); }
    b=strings_split_number(strings_at(a,1));
    if(strings_size(b)>=1 && strings_at(b,0)!=NULL){ strings_item_set(c,2,strings_at(b,0),"()[] \n\t"); }
    if(strings_size(b)>=2 && strings_at(b,1)!=NULL){ strings_item_set(c,3,strings_at(b,1),"()[] \n\t"); }
    b=strings_del(b);
    if(str_has_any_char(strings_at(a,1),"iIjJ") && strings_at(c,2)!=NULL && strings_at(c,3)==NULL){
      strings_item_set(c,3,strings_at(c,2),"()[] \n\t");
      strings_item_del(c,2);
    }    
  }
  for(i=0; i<4; i++){
    if(strings_at(c,i)!=NULL && STR_EQ(strings_at(c,i),"+")){ strings_item_set(c,i,"1","()[] \n\t"); }
    if(strings_at(c,i)!=NULL && STR_EQ(strings_at(c,i),"-")){ strings_item_set(c,i,"-1","()[] \n\t"); }
  }
  if(strings_at(c,0)!=NULL && strings_at(c,1)==NULL){ strings_item_set(c,1,"0","()[] \n\t"); }
  if(strings_at(c,2)!=NULL && strings_at(c,3)==NULL){ strings_item_set(c,3,"0","()[] \n\t"); }
  if(strings_at(c,0)==NULL && strings_at(c,1)!=NULL){ strings_item_set(c,0,"0","()[] \n\t"); }
  if(strings_at(c,2)==NULL && strings_at(c,3)!=NULL){ strings_item_set(c,2,"0","()[] \n\t"); }
  a=strings_del(a);
  return c;
}


//EOF
