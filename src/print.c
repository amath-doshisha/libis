#include<stdio.h>
#include<stdlib.h>
#include<stdarg.h>
#include"is_print.h"

void print_reset()        { fprint_reset(stdout); }
void print_black()        { fprint_black(stdout); }
void print_red()          { fprint_red(stdout); }
void print_green()        { fprint_green(stdout); }
void print_brown()        { fprint_brown(stdout); }
void print_blue()         { fprint_blue(stdout); }
void print_purple()       { fprint_purple(stdout); }
void print_cyan()         { fprint_cyan(stdout); }
void print_light_gray()   { fprint_light_gray(stdout); }
void print_dark_gray()    { fprint_dark_gray(stdout); }
void print_light_red()    { fprint_light_red(stdout); }
void print_light_green()  { fprint_light_green(stdout); }
void print_yellow()       { fprint_yellow(stdout); }
void print_light_blue()   { fprint_light_blue(stdout); }
void print_light_purple() { fprint_light_purple(stdout); }
void print_light_cyan()   { fprint_light_cyan(stdout); }
void print_white()        { fprint_white(stdout); }
void print_color(int color){ fprint_color(stdout,color); }

void fprint_ansi(FILE *fp, int d)  { fprintf(fp,"\033[%dm",d); }
void fprint_ansi2(FILE *fp, int d) { fprintf(fp,"\033[1;%dm",d); }
void fprint_reset(FILE *fp)        { fprint_ansi(fp,0);  }
void fprint_black(FILE *fp)        { fprint_ansi(fp,30); }
void fprint_red(FILE *fp)          { fprint_ansi(fp,31); }
void fprint_green(FILE *fp)        { fprint_ansi(fp,32); }
void fprint_brown(FILE *fp)        { fprint_ansi(fp,33); }
void fprint_blue(FILE *fp)         { fprint_ansi(fp,34); }
void fprint_purple(FILE *fp)       { fprint_ansi(fp,35); }
void fprint_cyan(FILE *fp)         { fprint_ansi(fp,36); }
void fprint_light_gray(FILE *fp)   { fprint_ansi(fp,37); }
void fprint_dark_gray(FILE *fp)    { fprint_ansi2(fp,30); }
void fprint_light_red(FILE *fp)    { fprint_ansi2(fp,31); }
void fprint_light_green(FILE *fp)  { fprint_ansi2(fp,32); }
void fprint_yellow(FILE *fp)       { fprint_ansi2(fp,33); }
void fprint_light_blue(FILE *fp)   { fprint_ansi2(fp,34); }
void fprint_light_purple(FILE *fp) { fprint_ansi2(fp,35); }
void fprint_light_cyan(FILE *fp)   { fprint_ansi2(fp,36); }
void fprint_white(FILE *fp)        { fprint_ansi2(fp,37); }

void fprint_color(FILE *fp, int color)
{
  switch(color){
  case COLOR_RESET:             fprint_reset(fp);        break;
  case COLOR_BLACK:             fprint_black(fp);        break;
  case COLOR_RED:               fprint_red(fp);          break;
  case COLOR_GREEN:             fprint_green(fp);        break;
  case COLOR_BROWN:             fprint_brown(fp);        break;
  case COLOR_BLUE:              fprint_blue(fp);         break;
  case COLOR_PURPLE:            fprint_purple(fp);       break;
  case COLOR_CYAN:              fprint_cyan(fp);         break;
  case COLOR_LIGHT_GRAY:        fprint_light_gray(fp);   break;
  case COLOR_DARK_GRAY:         fprint_dark_gray(fp);    break;
  case COLOR_LIGHT_RED:         fprint_light_red(fp);    break;
  case COLOR_LIGHT_GREEN:       fprint_light_green(fp);  break;
  case COLOR_YELLOW:            fprint_yellow(fp);       break;
  case COLOR_LIGHT_BLUE:        fprint_light_blue(fp);   break;
  case COLOR_LIGHT_PURPLE:      fprint_light_purple(fp); break;
  case COLOR_LIGHT_CYAN:        fprint_light_cyan(fp);   break;
  case COLOR_WHITE:             fprint_white(fp);        break;
  default:                                               break;
  }
}

///////////////////////////////////////////////////////////////////////////

int file_can_read(char* fmt, ...)
{
  int value=0;
  char fname[1025];
  va_list argp;
  FILE *fid;
  // file name
  va_start(argp,fmt);
  vsprintf(fname,fmt,argp);
  // open file
  if((fid=fopen(fname,"r"))==0){ value=0; }
  else{ value=1; }
  // close file
  fclose(fid);
  // done
  return value;
}

//EOF
