#ifndef IS_PRINT_H
#define IS_PRINT_H

enum { COLOR_RESET=0, COLOR_BLACK, COLOR_DARK_GRAY, COLOR_LIGHT_GRAY, COLOR_WHITE,
       COLOR_RED,    COLOR_LIGHT_RED,  COLOR_GREEN,  COLOR_LIGHT_GREEN, COLOR_BROWN,
       COLOR_BLUE,   COLOR_LIGHT_BLUE, COLOR_PURPLE, COLOR_LIGHT_PURPLE,
       COLOR_CYAN,   COLOR_LIGHT_CYAN, COLOR_YELLOW
};

void print_reset(void);
void print_black(void);
void print_red(void);
void print_green(void);
void print_brown(void);
void print_blue(void);
void print_purple(void);
void print_cyan(void);
void print_light_gray(void);
void print_dark_gray(void);
void print_light_red(void);
void print_light_green(void);
void print_yellow(void);
void print_light_blue(void);
void print_light_purple(void);
void print_light_cyan(void);
void print_white(void);
void print_color(int color);

void fprint_ansi(FILE *fp, int d);
void fprint_ansi2(FILE *fp, int d);
void fprint_reset(FILE *fp);
void fprint_black(FILE *fp);
void fprint_red(FILE *fp);
void fprint_green(FILE *fp);
void fprint_brown(FILE *fp);
void fprint_blue(FILE *fp);
void fprint_purple(FILE *fp);
void fprint_cyan(FILE *fp);
void fprint_light_gray(FILE *fp);
void fprint_dark_gray(FILE *fp);
void fprint_light_red(FILE *fp);
void fprint_light_green(FILE *fp);
void fprint_yellow(FILE *fp);
void fprint_light_blue(FILE *fp);
void fprint_light_purple(FILE *fp);
void fprint_light_cyan(FILE *fp);
void fprint_white(FILE *fp);
void fprint_color(FILE *fp, int color);

int file_can_read(char* fmt, ...);

#endif
