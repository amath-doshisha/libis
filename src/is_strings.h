#ifndef IS_STRINGS_H
#define IS_STRINGS_H

#include<stdarg.h>

typedef struct {
  int id;
  char **str;
  int n;
} strings;

strings* strings_new(int n);
strings* strings_new_str(char *str[], char *skip);
strings* strings_new_clone(strings *src);
int strings_size(strings *str);
char *strings_at(strings *str, int i);
void strings_del_str(strings *list);
strings* strings_del(strings *list);
void strings_reset(strings *list, int n, char **str);
void strings_resize(strings *list, int n);
int strings_index(strings *list, char *str);
void strings_print(strings *list);
void strings_put(strings *list, char *name);
void strings_item_del(strings *list, int n);
void strings_item_set(strings *list, int n, char *str, char *skip);
void strings_item_replace(strings *list, int k, strings *list2);
strings* strings_split(char *str, char *sep, char *mask_start, char *mask_end, char *skip);
strings* strings_split_mask(char *str, char *mask_begin, char *mask_end, char *skip);
strings* strings_split_path(char *str, char *suffix);
strings* strings_split_number(char *s);
int strings_cmp(strings *f, strings *g);
char *str_create_mask(char *s, char *mask_begin, char *mask_end);
int str_has_any_char(char *s, char *c);

char *char_new(char *str, char *skip);
char *char_del(char *p);
char *char_clone(char *str);
char *char_renew(char *old_str, char *str, char *skip);
char *char_renew_sprintf(char *old, char *skip, char* fmt, ...);
char *char_cat(char *str_old, char *str_append);
int char_cmp(char *a, char *b);
int char_eq(char *a, char *b);
int char_match_any(char c, char *str, int *position);
int str_match(char *str, char *pattern, char *mask_start, char *mask_end);

strings *strings_split_str_to_irmulti(char *str);

#endif
