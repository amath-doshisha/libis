#ifndef IS_STRINGS_H
#define IS_STRINGS_H

#include<stdarg.h>

typedef struct {
  int id;
  char **str;
  int n;
} strings;

strings* strings_new(int n);
strings* strings_new_str(const char *str[], const char *skip);
strings* strings_new_clone(strings *src);
int strings_size(strings *str);
char *strings_at(strings *str, int i);
void strings_del_str(strings *list);
strings* strings_del(strings *list);
void strings_reset(strings *list, int n, char **str);
void strings_resize(strings *list, int n);
int strings_index(strings *list, const char *str);
void strings_print(strings *list);
void strings_put(strings *list, char *name);
void strings_item_del(strings *list, int n);
void strings_item_set(strings *list, int n, const char *str, const char *skip);
void strings_item_replace(strings *list, int k, strings *list2);
strings* strings_split(const char *str, const char *sep, const char *mask_start, const char *mask_end, const char *skip);
strings* strings_split_mask(const char *str, const char *mask_begin, const char *mask_end, const char *skip);
strings* strings_split_path(const char *str, const char *suffix);
int strings_cmp(strings *f, strings *g);
char *str_create_mask(const char *s, const char *mask_begin, const char *mask_end);

char *char_new(const char *str, const char *skip);
char *char_del(char *p);
char *char_clone(const char *str);
char *char_renew(char *old_str, const char *str, const char *skip);
char *char_renew_sprintf(char *old, char *skip, char* fmt, ...);
char *char_cat(char *str_old, const char *str_append);
int char_cmp(const char *a, const char *b);
int char_eq(const char *a, const char *b);
int char_match_any(char c, const char *str, int *position);
int str_match(const char *str, const char *pattern, const char *mask_start, const char *mask_end);

#endif
