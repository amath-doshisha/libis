#ifndef IS_SVEC_H
#define IS_SVEC_H


char **svec_allocate(int n);
char **svec_free(int n, char **x);
void svec_copy(int n, char **y, char **x); // y=x
void svec_set_si(int n, char **y, int *x); // y=x
void svec_set_all_d(int n, char **y, double x, char *format); // y=ones(n,1)*x
void svec_set_all(int n, char **y, char *x); // y=ones(n,1)*x

int svec_has(int n, char **x, char *match);

  int svec_max_length(int n, char **x); // y=max(length(x))
int smat_max_length_row(int n, char **A, int LDA, int k); // y=max(length(A(k,:))
int smat_max_length_column(int m, char **A, int LDA, int k); // y=max(length(A(:,k))

void svec_put(int n, char **x, char *sep);
void svec_print(int n, char **x, const char *name);
void smat_print(int m, int n, char **A, int LDA, const char *name);


#endif
