#ifndef ISYS_FUNC_H
#define ISYS_FUNC_H

//#define DEBUG

#include<is_strings.h>
#include<is_rmulti.h>
#include<is_cmulti.h>
#include<is_bigint.h>


typedef struct { int n,*num,*pow; } func_var_struct;
typedef struct { int n,*x; } func_ivec_struct;
typedef struct { int n; rmulti **x; } func_rvec_struct;
typedef struct { int n; cmulti **x; } func_cvec_struct;
typedef struct { int LD,m,n; rmulti **A; } func_rmat_struct;
typedef struct { int LD,m,n; cmulti **A; } func_cmat_struct;

enum { FUNC_P_NULL=0, FUNC_P_BUILTIN, FUNC_P_DEF, FUNC_P_POWER,
       FUNC_P_BIGINT, FUNC_P_REAL, FUNC_P_COMPLEX, FUNC_P_STRINGS,
       FUNC_P_VAR, FUNC_P_IVEC, FUNC_P_RVEC, FUNC_P_CVEC, FUNC_P_RMAT, FUNC_P_CMAT,
       FUNC_P_SIZE };

enum { FUNC_NO=0, FUNC_YES=1, FUNC_NOLIMIT=-1 };
enum { FUNC_SCALAR_NO=0, FUNC_SCALAR=1, FUNC_SCALAR_DEPEND, FUNC_COMMAND };
enum { FUNC_ORDER_DEFAULT=-1, FUNC_ORDER_NULL=0, FUNC_ORDER_null, FUNC_ORDER_NAN, FUNC_ORDER_INF,
       FUNC_ORDER_ZERO, FUNC_ORDER_ONE, FUNC_ORDER_BIGINT, FUNC_ORDER_REAL, FUNC_ORDER_COMPLEX,
       FUNC_ORDER_VAR, FUNC_ORDER_ADD, FUNC_ORDER_MUL,
       FUNC_ORDER_SQRT, FUNC_ORDER_EXP, FUNC_ORDER_LOG, FUNC_ORDER_POW,
       FUNC_ORDER_SIN, FUNC_ORDER_COS, FUNC_ORDER_TAN,
       FUNC_ORDER_ASIN, FUNC_ORDER_ACOS, FUNC_ORDER_ATAN,
       FUNC_ORDER_SINH, FUNC_ORDER_COSH, FUNC_ORDER_TANH,
       FUNC_ORDER_ASINH, FUNC_ORDER_ACOSH, FUNC_ORDER_ATANH,
       FUNC_ORDER_IVEC, FUNC_ORDER_RVEC, FUNC_ORDER_CVEC, FUNC_ORDER_RMAT, FUNC_ORDER_CMAT,
       FUNC_ORDER_LIST, FUNC_ORDER_TABLE, FUNC_ORDER_STRINGS,
       FUNC_ORDER_SCOPE, FUNC_ORDER_DEF, FUNC_ORDER_BUILTIN };

union func_p_struct {
  void *mem;
  struct func_builtin_struct *builtin;
  struct func_def_struct *def;
  bigint *bi;
  rmulti *rm;
  cmulti *cm;
  int *i;
  strings *s;
  func_var_struct *var;
  func_ivec_struct *ivec;
  func_rvec_struct *rvec;
  func_cvec_struct *cvec;
  func_rmat_struct *rmat;
  func_cmat_struct *cmat;
};
typedef union func_p_struct func_p_t;

struct func_struct {
  int count;
  char *op;
  int n;
  struct func_struct **a;
  int ptype;
  union func_p_struct p;
};
typedef struct func_struct func_t;

typedef void (func_set_t)(func_t *f);
typedef void (func_set2_t)(func_t *f, func_t *g);
typedef func_t *(func_arg1_t)(func_t *f);
typedef int (func_cmp_t)(func_t *f, func_t *g);
typedef int (func_is_t)(func_t *f);

struct func_builtin_struct {
  char *name;
  int order;
  int amin;
  int amax;
  int type;
  int ptype;
  func_set_t *p_new,*p_del;
  func_set2_t *p_clone;
  func_cmp_t *p_cmp;
  func_arg1_t *eval;
  func_set_t *print;
};
typedef struct func_builtin_struct func_builtin_t;

struct func_def_struct {
  char *name;
  int amin;
  int amax;
};
typedef struct func_def_struct func_def_t;

// general command
int func_size(func_t *f);
int func_rows(func_t *f);
int func_cols(func_t *f);

// built-in
int func_is_builtin(func_t *f);
func_t *func_builtin_new(const char *name);
void func_builtin_p_new(func_t *f);
void func_builtin_p_del(func_t *f);
void func_builtin_p_clone(func_t *f, func_t *g);
int func_builtin_p_cmp(func_t *f, func_t *g);
func_builtin_t *func_builtin_struct_new(void);
func_builtin_t *func_builtin_struct_del(func_builtin_t *op);
func_builtin_t *func_builtin_p(func_t *f);
char *func_builtin_name(func_t *f);
int func_builtin_order(func_t *f);
int func_builtin_amin(func_t *f);
int func_builtin_amax(func_t *f);
int func_builtin_type(func_t *f);
int func_builtin_ptype(func_t *f);
func_set_t *func_builtin_p_new_call(func_t *f);
func_set_t *func_builtin_p_del_call(func_t *f);
func_set2_t *func_builtin_p_clone_call(func_t *f);
func_arg1_t *func_builtin_eval_call(func_t *f);
func_set_t *func_builtin_print_call(func_t *f);
func_cmp_t *func_builtin_cmp_call(func_t *f);
void func_builtin_print(func_t *f);
func_t *func_builtin_script(func_t *f);

// op-list
func_t *func_op_builtin_new(void);
func_t *func_op_list_new(void);
func_t *func_op_table_new(void);
func_t *func_op_strings_new(void);
func_t *func_op_scope_new(void);
func_t *func_op_nan_new(void);
func_t *func_op_inf_new(void);
func_t *func_op_zero_new(void);
func_t *func_op_one_new(void);
func_t *func_op_bigint_new(void);
func_t *func_op_real_new(void);
func_t *func_op_complex_new(void);
func_t *func_op_var_new(void);
func_t *func_op_add_new(void);
func_t *func_op_mul_new(void);
func_t *func_op_sqrt_new(void);
func_t *func_op_exp_new(void);
func_t *func_op_log_new(void);
func_t *func_op_pow_new(void);
func_t *func_op_sin_new(void);
func_t *func_op_cos_new(void);
func_t *func_op_tan_new(void);
func_t *func_op_asin_new(void);
func_t *func_op_acos_new(void);
func_t *func_op_atan_new(void);
func_t *func_op_sinh_new(void);
func_t *func_op_cosh_new(void);
func_t *func_op_tanh_new(void);
func_t *func_op_asinh_new(void);
func_t *func_op_acosh_new(void);
func_t *func_op_atanh_new(void);
func_t *func_op_ivec_new(void);
func_t *func_op_rvec_new(void);
func_t *func_op_cvec_new(void);
func_t *func_op_rmat_new(void);
func_t *func_op_cmat_new(void);
func_t *func_op_def_new(void);
func_t *func_op_begin_new(void);
func_t *func_op_end_new(void);
func_t *func_op_print_new(void);
func_t *func_op_eval_new(void);
func_t *func_op_evalf_new(void);
func_t *func_op_expand_new(void);
func_t *func_op_diff_new(void);
func_t *func_op_grad_new(void);
func_t *func_op_gbasis_new(void);

// op
int func_ptype(func_t *f);
const char *func_op(func_t *f);
int func_is(func_t *f, const char *op);

// find-op
int func_find_ptype(func_t *f);
int func_find_amin(func_t *f);
int func_find_amax(func_t *f);
int func_find_order(func_t *f);
int func_find_type(func_t *f);
func_set_t *func_find_p_new(func_t *f);
func_set_t *func_find_p_del(func_t *f);
func_set2_t *func_find_p_clone(func_t *f);
func_cmp_t *func_find_p_cmp(func_t *f);
func_arg1_t *func_find_eval(func_t *f);
func_set_t *func_find_print(func_t *f);

// init, import, and clear
void func_init(void);
void func_clear(void);
void func_import_basic(int n);
func_t *func_begin(func_t *f);
func_t *func_begin_eval(func_t *f);
func_t *func_end(func_t *f);
func_t *func_end_eval(func_t *f);
func_t *func_set(func_t *f);
func_t *func_set_eval(func_t *f);
func_t *func_set_script(func_t *f);
func_t *func_op_set_new(void);

// scope
int func_is_scope(func_t *f);
void func_scope_begin(func_t *vlist);
void func_scope_end(void);
func_t *func_scope(int n);
func_t *func_scope_previous(int n);
func_t *func_scope_vlist(int n);
func_t *func_scope_table(int n);
void func_scope_set(int n, func_t *f);
func_t *func_scope_find(int n, const char *name);

// scope (private)
func_t *func_scope_new(func_t *previous, func_t *vlist);
void func_scope_init_vlist_to_table(int n);

// monomial order
enum { FUNC_MONO_LEX=0, FUNC_MONO_GRLEX, FUNC_MONO_GREVLEX, FUNC_MONO_ORDER_SIZE };
void func_set_mono_order(int order);
int func_get_mono_order(void);

// allocation
int func_new_del_check_sum();
func_t *func_new(const char *op);
func_t *func_retain(func_t *f);
func_t *func_del(func_t *f);
func_t *func_clone(func_t *f);
func_t *func_replace(func_t *f, func_t *g);

// allocation (private)
func_t *func_struct_new(void);
func_t *func_struct_del(func_t *f);

// allocation (private)
void func_p_new(func_t *f);
void func_p_del(func_t *f);
void func_p_clone(func_t *f, func_t *g);
func_t *func_arg1_new(const char *op, func_t *g);
func_t *func_arg2_new(const char *op, func_t *g0, func_t *g1);

// script
func_t *func_script(const char *str);
func_t *func_scriptf(char* fmt, ...);
func_t *func_scriptf_echo(char* fmt, ...);
func_t *func_null_script(func_t *f);
func_t *func_bracket_script(func_t *f);

// power
int func_has_power(func_t *f);
int func_power(func_t *f);
void func_set_power(func_t *f, int pow);
int func_power_cmp(func_t *f, func_t *g);

// power (private)
void func_pow_p_new(func_t *f);
void func_pow_p_del(func_t *f);
void func_pow_p_clone(func_t *f, func_t *g);

// nan
func_t *func_nan(void);
int func_is_nan(func_t *f);
void func_nan_print(func_t *f);
func_t *func_nan_script(func_t *f);

// inf
func_t *func_inf(void);
int func_is_inf(func_t *f);
void func_inf_print(func_t *f);
func_t *func_inf_script(func_t *f);

// zero
func_t *func_zero(void);
int func_is_zero(func_t *f);
void func_zero_print(func_t *f);
func_t *func_zero_script(func_t *f);

// one
func_t *func_one(void);
int func_is_one(func_t *f);
void func_one_print(func_t *f);
func_t *func_one_script(func_t *f);

// big integer
void func_bigint_print(func_t *f);
func_t *func_bigint_script(func_t *f);
func_t *func_bigint(void);
func_t *func_bigint_int(int num, int den);
func_t *func_bigint_str(const char *num, const char *den);
int func_bigint_is_integer(func_t *f);
int func_bigint_sgn(func_t *f);
int func_bigint_get_si(func_t *f);
func_t *func_bigint_add(func_t *f1, func_t *f2);
func_t *func_bigint_sub(func_t *f1, func_t *f2);
func_t *func_bigint_mul(func_t *f1, func_t *f2);
func_t *func_bigint_div(func_t *f1, func_t *f2);
func_t *func_bigint_inv(func_t *f1);
func_t *func_bigint_pow_n(func_t *f, int p);
int func_bigint_cmp(func_t *f1, func_t *f2);
int func_is_bigint(func_t *f);
int func_in_bigint(func_t *f);
func_t *func_bigint_cast(func_t *f);

// big integer (private)
void func_bigint_p_new(func_t *f);
void func_bigint_p_del(func_t *f);
void func_bigint_p_clone(func_t *f, func_t *g);
func_t *func_bigint_eval(func_t *f);
bigint *func_bigint_p(func_t *f);

// convertor
func_t *func_evalf(func_t *f);
func_t *func_evalf_eval(func_t *g);
func_t *func_eval(func_t *f);
func_t *func_eval_eval(func_t *g);

// real number
func_t *func_real(void);
func_t *func_real_rmulti(rmulti *x);
func_t *func_real_str(const char *str);
int func_real_cmp(func_t *f, func_t *g);
func_t *func_real_get_bigint(func_t *f);
int func_is_real(func_t *f);
int func_in_real(func_t *f);
func_t *func_real_cast(func_t *g);
void func_real_print(func_t *f);
func_t *func_real_script(func_t *f);

// real number (private)
void func_real_p_new(func_t *f);
void func_real_p_del(func_t *f);
void func_real_p_clone(func_t *f, func_t *g);
rmulti *func_real_p(func_t *f);
func_t *func_real_eval(func_t *f);
func_t *func_real_add(func_t *f1, func_t *f2);
func_t *func_real_sub(func_t *f1, func_t *f2);
func_t *func_real_mul(func_t *f1, func_t *f2);
func_t *func_real_inv(func_t *f1);
func_t *func_real_div(func_t *f1, func_t *f2);
func_t *func_real_pow_n(func_t *f, int p);
func_t *func_real_pow(func_t *f1, func_t *f2);

// complex nubmer
func_t *func_complex(void);
func_t *func_complex_cmulti(cmulti *x);
int func_complex_cmp(func_t *f, func_t *g);
func_t *func_complex_i(void);
int func_in_complex(func_t *f);
int func_is_complex(func_t *f);
func_t *func_complex_cast(func_t *g);

// complex nubmer (private)
void func_complex_p_new(func_t *f);
void func_complex_p_del(func_t *f);
void func_complex_p_clone(func_t *f, func_t *g);
cmulti *func_complex_p(func_t *f);
func_t *func_complex_eval(func_t *f);
func_t *func_complex_add(func_t *f1, func_t *f2);
func_t *func_complex_sub(func_t *f1, func_t *f2);
func_t *func_complex_mul(func_t *f1, func_t *f2);
func_t *func_complex_inv(func_t *f1);
func_t *func_complex_div(func_t *f1, func_t *f2);
func_t *func_complex_pow_n(func_t *f1, int p);
func_t *func_complex_pow(func_t *f1, func_t *f2);

// number
int func_is_number(func_t *f);
int func_number_cmp(func_t *f, func_t *g);
func_t *func_number_pull(func_t *f);

// var
func_t *func_var(int n);
func_t *func_var1(int var, int pow);
int func_var_size(func_t *f);
int func_var_num(func_t *f, int i);
int func_var_pow(func_t *f, int i);
int func_is_1var(func_t *f);
int func_is_var(func_t *f);
void func_var_print(func_t *f);
func_t *func_var_script(func_t *f);

// var (private)
void func_var_p_new(func_t *f);
void func_var_p_del(func_t *f);
void func_var_p_clone(func_t *f, func_t *g);
void func_var_p_resize(func_t *f, int n);
void func_var_replace(func_t *f, int n, int *num, int *pow);
 void func_var_sort(func_t *f);
 void func_var_mv_pow_mul(func_t *f);
 void func_var_mv_pow_lcm(func_t *f);
 void func_var_rm_pow0(func_t *f);
func_t *func_var_check_mul(func_t *f);
func_t *func_var_check_lcm(func_t *f);
func_t *func_var_mul(func_t *f1, func_t *f2);
func_t *func_var_lcm(func_t *f1, func_t *f2);
func_t *func_var_div(func_t *f1, func_t *f2);
func_t *func_var_pow_n(func_t *f, int p);
  int func_var_degree(func_t *f);
  int func_poly_degree_max(func_t *f);
  int func_degree(func_t *f);
  int func_var_cmp_lex(func_t *f, func_t *g);
  int func_var_cmp_grlex(func_t *f, func_t *g);
  int func_var_cmp_grevlex(func_t *f, func_t *g);
  int func_var_cmp(func_t *f, func_t *g);
  int func_var_can_div(func_t *f, func_t *g);
int func_var_n_var(func_t *f);
int func_var_exist_varn(func_t *f, int n);
int func_var_get_index(func_t *f, int n);
int func_var1n(func_t *f);
int func_var_var1n(func_t *f);

// mono
int func_is_mono_not_num(func_t *f);
int func_is_mono(func_t *f);
int func_mono_degree(func_t *f);

// mono (private)
int func_mono_n_var(func_t *f);
int func_mono_exist_varn(func_t *f, int n);
int func_mono_var1n(func_t *f);
func_t *func_mono_get_var(func_t *f);
func_t *func_mono_get_coeff(func_t *f);
int func_mono_can_div(func_t *f, func_t *g);

// table
int func_is_table(func_t *f);
char *func_key_name(func_t *f);
func_t *func_table(void);
void func_table_set(func_t *f, func_t *g);
func_t *func_table_find(func_t *f, const char *name);

// table (private)
void func_table_grow(func_t *f);
func_t *func_table_bigger(func_t *g);
int func_table_bigger_size(func_t *f);
int func_table_key_index(func_t *f, const char *name);

// list
int func_is_list(func_t *f);
func_t *func_list(int n);
func_t *func_list_zeros(int n);
func_t *func_list_zeros2(int m, int n);
func_t *func_list_concat(func_t *f,func_t *g);
func_t *func_var1_list(int n);
void func_list_print(func_t *f);
func_t *func_list_script(func_t *f);

// list (private)
func_t *func_list_eval(func_t *g);
func_t *func_list_concat(func_t *f,func_t *g);
func_t *func_list_sol_convert_tree(func_t *f);
int func_list_degree_max(func_t *f);
func_t *func_list_find_var1n(func_t *f);
int func_list_real_array_size(func_t *f);
int func_list_complex_array_size(func_t *f);

// poly
int func_is_poly(func_t *f);
int func_is_poly_list(func_t *f);
func_t *func_poly_div_r(func_t *f1, func_t *f2);
func_t *func_poly_div_q(func_t *f1, func_t *f2);
func_t *func_poly_list_div_r(func_t *f, func_t *g);
func_t *func_poly_monic(func_t *f);
func_t *func_poly_list_monic(func_t *f);
func_t *func_poly_list_reduced(func_t *f);
// poly (private)
func_t *func_poly_get_lt(func_t *f);
func_t *func_poly_get_lt_coeff(func_t *f);
func_t *func_poly_get_lc(func_t *f);
func_t *func_poly_get_lm(func_t *f);
func_t *func_poly_get_mono_n(func_t *f,int n);
func_t *func_poly_get_var_n(func_t *f,int n);
func_t *func_poly_get_coeff_n(func_t *f,int n);
func_t *func_poly_get_mono_ntarm(func_t *f,int n);
func_t *func_poly_get_var_ntarm(func_t *f,int n);
func_t *func_poly_get_coeff_ntarm(func_t *f,int n);
int func_poly_n_var(func_t *f);
//////int func_poly_list_n_var(func_t *f);
int func_poly_exist_varn(func_t *f, int n);
int func_poly_var1n(func_t *f);
func_t *func_poly_calculate_soltion(func_t *a,func_t *b);
func_t *func_poly_quadratic_formula(func_t *a,func_t *b,func_t *c);
int func_poly_can_div(func_t *f, func_t *g);
func_t *func_poly_rm_lt(func_t *f);
func_t *func_poly_lm_lcm(func_t *f,func_t *g);
func_t *func_poly_s(func_t *f, func_t *g);
func_t *func_poly_div_r_and_q(func_t *f, func_t *g);
func_t *func_poly_list_buchberger(func_t *f, int debug);
func_t *func_poly_list_buchberger_simple(func_t *f, int debug);
func_t *func_poly_list_buchberger_selected(func_t *f, int debug);
func_t *func_poly_list_minimal_groebner(func_t *f, int debug);
func_t *func_poly_list_reduced_groebner(func_t *f, int debug);
func_t *func_poly_clone_coeff_ntarm(func_t *f, int n);

// groebner basis
func_t *func_gbasis_eval(func_t *g);

// ivec
func_t *func_ivec(int n);
int func_ivec_cmp(func_t *f, func_t *g);
int func_ivec_size(func_t *f);
int *func_ivec_p(func_t *f);
int func_ivec_at(func_t *f, int i);

// ivec (private)
void func_ivec_p_new(func_t *f);
void func_ivec_p_del(func_t *f);
void func_ivec_p_clone(func_t *f, func_t *g);
void func_ivec_p_resize(func_t *f, int n);

// vector, matrix
int func_in_vec(func_t *f);
int func_in_mat(func_t *f);
func_t *func_vec_add(func_t *f1, func_t *f2);
func_t *func_vec_sub(func_t *f1, func_t *f2);
func_t *func_mat_add(func_t *f1, func_t *f2);
func_t *func_mat_sub(func_t *f1, func_t *f2);

// rvec
int func_is_rvec(func_t *f);
func_t *func_rvec(int n);
int func_rvec_cmp(func_t *f, func_t *g);
int func_rvec_size(func_t *f);
rmulti **func_rvec_p(func_t *f);
rmulti *func_rvec_at(func_t *f, int i);
func_t *func_rvec_get(func_t *f, int i);
void func_rvec_set(func_t *f, int i, func_t *g);
func_t *func_rvec_get_cvec(func_t *g);

// rvec (private)
void func_rvec_p_new(func_t *f);
void func_rvec_p_del(func_t *f);
void func_rvec_p_clone(func_t *f, func_t *g);
void func_rvec_p_resize(func_t *f, int n);

// cvec
int func_is_cvec(func_t *f);
func_t *func_cvec(int n);
int func_cvec_cmp(func_t *f, func_t *g);
int func_cvec_size(func_t *f);
cmulti **func_cvec_p(func_t *f);
cmulti *func_cvec_at(func_t *f, int i);
func_t *func_cvec_get(func_t *f, int i);
void func_cvec_set(func_t *f, int i, func_t *g);

// cvec (private)
void func_cvec_p_new(func_t *f);
void func_cvec_p_del(func_t *f);
void func_cvec_p_clone(func_t *f, func_t *g);
void func_cvec_p_resize(func_t *f, int n);

// rmat
int func_is_rmat(func_t *f);
func_t *func_rmat(int m, int n);
int func_rmat_cmp(func_t *f, func_t *g);
int func_rmat_rows(func_t *f);
int func_rmat_cols(func_t *f);
int func_rmat_ld(func_t *f);
rmulti *func_rmat_at(func_t *f, int i, int j);
rmulti **func_rmat_p(func_t *f);
void func_rmat_set_row(func_t *f, int i, func_t *g);
func_t *func_rmat_get_cmat(func_t *g);

// rmat (private)
void func_rmat_p_new(func_t *f);
void func_rmat_p_del(func_t *f);
void func_rmat_p_clone(func_t *f, func_t *g);
void func_rmat_p_resize(func_t *f, int m, int n);

// cmat
int func_is_cmat(func_t *f);
func_t *func_cmat(int m, int n);
int func_cmat_cmp(func_t *f, func_t *g);
int func_cmat_rows(func_t *f);
int func_cmat_cols(func_t *f);
int func_cmat_ld(func_t *f);
cmulti *func_cmat_at(func_t *f, int i, int j);
cmulti **func_cmat_p(func_t *f);
void func_cmat_set_row(func_t *f, int i, func_t *g);

// cmat (private)
void func_cmat_p_new(func_t *f);
void func_cmat_p_del(func_t *f);
void func_cmat_p_clone(func_t *f, func_t *g);
void func_cmat_p_resize(func_t *f, int m, int n);

// add
func_t *func_add_new(int n);
func_t *func_add(func_t *f1, func_t *f2);
func_t *func_add_eval(func_t *f);
int func_is_add(func_t *f);
void func_add_print(func_t *f);
func_t *func_add_script(func_t *f);

// add (private)
void func_add_args(func_t *f, func_is_t *fis, func_is_t *fin);
void func_add_args_collect(func_t *f, func_is_t *fis);

// mul
func_t *func_mul_new(int n);
func_t *func_mul(func_t *f1, func_t *f2);
func_t *func_mul_eval(func_t *f);
func_t *func_mul_eval_size(func_t *f);
func_t *func_mul_split(func_t *f);
func_t *func_mul_pow_n(func_t *f, int pow);
int func_is_mul(func_t *f);
void func_mul_print(func_t *f);
func_t *func_mul_script(func_t *f);

// mul (private)
void func_mul_args(func_t *f, func_is_t *fis, func_is_t *fin);
void func_mul_args_pow(func_t *f);
func_t *func_mul_split_list(func_t *f, func_is_t *fis);
func_t *func_mul_split_pow_list(func_t *f);

// pow
func_t *func_pow_n(func_t *f, int pow);
func_t *func_pow(func_t *f1, func_t *f2);
func_t *func_pow_eval(func_t *f);
void func_pow_print(func_t *f);
func_t *func_pow_script(func_t *f);

// inv, div
func_t *func_inv(func_t *f1);
func_t *func_div(func_t *f1, func_t *f2);
func_t *func_div_script(func_t *f);

// sub
func_t *func_sub(func_t *f1, func_t *f2);
func_t *func_sub_script(func_t *f);

// math
func_t *func_sqrt(func_t *g);
func_t *func_sqrt_eval(func_t *f);
func_t *func_exp(func_t *g);
func_t *func_exp_eval(func_t *f);
func_t *func_log(func_t *g);
func_t *func_log_eval(func_t *f);
func_t *func_sin(func_t *g);
func_t *func_sin_eval(func_t *f);
func_t *func_cos(func_t *g);
func_t *func_cos_eval(func_t *f);
func_t *func_tan(func_t *g);
func_t *func_tan_eval(func_t *f);
func_t *func_asin(func_t *g);
func_t *func_asin_eval(func_t *f);
func_t *func_acos(func_t *g);
func_t *func_acos_eval(func_t *f);
func_t *func_atan(func_t *g);
func_t *func_atan_eval(func_t *f);
func_t *func_sinh(func_t *g);
func_t *func_sinh_eval(func_t *f);
func_t *func_cosh(func_t *g);
func_t *func_cosh_eval(func_t *f);
func_t *func_tanh(func_t *g);
func_t *func_tanh_eval(func_t *f);
func_t *func_asinh(func_t *g);
func_t *func_asinh_eval(func_t *f);
func_t *func_acosh(func_t *g);
func_t *func_acosh_eval(func_t *f);
func_t *func_atanh(func_t *g);
func_t *func_atanh_eval(func_t *f);
// operator (private)
func_t *func_flatten(func_t *f, const char *op);
int func_is_split_mul(func_t *f);
int func_is_split(func_t *f);
func_t *func_get_number_part(func_t *f);
func_t *func_get_split_part(func_t *f);
int func_is_coeff(func_t *f);
func_t *func_sqrt_pow_n(func_t *f, int pow);

// expand
func_t *func_expand(func_t *f);
func_t *func_expand_eval(func_t *f);

// expand (private)
func_t *func_expand_mul_pow_n(func_t *f);
func_t *func_expand_mul_pow1(func_t *f);
func_t *func_expand_add_pow_n(func_t *f);
func_t *func_expand_add_pow1(func_t *f);
func_t *func_expand_mul1_add1(func_t *f);
func_t *func_expand_list(func_t *f);

// diff
func_t *func_diff(func_t *f, func_t *x);
func_t *func_diff_eval(func_t *f);
func_t *func_grad(func_t *f, func_t *x);
func_t *func_grad_eval(func_t *f);
// diff (private)
func_t *func_diff_var(func_t *f, int var);
func_t *func_diff_pow_n(func_t *f, int var);
func_t *func_diff_pow(func_t *f, int var);
func_t *func_diff_list(func_t *f, int var);
func_t *func_diff_add_pow1(func_t *f, int var);
func_t *func_diff_mul_pow1(func_t *f, int var);
func_t *func_diff_sqrt_pow1(func_t *f, int var);
func_t *func_diff_exp_pow1(func_t *f, int var);
func_t *func_diff_log_pow1(func_t *f, int var);
func_t *func_diff_sin_pow1(func_t *f, int var);
func_t *func_diff_cos_pow1(func_t *f, int var);
func_t *func_diff_tan_pow1(func_t *f, int var);
func_t *func_diff_sinh_pow1(func_t *f, int var);
func_t *func_diff_cosh_pow1(func_t *f, int var);
func_t *func_diff_tanh_pow1(func_t *f, int var);
func_t *func_diff_asin_pow1(func_t *f, int var);
func_t *func_diff_acos_pow1(func_t *f, int var);
func_t *func_diff_atan_pow1(func_t *f, int var);
func_t *func_diff_asinh_pow1(func_t *f, int var);
func_t *func_diff_acosh_pow1(func_t *f, int var);
func_t *func_diff_atanh_pow1(func_t *f, int var);

// diff for vec, mat
func_t *func_grad(func_t *f, func_t *x);

// mapping
func_t *func_maps(func_t *f, int n0, func_t *x);

// mapping (private)
func_t *func_maps_var(func_t *f, int n0, func_t *x);

// cmp
int func_cmp(func_t *f, func_t *g);
int func_lt(func_t *f, func_t *g);
int func_le(func_t *f, func_t *g);
int func_gt(func_t *f, func_t *g);
int func_ge(func_t *f, func_t *g);

// args
int func_asize(func_t *f);
func_t *func_aget(func_t *f, int i);
void func_aset(func_t *f, int i, func_t *g);
void func_adel(func_t *f, int i);

// args (private)
void func_a_new(func_t *f, int n);
void func_a_del(func_t *f);
void func_a_clone(func_t *f, func_t *g);
void func_a_replace(func_t *f, int n, func_t **arg);
void func_a_swap(func_t *f, func_t *g);
void func_a_resize(func_t *f, int n);
void func_a_insert(func_t *f, int n, func_t *g);
void func_a_append(func_t *f, func_t *g);
void func_a_rm(func_t *f, int i);
void func_a_rm_end(func_t *f);
void func_a_rm_null(func_t *f);
void func_a_rm_op(func_t *f, func_is_t *fis);
void func_a_del_op(func_t *f, func_is_t *fis);
void func_a_rm_not_op(func_t *f, func_is_t *fis);
void func_a_del_not_op(func_t *f, func_is_t *fis);
void func_a_del_not_op_op(func_t *f, func_is_t *fis1, func_is_t *fis2);
int func_a_has_op(func_t *f, func_is_t *fis);
int func_a_has_op_pow1(func_t *f, func_is_t *fis);

void func_args_arrange(func_t *f, int *I);
void func_args_rm_op(func_t *f, const char *op);
void func_args_rm_add_pow1(func_t *f);
void func_args_swap(func_t *f, int i, int j);
void func_args_quick_sort(func_t *f, int *I, int left, int right);
void func_args_sort(func_t *f);
void func_args_sort_index(func_t *f, int *I);
void func_args_reverse(func_t *f);
int func_args_count_non_null(func_t *f);
int func_args_have(func_t *f, func_t *g);
int func_args_count_op(func_t *f, const char *op);

// solve
func_t *func_mat_solve(func_t *A, func_t *b); // x=A\b
func_t *func_poly_solve(func_t *f);
func_t *func_poly_solve_var1(func_t *f);
func_t *func_poly_solve_var2(func_t *f);
func_t *func_poly_solve_var2_have_ac(func_t *f);
func_t *func_poly_solve_var2_have_ab(func_t *f);
func_t *func_poly_solve_var2_have_abc(func_t *f);
func_t *func_poly_calculate_soltion(func_t *a,func_t *b);
func_t *func_poly_quadratic_formula(func_t *a,func_t *b,func_t *c);
func_t *func_poly_solve_varn(func_t *f);
func_t *func_poly_list_solve(func_t *f);
func_t *func_poly_list_solve_and_arrange(func_t *f);
void func_change_cmulti(cmulti *x, func_t *f);
func_t *func_change_cmulti_list(int n, cmulti **x);
void func_ccopy_coeff(int n, cmulti **z, func_t *f);
void func_init_val_r_and_balance(rmulti *r, cmulti *g, int n, cmulti **a);
void func_init_val(int n, cmulti **z, cmulti **a);
int func_weierstrass(int n, cmulti **z1, cmulti **z0, cmulti **a);

// def
int func_is_def(func_t *f);
func_t *func_def(const char *name, func_t *g, int amin, int amax);
void func_def_print(func_t *f);
func_t *func_def_script(func_t *f);
func_t *func_def_any_eval(func_t *f);
void func_def_p_new(func_t *f);
void func_def_p_del(func_t *f);
void func_def_p_clone(func_t *f, func_t *g);
int func_def_p_cmp(func_t *f, func_t *g);
func_def_t *func_def_p(func_t *f);
char *func_def_name(func_t *f);
int func_def_amin(func_t *f);
int func_def_amax(func_t *f);
func_def_t *func_def_struct_new(void);
func_def_t *func_def_struct_del(func_def_t *def);

// strings
int func_is_strings(func_t *f);
func_t *func_strings(int n);
func_t *func_strings_str(const char *str[]);
func_t *func_strings_strings(strings *str);
func_t *func_strings_char(const char *str);
func_t *func_strings_split(const char *str, const char *sep, const char *mask_begin, const char *mask_end, const char *skip);
int func_strings_p_cmp(func_t *f, func_t *g);
int func_strings_size(func_t *f);
char *func_strings_at(func_t *f, int i);
void func_strings_set(func_t *f, int i, const char *str);
strings *func_strings_p(func_t *f);

// strings (private)
void func_strings_p_new(func_t *f);
void func_strings_p_del(func_t *f);
void func_strings_p_clone(func_t *f, func_t *g);
void func_strings_p_replace(func_t *f, strings *s);
void func_strings_print(func_t *f);

// output
void func_print(func_t *f);
void func_put(func_t *f);

// output (private)
func_t *func_print_eval(func_t *f);
void func_print_hat(int pow);
void func_print_complex(func_t *f);
void func_print_args(func_t *f);
void func_print_ivec(func_t *f);
void func_print_rvec(func_t *f);
void func_print_cvec(func_t *f);
void func_print_rmat(func_t *f);
void func_print_cmat(func_t *f);
void func_print_table(func_t *f);
void func_print_scope(func_t *f);
void func_print_scope(func_t *f);

// macros
#define FUNC_ERROR_ARG1(F,F1)      { printf("Error in %s(%s)\n",F,func_op(F1)); printf("#0="); func_print(F1); printf("\n");  exit(0); }
#define FUNC_ERROR_ARG2(F,F1,F2)   { printf("Error in %s(%s, %s)\n",F,func_op(F1),func_op(F2)); printf("#0="); func_print(F1); printf("\n#1=");; func_print(F2); printf("\n"); exit(0); }
#define FUNC_ERROR_ARG3(F,F1,F2,F3){ printf("Error in %s(%s, %s, %s)\n",F,func_op(F1),func_op(F2),func_op(F3)); exit(0); }


#endif
