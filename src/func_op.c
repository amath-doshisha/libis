#include<stdio.h>
#include<stdlib.h>
#include"is_macros.h"
#include"is_func.h"

////////////////////////////////////////////////////////////////////

func_t *func_op_print_new(void)
{
  func_t *f=func_builtin_new("print");
  func_builtin_p(f)->amin=1;
  func_builtin_p(f)->amax=1;
  func_builtin_p(f)->type=FUNC_COMMAND;
  func_builtin_p(f)->eval=func_print_eval;
  return f;
}

///////////////////////////////////////////////////////////////////////

func_t *func_op_evalf_new(void)
{
  func_t *f=func_builtin_new("evalf");
  func_builtin_p(f)->amin=1;
  func_builtin_p(f)->amax=1;
  func_builtin_p(f)->type=FUNC_SCALAR_DEPEND;
  func_builtin_p(f)->eval=func_evalf_eval;
  return f;
}

///////////////////////////////////////////////////////////////////////

func_t *func_op_expand_new(void)
{
  func_t *f=func_builtin_new("expand");
  func_builtin_p(f)->amin=1;
  func_builtin_p(f)->amax=1;
  func_builtin_p(f)->type=FUNC_SCALAR_DEPEND;
  func_builtin_p(f)->eval=func_expand_eval;
  return f;
}

////////////////////////////////////////////////////////////////////

func_t *func_op_gbasis_new(void)
{
  func_t *f=func_builtin_new("gbasis");
  func_builtin_p(f)->amin=1;
  func_builtin_p(f)->amax=1;
  func_builtin_p(f)->eval=func_gbasis_eval;
  return f;
}

////////////////////////////////////////////////////////////////////

func_t *func_op_sqrt_new(void)
{
  func_t *f=func_builtin_new("sqrt");
  func_builtin_p(f)->order=FUNC_ORDER_SQRT;
  func_builtin_p(f)->amin=1;
  func_builtin_p(f)->amax=1;
  func_builtin_p(f)->type=FUNC_SCALAR_DEPEND;
  func_builtin_p(f)->ptype=FUNC_P_POWER;
  func_builtin_p(f)->p_new=func_pow_p_new;
  func_builtin_p(f)->p_del=func_pow_p_del;
  func_builtin_p(f)->p_clone=func_pow_p_clone;
  func_builtin_p(f)->p_cmp=func_power_cmp;
  func_builtin_p(f)->eval=func_sqrt_eval;
  func_builtin_p(f)->print=func_print_args;
  return f;
}

func_t *func_op_exp_new(void)
{
  func_t *f=func_builtin_new("exp");
  func_builtin_p(f)->order=FUNC_ORDER_EXP;
  func_builtin_p(f)->amin=1;
  func_builtin_p(f)->amax=1;
  func_builtin_p(f)->type=FUNC_SCALAR_DEPEND;
  func_builtin_p(f)->ptype=FUNC_P_POWER;
  func_builtin_p(f)->p_new=func_pow_p_new;
  func_builtin_p(f)->p_del=func_pow_p_del;
  func_builtin_p(f)->p_clone=func_pow_p_clone;
  func_builtin_p(f)->p_cmp=func_power_cmp;
  func_builtin_p(f)->eval=func_exp_eval;
  func_builtin_p(f)->print=func_print_args;
  return f;
}

func_t *func_op_log_new(void)
{
  func_t *f=func_builtin_new("log");
  func_builtin_p(f)->order=FUNC_ORDER_LOG;
  func_builtin_p(f)->amin=1;
  func_builtin_p(f)->amax=1;
  func_builtin_p(f)->type=FUNC_SCALAR_DEPEND;
  func_builtin_p(f)->ptype=FUNC_P_POWER;
  func_builtin_p(f)->p_new=func_pow_p_new;
  func_builtin_p(f)->p_del=func_pow_p_del;
  func_builtin_p(f)->p_clone=func_pow_p_clone;
  func_builtin_p(f)->p_cmp=func_power_cmp;
  func_builtin_p(f)->eval=func_log_eval;
  func_builtin_p(f)->print=func_print_args;
  return f;
}

func_t *func_op_sin_new(void)
{
  func_t *f=func_builtin_new("sin");
  func_builtin_p(f)->order=FUNC_ORDER_SIN;
  func_builtin_p(f)->amin=1;
  func_builtin_p(f)->amax=1;
  func_builtin_p(f)->type=FUNC_SCALAR_DEPEND;
  func_builtin_p(f)->ptype=FUNC_P_POWER;
  func_builtin_p(f)->p_new=func_pow_p_new;
  func_builtin_p(f)->p_del=func_pow_p_del;
  func_builtin_p(f)->p_clone=func_pow_p_clone;
  func_builtin_p(f)->p_cmp=func_power_cmp;
  func_builtin_p(f)->eval=func_sin_eval;
  func_builtin_p(f)->print=func_print_args;
  return f;
}

func_t *func_op_cos_new(void)
{
  func_t *f=func_builtin_new("cos");
  func_builtin_p(f)->order=FUNC_ORDER_COS;
  func_builtin_p(f)->amin=1;
  func_builtin_p(f)->amax=1;
  func_builtin_p(f)->type=FUNC_SCALAR_DEPEND;
  func_builtin_p(f)->ptype=FUNC_P_POWER;
  func_builtin_p(f)->p_new=func_pow_p_new;
  func_builtin_p(f)->p_del=func_pow_p_del;
  func_builtin_p(f)->p_clone=func_pow_p_clone;
  func_builtin_p(f)->p_cmp=func_power_cmp;
  func_builtin_p(f)->eval=func_cos_eval;
  func_builtin_p(f)->print=func_print_args;
  return f;
}

func_t *func_op_tan_new(void)
{
  func_t *f=func_builtin_new("tan");
  func_builtin_p(f)->order=FUNC_ORDER_TAN;
  func_builtin_p(f)->amin=1;
  func_builtin_p(f)->amax=1;
  func_builtin_p(f)->type=FUNC_SCALAR_DEPEND;
  func_builtin_p(f)->ptype=FUNC_P_POWER;
  func_builtin_p(f)->p_new=func_pow_p_new;
  func_builtin_p(f)->p_del=func_pow_p_del;
  func_builtin_p(f)->p_clone=func_pow_p_clone;
  func_builtin_p(f)->p_cmp=func_power_cmp;
  func_builtin_p(f)->eval=func_tan_eval;
  func_builtin_p(f)->print=func_print_args;
  return f;
}

func_t *func_op_asin_new(void)
{
  func_t *f=func_builtin_new("asin");
  func_builtin_p(f)->order=FUNC_ORDER_ASIN;
  func_builtin_p(f)->amin=1;
  func_builtin_p(f)->amax=1;
  func_builtin_p(f)->type=FUNC_SCALAR_DEPEND;
  func_builtin_p(f)->ptype=FUNC_P_POWER;
  func_builtin_p(f)->p_new=func_pow_p_new;
  func_builtin_p(f)->p_del=func_pow_p_del;
  func_builtin_p(f)->p_clone=func_pow_p_clone;
  func_builtin_p(f)->p_cmp=func_power_cmp;
  func_builtin_p(f)->eval=func_asin_eval;
  func_builtin_p(f)->print=func_print_args;
  return f;
}

func_t *func_op_acos_new(void)
{
  func_t *f=func_builtin_new("acos");
  func_builtin_p(f)->order=FUNC_ORDER_ACOS;
  func_builtin_p(f)->amin=1;
  func_builtin_p(f)->amax=1;
  func_builtin_p(f)->type=FUNC_SCALAR_DEPEND;
  func_builtin_p(f)->ptype=FUNC_P_POWER;
  func_builtin_p(f)->p_new=func_pow_p_new;
  func_builtin_p(f)->p_del=func_pow_p_del;
  func_builtin_p(f)->p_clone=func_pow_p_clone;
  func_builtin_p(f)->p_cmp=func_power_cmp;
  func_builtin_p(f)->eval=func_acos_eval;
  func_builtin_p(f)->print=func_print_args;
  return f;
}

func_t *func_op_atan_new(void)
{
  func_t *f=func_builtin_new("atan");
  func_builtin_p(f)->order=FUNC_ORDER_ATAN;
  func_builtin_p(f)->amin=1;
  func_builtin_p(f)->amax=1;
  func_builtin_p(f)->type=FUNC_SCALAR_DEPEND;
  func_builtin_p(f)->ptype=FUNC_P_POWER;
  func_builtin_p(f)->p_new=func_pow_p_new;
  func_builtin_p(f)->p_del=func_pow_p_del;
  func_builtin_p(f)->p_clone=func_pow_p_clone;
  func_builtin_p(f)->p_cmp=func_power_cmp;
  func_builtin_p(f)->eval=func_atan_eval;
  func_builtin_p(f)->print=func_print_args;
  return f;
}

func_t *func_op_sinh_new(void)
{
  func_t *f=func_builtin_new("sinh");
  func_builtin_p(f)->order=FUNC_ORDER_SINH;
  func_builtin_p(f)->amin=1;
  func_builtin_p(f)->amax=1;
  func_builtin_p(f)->type=FUNC_SCALAR_DEPEND;
  func_builtin_p(f)->ptype=FUNC_P_POWER;
  func_builtin_p(f)->p_new=func_pow_p_new;
  func_builtin_p(f)->p_del=func_pow_p_del;
  func_builtin_p(f)->p_clone=func_pow_p_clone;
  func_builtin_p(f)->p_cmp=func_power_cmp;
  func_builtin_p(f)->eval=func_sinh_eval;
  func_builtin_p(f)->print=func_print_args;
  return f;
}

func_t *func_op_cosh_new(void)
{
  func_t *f=func_builtin_new("cosh");
  func_builtin_p(f)->order=FUNC_ORDER_COSH;
  func_builtin_p(f)->amin=1;
  func_builtin_p(f)->amax=1;
  func_builtin_p(f)->type=FUNC_SCALAR_DEPEND;
  func_builtin_p(f)->ptype=FUNC_P_POWER;
  func_builtin_p(f)->p_new=func_pow_p_new;
  func_builtin_p(f)->p_del=func_pow_p_del;
  func_builtin_p(f)->p_clone=func_pow_p_clone;
  func_builtin_p(f)->p_cmp=func_power_cmp;
  func_builtin_p(f)->eval=func_cosh_eval;
  func_builtin_p(f)->print=func_print_args;
  return f;
}

func_t *func_op_tanh_new(void)
{
  func_t *f=func_builtin_new("tanh");
  func_builtin_p(f)->order=FUNC_ORDER_TANH;
  func_builtin_p(f)->amin=1;
  func_builtin_p(f)->amax=1;
  func_builtin_p(f)->type=FUNC_SCALAR_DEPEND;
  func_builtin_p(f)->ptype=FUNC_P_POWER;
  func_builtin_p(f)->p_new=func_pow_p_new;
  func_builtin_p(f)->p_del=func_pow_p_del;
  func_builtin_p(f)->p_clone=func_pow_p_clone;
  func_builtin_p(f)->p_cmp=func_power_cmp;
  func_builtin_p(f)->eval=func_tanh_eval;
  func_builtin_p(f)->print=func_print_args;
  return f;
}

func_t *func_op_asinh_new(void)
{
  func_t *f=func_builtin_new("asinh");
  func_builtin_p(f)->order=FUNC_ORDER_ASINH;
  func_builtin_p(f)->amin=1;
  func_builtin_p(f)->amax=1;
  func_builtin_p(f)->type=FUNC_SCALAR_DEPEND;
  func_builtin_p(f)->ptype=FUNC_P_POWER;
  func_builtin_p(f)->p_new=func_pow_p_new;
  func_builtin_p(f)->p_del=func_pow_p_del;
  func_builtin_p(f)->p_clone=func_pow_p_clone;
  func_builtin_p(f)->p_cmp=func_power_cmp;
  func_builtin_p(f)->eval=func_asinh_eval;
  func_builtin_p(f)->print=func_print_args;
  return f;
}

func_t *func_op_acosh_new(void)
{
  func_t *f=func_builtin_new("acosh");
  func_builtin_p(f)->order=FUNC_ORDER_ACOSH;
  func_builtin_p(f)->amin=1;
  func_builtin_p(f)->amax=1;
  func_builtin_p(f)->type=FUNC_SCALAR_DEPEND;
  func_builtin_p(f)->ptype=FUNC_P_POWER;
  func_builtin_p(f)->p_new=func_pow_p_new;
  func_builtin_p(f)->p_del=func_pow_p_del;
  func_builtin_p(f)->p_clone=func_pow_p_clone;
  func_builtin_p(f)->p_cmp=func_power_cmp;
  func_builtin_p(f)->eval=func_acosh_eval;
  func_builtin_p(f)->print=func_print_args;
  return f;
}

func_t *func_op_atanh_new(void)
{
  func_t *f=func_builtin_new("atanh");
  func_builtin_p(f)->order=FUNC_ORDER_ATANH;
  func_builtin_p(f)->amin=1;
  func_builtin_p(f)->amax=1;
  func_builtin_p(f)->type=FUNC_SCALAR_DEPEND;
  func_builtin_p(f)->ptype=FUNC_P_POWER;
  func_builtin_p(f)->p_new=func_pow_p_new;
  func_builtin_p(f)->p_del=func_pow_p_del;
  func_builtin_p(f)->p_clone=func_pow_p_clone;
  func_builtin_p(f)->p_cmp=func_power_cmp;
  func_builtin_p(f)->eval=func_atanh_eval;
  func_builtin_p(f)->print=func_print_args;
  return f;
}

func_t *func_op_ivec_new(void)
{
  func_t *f=func_builtin_new("ivec");
  func_builtin_p(f)->order=FUNC_ORDER_IVEC;
  func_builtin_p(f)->ptype=FUNC_P_IVEC;
  func_builtin_p(f)->p_new=func_ivec_p_new;
  func_builtin_p(f)->p_del=func_ivec_p_del;
  func_builtin_p(f)->p_clone=func_ivec_p_clone;
  func_builtin_p(f)->p_cmp=func_ivec_cmp;
  func_builtin_p(f)->print=func_print_ivec;
  return f;
}

func_t *func_op_rvec_new(void)
{
  func_t *f=func_builtin_new("rvec");
  func_builtin_p(f)->order=FUNC_ORDER_RVEC;
  func_builtin_p(f)->ptype=FUNC_P_RVEC;
  func_builtin_p(f)->p_new=func_rvec_p_new;
  func_builtin_p(f)->p_del=func_rvec_p_del;
  func_builtin_p(f)->p_clone=func_rvec_p_clone;
  func_builtin_p(f)->p_cmp=func_rvec_cmp;
  func_builtin_p(f)->print=func_print_rvec;
  return f;
}

func_t *func_op_cvec_new(void)
{
  func_t *f=func_builtin_new("cvec");
  func_builtin_p(f)->order=FUNC_ORDER_CVEC;
  func_builtin_p(f)->ptype=FUNC_P_CVEC;
  func_builtin_p(f)->p_new=func_cvec_p_new;
  func_builtin_p(f)->p_del=func_cvec_p_del;
  func_builtin_p(f)->p_clone=func_cvec_p_clone;
  func_builtin_p(f)->p_cmp=func_cvec_cmp;
  func_builtin_p(f)->print=func_print_cvec;
  return f;
}

func_t *func_op_rmat_new(void)
{
  func_t *f=func_builtin_new("rmat");
  func_builtin_p(f)->order=FUNC_ORDER_RMAT;
  func_builtin_p(f)->ptype=FUNC_P_RMAT;
  func_builtin_p(f)->p_new=func_rmat_p_new;
  func_builtin_p(f)->p_del=func_rmat_p_del;
  func_builtin_p(f)->p_clone=func_rmat_p_clone;
  func_builtin_p(f)->p_cmp=func_rmat_cmp;
  func_builtin_p(f)->print=func_print_rmat;
  return f;
}

func_t *func_op_cmat_new(void)
{
  func_t *f=func_builtin_new("cmat");
  func_builtin_p(f)->order=FUNC_ORDER_CMAT;
  func_builtin_p(f)->ptype=FUNC_P_CMAT;
  func_builtin_p(f)->p_new=func_cmat_p_new;
  func_builtin_p(f)->p_del=func_cmat_p_del;
  func_builtin_p(f)->p_clone=func_cmat_p_clone;
  func_builtin_p(f)->p_cmp=func_cmat_cmp;
  func_builtin_p(f)->print=func_print_cmat;
  return f;
}

//private//////////////////////////////////////////////



//EOF
