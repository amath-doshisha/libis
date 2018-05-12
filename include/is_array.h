#ifndef IS_ARRAY_H
#define IS_ARRAY_H

#include<is_rmulti.h>
#include<is_cmulti.h>

/**
 @file  is_array.h
 @brief 配列型arrayの宣言と関数．
 */

/*
 * struct
 */
typedef struct {
  int type;   /**< 型：type='i','d','z','r','c','R','C' */
  int size;   /**< サイズ */
  int ndim;   /**< 次元数 */
  int *dim;   /**< 次元配列 */
  int *LD;    /**< LD */
  void *p0;   /**< 配列１ */
  void *p1;   /**< 配列２ */
} array_struct;
typedef array_struct array;

/*
 * macros
 */
#define ARRAY_TYPE(A)   (A)->type
#define ARRAY_NDIM(A)   (A)->ndim
#define ARRAY_DIM_P(A)  (A)->dim
#define ARRAY_DIM(A,I)  ((I)<0?0:(((I)>=(A)->ndim)?1:(A)->dim[I]))
#define ARRAY_M(A)      (A)->dim[0]
#define ARRAY_N(A)      (A)->dim[1]
#define ARRAY_L(A)      (A)->dim[2]
#define ARRAY_LD_P(A)   (A)->LD
#define ARRAY_LD(A,I)   ((I)<0?0:(((I)<(A)->ndim)?(A)->LD[I]:(A)->LD[((A)->ndim)-1]))
#define ARRAY_LD0(A)    (A)->LD[0]
#define ARRAY_LD1(A)    (A)->LD[1]
#define ARRAY_LD2(A)    (A)->LD[2]
#define ARRAY_LD3(A)    (A)->LD[3]
#define ARRAY_LD_END(A) (A)->LD[((A)->ndim)-1]
#define ARRAY_SIZE(A)   (A)->size
#define ARRAY_P0(A)     (A)->p0
#define ARRAY_P0_A(A)   ((array**)((A)->p0))
#define ARRAY_P0_S(A)   ((char**)((A)->p0))
#define ARRAY_P0_I(A)   ((int*)((A)->p0))
#define ARRAY_P0_D(A)   ((double*)((A)->p0))
#define ARRAY_P0_Z(A)   ((dcomplex*)((A)->p0))
#define ARRAY_P0_R(A)   ((rmulti**)((A)->p0))
#define ARRAY_P0_C(A)   ((cmulti**)((A)->p0))
#define ARRAY_P1(A)     (A)->p1
#define ARRAY_P1_I(A)   ((int*)((A)->p1))
#define ARRAY_P1_D(A)   ((double*)((A)->p1))
#define ARRAY_P1_Z(A)   ((dcomplex*)((A)->p1))
#define ARRAY_P1_R(A)   ((rmulti**)((A)->p1))
#define ARRAY_P1_C(A)   ((cmulti**)((A)->p1))
#define ARRAY_AVEC(A,I) ((array**)((A)->p0))[I]
#define ARRAY_SVEC(A,I) ((char**)((A)->p0))[I]
#define ARRAY_IVEC(A,I) ((int*)((A)->p0))[I]
#define ARRAY_DVEC(A,I) ((double*)((A)->p0))[I]
#define ARRAY_ZVEC(A,I) ((dcomplex*)((A)->p0))[I]
#define ARRAY_RVEC(A,I) ((rmulti**)((A)->p0))[I]
#define ARRAY_CVEC(A,I) ((cmulti**)((A)->p0))[I]
#define ARRAY_DVEC0(A,I) ((double*)((A)->p0))[I]
#define ARRAY_ZVEC0(A,I) ((dcomplex*)((A)->p0))[I]
#define ARRAY_RVEC0(A,I) ((rmulti**)((A)->p0))[I]
#define ARRAY_CVEC0(A,I) ((cmulti**)((A)->p0))[I]
#define ARRAY_RVEC1(A,I) ((rmulti**)((A)->p1))[I]
#define ARRAY_CVEC1(A,I) ((cmulti**)((A)->p1))[I]
#define ARRAY_SMAT(A,I,J) ((char**)((A)->p0))[I+J*((A)->LD[0])]
#define ARRAY_IMAT(A,I,J) ((int*)((A)->p0))[I+J*((A)->LD[0])]
#define ARRAY_DMAT(A,I,J) ((double*)((A)->p0))[I+J*((A)->LD[0])]
#define ARRAY_ZMAT(A,I,J) ((dcomplex*)((A)->p0))[I+J*((A)->LD[0])]
#define ARRAY_RMAT(A,I,J) ((rmulti**)((A)->p0))[I+J*((A)->LD[0])]
#define ARRAY_CMAT(A,I,J) ((cmulti**)((A)->p0))[I+J*((A)->LD[0])]
#define ARRAY_RMAT0(A,I,J) ((rmulti**)((A)->p0))[I+J*((A)->LD[0])]
#define ARRAY_CMAT0(A,I,J) ((cmulti**)((A)->p0))[I+J*((A)->LD[0])]
#define ARRAY_RMAT1(A,I,J) ((rmulti**)((A)->p1))[I+J*((A)->LD[0])]
#define ARRAY_CMAT1(A,I,J) ((cmulti**)((A)->p1))[I+J*((A)->LD[0])]

/*
 * constructions and destruction
 */
array *array_allocate(int type, int ndim, int *dim);
array *array_free(array *x);
array *array_copy(array *x);  // y=x
array *array_clone(array *x); // y=x
array **avec_allocate(int n);
array **avec_free(int n, array **x);


/*
 * member variables
 */
void array_get_index(int *index, array *x, int t);
int array_get_position(int *index, array *x);

/*
 * query
 */
int array_type_check(int type);
int array_is_empty(array *x);
int array_same_dim_check(array *x, array *y);
int array_get_ndim(array *x);
int array_get_subdim(array *x, array *y);
int array_compatible_dim_check(array * x, array *y);
int *array_get_inclusion_dim(array * x, array *y);
char array_get_inclusion_type(char x, char y);
int array_is_scalar(array *x);

/*
 * I/O
 */
void array_put(array *x);
void array_print(array *x, char *name, char format, int digits);

/*
 * setting
 */
void array_set_zeros(array *x);
void array_set_ones(array *x);
void array_set_nan(array *x);
void array_set_inf(array *x);
void array_set_all_d(array *x, double y);
void array_set_rand(array *x, double a, double b);
void array_set_grid(array *x);

/*
 * casting
 */
array *array_get(int type, array *x); // y=cast(x)
array *array_get_char(array *x, char format, int digits); // y=char(x)
array *array_get_int(array *x);     // y=int(x)
array *array_get_real(array *x);    // y=real(x)
array *array_get_complex(array *x); // y=complex(x)
array *array_get_double(array *x);  // y=double(x)
array *array_get_multi(array *x);   // y=multi(x)
array *array_get_imulti(array *x);  // y=imulti(x)
array *array_get_imulti2(array *x0, array *x1); // y=imulti2(x0,x1)
array *array_get_complex2(array *xr, array *xi); // y=complex2(xr,xi)

/*
 * operatior for one argumentr
 */
array *array_conj(array *x); // y=conj(x)

/*
 * operatior for two arguments
 */
array *array_add(array *x, array *y); // z=x+y
array *array_sub(array *x, array *y); // z=x-y
array *array_mul(array *x, array *y); // z=x.*y
array *array_div(array *x, array *y); // z=x./y


/*
 * operatior for tree arguments
 */

/*
 * comparision
 */

#endif
