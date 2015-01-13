
OS=MACOSX
GCC=gcc
INST_DIR_INC=../include
INST_DIR_LIB=../$(OS)
CFLAGS=-I/usr/local/include -I. -Wall -Wextra -Wcast-align -O3 -funroll-loops -D__$(OS)__

##################################################################

SHARE_HDRS = mt19937ar.h GeneralHashFunctions.h \
		isys.h \
		is_macros.h is_rounding.h is_strings.h \
		is_bigint.h \
		is_print.h is_dcomplex.h \
		is_rmulti.h is_cmulti.h \
		is_irmulti.h is_irvec.h is_irmat.h is_ireig.h \
		is_icmulti.h is_icvec.h is_icmat.h is_iceig.h \
		is_ivec.h is_dvec.h is_zvec.h is_rvec.h is_cvec.h \
		is_imat.h is_dmat.h is_zmat.h is_rmat.h is_cmat.h \
		is_dsolve.h is_zsolve.h is_rsolve.h is_csolve.h \
		is_rsolve_newton.h is_csolve_newton.h \
		is_deig.h is_zeig.h is_reig.h is_ceig.h \
		is_reig_hqr.h is_ceig_hqr.h is_reig_ii.h is_ceig_ii.h \
		is_dsvd.h \
		is_dtestmat.h is_ztestmat.h is_rtestmat.h is_ctestmat.h \
		is_dhshldr.h is_zhshldr.h is_rhshldr.h is_chshldr.h \
		is_rhsnbrg.h is_chsnbrg.h \
		is_dhpeig.h is_zhpeig.h is_rhpeig.h is_chpeig.h \
		is_dhpsvd.h \
		is_func.h

SHARE_SRCS = mt19937ar.c GeneralHashFunctions.c \
		strings.c \
		bigint.c \
		print.c dcomplex.c \
		rmulti.c cmulti.c \
		irmulti.c irvec.c irmat.c ireig.c \
		icmulti.c icvec.c icvec.c iceig.c \
		ivec.c dvec.c zvec.c rvec.c cvec.c \
		imat.c dmat.c zmat.c rmat.c cmat.c \
		dsolve.c zsolve.c rsolve.c csolve.c \
		rsolve_newton.c csolve_newton.c \
		deig.c zeig.c reig.c ceig.c \
		reig_hqr.c ceig_hqr.c reig_ii.c ceig_ii.c \
		dsvd.c \
		dtestmat.c ztestmat.c rtestmat.c ctestmat.c \
		dhshldr.c zhshldr.c rhshldr.c chshldr.c \
		rhsnbrg.c chsnbrg.c \
		dhpeig.c zhpeig.c rhpeig.c chpeig.c \
		dhpsvd.c \
		func_builtin.c func_find.c func_op.c func_table.c func.c func_strings.c func_def.c \
		func_scope.c func_import.c func_power.c func_cmp.c \
		func_print.c func_script.c func_args.c \
		func_add.c func_mul.c func_math.c \
		func_var.c func_mono.c func_poly.c func_gb.c func_list.c \
		func_solve.c func_expand.c \
		func_ivec.c func_rvec.c func_cvec.c func_rmat.c func_cmat.c func_vec.c \
		func_maps.c func_bigint.c func_real.c func_complex.c \
		func_number.c func_diff.c

SHARE_OBJS = mt19937ar.o GeneralHashFunctions.o \
		strings.o \
		bigint.o \
		print.o dcomplex.o \
		rmulti.o cmulti.o \
		irmulti.o irvec.o irmat.o ireig.o \
		icmulti.o icvec.o icmat.o iceig.o \
		ivec.o dvec.o zvec.o rvec.o cvec.o \
		imat.o dmat.o zmat.o rmat.o cmat.o \
		dsolve.o zsolve.o rsolve.o csolve.o \
		rsolve_newton.o csolve_newton.o \
		deig.o zeig.o reig.o ceig.o \
		reig_hqr.o ceig_hqr.o reig_ii.o ceig_ii.o \
		dsvd.o \
		dtestmat.o ztestmat.o rtestmat.o ctestmat.o \
		dhshldr.o zhshldr.o rhshldr.o chshldr.o \
		rhsnbrg.o chsnbrg.o \
		dhpeig.o zhpeig.o rhpeig.o chpeig.o \
		dhpsvd.o \
		func_builtin.o func_find.o func_op.o func_table.o func.o func_strings.o func_def.o \
		func_scope.o func_import.o func_power.o func_cmp.o \
		func_print.o func_script.o func_args.o \
		func_add.o func_mul.o func_math.o \
		func_var.o func_mono.o func_poly.o func_gb.o func_list.o \
		func_solve.o func_expand.o \
		func_ivec.o func_rvec.o func_cvec.o func_rmat.o func_cmat.o func_vec.o \
		func_maps.o func_bigint.o func_real.o func_complex.o \
		func_number.o func_diff.o

##################################################################

HDRS =	$(SHARE_HDRS)

SRCS =  $(SHARE_SRCS)

OBJS =  $(SHARE_OBJS)

PROG = libis.a

##################################################################

all:	$(PROG)

install: $(PROG) $(HDRS)
	mkdir -p $(INST_DIR_INC)
	rm -f $(INST_DIR_INC)/*
	cp -vf $(SHARE_HDRS) $(INST_DIR_INC)
	mkdir -p $(INST_DIR_LIB)
	cp -vf $(PROG) $(INST_DIR_LIB)

##################################################################

$(PROG):	$(SHARE_OBJS) $(SHARE_HDRS)
	rm -f $(PROG)
	ar rc $(PROG) $(SHARE_OBJS)

.c.o:	$(SRCS) $(HDRS)
	$(GCC) -c $< $(CFLAGS)

clean:
	rm -f *.o $(PROG) *.a a.out *.out foo foo.* \#* __*__ *~

#EOF
