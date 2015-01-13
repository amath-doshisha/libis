
SUFFIX=mexmaci64

MATLAB_BIN=/Applications/MATLAB_R2014b.app/bin

MEX=$(MATLAB_BIN)/mex

CFLAGS=-I../../../include -I/usr/local/include

LIBS=-L../../../MACOSX -I/usr/local/lib -lis -lmpfr -lgmp

SRCS=rmulti_allocate.c rmulti_double.c rmulti_two.c rmulti_one.c

all:
	for x in $(SRCS) ; do \
		$(MEX) $(CFLAGS) $(LIBS) $$x ; \
	done

clean:
	rm -fv \#* *~ *.mex*
