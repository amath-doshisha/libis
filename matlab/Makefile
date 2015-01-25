
SUFFIX=mexmaci64

MATLAB_BIN=/Applications/MATLAB_R2014b.app/bin

MEX=$(MATLAB_BIN)/mex

CFLAGS=-I../include -I/usr/local/include

LIBS=-L../MACOSX -I/usr/local/lib -lis -lmpfr -lgmp

SRCS=rmulti_data.c cmulti_data.c

all:
	for x in $(SRCS) ; do \
		echo $(MEX) $(CFLAGS) $(LIBS) $$x ; \
		$(MEX) $(CFLAGS) $(LIBS) $$x ; \
	done

clean:
	( cd @rmulti ; make clean )
	( cd @cmulti ; make clean )
	rm -fv \#* *~ *.mex*
