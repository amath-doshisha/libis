
SUFFIX=mexmaci64

#MATLAB_BIN=/Applications/MATLAB_R2015b.app/bin
MATLAB_BIN=/Applications/MATLAB_R2017a.app/bin

MEX=$(MATLAB_BIN)/mex

CFLAGS=-I../include -I/usr/local/include

LIBS=-L../MACOSX -I/usr/local/lib -lis -lmpfr -lgmp

SRCS=multi_mex.c

all:
	for x in $(SRCS) ; do \
		echo $(MEX) $(CFLAGS) $(LIBS) $$x ; \
		$(MEX) $(CFLAGS) $(LIBS) $$x ; \
	done

clean:
	( cd @multi ; make clean )
	rm -fv \#* *~ *.mex*
