
OS = MACOSX

CC = gcc

CFLAGS = -I/usr/local/include -I../include -Wall -Wextra -Werror

LIBS = -L/usr/local/lib -L../$(OS) -lis -lmpfr -lgmp -lm

SRCS = HashTest.c \
	test_dsolve.c test_zsolve.c test_rsolve.c test_csolve.c \
	test_deig.c test_zeig.c test_reig.c \
	test_rmulti.c test_bigint.c \
	test_strings.c \
	test_func.c test_fvec.c test_gb.c \
	test_cnewton.c test_rnewton.c

PROGS = HashTest \
	test_dsolve test_zsolve test_rsolve test_csolve \
	test_deig test_zeig test_reig \
	test_rmulti test_bigint \
	test_strings \
	test_func test_fvec test_gb \
	test_cnewton test_rnewton

all:	$(PROGS)

HashTest:	HashTest.o
	$(CC) -o $@ $< $(CFLAGS) $(LIBS)

test_dsolve:	test_dsolve.o
	$(CC) -o $@ $< $(CFLAGS) $(LIBS)

test_zsolve:	test_zsolve.o
	$(CC) -o $@ $< $(CFLAGS) $(LIBS)

test_rsolve:	test_rsolve.o
	$(CC) -o $@ $< $(CFLAGS) $(LIBS)

test_csolve:	test_csolve.o
	$(CC) -o $@ $< $(CFLAGS) $(LIBS)

test_deig:	test_deig.o
	$(CC) -o $@ $< $(CFLAGS) $(LIBS)

test_zeig:	test_zeig.o
	$(CC) -o $@ $< $(CFLAGS) $(LIBS)

test_reig:	test_reig.o
	$(CC) -o $@ $< $(CFLAGS) $(LIBS)

test_rmulti:	test_rmulti.o
	$(CC) -o $@ $< $(CFLAGS) $(LIBS)

test_bigint:	test_bigint.o
	$(CC) -o $@ $< $(CFLAGS) $(LIBS)

test_strings:	test_strings.o
	$(CC) -o $@ $< $(CFLAGS) $(LIBS)

test_func:	test_func.o
	$(CC) -o $@ $< $(CFLAGS) $(LIBS)

test_fvec:	test_fvec.o
	$(CC) -o $@ $< $(CFLAGS) $(LIBS)

test_gb:	test_gb.o
	$(CC) -o $@ $< $(CFLAGS) $(LIBS)

test_cnewton:	test_cnewton.o
	$(CC) -o $@ $< $(CFLAGS) $(LIBS)

test_rnewton:	test_rnewton.o
	$(CC) -o $@ $< $(CFLAGS) $(LIBS)

.c.o:	$(SRCS) $(HDRS)
	$(CC) -c $< $(CFLAGS)

clean:
	rm -f $(PROGS) *.o *~ foo foo.* bar bar.* *.dat *.txt
