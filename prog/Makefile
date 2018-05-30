
CC = gcc

CFLAGS = -I/usr/local/include -I../include -Wall -Wextra -Werror

LIBS = -L/usr/local/lib -L../lib -lis -lmpfr -lgmp -lm

SRCS = 	matshow.c \
	matshow2.c \
	mattrans.c \
	mateig.c \
	matgen_dhToda.c \
	mateig_krawczyk.c \
	imatshow2.c \
	sovle_newton.c

PROGS =	matshow \
	matshow2 \
	mattrans \
	mateig \
	matgen \
	matgen_dhToda \
	mateig_krawczyk \
	imatshow2 \
	solve_newton

all:	$(PROGS)

matshow:	matshow.o
	$(CC) -o $@ $< $(CFLAGS) $(LIBS)

matshow2:	matshow2.o
	$(CC) -o $@ $< $(CFLAGS) $(LIBS)

mattrans:	mattrans.o
	$(CC) -o $@ $< $(CFLAGS) $(LIBS)

mateig:		mateig.o
	$(CC) -o $@ $< $(CFLAGS) $(LIBS)

matgen:	matgen.o
	$(CC) -o $@ $< $(CFLAGS) $(LIBS)

matgen_dhToda:	matgen_dhToda.o
	$(CC) -o $@ $< $(CFLAGS) $(LIBS)

mateig_krawczyk:	mateig_krawczyk.o
	$(CC) -o $@ $< $(CFLAGS) $(LIBS)

imatshow2:	imatshow2.o
	$(CC) -o $@ $< $(CFLAGS) $(LIBS)

solve_newton:	solve_newton.o
	$(CC) -o $@ $< $(CFLAGS) $(LIBS)

.c.o:	$(SRCS) $(HDRS)
	$(CC) -c $< $(CFLAGS)

clean:
	rm -f $(PROGS) *.o *~ \#* foo foo.* bar bar.* *.dat *.txt
