clear all
close all

mex -I/usr/local/include -I../../../include -L/usr/local/lib -L../../../MACOSX -lis -lmpfr -lgmp rmulti_allocate.c
mex -I/usr/local/include -I../../../include -L/usr/local/lib -L../../../MACOSX -lis -lmpfr -lgmp rmulti_double.c
mex -I/usr/local/include -I../../../include -L/usr/local/lib -L../../../MACOSX -lis -lmpfr -lgmp rmulti_two.c

