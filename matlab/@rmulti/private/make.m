clear all
close all

mex -I/usr/local/include -I../../../include -L/usr/local/lib -L../../../MACOSX -lis -lmpfr -lgmp rmulti_data.c
