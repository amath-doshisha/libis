clear all
close all
set_default_prec(128);

A=[1+i; 3+4*i; 1-i; -3+4*i; -1+i; -1-i; -3-4*i]
A=multi(A);
%B=sort(A)
max(A)
min(A)
