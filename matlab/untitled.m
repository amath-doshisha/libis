clear all
close all
set_default_prec(128);

A=[1 2; 3 4]
A=multi(A)
B=A^-4
A^4*B-eye(2)
inv(A)
inv(A)*A-eye(2)