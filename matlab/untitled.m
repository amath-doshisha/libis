clear all
close all
set_default_prec(256);

A=[123.4+i*0.5678 1-i; 0.12345e-8 -i*0.12345e-8]
A=multi(A)
%A=multi('0.1234 1; 0.12345 0.12345');
S=get_s(A,'%.7e')
B=multi(S)
