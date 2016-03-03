clear all
close all
set_default_prec(256);

S='1+0i  -1+0i; 0+1i 0-1i'
%S='1  1e-3; -1e+2 +3e2'
A=multi(S)
%s=get_s(A,'%.8e     ')
s=num2str(A,'%.8e     ')
B=multi(s)