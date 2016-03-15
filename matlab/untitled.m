clear all
close all
set_default_prec(256);

%A=[1+1i -1-1i; -1+1i 1-1i]
A=[1 -1; -1 1]
A=multi(A);
% abs(A)
% angle(A)*(180/pi)
% real(A)
% imag(A)
conj(A)
