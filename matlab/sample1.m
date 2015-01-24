clear all
close all

set_prec(512);

A=cmulti(rand(5,5)*2-1+i*(rand(5,5)*2-1))
b=cmulti([1 2; 3 -1; 4 5; -1 2; 2 2])
x=A\b
r=A*x-b
r.abs.max.max


%A=cmulti([1-i 2-2*i 3-3*i; 4-4*i 5-5*i 6-6*i; 7-7*i 8-8*i 9-9*i])
%A=cmulti(rmulti([1 2; 3 -1; 4 5; -1 2; 2 2]))
% B=A.real
% B=A.imag
% B=A.abs

%A=cmulti.zeros(5,5)
%A=cmulti.ones(5,5)
%A=cmulti.rand(5,5)
%B=A.round(53)