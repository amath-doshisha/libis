clear all
close all
set_default_prec(128);

L=2;
A=ones(3,3,L)
B=ones(2,3,L)*2
C=ones(1,3,L)*3
A=multi(A);
B=multi(B);
C=multi(C);

D=[A; B; C]