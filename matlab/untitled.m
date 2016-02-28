clear all
close all
set_default_prec(128);

M=3;N=3;L=2;
B=zeros(M,N,L);
A=multi(B)

A(2:4,2:4,1:3)=ones(3,3,3)