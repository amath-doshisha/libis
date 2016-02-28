clear all
close all
set_default_prec(128);

M=3;N=4;L=2;
B=zeros(M,N,L);
B(:)=1:(M*N*L);
A=multi(B)

B(1:3:20)
A(1:3:20)