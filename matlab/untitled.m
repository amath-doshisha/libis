clear all
close all
set_default_prec(128);

M=3;N=3;L=2;
A=zeros(M,N,L);
A=multi(A)
A(:)=1:(M*N*L)
A(1:2:M,1:2:N,2)
