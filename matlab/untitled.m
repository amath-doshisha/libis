clear all
close all
set_default_prec(128);

M=5;N=3;L=2;
B=zeros(M,N,L);
B(:)=1:(M*N*L);
A=multi(B)

B(:,2:3,:)
A(:,2:3,:)