clear all
close all
set_default_prec(128);

M=3;N=4;L=1;
B=zeros(M,N,L);
A=multi(B)

B(1:2:12)=ones(2,1,3)
A(1:2:12)=ones(2,1,3)

