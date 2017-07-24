%function sample5

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
%n=10;A=rand(n)
A=toeplitz([1 2 0 1 0],[1 2 0 0 0])
lambda=sort(eig(A))

set_default_prec(512)
A=multi(A);
lambda2=sort(double(eig(A)))

plot(real(lambda),imag(lambda),'o',real(lambda2),imag(lambda2),'x')
axis square
grid on




%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clear all
% close all
% 
% set_default_prec(512)
% A=zeros(10,10,'multi');
% A(:)=1:100;
% A
% A(44)
% A(7,:)
% A(:,5)
% A(1:2:10,2:2:10)
% i=1;
% while i<=100
%     A(i)=-1*A(i)
%     i=i+1;
% end