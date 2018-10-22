clear all

% m=5;
% x=(1:m)'
% %x=multi(x)
% A=double(diag(x))
% A=double(diag(x,-1))
% A=double(diag(x,1))
% A=double(diag(x,-2))
% A=double(diag(x,2))
% 

A=rand(3,5)
A=multi(A)
diag(A)
diag(A,1)
diag(A,-1)
diag(A,2)
diag(A,-2)
diag(A,3)
diag(A,-3)
diag(A,4)
diag(A,-4)

