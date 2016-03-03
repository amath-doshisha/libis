%function sample1
addpath('../');
clear all;

% Initialization by double
n=10;
b=ones(n,1);
A=rand(n,n);

% Cast double to multi
set_default_prec(1024)
b=multi(b);
A=multi(A);

% Calculations
x=A\b;
r=A*x-b;
rmax=max(abs(r));
disp('x=inv(A)*b=');
disp(num2str(x,'   %+.30e '));
disp(sprintf('||A*x-b||_inf=%s',num2str(rmax,'%.1e')));
