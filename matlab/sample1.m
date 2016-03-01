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

% Done
disp(sprintf('||A*x-b||_inf=%.1e',double(max(abs(r)))));