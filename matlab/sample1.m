clear all;

n=10;
b=ones(n,1);
A=rand(n,n);

set_default_prec(1024)
b=multi(b);
A=multi(A);

x=A\b;
r=A*x-b;
disp(sprintf('||A*x-b||_inf=%.1e',double(max(abs(r)))));