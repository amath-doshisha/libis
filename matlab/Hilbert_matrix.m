clear all
close all

set_default_prec(1000);
N=5;

A=zeros(N);
A=multi(A);
x=0;
x=multi(x);
j=1;k=1;
while(j<=N)
    while(k<=N)
        x=1/(j+k-1);
        A(j,k)=x;
        k=k+1;
    end
    k=1;j=j+1;
end

% H=hilb(N);
% H=imulti(H);

B=inv(A);
E=A*B;

disp('A=');disp(num2str(A,'%+.7e'))
disp('E=');disp(num2str(E,'%+.7e'))