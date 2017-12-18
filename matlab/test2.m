clear all
close all

set_default_prec(4000);
N=100;
%A=2*rand(3,3)-1;  
%A=2*rand(4,4)-1+i*(2*rand(4,4)-1); 
%a=2*rand(1,4)-1%+i*(2*rand(1,4)-1); 
%A=15+1i*15;
%a=15+1i*15; 
%A=floor(A*10);

%B=2*rand(3,3)-1;
%B=2*rand(4,4)-1%+i*(2*rand(4,4)-1);
%b=2*rand(4,1)-1%+i*(2*rand(4,1)-1);
% B=floor(B*10);
A=rand(N);
B=rand(N,1);
%B=rand()+i*rand();

x=A\B;
A=imulti(A);
B=imulti(B);
%X=A\B;
%r=A*x-B;
%R=A*X-B;

%disp('A=');disp(num2str(A,'%+.7e'))
%disp('B=');disp(num2str(B,'%+.7e'))
%disp('r=');disp(num2str(r,'%+.7e'))
%disp('R=');disp(num2str(R,'%+.7e'))

%get_prec(A)  %A‚Ìprec