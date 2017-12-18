clear all
close all

set_default_prec(26);
A=2*rand(4,4)-1;  
cA=2*rand(4,4)-1+i*(2*rand(4,4)-1); 
x=rand();
cx=rand()+i*rand();

B=2*rand(4,4)-1;
cB=2*rand(4,4)-1+i*(2*rand(4,4)-1);
y=rand();
cy=rand()+i*rand();


A=imulti(A)
cA=imulti(cA)
x=imulti(x)
cx=imulti(cx)

B=imulti(B)
cB=imulti(cB)
y=imulti(y)
cy=imulti(cy)

A+B
A+cB
A+y
A+cy
cA+B
cA+cB
cA+y
cA+cy
x+B
x+cB
x+y
x+cy
cx+B
cx+cB
cx+y
cx+cy


%a=imulti(a)
%b=imulti(b)
%c=a+b
C=A+B
%C=A+a
%c=a+b
%C=A+B;
% C=C+B;C=C+B;C=C+B;C=C+B;
% C=C+A;C=C+A;C=C+A;C=C+A;

disp('A=');disp(num2str(A,'%+.7e'))
disp('B=');disp(num2str(B,'%+.7e'))
disp('C=');disp(num2str(C,'%+.7e'))