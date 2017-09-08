clear all
close all

set_default_prec(26);
A=2*rand(4,4)-1;
B=2*rand(4,4)-1;
%A=2*rand(4,4)-1+i*(2*rand(4,4)-1);
A=imulti(A)
B=imulti(B)
C=A+B;
C=C+B;C=C+B;C=C+B;C=C+B;
C=C+A;C=C+A;C=C+A;C=C+A;
disp('A=');disp(num2str(A,'%+.7e'))
disp('B=');disp(num2str(B,'%+.7e'))
disp('C=');disp(num2str(C,'%+.7e'))

mid(A)
mid(B)
mid(C)

rad(A)
rad(B)
rad(C)
(up(C)-down(C))/2


