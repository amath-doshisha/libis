clear all
close all

set_default_prec(128);
%auto_prec_enabled;
%auto_prec_disabled;


% %% eq
A=multi([2 3; 4 5; -2 -3])
B=multi([2 -3; 4 -5; -2 3])
C=multi([2+i 3; 4 5-i; -2*i -3])
A==B
A==4
5==A
C==A
C==2+i
C==4
C==-2*i
3==C
5-i==C
-3==C


% %% plus
A=multi([2 3; 4 5; -2 -3])
B=multi([2+i 3-i; 4*i 5-i; -2*i -3+i])
A+B
B+A
A+A
B+B
1+A
i+A
1+B
i+B
A+1
A+i
B+1
B+i

%% neg
A=multi([2 3; 4 5; -2 -3])
B=multi([2+i 3*i; 4 5; -2 -3])
-A
-B

%% copy
A=multi([2 3; 4 5; -2 -3])
B=multi([2+i 3*i; 4 5; -2 -3])
multi(A)
multi(B)

%% zeros
A=zeros(3,2,'multi')
size(A)
B=zeros([3,2,2],'multi')
size(B)

%% ones
A=ones(3,2,'multi')
B=ones([3,2,2],'multi')

%% rand
% A=rand(3,2,'multi')
% B=rand([3,2,2],'multi')

%% disp
disp(A);

%% double
double(A)
