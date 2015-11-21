clear all
close all

%% eq
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
% A=multi([2 3; 4 5; -2 -3])
% B=multi([2+i 3-i; 4*i 5-i; -2*i -3+i])
% A+B
% B+A
% A+A
% B+B
% 1+A
% i+A
% 1+B
% i+B
% A+1
% A+i
% B+1
% B+i

%% neg, type='r'
% A=multi([2 3; 4 5; -2 -3])
% B=-A

%% neg, type='c'
% A=multi([2+i 3*i; 4 5; -2 -3])
% B=-A

%% copy, type='r'
%A=multi([2 3; 4 5; -2 -3])
%B=multi(A)

%% copy, type='c'
% A=multi([2+i 3*i; 4 5; -2 -3])
% B=multi(A)

%% zeros
% Ar=multi.zeros('r',3,2)
% Ac=multi.zeros('c',3,2)
% AR=multi.zeros('R',3,2)
% AC=multi.zeros('C',3,2)

%% ones
% Ar=multi.ones('r',3,2)
% Ac=multi.ones('c',3,2)
% AR=multi.ones('R',3,2)
% AC=multi.ones('C',3,2)

%% disp
% disp(Ar);
% disp(Ac);
% disp(AR);
% disp(AC);

%% double
% double(Ar)
% double(Ac)
% double(AR)
% double(AC)
