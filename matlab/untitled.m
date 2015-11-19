clear all
close all

%% plus(r,r)
A=multi([2 3; 4 5; -2 -3])
B=multi([-1 2; 3 1; 1 2])
A+B

%% plus(c,c)
A=multi([2+i 3-i; 4*i 5-i; -2*i -3+i])
B=multi([-1+i 2-i; 3*i 1; 1+i 2-i])
A+B




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
