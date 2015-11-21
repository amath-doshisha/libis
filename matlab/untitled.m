clear all
close all

% %% eq
% A=multi([2 3; 4 5; -2 -3])
% B=multi([2 -3; 4 -5; -2 3])
% C=multi([2+i 3; 4 5-i; -2*i -3])
% A==B
% A==4
% 5==A
% C==A
% C==2+i
% C==4
% C==-2*i
% 3==C
% 5-i==C
% -3==C


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
A=zeros(3,2,'multi')
B=zeros([3,2,2],'multi')

%% ones
% A=ones(3,2,'multi')
% B=ones([3,2,2],'multi')

%% rand
% A=rand(3,2,'multi')
% B=rand([3,2,2],'multi')

%% disp
%disp(A);

%% double
%double(A)
