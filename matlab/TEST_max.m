clear all
clc

a=[1 2 3;4 5 6;7 8 9];
a(:,:,2)=[10 20 30;40 50 60;70 80 90+i];
A=multi(a)

disp('=================');

max(A)
max(max(A))
max(max(max(A)))

max(max(max(max(A))))


%{
angle([-1-i,1-i,1+i,-1+i])%���f���Ő�Βl���������ꍇ�́A�ʑ��p�̑傫����(����̑傫����)
max(multi([-1-i,1-i,1+i,-1+i]))
max(multi([-1-i,1-i,1+i]))
max(multi([-1-i,1-i]))
%}

