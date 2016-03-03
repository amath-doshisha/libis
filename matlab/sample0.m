%function sample0
clear all
set_default_prec(128);
disp('Initialize double matrix D, then');
D=zeros(3,3,'double')
disp('Initialize multi matrix M, then');
M=zeros(3,3,'multi')
disp('Copy the multi matrix M to');
disp('a multi matrix A, then');
A=M
disp('The entries of the multi matrix A can');
disp('be overwritten by a double matrix as');
A(1:3,1:3)=[1:3; 4:6; 7:9]
disp('The double matrix D, whose entries are overwritten by those');
disp('of the multi matrix A, dose NOT become a multi matrix.');
disp('Since the lhs is double, the rhs. is casted to doube before copying.');
disp('Note that the matrix D is still double.');
D(1:3,1:3)=A(1:3,1:3)
disp('The multi array can be allocated as 1, 2 and 3 dimensions.');
A=zeros(1,3,'multi')
A=zeros(3,1,'multi')
A=zeros(3,3,'multi')
A=zeros(3,3,2,'multi')
disp('The multi array can be initialized by functions zeros(), ones(), eye() and rand().');
A=zeros(1,3,'multi')
A=ones(1,3,'multi')
%A=eye(3,3,'multi')
A=rand(3,3,'multi')
disp('The multi array can be obtained by function multi() from a double array.');
D=[1.245678901234567890123456789 -1; -2 3]
A=multi(D)
disp('The multi vector can be initialized by string.');
S='1 2 3 4 5 6 7 8 9'
A=multi(S)
disp('The multi matrix can be initialized by string.');
S='1.2456789012345678901234567890 -1; -2 3'
A=multi(S)
disp('The multi array can be transformed to double array by function double().');
D=double(A)
disp('The multi array can be transformed to string.');
S=num2str(A,'%12.7f ')
disp('The multi array can be transformed to cell of string.');
s=get_s(A,'%.30e')





