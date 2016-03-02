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
disp('You can use the class multi as 1, 2 and 3 dimensional arrays.');
A=zeros(1,3,'multi')
A=zeros(3,1,'multi')
A=zeros(3,3,'multi')
A=zeros(3,3,2,'multi')

