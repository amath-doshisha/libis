%function sample0
clear all
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
disp('The multi array can be allocated as any dimensions, e.g.,');
A=zeros(1,3,'multi')
A=zeros(3,1,'multi')
A=zeros(3,3,'multi')
A=zeros(3,3,2,'multi')
disp('The multi array can be initialized by functions zeros(), ones(), eye() and rand().');
A=zeros(1,3,'multi')
A=ones(1,3,'multi')
A=eye(3,3,'multi')
A=rand(3,3,'multi')
disp('The multi array can be obtained by function multi() from a double array.');
A=multi([1.2456789 -1; -2 3])
A=multi([1+i -1-i; +i -i])
disp('The multi vector can be initialized by string.');
S='1 2 3 4'
A=multi(S)
disp('The multi matrix can be initialized by string.');
S='1.2456789 -1; -2 3'
A=multi(S)
disp('The multi matrix whose entries are complex can be initialized by string.');
S='1+1i  -1-1i; 0+1i 0-1i'
A=multi(S)
disp('[Important] The following is NOT correct.');
S='1+i -1-i; +i -i'
A=multi(S)
disp('The multi array can be transformed to double array by function double().');
A=multi(rand(3))
D=double(A)
disp('The multi array can be transformed to string.');
S=num2str(A)
S=num2str(A,'%.15e ')
disp('The multi array can be transformed to cell of string.');
s=get_s(A)
s=get_s(A,'e',15)
disp('The multi array can be displayed by long format mode.');
disp('Also, the displaid digits can be set with set_multi_disp_digits(digits).');
set_multi_disp_digits(20);
format long
A
disp('Then it is turned back to the short format mode.');
format short
A
disp('At the start, the default precision of multi class is 64-bit.');
disp('The default precision can be changed by set_default_prec(prec).');
set_default_prec(128);
A=rand(1,2,'multi')
disp('The precisions of values which stored to the variable can be given by get_prec().');
get_prec(A)
disp('If set_default_prec(prec) is called again, then the default precision is changed.');
set_default_prec(256);
B=rand(1,2,'multi')
get_prec(B)
disp('But, the precision of value which stored to variable never be changed.');
A
get_prec(A)
B
get_prec(B)
disp('Even in the case where the precisions of variables are different,');
disp('any calculations are executed by the default precision of set_default_prec().');
set_default_prec(512);
C=A+B
get_prec(C)

disp('')
disp('The values can be rounded with naer/up/down rouding mode.');
set_default_prec(64);
A=rand(1,2,'multi');
disp(['64-bit with near mode: ' num2str(A,'%+.10e ')])
set_default_prec(26);
set_default_round_mode(0);
B=multi(A);
disp(['26-bit with near mode: ' num2str(B,'%+.10e ')]);
set_default_round_mode(1);
B=multi(A);
disp(['26-bit with up mode:   ' num2str(B,'%+.10e ')]);
set_default_round_mode(-1);
B=multi(A);
disp(['26-bit with down mode: ' num2str(B,'%+.10e ')])

disp('   ');
disp('The operations can be executed under naer/up/down rounding mode.');
set_default_prec(26);
set_default_round_mode(0);
a=multi(1.23456789012345);
b=multi(7.89012345678901);
disp(sprintf('%.10f + %.10f=',a,b));
set_default_round_mode(0);
c=a+b;
disp(['     ' num2str(c,'%+.10e') '  (near mode)' ]);
set_default_round_mode(1);
c=a+b;
disp(['     ' num2str(c,'%+.10e') '  (up mode)' ]);
set_default_round_mode(-1);
c=a+b;
disp(['     ' num2str(c,'%+.10e') '  (down mode)']);

disp('');
disp('Interval arithmetic can work.');
disp('Double precision matrix is');
A=rand(3)
disp('Intrval multi-precision matrix is');
B=imulti(A)
disp('Intrval multi-precision matrix is');
A0=rand(3)
A1=rand(3)
B=imulti(A0,A1)

