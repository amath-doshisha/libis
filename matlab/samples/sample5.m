clear all
close all
A = [1 2 3; 4 5 6]
b = [10; 100]
B=zeros(2,3,2,2);
size(B)
C=zeros(2,3,2,2,2);
size(C)
% A+b
% b+b'
% B+A
% B+b
D=B+C;
size(D)

% A=multi([1 2 3; 4 5 6])
% b=multi([10; 100])
% B=zeros(2,3,2,2)
% A+b
% b+b'
% B+A
% B+b

% A=zeros(4,3);A(:)=1:12
% B=zeros(1,3,3);B(:)=1:9
% A+B

% Ad=rand(3,3)
% Ar=multi(Ad)
% AR=imulti(Ad)
% Az=Ad+i*rand(3,3)
% Ac=multi(Az)
% AC=imulti(Az)
% 
% B=Ar;
% B(:,2)=Ar(:,3)
% B(:,2)=Ad(:,3)
% B(:,2)=Az(:,3)
% 
% 
% disp('----------------');
% B=Ad
% Bd=multi_cast(B,'d')
% Bz=multi_cast(B,'z')
% Br=multi_cast(B,'r')
% Bc=multi_cast(B,'c')
% BR=multi_cast(B,'R')
% BC=multi_cast(B,'C')
% 
% 
% disp('----------------');
% B=Az
% Bd=multi_cast(B,'d')
% Bz=multi_cast(B,'z')
% Br=multi_cast(B,'r')
% Bc=multi_cast(B,'c')
% BR=multi_cast(B,'R')
% BC=multi_cast(B,'C')
% 
% 
% disp('----------------');
% B=Ar
% Bd=multi_cast(B,'d')
% Bz=multi_cast(B,'z')
% Br=multi_cast(B,'r')
% Bc=multi_cast(B,'c')
% BR=multi_cast(B,'R')
% BC=multi_cast(B,'C')
% 
% 
% disp('----------------');
% B=Ac
% Bd=multi_cast(B,'d')
% Bz=multi_cast(B,'z')
% Br=multi_cast(B,'r')
% Bc=multi_cast(B,'c')
% BR=multi_cast(B,'R')
% BC=multi_cast(B,'C')
% 
% 
% disp('----------------');
% B=AR
% Bd=multi_cast(B,'d')
% Bz=multi_cast(B,'z')
% Br=multi_cast(B,'r')
% Bc=multi_cast(B,'c')
% BR=multi_cast(B,'R')
% BC=multi_cast(B,'C')
% 
% disp('----------------');
% B=AC
% Bd=multi_cast(B,'d')
% Bz=multi_cast(B,'z')
% Br=multi_cast(B,'r')
% Bc=multi_cast(B,'c')
% BR=multi_cast(B,'R')
% BC=multi_cast(B,'C')

% disp('----------------');
% A=rand(3,1)
% disp('----------------');
% B=multi(A)
% disp('----------------');
% A=rand(1,3)
% disp('----------------');
% B=multi(A)
% disp('----------------');
% 
% format short
% A=rand(3,3)
% B=multi(A)
% C=imulti(A)
% 
% 
% format short
% A=rand(3,3)+i*rand(3,3)
% B=multi(A)
% C=imulti(A)
% 
% set_multi_disp_digits(15);
% format long
% A=rand(2,2)
% B=multi(A)
% C=imulti(A)
% 
% set_multi_disp_digits(15);
% format long
% A=rand(2,2)+i*rand(2,2)
% B=multi(A)
% C=imulti(A)
% 
% format short
% disp('----------------');
% A=rand(3,3,3)
% disp('----------------');
% B=multi(A)
% disp('----------------');
% C=imulti(A)
% disp('----------------');
% 
% 
% format short
% disp('----------------');
% A=rand(3,3,3,2)
% disp('----------------');
% B=multi(A)
% disp('----------------');
% C=imulti(A)
% disp('----------------');
% 
% 
% 
% disp('----------------');
% A=rand(3,1);
% disp(A)
% disp('----------------');
% B=multi(A);
% disp(B)
% disp('----------------');
% A=rand(1,3);
% disp(A);
% disp('----------------');
% B=multi(A);
% disp(B)
% disp('----------------');
% 
% 
% format short
% disp('----------------');
% A=rand(3,3,3,2);
% disp(A);
% disp('----------------');
% B=multi(A);
% disp(B);
% disp('----------------');
% C=imulti(A);
% disp(C);
% disp('----------------');
% 


% format short
% A=rand(3,3,3)
% B=multi(A)
%C=imulti(A)


% A=zeros(3);
% A(:)=1:9;
% B=zeros(3);
% B(:)=10:18;
% % 
% % complex(multi(A))
% % complex(multi(B))
% % 
% % complex(multi(A),multi(B))
% complex(imulti(A,B),imulti(B,A))
% 
% 
% 
% C=multi(A)
% get_type(C)
% Cr=multi_cast(C,'r')
% Cc=multi_cast(C,'c')
% Cd=multi_cast(C,'d')
% Cz=multi_cast(C,'z')
% CR=multi_cast(C,'R')
% CC=multi_cast(C,'C')
% 
% disp('------------')
% 
% C=CR
% get_type(C)
% Cr=multi_cast(C,'r')
% Cc=multi_cast(C,'c')
% Cd=multi_cast(C,'d')
% Cz=multi_cast(C,'z')
% CR=multi_cast(C,'R')
% CC=multi_cast(C,'C')
% 
% disp('------------')
% 
% C=multi(A+i*B)
% get_type(C)
% Cr=multi_cast(C,'r')
% Cc=multi_cast(C,'c')
% Cd=multi_cast(C,'d')
% Cz=multi_cast(C,'z')
% CR=multi_cast(C,'R')
% CC=multi_cast(C,'C')
% 
% disp('------------')
% 
% C=imulti(A+i*B)
% get_type(C)
% Cr=multi_cast(C,'r')
% Cc=multi_cast(C,'c')
% Cd=multi_cast(C,'d')
% Cz=multi_cast(C,'z')
% CR=multi_cast(C,'R')
% CC=multi_cast(C,'C')
% 
% disp('------------')
% 
% C=A
% get_type(C)
% Cr=multi_cast(C,'r')
% Cc=multi_cast(C,'c')
% Cd=multi_cast(C,'d')
% Cz=multi_cast(C,'z')
% CR=multi_cast(C,'R')
% CC=multi_cast(C,'C')
% 
% disp('------------')
% 
% C=A+i*B
% get_type(C)
% Cr=multi_cast(C,'r')
% Cc=multi_cast(C,'c')
% Cd=multi_cast(C,'d')
% Cz=multi_cast(C,'z')
% CR=multi_cast(C,'R')
% CC=multi_cast(C,'C')
% 
% disp('------------')
% 
% A=zeros(3)
% A(:)=1:9
% B=zeros(3)
% B(:)=10:18
% 
% disp('imulti(d)')
% imulti(A)
% 
% disp('imulti(r)')
% imulti(multi(A))
% 
% disp('imulti(z)')
% imulti(A+i*B)
% 
% disp('imulti(c)')
% imulti(multi(A+i*B))
% 
% disp('imulti(z,d)')
% imulti(A-i*B,B)
% 
% disp('imulti(d,z)')
% imulti(A,B-i*A)
% 
% disp('imulti(z,z)')
% imulti(A-i*B,B-i*A)
% 
% disp('imulti(c,c)')
% imulti(multi(A-i*B),multi(B-i*A))
% 
% disp('imulti(d,d)')
% imulti(A,B)
% 
% disp('imulti(r,r)')
% imulti(multi(A),multi(B))
% 
% disp('imulti(r,c)')
% imulti(multi(A),multi(B-i*A))
% 
% disp('imulti(c,r)')
% imulti(multi(A-i*B),multi(B))
% 
% disp('------------')
% 
% 
