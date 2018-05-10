clear all
close all

A=zeros(3,1);
A(:)=1:length(A)
B=zeros(3,3,3)
A=multi(A);
B=multi(B);
A+B
B+A

% A = [1 2 3; 4 5 6]
% b = [10; 100]
% B=zeros(2,3,2,2);
% size(B)
% C=zeros(2,3,2,2,2);
% size(C)
% % A+b
% % b+b'
% % B+A
% % B+b
% D=B+C;
% size(D)

% A=multi([1 2 3; 4 5 6])
% b=multi([10; 100])
% B=zeros(2,3,2,2)
% A+b
% b+b'
% B+A
% B+b

% Az=rand(3,3);
% Az(:)=(1:9)-i*(9:-1:1)
% Ad=real(Az)
% Ar=multi(Ad)
% AR=imulti(Ad)
% Ac=multi(Az)
% AC=imulti(Az)
% 
% disp('----------------');
% % Bd=multi_cast(Ad,'d')
% % Bz=multi_cast(Ad,'z')
% % Br=multi_cast(Ad,'r')
% % Bc=multi_cast(Ad,'c')
% % BR=multi_cast(Ad,'R')
% % BC=multi_cast(Ad,'C')
% % disp('----------------');
% % Bd=multi_cast(Ar,'d')
% % Bz=multi_cast(Ar,'z')
% % Br=multi_cast(Ar,'r')
% % Bc=multi_cast(Ar,'c')
% % BR=multi_cast(Ar,'R')
% % BC=multi_cast(Ar,'C')
% % disp('----------------');
% % Bd=multi_cast(AR,'d')
% % Bz=multi_cast(AR,'z')
% % Br=multi_cast(AR,'r')
% % Bc=multi_cast(AR,'c')
% % BR=multi_cast(AR,'R')
% % BC=multi_cast(AR,'C')
% % disp('----------------');
% % Bd=multi_cast(Az,'d')
% % Bz=multi_cast(Az,'z')
% % Br=multi_cast(Az,'r')
% % Bc=multi_cast(Az,'c')
% % BR=multi_cast(Az,'R')
% % BC=multi_cast(Az,'C')
% % disp('----------------');
% % Bd=multi_cast(Ac,'d')
% % Bz=multi_cast(Ac,'z')
% % Br=multi_cast(Ac,'r')
% % Bc=multi_cast(Ac,'c')
% % BR=multi_cast(Ac,'R')
% % BC=multi_cast(Ac,'C')
% disp('----------------');
% Bd=multi_cast(AC,'d')
% Bz=multi_cast(AC,'z')
% Br=multi_cast(AC,'r')
% Bc=multi_cast(AC,'c')
% BR=multi_cast(AC,'R')
% BC=multi_cast(AC,'C')
% % 
% % disp('----------------');
% % double(Ad+Bd)-2*Ad
% % double(Ad+Bz)-(Ad+Bz)
% % double(Ad+Br)-2*Ad
% % double(Ad+Bc)-(Ad+Bz)
% % double(Ad+BR)-2*Ad
% % double(Ad+BC)-(Ad+Bz)
% % disp('----------------');
% % double(Az+Bd)-(Az+Bd)
% % double(Az+Bz)-(Az+Bz)
% % double(Az+Br)-(Az+Bd)
% % double(Az+Bc)-(Az+Bz)
% % double(Az+BR)-(Az+Bd)
% % double(Az+BC)-(Az+Bz)
% % disp('----------------');
% % double(Ar+Bd)-(Ad+Bd)
% % double(Ar+Bz)-(Ad+Bz)
% % double(Ar+Br)-(Ad+Bd)
% % double(Ar+Bc)-(Ad+Bz)
% % double(Ar+BR)-(Ad+Bd)
% % double(Ar+BC)-(Ad+Bz)
% % disp('----------------');
% % double(Ac+Bd)-(Az+Bd)
% % double(Ac+Bz)-(Az+Bz)
% % double(Ac+Br)-(Az+Bd)
% % double(Ac+Bc)-(Az+Bz)
% % double(Ac+BR)-(Az+Bd)
% % double(Ac+BC)-(Az+Bz)
% % disp('----------------');
% % double(AR+Bd)-(Ad+Bd)
% % double(AR+Bz)-(Ad+Bz)
% % double(AR+Br)-(Ad+Bd)
% % double(AR+Bc)-(Ad+Bz)
% % double(AR+BR)-(Ad+Bd)
% % double(AR+BC)-(Ad+Bz)
% % disp('----------------');
% % double(AC+Bd)-(Az+Bd)
% % double(AC+Bz)-(Az+Bz)
% % double(AC+Br)-(Az+Bd)
% % double(AC+Bc)-(Az+Bz)
% % double(AC+BR)-(Az+Bd)
% % double(AC+BC)-(Az+Bz)
% % disp('----------------');
% 
% ad=1
% ar=multi(ad)
% aR=imulti(1,1.01)
% az=1+i
% ac=multi(az)
% aC=imulti(1+i,1.01+1.01*i)
% 
% disp('----------------');
% (Ad-ad)+(ad-Ad)
% (Ad-ar)+(ar-Ad)
% (Ad-aR)+(aR-Ad)
% (Ad-az)+(az-Ad)
% (Ad-ac)+(ac-Ad)
% (Ad-aC)+(aC-Ad)
% disp('----------------');
% (Az-ad)+(ad-Az)
% (Az-ar)+(ar-Az)
% (Az-aR)+(aR-Az)
% (Az-az)+(az-Az)
% (Az-ac)+(ac-Az)
% (Az-aC)+(aC-Az)
% disp('----------------');
% (Ar-ad)+(ad-Ar)
% (Ar-ar)+(ar-Ar)
% (Ar-aR)+(aR-Ar)
% (Ar-az)+(az-Ar)
% (Ar-ac)+(ac-Ar)
% (Ar-aC)+(aC-Ar)
% disp('----------------');
% (Ac-ad)+(ad-Ac)
% (Ac-ar)+(ar-Ac)
% (Ac-aR)+(aR-Ac)
% (Ac-az)+(az-Ac)
% (Ac-ac)+(ac-Ac)
% (Ac-aC)+(aC-Ac)
% disp('----------------');
% (AR-ad)+(ad-AR)
% (AR-ar)+(ar-AR)
% (AR-aR)+(aR-AR)
% (AR-az)+(az-AR)
% (AR-ac)+(ac-AR)
% (AR-aC)+(aC-AR)
% disp('----------------');
% (AC-ad)+(ad-AC)
% (AC-ar)+(ar-AC)
% (AC-aR)+(aR-AC)
% (AC-az)+(az-AC)
% (AC-ac)+(ac-AC)
% (AC-aC)+(aC-AC)
% disp('----------------');

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
