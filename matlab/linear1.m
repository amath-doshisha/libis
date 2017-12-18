clear all
close all

N=200;prec=53;

set_default_prec(prec);    
A=rand(N);
A=multi(A);
b=rand(N,1);
B=multi(b);
R=inv(A);
x=R*b;

%Gd=0;Gd=multi(Gd);
%rd=0;rd=multi(rd);
%Gu=0;Gu=multi(Gu);
%ru=0;ru=multi(ru);
%down();
set_default_round_mode(-1)
Gd=abs(R*A-eye(N));
rd=abs(A*x-b);
%up();
set_default_round_mode(1)
Gu=abs(R*A-eye(N));
ru=abs(A*x-b);

%Gu=max(Gd,Gu);
R_norm=max(sum(R));
G_norm=max(sum(Gu));
r_norm=max(sum(ru));
%down();
set_default_round_mode(-1)
D=1-G_norm;
if D>0
    %up();
    set_default_round_mode(1)
    A_inv=R_norm/D
    error=A_inv*r_norm
else
    disp('false')
end