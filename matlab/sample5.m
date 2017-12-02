clear all
close all

type='double';
%type='multi';
set_default_prec(128);
N=20;
rand('seed',0);
A=rand(N,type);
b=rand(N,1,type);
set_default_round_mode(0);
R=inv(A);
x=R*b;
set_default_round_mode(-1);
Gd=abs(R*A-eye(N,type));
rd=A*x-b
set_default_round_mode(1);
Gu=abs(R*A-eye(N,type));
ru=A*x-b
center=(ru+rd)/2;
radius=center-rd;
ru=max(rd,ru);
Gu=max(Gd,Gu);
R_norm=norm(R,inf);
G_norm=norm(Gu,inf);
Rru=abs(R*center+abs(R)*radius);
set_default_round_mode(-1);
D=1-G_norm;
if D>0
    set_default_round_mode(1);
    Rru=max(Rrd,Rru);
    Rr_norm=norm(Rru,inf);
    error=Rr_norm/D;
    kappa=Gu*ones(N,1,type);
    error2=Rru+(error*kappa);
    cwerror=max(error2);
else
    disp('false');
end