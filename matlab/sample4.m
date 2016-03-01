%function A=sample4(lambda,M)
clear all

%% initializations
M=3;
lambda=[1 2 3 4 5 6 7 8 9 10];
m=size(lambda,2);
c=ones(1,m);
q=[];
e=[];
R=[];
L=[];
A=[];
%% initializations for multi precision arithmetic
set_default_prec(128);
lambda=multi(lambda);
q=multi(q);
e=multi(e);
R=multi(R);
L=multi(L);
A=multi(A);
%% operations
sigma=lambda.^(1/M);
n=0;
while n<=(M+1)*(m-1)+M
    f(n+1)=c*(sigma.^n)';
    n=n+1;
end
n=0;
while n<=(M+1)*(m-1)+M-1
    e(1,n+1)=0;
    q(1,n+1)=f(n+2)/f(n+1);
    n=n+1;
end
k=1;
while k<=m-1
    n=0;
    while n<=(M+1)*(m-k-1)+M
        e(k+1,n+1)=q(k,n+M+1)-q(k,n+1)+e(k,n+2);
        n=n+1;
    end
    n=0;
    while n<=(M+1)*(m-k-1)+M
        q(k+1,n+1)=(e(k+1,n+2)/e(k+1,n+1))*q(k,n+M+1);
        n=n+1;
    end
    k=k+1;
end
L(:,:,1)=sample_bidiag_lower(e(2:end,1));
n=0;
while n<=M-1
    R(:,:,n+1)=sample_bidiag_upper(q(:,n+1));
    n=n+1;
end
A=L;
n=0;
while n<=M-1
    A=A*R(:,:,M-1-n+1);
    n=n+1;
end
format short
A
format long
sort(eig(double(A)))
