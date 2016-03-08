function A=sample_matTN(lambda,M)
%% initializations
m=length(lambda);
c=ones(1,m,'like',lambda);
q=zeros(0,0,'like',lambda);
e=zeros(0,0,'like',lambda);
R=zeros(0,0,'like',lambda);
L=zeros(0,0,'like',lambda);
A=zeros(0,0,'like',lambda);
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