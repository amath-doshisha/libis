function [A,L,R,q,e]=sample_matTN(lambda,N,M)
%% initializations
m=length(lambda);
nmax=(M+N)*(m-1)+M*N;
f=zeros(1,nmax+1,'like',lambda);
c=ones(1,m,'like',lambda);
q=zeros(m,nmax+1,'like',lambda);
e=zeros(m,nmax+1,'like',lambda);
R=zeros(m,m,M*(N-1)+1,'like',lambda);
L=zeros(m,m,N*(M-1)+1,'like',lambda);
A=zeros(m,m,'like',lambda);
%% operations
sigma=lambda.^(1/(N*M));
n=0;
while n<=(M+N)*(m-1)+M*N
    f(n+1)=c*(sigma.^n)';
    n=n+1;
end
n=0;
while n<=(M+N)*(m-1)+(N-1)*M
    e(1,n+1)=0;
    n=n+1;
end
n=0;
while n<=(M+N)*(m-1)+(M-1)*N
    q(1,n+1)=f(n+N+1)/f(n+1);
    n=n+1;
end
k=1;
while k<=m-1
    n=0;
    while n<=(M+N)*(m-k-1)+M*N
        e(k+1,n+1)=q(k,n+M+1)-q(k,n+1)+e(k,n+N+1);
        n=n+1;
    end
    n=0;
    while n<=(M+N)*(m-k-1)+(M-1)*N
        q(k+1,n+1)=(e(k+1,n+N+1)/e(k+1,n+1))*q(k,n+M+1);
        n=n+1;
    end
    k=k+1;
end
n=0;
while n<=M*(N-1)
    L(:,:,n+1)=sample_bidiag_lower(e(2:end,n+1));
    n=n+1;
end
n=0;
while n<=N*(M-1)
    R(:,:,n+1)=sample_bidiag_upper(q(:,n+1));
    n=n+1;
end
A=eye(m);
n=0;
while n<=N-1
    A=A*L(:,:,M*n+1);
    n=n+1;
end
n=0;
while n<=M-1
    A=A*R(:,:,N*(M-1-n)+1);
    n=n+1;
end