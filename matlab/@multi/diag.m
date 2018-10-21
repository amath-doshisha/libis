function y=diag(x,k)

if nargin<2
    k=0;
end

s=size(x);
if length(s)>=3
    error('Å‰‚Ì“ü—Í‚Í 2 ŽŸŒ³‚Å‚È‚¯‚ê‚Î‚È‚è‚Ü‚¹‚ñB');
end

% vector -> matrix
if s(1)==1 || s(2)==1
    n=length(x);
    m=n+abs(k);
    y=zeros(m,m,'like',x);
    i=1;
    while i<=n
        if k<0
            y.data(i-k,i)=x.data(i);
        else
            y.data(i,i+k)=x.data(i);
        end
        i=i+1;
    end
else
    if k<0
        i0=1-k;
        j0=1;
    else
        i0=1;
        j0=1+k;
    end
    n=0;
    while i0+n<=s(1) && j0+n<=s(2)
        n=n+1;
    end
    y=zeros(n,1,'like',x);
    i=1;
    while i<=n
        y.data(i)=x.data(i0+i-1,j0+i-1);
        i=i+1;
    end
end

