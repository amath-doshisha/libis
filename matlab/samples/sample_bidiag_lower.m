function L=sample_bidiag_lower(e)
m=length(e)+1;
L=zeros(m,m,'like',e);
i=1;
while i<=m
    L(i,i)=1;
    i=i+1;
end
i=1;
while i<=m-1
    L(i+1,i)=e(i);
    i=i+1;
end
