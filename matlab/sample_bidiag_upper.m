function R=sample_bidiag_upper(q)
m=length(q);
R=zeros(m,m,class(q)); % point!
i=1;
while i<=m
    R(i,i)=q(i);
    i=i+1;
end
i=1;
while i<=m-1
    R(i,i+1)=1;
    i=i+1;
end