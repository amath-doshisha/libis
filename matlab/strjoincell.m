function s=strjoincell(c)
S=size(c);
M=S(1);
if length(S)>=2
    N=S(2);
else
    N=1;
end
if length(S)>=3
    L=S(3);
else
    L=1;
end
s={};
i=1;
while i<=M
    s{i}=c{i,1};
    j=2;
    while j<=N*L
        s{i}=[s{i} '  ' c{i,j}];
        j=j+1;
    end
    i=i+1;
end
s=char(s);