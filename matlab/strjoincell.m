function s=strjoincell(c)
s={};
i=1;
while i<=size(c,1)
    s{i}=c{i,1};
    j=2;
    while j<=size(c,2)
        s{i}=[s{i} c{i,j}];
        j=j+1;
    end
    i=i+1;
end
s=char(s);