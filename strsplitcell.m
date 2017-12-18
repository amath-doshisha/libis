function c=strsplitcell(s)
if size(s,1)==1
    a=strsplit(s,';');
else
    a={};
    i=1;
    while i<=size(s,1)
        a{i}=s(i,:);
        i=i+1;
    end
end
c={};
i=1;
while i<=length(a)
    s=a{i};
    while s(1)==' ' || s(1)=='['
        s=s(2:end);
    end
    while s(end)==' ' || s(end)==']'
        s=s(1:end-1);
    end
    b=strsplit(s,{' ',','});
    j=1;
    while j<=length(b)
        c{i,j}=b{j};
        j=j+1;
    end
    i=i+1;
end