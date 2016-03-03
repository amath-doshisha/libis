function c=strcellchomp(c)
i=1;
while i<=size(c,1)
    j=1;
    while j<=size(c,2)
        s=c{i,j};
        while s(1)==' '
            s=s(2:end);
        end
        while s(end)==' '
            s=s(1:end-1);
        end
        c{i,j}=s;
        j=j+1;
    end
    i=i+1;
end

