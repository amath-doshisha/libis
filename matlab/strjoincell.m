function s=strjoincell(c,sp)
S=size(c);
M=S(1);
N=prod(S(2:end));
c=reshape(c,[M,N]);
j=1;
while j<=N
    l=0;
    i=1;
    while i<=M        
        a=length(c{i,j});
        if a>l
            l=a;
        end
        i=i+1;
    end
    i=1;
    while i<=M
        a=length(c{i,j});
        k=1;
        while k<=l-a
            c{i,j}=[' ' c{i,j}];
            k=k+1;
        end        
        i=i+1;
    end    
    j=j+1;
end
s={};
i=1;
while i<=M
    s{i}=c{i,1};
    j=2;
    while j<=N
        s{i}=[s{i} sp c{i,j}];
        j=j+1;
    end
    i=i+1;
end
s=char(s);