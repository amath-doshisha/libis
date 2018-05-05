function disp_display(x,name)
% display format
if strcmp(get(0,'Format'),'long')
    f='f';
    digits=get_multi_disp_digits();
    if get_type(x)=='R' || get_type(x)=='c' || get_type(x)=='z'
        spc='  ';
    elseif get_type(x)=='C'
        spc=' ';
    else
        spc='   ';
    end
else
    f='f';
    digits=4;
    if get_type(x)=='d' || get_type(x)=='r'
        spc='    ';
    elseif get_type(x)=='R' || get_type(x)=='z' || get_type(x)=='c'
        spc='   ';
    elseif get_type(x)=='C'
        spc='  ';
    else
        spc='    ';
    end
end
% c=get_s(x)
c=get_s(x,f,digits);
% get size
S=size(c);
M=S(1);
N=S(2);
if length(S)>=3
    L=prod(S(3:end));
else
    L=1;
end
% reshap
c=reshape(c,[M,N*L]);
% output
if L<=1 && nargin>=2
    disp(' ');
    disp([name,' = '])
    disp(' ');
end
% loop for 3d-array
t=1;
while t<=L
    % loop for 2d-array
    % get width of column
    j=1;
    while j<=N
        l=0;
        i=1;
        while i<=M
            a=length(c{i,j+(t-1)*N});
            if a>l
                l=a;
            end
            i=i+1;
        end
        i=1;
        while i<=M
            a=length(c{i,j+(t-1)*N});
            k=1;
            while k<=l-a
                c{i,j+(t-1)*N}=[' ' c{i,j+(t-1)*N}];
                k=k+1;
            end
            i=i+1;
        end
        j=j+1;
    end
    % add spaces to column
    s={};
    i=1;
    while i<=M
        s{i}='';
        j=1;
        while j<=N
            s{i}=[s{i} spc c{i,j+(t-1)*N}];
            j=j+1;
        end
        i=i+1;
    end
    % output column
    if L>=2
        index=[];
        p=t-1;
        r=3;
        while r<=length(S)
            index(r-2)=mod(p,S(r))+1;
            p=floor(p/S(r));
            r=r+1;
        end
        p=[',',num2str(index(1))];
        r=2;
        while r<=length(index)
            p=[p,',',num2str(index(r))];
            r=r+1;
        end
        disp(' ');
        if nargin>=2
            disp([name,'(:,:,',p,') = '])
        else
            disp(['(:,:,',p,') = '])
        end
        disp(' ');
    end
    s=char(s);
    disp(s);
    disp(' ');
    t=t+1;
end