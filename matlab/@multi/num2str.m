% s=num2str(A,opt)
function s=num2str(A,opt)
cmd='num2str';
digits=6;
f='g';
if isa(A,'multi')
    if nargin==2 && isa(opt,'char')
        if length(strfind(opt,'f'))>0
            f='f';
        end
        if length(strfind(opt,'e'))>0
            f='e';
        end
        if length(strfind(opt,'g'))>0
            f='g';
        end        
        a=strsplit(opt,'.');
        if length(a)>=2
            a=regexprep(a{2},'[^\d]+','');
            digits=str2num(a);
        end        
        sp='';
        a=length(strfind(opt,' '));
        i=1;
        while i<=a
            sp=[' ' sp];
            i=i+1;
        end
        c=get_s(A,f,digits);
        s=strjoincell(c,sp);
    elseif nargin==2
        sp='  ';
        f='g';        
        digits=sum(double(opt));
        c=get_s(A,f,digits);
        s=strjoincell(c,sp);
    else
        sp='  ';
        c=get_s(A,f,digits);
        s=strjoincell(c,sp);
    end
else
    if nargin==1
        s=num2str(A);
    else
        s=num2str(A,opt);
    end
end
