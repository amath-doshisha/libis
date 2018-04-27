function y=imulti(x,u)
if nargin==1
    if isa(x,'double')
        y=multi('iset_d',x);
    elseif isnumeric(x) || islogical(x)
        y=multi('iset_d',double(x));
    elseif isa(x,'multi')
        y=multi('icopy',x.data);
    else
        y=multi('icopy',multi(x).data);
    end
    %’Ç‰Á
elseif nargin==2
    if isa(x,'double') && isa(u,'double') 
        y=multi('iset_dd',x,u);
    elseif (isnumeric(x) || islogical(x)) &&  (isnumeric(u) || islogical(u))
        y=multi('iset_dd',double(x),double(u));
    elseif (isa(x,'multi') && isa(u,'multi'))
        y=multi('icopy2',x.data,u.data);
    else
        error('imulti:not supported');
    end
else
    error('imulti:nargin');
end