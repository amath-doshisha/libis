function y=complex(a,b)
cmd='complex';
if nargin==1
    if isa(a,'double')
        y=multi(cmd,a);
    elseif isa(a,'multi')
        y=multi(cmd,a.data);
    else
        error('complex(): Not supported');
    end
elseif nargin==2
    if isa(a,'double') && isa(b,'double')
        y=multi(cmd,a,b);
    elseif (isa(a,'multi') && isa(b,'multi'))
        y=multi(cmd,a.data,b.data);
    else
        error('complex(): Not supported');
    end
else
    error('complex()');
end