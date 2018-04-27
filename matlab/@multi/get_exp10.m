%% z=get_exp10(x,offset)
function z=get_exp10(x,offset)
cmd='get_exp10';
if nargin<=1
    if isa(x,'multi')
        z=multi(cmd,x.data,0);
    else
        z=multi(cmd,multi(x).data,0);
    end
elseif nargin==2
    if isa(x,'multi')
        z=multi(cmd,x.data,offset);
    else
        z=multi(cmd,multi(x).data,offset);
    end
else
    error('get_exp10:nargin');
end
