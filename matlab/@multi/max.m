%% z=max(x)
function z=max(x,y)
cmd='max';
if nargin==1
    if isa(x,'multi')
        z=multi(cmd,x.data);
    else
        z=multi(cmd,multi(x).data);
    end
elseif nargin==2
    if isa(x,'multi') && isa(y,'multi')
        z=multi(cmd,x.data,y.data);
    elseif isa(x,'multi')
        z=multi(cmd,x.data,multi(y).data);
    elseif isa(y,'multi')
        z=multi(cmd,multi(x).data,y.data);
    else
        z=multi(cmd,multi(x).data,multi(y).data);
    end
else
    error('max:nargin');
end
