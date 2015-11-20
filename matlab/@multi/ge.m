%% z=(x>=y)
function z=ge(x,y)
cmd='ge';
if isa(x,'multi') && isa(y,'multi')
    z=multi(cmd,x.data,y.data).data;
elseif isa(x,'multi')
    z=multi(cmd,x.data,multi(y).data).data;
elseif isa(y,'multi')
    z=multi(cmd,multi(x).data,y.data).data;
else
    z=multi(cmd,multi(x).data,multi(y).data).data;
end