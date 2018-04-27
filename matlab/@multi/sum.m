%% z=sum(x)
function z=sum(x)
cmd='sum';
if isa(x,'multi')
    z=multi(cmd,x.data);
else
    z=multi(cmd,multi(x).data);
end
