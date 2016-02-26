%% z=max(x)
function z=max(x)
cmd='max';
if isa(x,'multi')
    z=multi(cmd,x.data);
else
    z=multi(cmd,multi(x).data);
end
