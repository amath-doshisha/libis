%% z=min(x)
function z=min(x)
cmd='min';
if isa(x,'multi')
    z=multi(cmd,x.data);
else
    z=multi(cmd,multi(x).data);
end
