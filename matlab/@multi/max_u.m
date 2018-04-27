%% z=max_u(x)
function z=max_u(x)
cmd='max_u';
if isa(x,'multi')
    z=multi(cmd,x.data);
else
    z=multi(cmd,multi(x).data);
end
