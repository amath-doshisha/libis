%% z=abs(x)
function z=abs(x)
cmd='abs';
if isa(x,'multi')
    z=multi(cmd,x.data);
else
    z=multi(cmd,multi(x).data);
end
