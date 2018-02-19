%% z=norm_Inf(x)
function z=norm_Inf(x)
cmd='norm_Inf';
if isa(x,'multi')
    z=multi(cmd,x.data);
else
    z=multi(cmd,multi(x).data);
end
