%% y=get_exp(x)
function y=get_exp(x)
cmd='get_exp';
if isa(x,'multi')
    y=multi(cmd,x.data).data;
else
    y=multi(cmd,multi(x).data).data;
end