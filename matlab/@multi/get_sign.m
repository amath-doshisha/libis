%% y=get_sign(x)
function y=get_sign(x)
cmd='get_sign';
if isa(x,'multi')
    y=multi(cmd,x.data).data;
else
    y=multi(cmd,multi(x).data).data;
end