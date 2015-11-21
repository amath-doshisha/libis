%% y=get_prec(x)
function y=get_prec(x)
cmd='get_prec';
if isa(x,'multi')
    y=multi(cmd,x.data).data;
else
    y=multi(cmd,multi(x).data).data;
end