% y=log10(x)
function y=log10(x)
cmd='log10';
if isa(x,'multi')
    y=multi(cmd,x.data);
else
    y=multi(cmd,multi(x).data);
end