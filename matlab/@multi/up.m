% y=up(x)
function y=up(x)
cmd='up';
if isa(x,'multi')
    y=multi(cmd,x.data);
else
    y=multi(cmd,multi(x).data);
end