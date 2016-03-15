% y=angle(x)
function y=angle(x)
cmd='angle';
if isa(x,'multi')
    y=multi(cmd,x.data);
else
    y=multi(cmd,multi(x).data);
end