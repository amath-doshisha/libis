% y=mid(x)
function y=mid(x)
cmd='mid';
if isa(x,'multi')
    y=multi(cmd,x.data);
else
    y=multi(cmd,multi(x).data);
end