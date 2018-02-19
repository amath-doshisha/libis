% y=mag(x)
function y=mag(x)
cmd='mag';
if isa(x,'multi')
    y=multi(cmd,x.data);
else
    y=multi(cmd,multi(x).data);
end