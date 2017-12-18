% y=real(x)
function y=real(x)
cmd='real';
if isa(x,'multi')
    y=multi(cmd,x.data);
else
    y=multi(cmd,multi(x).data);
end