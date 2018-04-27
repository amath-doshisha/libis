% y=imag(x)
function y=imag(x)
cmd='imag';
if isa(x,'multi')
    y=multi(cmd,x.data);
else
    y=multi(cmd,multi(x).data);
end