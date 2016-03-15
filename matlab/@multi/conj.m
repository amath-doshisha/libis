% y=conj(x)
function y=conj(x)
cmd='conj';
if isa(x,'multi')
    y=multi(cmd,x.data);
else
    y=multi(cmd,multi(x).data);
end