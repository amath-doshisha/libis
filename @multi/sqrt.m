% y=sqrt(x)
function y=sqrt(x)
cmd='sqrt';
if isa(x,'multi')
    y=multi(cmd,x.data);
else
    y=multi(cmd,multi(x).data);
end