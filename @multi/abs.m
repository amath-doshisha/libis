% y=abs(x)
function y=abs(x)
cmd='abs';
if isa(x,'multi')
    y=multi(cmd,x.data);
else
    y=multi(cmd,multi(x).data);
end