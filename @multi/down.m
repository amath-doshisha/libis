% y=down(x)
function y=down(x)
cmd='down';
if isa(x,'multi')
    y=multi(cmd,x.data);
else
    y=multi(cmd,multi(x).data);
end