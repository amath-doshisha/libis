% y=-x
function y=uminus(x)
cmd='uminus';
if isa(x,'multi')
    y=multi(cmd,x.data);
else
    y=multi(cmd,multi(x).data);
end