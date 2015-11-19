% y=-x
function y=uminus(x)
if isa(x,'multi')
    y=multi('uminus',x.data);
else
    y=multi('uminus',multi(x).data);
end
