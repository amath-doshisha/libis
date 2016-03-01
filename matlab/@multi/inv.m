%% y=inv(x)
function y=inv(x)
cmd='inv';
if isa(x,'multi')
    y=multi(cmd,x.data);
else
    y=multi(cmd,multi(x).data);
end
