% y=mig(x)
function y=mig(x)
cmd='mig';
if isa(x,'multi')
    y=multi(cmd,x.data);
else
    y=multi(cmd,multi(x).data);
end
