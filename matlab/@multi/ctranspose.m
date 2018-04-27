%% y=x'
function y=ctranspose(x)
cmd='ctranspose';
if isa(x,'multi')
    y=multi(cmd,x.data);
else
    y=multi(cmd,multi(x).data);
end