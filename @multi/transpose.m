%% y=x.'
function y=transpose(x)
cmd='transpose';
if isa(x,'multi')
    y=multi(cmd,x.data);
else
    y=multi(cmd,multi(x).data);
end