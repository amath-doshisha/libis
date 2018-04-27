%% y=~x
function y=not(x)
cmd='not';
if isa(x,'multi')
    y=multi(cmd,x.data).data;
else
    y=multi(cmd,multi(x).data).data;
end