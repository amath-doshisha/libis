%% y=+x
function y=uplus(x)
cmd='uplus';
if isa(x,'multi')
    y=multi(cmd,x.data);
else
    y=multi(cmd,multi(x).data);
end