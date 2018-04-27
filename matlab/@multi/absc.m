% y=absc(x)
function y=absc(x)
cmd='absc';
if isa(x,'multi')
    y=multi(cmd,x.data);
else
    y=multi(cmd,multi(x).data);
end