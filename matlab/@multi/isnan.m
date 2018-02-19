function y=isnan(x)
cmd='isnan';
if isa(x,'multi')
    y=logical(multi(cmd,x.data).data);
else
    y=logical(multi(cmd,multi(x).data).data);
end