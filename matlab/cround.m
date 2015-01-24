function y=cround(prec,x)
if isa(x,'double')
    y=cmulti('set',prec,x);
elseif isa(x,'rmulti')
    y=cmulti('set_r',prec,x.data);
else
end