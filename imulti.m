function y=imulti(x)
if isa(x,'double');
    y=multi('iset_d',x);
elseif isnumeric(x) || islogical(x);
    y=multi('iset_d',double(x));
elseif isa(x,'multi')
    y=multi('icopy',x.data);
else
    y=multi('icopy',multi(x).data);
end