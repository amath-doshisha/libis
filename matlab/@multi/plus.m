% z=x+y
function z=plus(x,y)
if isa(x,'multi') && isa(y,'multi')
    z=multi('plus',x.data,y.data);
elseif isa(x,'multi')
    z=multi('plus',x.data,multi(y).data);
elseif isa(y,'multi')
    z=multi('plus',multi(x).data,y.data);
else
    z=multi('plus',multi(x).data,multi(y).data);
end