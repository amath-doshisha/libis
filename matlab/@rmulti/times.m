function z=times(x,y)
% z=x.*y
if (isa(x,'double') && ~isreal(x)) || (isa(y,'double') && ~isreal(y))
    x=cmulti(x);
    y=cmulti(y);
    z=x.*y;
else
    z=rmulti(2,'times',x,y);
end