function z=mtimes(x,y)
% z=x*y
z=rmulti;
if isreal(x)
    x=rround(53,x);
end
if isreal(y)
    y=rround(53,y);
end
z.data=rmulti_two('mtimes',x.data,y.data);
