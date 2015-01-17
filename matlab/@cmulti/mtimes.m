function z=mtimes(x,y)
% z=x*y
z=cmulti;
if isa(x,'double')
    x=cround(53,x);
end
if isa(y,'double')
    y=cround(53,y);
end
z.data=cmulti_two('mtimes',x.data,y.data);
