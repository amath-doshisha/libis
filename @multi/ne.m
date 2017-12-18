% z=(x~=y)
function z=ne(x,y)
cmd='ne';
if isa(x,'multi') && isa(y,'multi')
    z=logical(multi(cmd,x.data,y.data).data);
elseif isa(x,'multi')
    z=logical(multi(cmd,x.data,multi(y).data).data);
elseif isa(y,'multi')
    z=logical(multi(cmd,multi(x).data,y.data).data);
else
    z=logical(multi(cmd,multi(x).data,multi(y).data).data);
end