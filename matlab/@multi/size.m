% d = size(X)
%    [m,n] = size(X)
%   m = size(X,dim)
%    [d1,d2,d3,...,dn] = size(X)
function [d1,d2,d3]=size(x,dim)
s=size(x.data);
if nargin==2
    d1=s(dim);
elseif nargin==1
    if nargout<=1
        d1=s;
    elseif nargout<=2
        d1=s(1);
        d2=s(2);
    elseif nargout<=3
        d1=s(1);
        d2=s(2);
        d3=s(3);
    end
end