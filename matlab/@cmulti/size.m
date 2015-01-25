function [out1,out2]=size(X,dim)
if nargin==1
    if nargout<=1
        %d = size(X)
        out1=size(X.data);
    elseif nargout==2
        %[m,n] = size(X)
        [out1,out2]=size(X.data);
    else
        error('rmulti.size: Too much outputs.');
    end
elseif nargin==2
    if nargout<=1
        %m = size(X,dim)
        out1=size(X.data,dim);
    else
        error('rmulti.size: Too much outputs.');
    end
else
    error('rmulti.size: Too much inputs.');
end
