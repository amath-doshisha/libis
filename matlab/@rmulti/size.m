function [M,N]=size(A)
if nargout==2
    [M,N]=size(A.data);
else
    M=size(A.data);
end