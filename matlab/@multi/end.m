%% ind=end(A,k,n)
function ind=end(A,k,n)
if n==1
    ind=numel(A);
elseif k<n
    s=size(A.data);
    ind=s(k);
else
    s=size(A.data);
    ind=prod(s(k:n));
end