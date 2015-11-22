%% ind=end(A,k,n)
function ind=end(A,k,n)
if n==1
    ind=numel(A);
elseif n<=3
    s=size(A.data);
    ind=s(k);
else
    error('multi.end(A,k=%d,n=%d)\n',k,n);
end