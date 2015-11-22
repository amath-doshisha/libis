%% ind=end(A,k,n)
function ind=end(obj,k,n)
szd=size(obj.data);
if k<n
    ind=szd(k);
else
    ind=prod(szd(k:n));
end