function [A,H,alpha]=eig_mat_hessenberg(A)
[m,n]=size(A);
if m~=n
    disp('Error');
end

H=zeros(n,n-2);
alpha=zeros(n-2,1);

k=1;
while k<=n-2
    [H(:,k),alpha(k)]=eig_vec_householder(A(:,k),k+1);
    A=eig_vec_householder_left(A,H(:,k),alpha(k));
    A=eig_vec_householder_right(A,H(:,k),alpha(k));
    k=k+1;
end