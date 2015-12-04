function B=eig_vec_householder_left(A,u,alpha)
[m,n]=size(A);
B=zeros(m,n);
p=A'*u;
B=A-alpha*u*p';
