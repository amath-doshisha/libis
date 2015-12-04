function B=eig_vec_householder_right(A,u,alpha)
[m,n]=size(A);
B=zeros(m,n);
q=A*u;
B=A-alpha*q*u';
