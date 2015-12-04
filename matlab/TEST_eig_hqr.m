clear all
clear close

n=5;
A=rand(n,n);
disp('Input matrix:');
disp(A);

% QR•ª‰ð
[C,wr,wi]=eig_hqr(A);

disp('Upper triangular matrix:');
disp(C);

disp('Eigenvarue');
lambda=wr+i*wi;
disp(lambda.');