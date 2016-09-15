%function sample4
clear all

%% initializations
N=2;
M=4;
disp(sprintf('Band width: N=%d M=%d',N,M));
lambda=[1 2 3 4 5 6];
disp(sprintf('Prescribed eigenvalues: [%s]',num2str(lambda)));
%prec=53;lambda=double(lambda);
%prec=64;set_default_prec(prec);lambda=multi(lambda);
prec=128;set_default_prec(prec);lambda=multi(lambda);
%prec=256;set_default_prec(prec);lambda=multi(lambda);
disp(sprintf('Precision: %d',prec));
A=sample_matTN(lambda,N,M);
disp('Generated matrix: ');
disp(num2str(A,'%.2f     '));
disp('Compute eigenvalues of generated matrix A on double precision: ');
l=sort(eig(double(A)));
disp(num2str(l,'%.40f'))
disp(sprintf('Compute eigenvalues of A on %d-bit:',prec*64));
set_default_prec(prec*64);
D=eig(multi(A));
disp(num2str(D,'%.40f'))
disp(sprintf('Compute eigenvalues of A, which is rounded to double, on %d-bit:',prec*64));
A=double(A);
set_default_prec(prec*64);
D=eig(multi(A));
disp(num2str(D,'%.40f'))
