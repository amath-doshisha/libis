%function sample4
clear all
close all

N=100;
prec=1024;
set_default_prec(prec);
%A=eye(N,N);
A=rand(N,N)>0.5

% i=1;
% while i<=N
%     j=1;
%     while j<=N
%         if i==j-1
%             A(i,j)=1;
%         elseif i-1==j
%             A(i,j)=-3;
%         end
%         j=j+1;
%     end
%     i=i+1;
% end



disp('Generated matrix: ');
disp(num2str(A,'%.2f     '));
disp('Compute eigenvalues of generated matrix A on double precision: ');
l=sort(eig(double(A)));
disp(num2str(l,'%.40f'))
figure(1);plot(real(l),imag(l),'o')

disp(sprintf('Compute eigenvalues of A on %d-bit:',prec));
D=eig(multi(A));
disp(num2str(D,'%.40f'))
D=double(D);
figure(2);plot(real(D),imag(D),'o',real(l),imag(l),'x')
